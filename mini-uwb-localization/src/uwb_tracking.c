/**
 * mini-uwb-localization: UWB Tracking Filters Implementation
 * Extended Kalman Filter for UWB tracking with CV/CA/CT motion models.
 * Includes RTS smoother for offline trajectory optimization.
 * L5 Algorithms (EKF), L8 Advanced (RTS smoother)
 */
#include "uwb_tracking.h"
#include <math.h>
#include <string.h>

void ekf_init(ekf_state_t *ekf, const ekf_config_t *config, int dim,
              const uwb_pos3d_t *initial_pos, double dt) {
    int i;
    if (!ekf||!config||!initial_pos) return;
    memset(ekf,0,sizeof(*ekf));
    ekf->state_dim=dim; ekf->dt=dt;
    ekf->motion_model=MOTION_CONSTANT_VELOCITY;
    ekf->turn_rate=0.0;
    ekf->state[0]=initial_pos->x;
    ekf->state[1]=initial_pos->y;
    if (dim>=EKF_STATE_DIM_3D) ekf->state[2]=initial_pos->z;
    for (i=0;i<dim;i++) {
        int pd=(dim==EKF_STATE_DIM_2D)?2:3;
        ekf->P[i*dim+i]=(i<pd)?config->initial_pos_uncertainty:config->initial_vel_uncertainty;
    }
    motion_model_process_noise(dim,dt,config->process_noise_pos,config->process_noise_vel,ekf->Q);
    {int m;for(m=0;m<8;m++)ekf->R[m*8+m]=config->measurement_noise_range;ekf->meas_dim=0;}
    ekf->is_initialized=1;
}

void ekf_predict(ekf_state_t *ekf, double dt) {
    double F[144],x_pred[12],FP[144];
    int i,j,k,n;
    if (!ekf||!ekf->is_initialized) return;
    n=ekf->state_dim;
    motion_model_transition_matrix(n,dt,ekf->motion_model,ekf->turn_rate,F);
    memset(x_pred,0,n*sizeof(double));
    for (i=0;i<n;i++) for (j=0;j<n;j++) x_pred[i]+=F[i*n+j]*ekf->state[j];
    memcpy(ekf->state,x_pred,n*sizeof(double));
    memset(FP,0,sizeof(FP));
    for (i=0;i<n;i++) for (j=0;j<n;j++) for (k=0;k<n;k++) FP[i*n+j]+=F[i*n+k]*ekf->P[k*n+j];
    memset(ekf->P,0,sizeof(ekf->P));
    for (i=0;i<n;i++) for (j=0;j<n;j++) {
        for (k=0;k<n;k++) ekf->P[i*n+j]+=FP[i*n+k]*F[j*n+k];
        ekf->P[i*n+j]+=ekf->Q[i*n+j];
    }
    ekf->dt=dt;
}

void ekf_update_range(ekf_state_t *ekf, const uwb_pos3d_t *anchors,
                      const double *ranges, int num_meas) {
    double H[96],hx[8],nu[8],S[64],Sinv[64],HP[96];
    double K[96],IKH[144],Pnew[144];
    int n,m,i,j,k,pd;
    if (!ekf||!anchors||!ranges||num_meas<1||!ekf->is_initialized) return;
    n=ekf->state_dim; m=(num_meas<8)?num_meas:8;
    pd=(n>=EKF_STATE_DIM_3D)?3:2;
    ekf->meas_dim=m;
    memset(H,0,sizeof(H));
    for (i=0;i<m;i++) {
        double dx=ekf->state[0]-anchors[i].x,dy=ekf->state[1]-anchors[i].y;
        double dz=(pd>=3)?(ekf->state[2]-anchors[i].z):0.0;
        double dist=sqrt(dx*dx+dy*dy+dz*dz);
        if (dist<1e-10) dist=1e-10;
        hx[i]=dist;
        H[i*n+0]=dx/dist; H[i*n+1]=dy/dist;
        if (pd>=3) H[i*n+2]=dz/dist;
    }
    for (i=0;i<m;i++) nu[i]=ranges[i]-hx[i];
    memset(HP,0,sizeof(HP));
    for (i=0;i<m;i++) for (j=0;j<n;j++) for (k=0;k<n;k++) HP[i*n+j]+=H[i*n+k]*ekf->P[k*n+j];
    memset(S,0,sizeof(S));
    for (i=0;i<m;i++) for (j=0;j<m;j++) {
        for (k=0;k<n;k++) S[i*m+j]+=HP[i*n+k]*H[j*n+k];
        S[i*m+j]+=ekf->R[i*8+j];
    }
    /* Invert S */
    memset(Sinv,0,sizeof(Sinv));
    for (i=0;i<m;i++) Sinv[i*m+i]=1.0;
    {
        double aug[144];
        for (i=0;i<m;i++) {
            for (j=0;j<m;j++) aug[i*2*m+j]=S[i*m+j];
            for (j=0;j<m;j++) aug[i*2*m+m+j]=(i==j)?1.0:0.0;
        }
        for (i=0;i<m;i++) {
            int pivot=i; double maxv=fabs(aug[i*2*m+i]);
            for (j=i+1;j<m;j++) if(fabs(aug[j*2*m+i])>maxv){maxv=fabs(aug[j*2*m+i]);pivot=j;}
            if (maxv<1e-12) continue;
            if (pivot!=i) for(j=0;j<2*m;j++){double t=aug[i*2*m+j];aug[i*2*m+j]=aug[pivot*2*m+j];aug[pivot*2*m+j]=t;}
            double piv=aug[i*2*m+i];
            for(j=0;j<2*m;j++) aug[i*2*m+j]/=piv;
            for(j=0;j<m;j++) if(j!=i){double f=aug[j*2*m+i];for(k=0;k<2*m;k++)aug[j*2*m+k]-=f*aug[i*2*m+k];}
        }
        for(i=0;i<m;i++) for(j=0;j<m;j++) Sinv[i*m+j]=aug[i*2*m+m+j];
    }
    /* K = P*H^T*S^{-1} = HP^T * S^{-1} */
    memset(K,0,sizeof(K));
    for(i=0;i<n;i++) for(j=0;j<m;j++) for(k=0;k<m;k++) K[i*n+j]+=HP[k*n+i]*Sinv[k*m+j];
    for(i=0;i<n;i++) {double sum=0.0;for(j=0;j<m;j++)sum+=K[i*n+j]*nu[j];ekf->state[i]+=sum;}
    for(i=0;i<n;i++) for(j=0;j<n;j++) {double kh=0.0;for(k=0;k<m;k++)kh+=K[i*n+k]*H[k*n+j];IKH[i*n+j]=(i==j)?1.0-kh:-kh;}
    memset(Pnew,0,sizeof(Pnew));
    for(i=0;i<n;i++) for(j=0;j<n;j++) for(k=0;k<n;k++) Pnew[i*n+j]+=IKH[i*n+k]*ekf->P[k*n+j];
    memcpy(ekf->P,Pnew,n*n*sizeof(double));
}

void ekf_update_tdoa(ekf_state_t *ekf, const uwb_pos3d_t *anchors,
                     const double *tdoa_values, int num_tdoa) {
    int n,i,j,k,pd;
    if (!ekf||!anchors||!tdoa_values||num_tdoa<1||!ekf->is_initialized) return;
    n=ekf->state_dim; pd=(n>=EKF_STATE_DIM_3D)?3:2;
    for(i=0;i<num_tdoa;i++){
        double dx0=ekf->state[0]-anchors[0].x,dy0=ekf->state[1]-anchors[0].y;
        double dz0=(pd>=3)?(ekf->state[2]-anchors[0].z):0.0;
        double dist0=sqrt(dx0*dx0+dy0*dy0+dz0*dz0);
        double dxi=ekf->state[0]-anchors[i+1].x,dyi=ekf->state[1]-anchors[i+1].y;
        double dzi=(pd>=3)?(ekf->state[2]-anchors[i+1].z):0.0;
        double disti=sqrt(dxi*dxi+dyi*dyi+dzi*dzi);
        if(dist0<1e-10)dist0=1e-10;
         if(disti<1e-10)disti=1e-10;
        double hx=disti-dist0, nu=tdoa_values[i]-hx;
        double H[12]={0};
        H[0]=dxi/disti-dx0/dist0; H[1]=dyi/disti-dy0/dist0;
        if(pd>=3) H[2]=dzi/disti-dz0/dist0;
        double S=ekf->R[i*8+i];
        for(j=0;j<n;j++){double ps=0.0;for(k=0;k<n;k++)ps+=ekf->P[j*n+k]*H[k];S+=H[j]*ps;}
        if(S<1e-15)continue;
        double K[12];
        for(j=0;j<n;j++){double sum=0.0;for(k=0;k<n;k++)sum+=ekf->P[j*n+k]*H[k];K[j]=sum/S;}
        for(j=0;j<n;j++)ekf->state[j]+=K[j]*nu;
        for(j=0;j<n;j++)for(k=0;k<n;k++)ekf->P[j*n+k]-=K[j]*S*K[k];
    }
}

void ekf_get_position(const ekf_state_t *ekf, uwb_pos3d_t *pos) {
    if (!ekf||!pos) return;
    pos->x=ekf->state[0]; pos->y=ekf->state[1];
    pos->z=(ekf->state_dim>=EKF_STATE_DIM_3D)?ekf->state[2]:0.0;
}

void ekf_get_velocity(const ekf_state_t *ekf, uwb_pos3d_t *vel) {
    if (!ekf||!vel) return;
    int off=(ekf->state_dim==EKF_STATE_DIM_2D)?2:3;
    vel->x=ekf->state[off]; vel->y=ekf->state[off+1];
    vel->z=(ekf->state_dim>=EKF_STATE_DIM_3D)?ekf->state[off+2]:0.0;
}

void ekf_get_covariance(const ekf_state_t *ekf, uwb_covariance_t *cov) {
    if (!ekf||!cov) return;
    memset(cov,0,sizeof(*cov));
    int n=ekf->state_dim;
    cov->var_x=ekf->P[0*n+0]; cov->var_y=ekf->P[1*n+1];
    cov->cov_xy=ekf->P[0*n+1];
    if(n>=EKF_STATE_DIM_3D){cov->var_z=ekf->P[2*n+2];cov->cov_xz=ekf->P[0*n+2];cov->cov_yz=ekf->P[1*n+2];}
}

double ekf_normalized_innovation_squared(const ekf_state_t *ekf) {
    if (!ekf) return 0.0;
    double tP=0.0,tR=0.0; int i;
    for(i=0;i<ekf->state_dim;i++) tP+=ekf->P[i*ekf->state_dim+i];
    for(i=0;i<ekf->meas_dim&&i<8;i++) tR+=ekf->R[i*8+i];
    return (tR<1e-15)?0.0:tP/tR;
}

void motion_model_transition_matrix(int state_dim, double dt,
                                    motion_model_t motion, double turn_rate,
                                    double *F_out) {
    int i;
    if(!F_out)return;
    memset(F_out,0,state_dim*state_dim*sizeof(double));
    for(i=0;i<state_dim;i++)F_out[i*state_dim+i]=1.0;
    if(motion==MOTION_STATIC)return;
    if(state_dim==EKF_STATE_DIM_2D){F_out[0*4+2]=dt;F_out[1*4+3]=dt;}
    else if(state_dim==EKF_STATE_DIM_3D){F_out[0*6+3]=dt;F_out[1*6+4]=dt;F_out[2*6+5]=dt;}
    else if(state_dim==EKF_STATE_DIM_3D_ACC){
        double dt2=0.5*dt*dt;
        F_out[0*9+3]=dt;F_out[0*9+6]=dt2;F_out[1*9+4]=dt;F_out[1*9+7]=dt2;
        F_out[2*9+5]=dt;F_out[2*9+8]=dt2;F_out[3*9+6]=dt;F_out[4*9+7]=dt;F_out[5*9+8]=dt;
    }
    if(motion==MOTION_COORDINATED_TURN&&state_dim==EKF_STATE_DIM_2D&&fabs(turn_rate)>1e-10){
        double w=turn_rate,sw=sin(w*dt),cw=cos(w*dt);
        F_out[0*4+0]=cw;F_out[0*4+1]=-sw;F_out[1*4+0]=sw;F_out[1*4+1]=cw;
    }
}

void motion_model_process_noise(int state_dim, double dt,
                                double pos_noise, double vel_noise,
                                double *Q_out) {
    int i,pdim;
    if(!Q_out)return;
    memset(Q_out,0,state_dim*state_dim*sizeof(double));
    pdim=(state_dim==EKF_STATE_DIM_2D)?2:3;
    double dt2=dt*dt,dt3=dt*dt*dt;
    for(i=0;i<pdim;i++){
        int pi=i,vi=i+pdim;
        Q_out[pi*state_dim+pi]=dt3*pos_noise/3.0;
        Q_out[pi*state_dim+vi]=dt2*pos_noise/2.0;
        Q_out[vi*state_dim+pi]=dt2*pos_noise/2.0;
        Q_out[vi*state_dim+vi]=dt*vel_noise;
    }
}

void dead_reckon_position(uwb_pos3d_t *pos, const uwb_pos3d_t *vel, double dt) {
    if(!pos||!vel)return;
    pos->x+=vel->x*dt;pos->y+=vel->y*dt;pos->z+=vel->z*dt;
}

/*
 * RTS Smoother (L8 Advanced): backward recursion
 * C_k = P_k*F_k^T*inv(P_{k+1|k})
 * x_k^s = x_k + C_k*(x_{k+1}^s - x_{k+1|k})
 * Reference: Rauch, Tung, Striebel (1965) AIAA Journal
 */
void rts_smoother(ekf_state_t *states, int num_states, double dt) {
    int k,n,i,j,l;
    double F[144],x_pred[12],P_pred[144],PFT[144],C[144];
    if(!states||num_states<2)return;
    n=states[0].state_dim;
    for(k=num_states-2;k>=0;k--){
        motion_model_transition_matrix(n,dt,states[k].motion_model,states[k].turn_rate,F);
        memset(x_pred,0,sizeof(x_pred));
        for(i=0;i<n;i++)for(j=0;j<n;j++)x_pred[i]+=F[i*n+j]*states[k].state[j];
        {
            double FP[144]; memset(FP,0,sizeof(FP));
            for(i=0;i<n;i++)for(j=0;j<n;j++)for(l=0;l<n;l++)FP[i*n+j]+=F[i*n+l]*states[k].P[l*n+j];
            memset(P_pred,0,sizeof(P_pred));
            for(i=0;i<n;i++)for(j=0;j<n;j++){
                for(l=0;l<n;l++)P_pred[i*n+j]+=FP[i*n+l]*F[j*n+l];
                P_pred[i*n+j]+=states[k].Q[i*n+j];
            }
        }
        memset(PFT,0,sizeof(PFT));
        for(i=0;i<n;i++)for(j=0;j<n;j++)for(l=0;l<n;l++)PFT[i*n+j]+=states[k].P[i*n+l]*F[j*n+l];
        double P_inv[144]; memset(P_inv,0,sizeof(P_inv));
        {
            double aug[288];
            for(i=0;i<n;i++){for(j=0;j<n;j++)aug[i*2*n+j]=P_pred[i*n+j];for(j=0;j<n;j++)aug[i*2*n+n+j]=(i==j)?1.0:0.0;}
            for(i=0;i<n;i++){
                int pivot=i;double maxv=fabs(aug[i*2*n+i]);
                for(j=i+1;j<n;j++)if(fabs(aug[j*2*n+i])>maxv){maxv=fabs(aug[j*2*n+i]);pivot=j;}
                if(maxv<1e-12)continue;
                if(pivot!=i)for(j=0;j<2*n;j++){double t=aug[i*2*n+j];aug[i*2*n+j]=aug[pivot*2*n+j];aug[pivot*2*n+j]=t;}
                double piv=aug[i*2*n+i];
                for(j=0;j<2*n;j++)aug[i*2*n+j]/=piv;
                for(j=0;j<n;j++)if(j!=i){double f=aug[j*2*n+i];for(l=0;l<2*n;l++)aug[j*2*n+l]-=f*aug[i*2*n+l];}
            }
            for(i=0;i<n;i++)for(j=0;j<n;j++)P_inv[i*n+j]=aug[i*2*n+n+j];
        }
        memset(C,0,sizeof(C));
        for(i=0;i<n;i++)for(j=0;j<n;j++)for(l=0;l<n;l++)C[i*n+j]+=PFT[i*n+l]*P_inv[l*n+j];
        {
            double diff[12];
            for(i=0;i<n;i++)diff[i]=states[k+1].state[i]-x_pred[i];
            for(i=0;i<n;i++){double sum=0.0;for(j=0;j<n;j++)sum+=C[i*n+j]*diff[j];states[k].state[i]+=sum;}
        }
        {
            double Pd[144],CPd[144],CPdCT[144];
            for(i=0;i<n;i++)for(j=0;j<n;j++)Pd[i*n+j]=states[k+1].P[i*n+j]-P_pred[i*n+j];
            memset(CPd,0,sizeof(CPd));
            for(i=0;i<n;i++)for(j=0;j<n;j++)for(l=0;l<n;l++)CPd[i*n+j]+=C[i*n+l]*Pd[l*n+j];
            memset(CPdCT,0,sizeof(CPdCT));
            for(i=0;i<n;i++)for(j=0;j<n;j++)for(l=0;l<n;l++)CPdCT[i*n+j]+=CPd[i*n+l]*C[j*n+l];
            for(i=0;i<n;i++)for(j=0;j<n;j++)states[k].P[i*n+j]+=CPdCT[i*n+j];
        }
    }
}
