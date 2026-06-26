/**
 * mini-uwb-localization: UWB Positioning Implementation
 *
 * Implements multilateration positioning algorithms.
 * L3 Math Structures, L4 Fundamental Laws (GDOP, CRLB, FIM)
 * L5 Algorithms (LS, WLS, GN, LM, Chan, Taylor)
 * L6 Canonical Problems (Trilateration, Multilateration, TDoA)
 */

#include "uwb_positioning.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

static int solve_2x2(const double A[4], const double b[2], double x[2]) {
    double det = A[0]*A[3] - A[1]*A[2];
    if (fabs(det) < 1e-15) return 0;
    x[0] = (b[0]*A[3] - A[1]*b[1]) / det;
    x[1] = (A[0]*b[1] - b[0]*A[2]) / det;
    return 1;
}

static int solve_3x3(double A[9], double b[3], double x[3]) {
    int i,j,k,max_row; double max_val,tmp,factor,Ab[12];
    for (i=0;i<3;i++) {
        Ab[i*4+0]=A[i*3+0]; Ab[i*4+1]=A[i*3+1];
        Ab[i*4+2]=A[i*3+2]; Ab[i*4+3]=b[i];
    }
    for (i=0;i<3;i++) {
        max_row=i; max_val=fabs(Ab[i*4+i]);
        for (j=i+1;j<3;j++) if (fabs(Ab[j*4+i])>max_val) {max_val=fabs(Ab[j*4+i]);max_row=j;}
        if (max_val<1e-15) return 0;
        if (max_row!=i) for (k=0;k<4;k++) {tmp=Ab[i*4+k];Ab[i*4+k]=Ab[max_row*4+k];Ab[max_row*4+k]=tmp;}
        for (j=i+1;j<3;j++) {factor=Ab[j*4+i]/Ab[i*4+i]; for (k=i;k<4;k++) Ab[j*4+k]-=factor*Ab[i*4+k];}
    }
    for (i=2;i>=0;i--) {x[i]=Ab[i*4+3];for(j=i+1;j<3;j++)x[i]-=Ab[i*4+j]*x[j]; x[i]/=Ab[i*4+i];}
    return 1;
}

static int invert_2x2(const double A[4], double A_inv[4]) {
    double det = A[0]*A[3]-A[1]*A[2];
    if (fabs(det)<1e-15) return 0;
    A_inv[0]=A[3]/det; A_inv[1]=-A[1]/det; A_inv[2]=-A[2]/det; A_inv[3]=A[0]/det;
    return 1;
}

/*
 * 2D Trilateration (L6): 3 anchors -> closed-form position
 */
int trilateration_2d(const uwb_pos3d_t *anchors, const double *ranges,
                     uwb_pos2d_t *result) {
    double A[4],b[2],sol[2],x1,y1,r1;
    if (!anchors||!ranges||!result) return 0;
    x1=anchors[0].x; y1=anchors[0].y; r1=ranges[0];
    A[0]=2.0*(anchors[1].x-x1); A[1]=2.0*(anchors[1].y-y1);
    A[2]=2.0*(anchors[2].x-x1); A[3]=2.0*(anchors[2].y-y1);
    b[0]=anchors[1].x*anchors[1].x+anchors[1].y*anchors[1].y-x1*x1-y1*y1-ranges[1]*ranges[1]+r1*r1;
    b[1]=anchors[2].x*anchors[2].x+anchors[2].y*anchors[2].y-x1*x1-y1*y1-ranges[2]*ranges[2]+r1*r1;
    if (!solve_2x2(A,b,sol)) return 0;
    result->x=sol[0]; result->y=sol[1];
    return 1;
}

/*
 * 3D Trilateration (L6): 4 anchors -> closed-form position
 */
int trilateration_3d(const uwb_pos3d_t *anchors, const double *ranges,
                     uwb_pos3d_t *result) {
    double A[9],b[3],sol[3],x1,y1,z1,r1;
    int i;
    if (!anchors||!ranges||!result) return 0;
    x1=anchors[0].x; y1=anchors[0].y; z1=anchors[0].z; r1=ranges[0];
    for (i=1;i<=3;i++) {
        int row=i-1;
        A[row*3+0]=2.0*(anchors[i].x-x1); A[row*3+1]=2.0*(anchors[i].y-y1); A[row*3+2]=2.0*(anchors[i].z-z1);
        b[row]=anchors[i].x*anchors[i].x+anchors[i].y*anchors[i].y+anchors[i].z*anchors[i].z
               -x1*x1-y1*y1-z1*z1-ranges[i]*ranges[i]+r1*r1;
    }
    if (!solve_3x3(A,b,sol)) return 0;
    result->x=sol[0]; result->y=sol[1]; result->z=sol[2];
    return 1;
}

/*
 * Linear Least Squares Multilateration (L5): min ||Ap-b||_2
 */
int multilateration_linear_ls(const uwb_pos3d_t *anchors, const double *ranges,
                              int num_anchors, uwb_positioning_dim_t dim,
                              uwb_positioning_result_t *result) {
    int D=(dim==POS_DIM_3D)?3:2;
    double AtA[9],Atb[3];
    int i,j,k;
    if (!anchors||!ranges||!result||num_anchors<D+1) return 0;
    memset(result,0,sizeof(*result));
    result->solver_used=POS_SOLVER_LINEAR_LS;
    memset(AtA,0,sizeof(AtA)); memset(Atb,0,sizeof(Atb));
    double x1=anchors[0].x,y1=anchors[0].y,z1=anchors[0].z,r1=ranges[0];
    for (i=1;i<num_anchors;i++) {
        double ai[3],bi;
        ai[0]=2.0*(anchors[i].x-x1); ai[1]=2.0*(anchors[i].y-y1);
        if (D==3) ai[2]=2.0*(anchors[i].z-z1);
        bi=r1*r1-ranges[i]*ranges[i]+anchors[i].x*anchors[i].x+anchors[i].y*anchors[i].y
           +anchors[i].z*anchors[i].z-x1*x1-y1*y1-z1*z1;
        for (j=0;j<D;j++) {for(k=0;k<D;k++) AtA[j*D+k]+=ai[j]*ai[k]; Atb[j]+=ai[j]*bi;}
    }
    if (D==2) {
        double ATA_inv[4];
        if (!invert_2x2(AtA,ATA_inv)) return 0;
        result->position.x=ATA_inv[0]*Atb[0]+ATA_inv[1]*Atb[1];
        result->position.y=ATA_inv[2]*Atb[0]+ATA_inv[3]*Atb[1];
        result->position.z=0.0;
    } else {
        double sol[3];
        if (!solve_3x3(AtA,Atb,sol)) return 0;
        result->position.x=sol[0]; result->position.y=sol[1]; result->position.z=sol[2];
    }
    double residual=0.0;
    for (i=0;i<num_anchors;i++) {
        double dx=result->position.x-anchors[i].x,dy=result->position.y-anchors[i].y,dz=result->position.z-anchors[i].z;
        double err=sqrt(dx*dx+dy*dy+dz*dz)-ranges[i]; residual+=err*err;
    }
    result->residual_norm=sqrt(residual);
    result->converged=1; result->num_iterations=1;
    compute_dop_metrics(anchors,num_anchors,&result->position,result);
    return 1;
}

/*
 * Weighted Least Squares (L5): p = (A^T W A)^-1 A^T W b
 */
int multilateration_weighted_ls(const uwb_pos3d_t *anchors, const double *ranges,
                                const double *weights, int num_anchors,
                                uwb_positioning_dim_t dim,
                                uwb_positioning_result_t *result) {
    int D=(dim==POS_DIM_3D)?3:2;
    double AtWA[9],AtWb[3];
    int i,j,k;
    if (!anchors||!ranges||!weights||!result||num_anchors<D+1) return 0;
    memset(result,0,sizeof(*result));
    result->solver_used=POS_SOLVER_WEIGHTED_LS;
    memset(AtWA,0,sizeof(AtWA)); memset(AtWb,0,sizeof(AtWb));
    double x1=anchors[0].x,y1=anchors[0].y,z1=anchors[0].z,r1=ranges[0];
    for (i=1;i<num_anchors;i++) {
        double ai[3],bi,w=weights[i];
        ai[0]=2.0*(anchors[i].x-x1); ai[1]=2.0*(anchors[i].y-y1);
        if (D==3) ai[2]=2.0*(anchors[i].z-z1);
        bi=r1*r1-ranges[i]*ranges[i]+anchors[i].x*anchors[i].x+anchors[i].y*anchors[i].y
           +anchors[i].z*anchors[i].z-x1*x1-y1*y1-z1*z1;
        for (j=0;j<D;j++) {for(k=0;k<D;k++) AtWA[j*D+k]+=w*ai[j]*ai[k]; AtWb[j]+=w*ai[j]*bi;}
    }
    if (D==2) {
        double AtWA_inv[4];
        if (!invert_2x2(AtWA,AtWA_inv)) return 0;
        result->position.x=AtWA_inv[0]*AtWb[0]+AtWA_inv[1]*AtWb[1];
        result->position.y=AtWA_inv[2]*AtWb[0]+AtWA_inv[3]*AtWb[1];
        result->position.z=0.0;
    } else {
        double sol[3];
        if (!solve_3x3(AtWA,AtWb,sol)) return 0;
        result->position.x=sol[0]; result->position.y=sol[1]; result->position.z=sol[2];
    }
    result->converged=1; result->num_iterations=1;
    compute_dop_metrics(anchors,num_anchors,&result->position,result);
    return 1;
}

/*
 * Gauss-Newton iteration (L5): J_k^T J_k * delta = -J_k^T * r
 * Reference: Nocedal & Wright (2006)
 */
int multilateration_gauss_newton(const uwb_pos3d_t *anchors, const double *ranges,
                                 int num_anchors, uwb_positioning_dim_t dim,
                                 uwb_positioning_result_t *result) {
    int D=(dim==POS_DIM_3D)?3:2;
    double p[3],J[60],r[20],JTJ[9],JTr[3];
    int iter,i,j,k;
    const int max_iter=50;
    const double tol=1e-6;
    if (!anchors||!ranges||!result||num_anchors<D+1) return 0;
    memset(result,0,sizeof(*result));
    result->solver_used=POS_SOLVER_GAUSS_NEWTON;
    p[0]=result->position.x; p[1]=result->position.y; p[2]=result->position.z;
    if (fabs(p[0])<1e-10&&fabs(p[1])<1e-10&&fabs(p[2])<1e-10) {
        p[0]=p[1]=p[2]=0.0;
        for (i=0;i<num_anchors;i++) {p[0]+=anchors[i].x; p[1]+=anchors[i].y; p[2]+=anchors[i].z;}
        p[0]/=num_anchors; p[1]/=num_anchors; p[2]/=num_anchors;
    }
    for (iter=0;iter<max_iter;iter++) {
        for (i=0;i<num_anchors;i++) {
            double dx=p[0]-anchors[i].x,dy=p[1]-anchors[i].y,dz=p[2]-anchors[i].z;
            double dist=sqrt(dx*dx+dy*dy+dz*dz);
            if (dist<1e-10) dist=1e-10;
            r[i]=dist-ranges[i];
            J[i*D+0]=dx/dist; J[i*D+1]=dy/dist;
            if (D==3) J[i*D+2]=dz/dist;
        }
        memset(JTJ,0,sizeof(JTJ)); memset(JTr,0,sizeof(JTr));
        for (i=0;i<num_anchors;i++) {
            for (j=0;j<D;j++) {for(k=0;k<D;k++) JTJ[j*D+k]+=J[i*D+j]*J[i*D+k]; JTr[j]+=J[i*D+j]*r[i];}
        }
        double delta[3]={0};
        if (D==2) {
            double JTJ_inv[4];
            if (!invert_2x2(JTJ,JTJ_inv)) break;
            delta[0]=-(JTJ_inv[0]*JTr[0]+JTJ_inv[1]*JTr[1]);
            delta[1]=-(JTJ_inv[2]*JTr[0]+JTJ_inv[3]*JTr[1]);
        } else {
            double neg_JTr[3]={-JTr[0],-JTr[1],-JTr[2]},sol[3];
            if (!solve_3x3(JTJ,neg_JTr,sol)) break;
            delta[0]=sol[0]; delta[1]=sol[1]; delta[2]=sol[2];
        }
        p[0]+=delta[0]; p[1]+=delta[1]; p[2]+=delta[2];
        if (sqrt(delta[0]*delta[0]+delta[1]*delta[1]+delta[2]*delta[2])<tol) {result->converged=1;break;}
    }
    result->position.x=p[0]; result->position.y=p[1]; result->position.z=p[2];
    result->num_iterations=iter+1;
    double residual=0.0;
    for (i=0;i<num_anchors;i++) {
        double dx=p[0]-anchors[i].x,dy=p[1]-anchors[i].y,dz=p[2]-anchors[i].z;
        double diff=sqrt(dx*dx+dy*dy+dz*dz)-ranges[i]; residual+=diff*diff;
    }
    result->residual_norm=sqrt(residual);
    compute_dop_metrics(anchors,num_anchors,&result->position,result);
    return result->converged;
}

/*
 * Levenberg-Marquardt (L5): (J^T J + lambda*I)*delta = -J^T r
 * Reference: Marquardt (1963) SIAM Journal
 */
int multilateration_levenberg_marquardt(const uwb_pos3d_t *anchors,
                                        const double *ranges, int num_anchors,
                                        uwb_positioning_dim_t dim,
                                        uwb_positioning_result_t *result) {
    int D=(dim==POS_DIM_3D)?3:2;
    double p[3],p_try[3],J[60],r[20],r_try[20],JTJ[9],JTr[3];
    double lambda=0.01,cost,cost_try;
    int iter,i,j,k;
    const int max_iter=100;
    const double tol=1e-8;
    if (!anchors||!ranges||!result||num_anchors<D+1) return 0;
    memset(result,0,sizeof(*result));
    result->solver_used=POS_SOLVER_LEVENBERG_MARQUARDT;
    p[0]=result->position.x; p[1]=result->position.y; p[2]=result->position.z;
    if (fabs(p[0])<1e-10&&fabs(p[1])<1e-10&&fabs(p[2])<1e-10) {
        p[0]=p[1]=p[2]=0.0;
        for (i=0;i<num_anchors;i++) {p[0]+=anchors[i].x;p[1]+=anchors[i].y;p[2]+=anchors[i].z;}
        p[0]/=num_anchors;p[1]/=num_anchors;p[2]/=num_anchors;
    }
    cost=0.0;
    for (i=0;i<num_anchors;i++) {
        double dx=p[0]-anchors[i].x,dy=p[1]-anchors[i].y,dz=p[2]-anchors[i].z;
        r[i]=sqrt(dx*dx+dy*dy+dz*dz)-ranges[i]; cost+=r[i]*r[i];
    }
    for (iter=0;iter<max_iter;iter++) {
        for (i=0;i<num_anchors;i++) {
            double dx=p[0]-anchors[i].x,dy=p[1]-anchors[i].y,dz=p[2]-anchors[i].z;
            double dist=sqrt(dx*dx+dy*dy+dz*dz);
            if (dist<1e-10) dist=1e-10;
            J[i*D+0]=dx/dist; J[i*D+1]=dy/dist;
            if (D==3) J[i*D+2]=dz/dist;
        }
        memset(JTJ,0,sizeof(JTJ)); memset(JTr,0,sizeof(JTr));
        for (i=0;i<num_anchors;i++) {
            for (j=0;j<D;j++) {for(k=0;k<D;k++) JTJ[j*D+k]+=J[i*D+j]*J[i*D+k]; JTr[j]+=J[i*D+j]*r[i];}
        }
        double lhs[9],rhs[3];
        memcpy(lhs,JTJ,sizeof(JTJ));
        for (j=0;j<D;j++) lhs[j*D+j]+=lambda*JTJ[j*D+j];
        for (j=0;j<D;j++) rhs[j]=-JTr[j];
        double delta[3]={0};
        int solved=0;
        if (D==2) {double lhs_inv[4]; if(invert_2x2(lhs,lhs_inv)){delta[0]=lhs_inv[0]*rhs[0]+lhs_inv[1]*rhs[1];delta[1]=lhs_inv[2]*rhs[0]+lhs_inv[3]*rhs[1];solved=1;}}
        else solved=solve_3x3(lhs,rhs,delta);
        if (!solved) {lambda*=4.0; continue;}
        for (j=0;j<3;j++) p_try[j]=p[j]+delta[j];
        cost_try=0.0;
        for (i=0;i<num_anchors;i++) {
            double dx=p_try[0]-anchors[i].x,dy=p_try[1]-anchors[i].y,dz=p_try[2]-anchors[i].z;
            r_try[i]=sqrt(dx*dx+dy*dy+dz*dz)-ranges[i]; cost_try+=r_try[i]*r_try[i];
        }
        double rho_num=cost-cost_try,rho_den=0.0;
        for (j=0;j<D;j++) rho_den+=delta[j]*(lambda*JTJ[j*D+j]*delta[j]-JTr[j]);
        double gain=(rho_den>1e-15)?rho_num/rho_den:0.0;
        if (gain>0.0) {
            for (j=0;j<3;j++) p[j]=p_try[j];
            for (i=0;i<num_anchors;i++) r[i]=r_try[i];
            cost=cost_try;
            lambda*=fmax(1.0/3.0,1.0-pow(2.0*gain-1.0,3.0));
            if (lambda<1e-7) lambda=1e-7;
            if (sqrt(delta[0]*delta[0]+delta[1]*delta[1]+delta[2]*delta[2])<tol) {result->converged=1;break;}
        } else {lambda*=2.0; if(lambda>1e10) break;}
    }
    result->position.x=p[0]; result->position.y=p[1]; result->position.z=p[2];
    result->residual_norm=sqrt(cost); result->num_iterations=iter+1;
    compute_dop_metrics(anchors,num_anchors,&result->position,result);
    return result->converged;
}

/*
 * Chan TDoA Algorithm (L5): Two-step WLS hyperbolic positioning.
 * Reference: Chan & Ho (1994) IEEE TSP
 */
int tdoa_positioning_chan(const uwb_pos3d_t *anchors, const double *tdoa_values,
                          int num_anchors, uwb_positioning_dim_t dim,
                          uwb_positioning_result_t *result) {
    int D=(dim==POS_DIM_3D)?3:2, M=num_anchors-1, i;
    double ri1[20],x1,y1,z1,Ga[60],h[20],result_pos[3]={0};
    if (!anchors||!tdoa_values||!result||num_anchors<D+2) return 0;
    memset(result,0,sizeof(*result));
    result->solver_used=POS_SOLVER_CHAN_TDOA;
    x1=anchors[0].x; y1=anchors[0].y; z1=anchors[0].z;
    for (i=0;i<M;i++) ri1[i]=UWB_C*tdoa_values[i];
    for (i=0;i<M;i++) {
        int idx=i+1;
        Ga[i*(D+1)+0]=anchors[idx].x-x1; Ga[i*(D+1)+1]=anchors[idx].y-y1;
        if (D==3) Ga[i*(D+1)+2]=anchors[idx].z-z1;
        Ga[i*(D+1)+D]=ri1[i];
        double Ki=anchors[idx].x*anchors[idx].x+anchors[idx].y*anchors[idx].y+anchors[idx].z*anchors[idx].z;
        double K1=x1*x1+y1*y1+z1*z1;
        h[i]=0.5*(ri1[i]*ri1[i]-Ki+K1);
    }
    double GaTGa[16],GaTh[4];
    int cols=D+1;
    memset(GaTGa,0,sizeof(GaTGa)); memset(GaTh,0,sizeof(GaTh));
    for (i=0;i<M;i++) {
        int j,k;
        for (j=0;j<cols;j++) {for(k=0;k<cols;k++) GaTGa[j*cols+k]+=Ga[i*cols+j]*Ga[i*cols+k]; GaTh[j]+=Ga[i*cols+j]*h[i];}
    }
    double augmented[20],za[4];
    int n=cols;
    for (i=0;i<n;i++) {int j;for(j=0;j<n;j++)augmented[i*(n+1)+j]=GaTGa[i*n+j]; augmented[i*(n+1)+n]=GaTh[i];}
    for (i=0;i<n;i++) {
        int j,pivot_row=i;
        double pivot_val=fabs(augmented[i*(n+1)+i]);
        for (j=i+1;j<n;j++) if(fabs(augmented[j*(n+1)+i])>pivot_val){pivot_val=fabs(augmented[j*(n+1)+i]);pivot_row=j;}
        if (pivot_val<1e-15) return 0;
        if (pivot_row!=i) for(j=0;j<=n;j++){double tmp=augmented[i*(n+1)+j];augmented[i*(n+1)+j]=augmented[pivot_row*(n+1)+j];augmented[pivot_row*(n+1)+j]=tmp;}
        for (j=i+1;j<n;j++) {double factor=augmented[j*(n+1)+i]/augmented[i*(n+1)+i];int k;for(k=i;k<=n;k++)augmented[j*(n+1)+k]-=factor*augmented[i*(n+1)+k];}
    }
    for(i=n-1;i>=0;i--){int j;za[i]=augmented[i*(n+1)+n];for(j=i+1;j<n;j++)za[i]-=augmented[i*(n+1)+j]*za[j];za[i]/=augmented[i*(n+1)+i];}
    result_pos[0]=za[0]+x1; result_pos[1]=za[1]+y1;
    if (D==3) result_pos[2]=za[2]+z1;
    result->position.x=result_pos[0]; result->position.y=result_pos[1];
    result->position.z=(D==3)?result_pos[2]:0.0;
    result->converged=1; result->num_iterations=2;
    compute_dop_metrics(anchors,num_anchors,&result->position,result);
    return 1;
}

/*
 * Taylor-series TDoA expansion (L5). Reference: Foy (1976) IEEE AES
 */
int tdoa_positioning_taylor(const uwb_pos3d_t *anchors, const double *tdoa_values,
                            int num_anchors, uwb_positioning_dim_t dim,
                            uwb_positioning_result_t *result) {
    int D=(dim==POS_DIM_3D)?3:2, M=num_anchors-1;
    double p[3],J[60],r_tdoa[20],JTJ[9],JTr[3];
    int iter,i,j,k;
    const int max_iter=30;
    const double tol=1e-6;
    if (!anchors||!tdoa_values||!result||num_anchors<D+2) return 0;
    memset(result,0,sizeof(*result));
    result->solver_used=POS_SOLVER_TAYLOR_TDOA;
    p[0]=result->position.x; p[1]=result->position.y; p[2]=result->position.z;
    if (fabs(p[0])<1e-10&&fabs(p[1])<1e-10) {
        for(i=0;i<num_anchors;i++){p[0]+=anchors[i].x;p[1]+=anchors[i].y;p[2]+=anchors[i].z;}
        p[0]/=num_anchors;p[1]/=num_anchors;p[2]/=num_anchors;
    }
    for (iter=0;iter<max_iter;iter++) {
        for (i=0;i<M;i++) {
            double dx0=p[0]-anchors[0].x,dy0=p[1]-anchors[0].y,dz0=p[2]-anchors[0].z;
            double dist0=sqrt(dx0*dx0+dy0*dy0+dz0*dz0);
            double dxi=p[0]-anchors[i+1].x,dyi=p[1]-anchors[i+1].y,dzi=p[2]-anchors[i+1].z;
            double disti=sqrt(dxi*dxi+dyi*dyi+dzi*dzi);
            double pred=disti-dist0, meas=UWB_C*tdoa_values[i];
            r_tdoa[i]=meas-pred;
            if (dist0<1e-10) dist0=1e-10;
            if (disti<1e-10) disti=1e-10;
            J[i*D+0]=dxi/disti-dx0/dist0; J[i*D+1]=dyi/disti-dy0/dist0;
            if (D==3) J[i*D+2]=dzi/disti-dz0/dist0;
        }
        memset(JTJ,0,sizeof(JTJ)); memset(JTr,0,sizeof(JTr));
        for (i=0;i<M;i++) {
            for(j=0;j<D;j++) {for(k=0;k<D;k++) JTJ[j*D+k]+=J[i*D+j]*J[i*D+k]; JTr[j]+=J[i*D+j]*r_tdoa[i];}
        }
        double delta[3]={0};
        int solved=0;
        if (D==2) {double JTJ_inv[4]; if(invert_2x2(JTJ,JTJ_inv)){delta[0]=JTJ_inv[0]*JTr[0]+JTJ_inv[1]*JTr[1];delta[1]=JTJ_inv[2]*JTr[0]+JTJ_inv[3]*JTr[1];solved=1;}}
        else solved=solve_3x3(JTJ,JTr,delta);
        if (!solved) break;
        p[0]+=delta[0]; p[1]+=delta[1]; p[2]+=delta[2];
        if (sqrt(delta[0]*delta[0]+delta[1]*delta[1]+delta[2]*delta[2])<tol) {result->converged=1;break;}
    }
    result->position.x=p[0]; result->position.y=p[1]; result->position.z=p[2];
    result->num_iterations=iter+1;
    compute_dop_metrics(anchors,num_anchors,&result->position,result);
    return result->converged;
}

/*
 * DOP Metrics (L4): G=(H^T H)^(-1)
 * GDOP=sqrt(trace(G)), PDOP=sqrt(Gxx+Gyy+Gzz)
 * HDOP=sqrt(Gxx+Gyy), VDOP=sqrt(Gzz)
 */
void compute_dop_metrics(const uwb_pos3d_t *anchors, int num_anchors,
                         const uwb_pos3d_t *tag_pos,
                         uwb_positioning_result_t *result) {
    double HTH[9],HTH_inv[9],det;
    int i;
    if (!anchors||!tag_pos||!result||num_anchors<3) return;
    memset(HTH,0,sizeof(HTH));
    for (i=0;i<num_anchors;i++) {
        double dx=tag_pos->x-anchors[i].x,dy=tag_pos->y-anchors[i].y,dz=tag_pos->z-anchors[i].z;
        double range=sqrt(dx*dx+dy*dy+dz*dz);
        if (range<1e-6) continue;
        double ux=dx/range,uy=dy/range,uz=dz/range;
        HTH[0]+=ux*ux; HTH[1]+=ux*uy; HTH[2]+=ux*uz;
        HTH[3]+=uy*ux; HTH[4]+=uy*uy; HTH[5]+=uy*uz;
        HTH[6]+=uz*ux; HTH[7]+=uz*uy; HTH[8]+=uz*uz;
    }
    det=HTH[0]*(HTH[4]*HTH[8]-HTH[5]*HTH[7])-HTH[1]*(HTH[3]*HTH[8]-HTH[5]*HTH[6])+HTH[2]*(HTH[3]*HTH[7]-HTH[4]*HTH[6]);
    if (fabs(det)<1e-15) {result->gdop=result->pdop=result->hdop=result->vdop=99.0; return;}
    HTH_inv[0]=(HTH[4]*HTH[8]-HTH[5]*HTH[7])/det; HTH_inv[1]=-(HTH[1]*HTH[8]-HTH[2]*HTH[7])/det;
    HTH_inv[2]=(HTH[1]*HTH[5]-HTH[2]*HTH[4])/det; HTH_inv[3]=-(HTH[3]*HTH[8]-HTH[5]*HTH[6])/det;
    HTH_inv[4]=(HTH[0]*HTH[8]-HTH[2]*HTH[6])/det; HTH_inv[5]=-(HTH[0]*HTH[5]-HTH[2]*HTH[3])/det;
    HTH_inv[6]=(HTH[3]*HTH[7]-HTH[4]*HTH[6])/det; HTH_inv[7]=-(HTH[0]*HTH[7]-HTH[1]*HTH[6])/det;
    HTH_inv[8]=(HTH[0]*HTH[4]-HTH[1]*HTH[3])/det;
    result->gdop=sqrt(HTH_inv[0]+HTH_inv[4]+HTH_inv[8]);
    result->pdop=sqrt(HTH_inv[0]+HTH_inv[4]+HTH_inv[8]);
    result->hdop=sqrt(HTH_inv[0]+HTH_inv[4]);
    result->vdop=sqrt(HTH_inv[8]);
    result->fisher_trace=HTH[0]+HTH[4]+HTH[8];
}

double evaluate_anchor_geometry(const uwb_pos3d_t *anchors, int num_anchors,
                                const uwb_pos3d_t *tag_pos) {
    uwb_positioning_result_t tmp;
    compute_dop_metrics(anchors,num_anchors,tag_pos,&tmp);
    if (tmp.gdop>10.0) return 0.0;
    if (tmp.gdop<1.0) return 1.0;
    return 1.0-(tmp.gdop-1.0)/9.0;
}

/*
 * CRLB for Position (L4): sqrt(trace(F^(-1)))
 * F = (1/sigma^2) * sum_i u_i * u_i^T
 */
double crlb_position(const uwb_pos3d_t *anchors, int num_anchors,
                     const uwb_pos3d_t *tag_pos, double range_variance) {
    double HTH[9],HTH_inv[9],det;
    int i;
    if (!anchors||!tag_pos||num_anchors<3||range_variance<=0.0) return 1e100;
    memset(HTH,0,sizeof(HTH));
    for (i=0;i<num_anchors;i++) {
        double dx=tag_pos->x-anchors[i].x,dy=tag_pos->y-anchors[i].y,dz=tag_pos->z-anchors[i].z;
        double range=sqrt(dx*dx+dy*dy+dz*dz);
        if (range<1e-6) continue;
        double ux=dx/range,uy=dy/range,uz=dz/range;
        HTH[0]+=ux*ux;HTH[1]+=ux*uy;HTH[2]+=ux*uz;
        HTH[3]+=uy*ux;HTH[4]+=uy*uy;HTH[5]+=uy*uz;
        HTH[6]+=uz*ux;HTH[7]+=uz*uy;HTH[8]+=uz*uz;
    }
    det=HTH[0]*(HTH[4]*HTH[8]-HTH[5]*HTH[7])-HTH[1]*(HTH[3]*HTH[8]-HTH[5]*HTH[6])+HTH[2]*(HTH[3]*HTH[7]-HTH[4]*HTH[6]);
    if (fabs(det)<1e-15) return 1e100;
    HTH_inv[0]=(HTH[4]*HTH[8]-HTH[5]*HTH[7])/det;
    HTH_inv[4]=(HTH[0]*HTH[8]-HTH[2]*HTH[6])/det;
    HTH_inv[8]=(HTH[0]*HTH[4]-HTH[1]*HTH[3])/det;
    double trace=HTH_inv[0]+HTH_inv[4]+HTH_inv[8];
    if (trace<0.0) trace=0.0;
    return sqrt(range_variance*trace);
}

/*
 * Error Metrics Update - Welford online algorithm (L5)
 * Reference: Welford (1962) Technometrics 4(3)
 */
void error_metrics_update(uwb_error_metrics_t *metrics,
                          const uwb_pos3d_t *estimated,
                          const uwb_pos3d_t *ground_truth,
                          uwb_positioning_dim_t dim) {
    (void)dim;
    double error_2d,error_3d,dx,dy,dz;
    if (!metrics||!estimated||!ground_truth) return;
    dx=estimated->x-ground_truth->x; dy=estimated->y-ground_truth->y;
    dz=estimated->z-ground_truth->z;
    error_2d=sqrt(dx*dx+dy*dy);
    error_3d=sqrt(dx*dx+dy*dy+dz*dz);
    metrics->num_samples++;
    double delta2=error_2d-metrics->mean_error_2d;
    metrics->mean_error_2d+=delta2/metrics->num_samples;
    double delta3=error_3d-metrics->mean_error_3d;
    metrics->mean_error_3d+=delta3/metrics->num_samples;
    double old_std=metrics->stddev_error;
    double new_var=((metrics->num_samples-1)*old_std*old_std+
        delta2*(error_2d-metrics->mean_error_2d))/metrics->num_samples;
    metrics->stddev_error=(new_var>0?sqrt(new_var):0.0);
    if (error_2d>metrics->max_error) metrics->max_error=error_2d;
    if (error_2d<metrics->min_error) metrics->min_error=error_2d;
    metrics->rmse_2d=sqrt((metrics->rmse_2d*metrics->rmse_2d*
        (metrics->num_samples-1)+error_2d*error_2d)/metrics->num_samples);
    metrics->rmse_3d=sqrt((metrics->rmse_3d*metrics->rmse_3d*
        (metrics->num_samples-1)+error_3d*error_3d)/metrics->num_samples);
    double threshold=3.0*metrics->stddev_error;
    if (threshold<1e-6) threshold=3.0;
    if (error_2d>metrics->mean_error_2d+threshold) metrics->num_outliers++;
}

void error_metrics_finalize(uwb_error_metrics_t *metrics) {
    if (!metrics||metrics->num_samples==0) return;
    double sigma_2d=metrics->stddev_error;
    if (sigma_2d<1e-10) sigma_2d=metrics->rmse_2d;
    metrics->cep50=0.6764*sigma_2d;
    metrics->cep90=1.8221*sigma_2d;
    metrics->sep50=0.823*sigma_2d;
    metrics->sep90=1.818*sigma_2d;
}
