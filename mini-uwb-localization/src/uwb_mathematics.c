#include "uwb_mathematics.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

void vec_scale(double *v, int n, double a) {int i;if(!v)return;for(i=0;i<n;i++)v[i]*=a;}
void vec_add(double *v, const double *v1, const double *v2, int n) {int i;if(!v||!v1||!v2)return;for(i=0;i<n;i++)v[i]=v1[i]+v2[i];}
void vec_sub(double *v, const double *v1, const double *v2, int n) {int i;if(!v||!v1||!v2)return;for(i=0;i<n;i++)v[i]=v1[i]-v2[i];}
double vec_dot(const double *v1, const double *v2, int n) {double s=0.0;int i;if(!v1||!v2)return 0.0;for(i=0;i<n;i++)s+=v1[i]*v2[i];return s;}
double vec_norm2(const double *v, int n) {return sqrt(vec_dot(v,v,n));}
void vec_normalize(double *v, int n) {double nm=vec_norm2(v,n);if(nm>1e-15)vec_scale(v,n,1.0/nm);}
void mat_vec_mul(double *r, const double *A, const double *x, int m, int n) {int i,j;if(!r||!A||!x)return;for(i=0;i<m;i++){r[i]=0.0;for(j=0;j<n;j++)r[i]+=A[i*n+j]*x[j];}}
void mat_mul(double *C, const double *A, const double *B, int m, int k, int n) {int i,j,l;if(!C||!A||!B)return;for(i=0;i<m;i++)for(j=0;j<n;j++){C[i*n+j]=0.0;for(l=0;l<k;l++)C[i*n+j]+=A[i*k+l]*B[l*n+j];}}
void mat_transpose(double *B, const double *A, int m, int n) {int i,j;if(!B||!A)return;for(i=0;i<m;i++)for(j=0;j<n;j++)B[j*m+i]=A[i*n+j];}
void mat_add_identity(double *A, int n, double alpha) {int i;if(!A)return;for(i=0;i<n;i++)A[i*n+i]+=alpha;}
void mat_set_identity(double *A, int n) {int i,j;if(!A)return;for(i=0;i<n;i++)for(j=0;j<n;j++)A[i*n+j]=(i==j)?1.0:0.0;}
void mat_copy(double *B, const double *A, int m, int n) {if(!B||!A)return;memcpy(B,A,m*n*sizeof(double));}
double mat_norm_frobenius(const double *A, int m, int n) {double s=0.0;int i,j;if(!A)return 0.0;for(i=0;i<m;i++)for(j=0;j<n;j++)s+=A[i*n+j]*A[i*n+j];return sqrt(s);}

int solve_linear_gauss(double *A, double *b, int n) {
    int i,j,k,pivot; double maxv,tmp,factor;
    if(!A||!b||n<1||n>UWB_MATH_MAX_DIM)return 0;
    for(i=0;i<n;i++){pivot=i;maxv=fabs(A[i*n+i]);for(j=i+1;j<n;j++)if(fabs(A[j*n+i])>maxv){maxv=fabs(A[j*n+i]);pivot=j;}
        if(maxv<1e-15)return 0;
        if(pivot!=i){for(j=0;j<n;j++){tmp=A[i*n+j];A[i*n+j]=A[pivot*n+j];A[pivot*n+j]=tmp;}tmp=b[i];b[i]=b[pivot];b[pivot]=tmp;}
        for(j=i+1;j<n;j++){factor=A[j*n+i]/A[i*n+i];for(k=i;k<n;k++)A[j*n+k]-=factor*A[i*n+k];b[j]-=factor*b[i];}}
    for(i=n-1;i>=0;i--){b[i]/=A[i*n+i];for(j=i-1;j>=0;j--)b[j]-=A[j*n+i]*b[i];}
    return 1;
}

int solve_linear_cholesky(double *A, double *b, int n) {
    int i,j,k; double sum;
    if(!A||!b||n<1||n>UWB_MATH_MAX_DIM)return 0;
    for(j=0;j<n;j++){sum=0.0;for(k=0;k<j;k++)sum+=A[j*n+k]*A[j*n+k];A[j*n+j]-=sum;if(A[j*n+j]<=1e-15)return 0;A[j*n+j]=sqrt(A[j*n+j]);for(i=j+1;i<n;i++){sum=0.0;for(k=0;k<j;k++)sum+=A[i*n+k]*A[j*n+k];A[i*n+j]=(A[i*n+j]-sum)/A[j*n+j];}}
    for(i=0;i<n;i++){sum=0.0;for(j=0;j<i;j++)sum+=A[i*n+j]*b[j];b[i]=(b[i]-sum)/A[i*n+i];}
    for(i=n-1;i>=0;i--){sum=0.0;for(j=i+1;j<n;j++)sum+=A[j*n+i]*b[j];b[i]=(b[i]-sum)/A[i*n+i];}
    return 1;
}

int mat_inverse(double *A_inv, double *A, int n) {
    int i,j,k,pivot; double maxv,tmp,factor;
    double aug[288];
    if(!A_inv||!A||n<1||n>UWB_MATH_MAX_DIM)return 0;
    for(i=0;i<n;i++){for(j=0;j<n;j++)aug[i*2*n+j]=A[i*n+j];for(j=0;j<n;j++)aug[i*2*n+n+j]=(i==j)?1.0:0.0;}
    for(i=0;i<n;i++){pivot=i;maxv=fabs(aug[i*2*n+i]);for(j=i+1;j<n;j++)if(fabs(aug[j*2*n+i])>maxv){maxv=fabs(aug[j*2*n+i]);pivot=j;}
        if(maxv<1e-15)return 0;
        if(pivot!=i)for(j=0;j<2*n;j++){tmp=aug[i*2*n+j];aug[i*2*n+j]=aug[pivot*2*n+j];aug[pivot*2*n+j]=tmp;}
        double piv=aug[i*2*n+i];for(j=0;j<2*n;j++)aug[i*2*n+j]/=piv;
        for(j=0;j<n;j++)if(j!=i){factor=aug[j*2*n+i];for(k=0;k<2*n;k++)aug[j*2*n+k]-=factor*aug[i*2*n+k];}}
    for(i=0;i<n;i++)for(j=0;j<n;j++)A_inv[i*n+j]=aug[i*2*n+n+j];
    return 1;
}

double solve_least_squares(const double *A, const double *b, int m, int n, double *x) {
    double AtA[144],Atb[12]; int i,j,k;
    if(!A||!b||!x||m<n||n>UWB_MATH_MAX_DIM)return 1e100;
    memset(AtA,0,sizeof(AtA));memset(Atb,0,sizeof(Atb));
    for(i=0;i<m;i++)for(j=0;j<n;j++){for(k=0;k<n;k++)AtA[j*n+k]+=A[i*n+j]*A[i*n+k];Atb[j]+=A[i*n+j]*b[i];}
    memcpy(x,Atb,n*sizeof(double));
    if(!solve_linear_gauss(AtA,x,n))return 1e100;
    double residual=0.0;
    for(i=0;i<m;i++){double pred=0.0;for(j=0;j<n;j++)pred+=A[i*n+j]*x[j];residual+=(b[i]-pred)*(b[i]-pred);}
    return sqrt(residual);
}

double solve_weighted_least_squares(const double *A, const double *b,
    const double *weights, int m, int n, double *x) {
    double AtWA[144],AtWb[12]; int i,j,k;
    if(!A||!b||!weights||!x||m<n||n>UWB_MATH_MAX_DIM)return 1e100;
    memset(AtWA,0,sizeof(AtWA));memset(AtWb,0,sizeof(AtWb));
    for(i=0;i<m;i++){double w=weights[i];for(j=0;j<n;j++){for(k=0;k<n;k++)AtWA[j*n+k]+=w*A[i*n+j]*A[i*n+k];AtWb[j]+=w*A[i*n+j]*b[i];}}
    memcpy(x,AtWb,n*sizeof(double));
    if(!solve_linear_gauss(AtWA,x,n))return 1e100;
    return 0.0;
}

int svd_decompose(double *A, int m, int n, double *U, double *S, double *V) {
    double AtA[144]; int i,j,k;
    (void)U;(void)V;
    if(!A||m<n||n>UWB_MATH_MAX_DIM)return 0;
    memset(S,0,n*sizeof(double)); memset(AtA,0,sizeof(AtA));
    for(i=0;i<m;i++)for(j=0;j<n;j++)for(k=0;k<n;k++)AtA[j*n+k]+=A[i*n+j]*A[i*n+k];
    if(n<=3){double evals[3];symm_eigenvalues_3x3(AtA,evals);for(i=0;i<n;i++)S[i]=(evals[i]>0)?sqrt(evals[i]):0.0;}
    return 1;
}

double mat_condition_number(const double *A, int m, int n) {
    double S[12],Acopy[144]; int i;
    if(!A)return 1e100;
    memcpy(Acopy,A,m*n*sizeof(double));
    svd_decompose(Acopy,m,n,0,S,0);
    if(S[0]<1e-15)return 1e100;
    double smin=S[0],smax=S[0];
    for(i=1;i<n;i++){if(S[i]<smin)smin=S[i];if(S[i]>smax)smax=S[i];}
    return (smin>1e-15)?smax/smin:1e100;
}

double mat_determinant(double *A, int n) {
    double det=1.0,Acopy[144]; int i,j,k,pivot; double maxv,tmp,factor;
    if(!A||n<1||n>UWB_MATH_MAX_DIM)return 0.0;
    memcpy(Acopy,A,n*n*sizeof(double));
    for(i=0;i<n;i++){pivot=i;maxv=fabs(Acopy[i*n+i]);for(j=i+1;j<n;j++)if(fabs(Acopy[j*n+i])>maxv){maxv=fabs(Acopy[j*n+i]);pivot=j;}
        if(maxv<1e-15)return 0.0;
        if(pivot!=i){det=-det;for(j=0;j<n;j++){tmp=Acopy[i*n+j];Acopy[i*n+j]=Acopy[pivot*n+j];Acopy[pivot*n+j]=tmp;}}
        det*=Acopy[i*n+i];
        for(j=i+1;j<n;j++){factor=Acopy[j*n+i]/Acopy[i*n+i];for(k=i;k<n;k++)Acopy[j*n+k]-=factor*Acopy[i*n+k];}}
    return det;
}

double mat_trace(const double *A, int n) {double t=0.0;int i;if(!A)return 0.0;for(i=0;i<n;i++)t+=A[i*n+i];return t;}

int mat_is_spd(const double *A, int n) {
    double L[144]; int i,j,k; double sum;
    if(!A||n<1||n>UWB_MATH_MAX_DIM)return 0;
    memcpy(L,A,n*n*sizeof(double));
    for(j=0;j<n;j++){sum=0.0;for(k=0;k<j;k++)sum+=L[j*n+k]*L[j*n+k];L[j*n+j]-=sum;if(L[j*n+j]<=1e-15)return 0;L[j*n+j]=sqrt(L[j*n+j]);for(i=j+1;i<n;i++){sum=0.0;for(k=0;k<j;k++)sum+=L[i*n+k]*L[j*n+k];L[i*n+j]=(L[i*n+j]-sum)/L[j*n+j];}}
    return 1;
}

void solve_upper_triangular(const double *U, double *b, int n) {
    int i,j; if(!U||!b)return;
    for(i=n-1;i>=0;i--){for(j=i+1;j<n;j++)b[i]-=U[i*n+j]*b[j];b[i]/=U[i*n+i];}
}

void solve_lower_triangular(const double *L, double *b, int n) {
    int i,j; if(!L||!b)return;
    for(i=0;i<n;i++){for(j=0;j<i;j++)b[i]-=L[i*n+j]*b[j];b[i]/=L[i*n+i];}
}

void qr_decompose_mgs(const double *A, int m, int n, double *Q, double *R) {
    int i,j,k; double At[144];
    if(!A||!Q||!R||m<n)return;
    memcpy(At,A,m*n*sizeof(double)); memset(R,0,n*n*sizeof(double));
    for(k=0;k<n;k++){double norm=0.0;for(i=0;i<m;i++)norm+=At[i*n+k]*At[i*n+k];R[k*n+k]=sqrt(norm);
        if(R[k*n+k]>1e-15){for(i=0;i<m;i++)Q[i*n+k]=At[i*n+k]/R[k*n+k];
            for(j=k+1;j<n;j++){R[k*n+j]=0.0;for(i=0;i<m;i++)R[k*n+j]+=Q[i*n+k]*At[i*n+j];for(i=0;i<m;i++)At[i*n+j]-=Q[i*n+k]*R[k*n+j];}}}
}

int symm_eigenvalues_2x2(const double *A, double *eigenvalues) {
    double a=A[0],b=A[1],c=A[2],d=A[3];
    double trace=a+d,det=a*d-b*c,disc=trace*trace-4.0*det;
    if(disc<0.0)disc=0.0;
    eigenvalues[0]=(trace+sqrt(disc))/2.0; eigenvalues[1]=(trace-sqrt(disc))/2.0;
    return 2;
}

int symm_eigenvalues_3x3(const double *A, double *eigenvalues) {
    double a=A[0],b=A[1],c=A[2],d=A[4],e_=A[5],f=A[8];
    double trace=a+d+f;
    double det_v=a*(d*f-e_*e_)-b*(b*f-e_*c)+c*(b*e_-d*c);
    double p=b*b+c*c+e_*e_;
    double q=(a*a+d*d+f*f+p)/3.0-trace*trace/9.0;
    double r=(det_v+(a+d+f)*(q-trace*trace/6.0)/3.0)/2.0;
    double phi;
    if(q*q*q-r*r>0)phi=acos(r/sqrt(q*q*q))/3.0; else phi=0.0;
    double sq=sqrt(q);
    eigenvalues[0]=trace/3.0+2.0*sq*cos(phi);
    eigenvalues[1]=trace/3.0-sq*(cos(phi)+sqrt(3.0)*sin(phi));
    eigenvalues[2]=trace/3.0-sq*(cos(phi)-sqrt(3.0)*sin(phi));
    return 3;
}

void householder_vector(const double *x, int n, double *v, double *beta) {
    double norm_x; int i;
    if(!x||!v||!beta||n<2)return;
    norm_x=0.0; for(i=1;i<n;i++)norm_x+=x[i]*x[i];
    if(norm_x<1e-15){*beta=0.0;return;}
    double mu=sqrt(x[0]*x[0]+norm_x);
    if(x[0]<=0)v[0]=x[0]-mu; else v[0]=-norm_x/(x[0]+mu);
    *beta=2.0*v[0]*v[0]/(norm_x+v[0]*v[0]);
    for(i=1;i<n;i++) v[i]=x[i];
    v[0]=1.0;
}

double wilkinson_shift(const double *diag, const double *offdiag, int n) {
    double d=(diag[n-2]-diag[n-1])/2.0;
    double sign=(d>=0)?1.0:-1.0;
    return diag[n-1]-sign*offdiag[n-2]*offdiag[n-2]/(fabs(d)+sqrt(d*d+offdiag[n-2]*offdiag[n-2]));
}

double hypot2(double a, double b) {
    double aa=fabs(a),ab=fabs(b);
    if(aa>ab)return aa*sqrt(1.0+(ab/aa)*(ab/aa));
    return (ab==0.0)?0.0:ab*sqrt(1.0+(aa/ab)*(aa/ab));
}

int mat_inverse_2x2(double *A_inv, const double *A) {
    double det=A[0]*A[3]-A[1]*A[2];
    if(fabs(det)<1e-15)return 0;
    A_inv[0]=A[3]/det;A_inv[1]=-A[1]/det;A_inv[2]=-A[2]/det;A_inv[3]=A[0]/det;
    return 1;
}

int mat_inverse_3x3(double *A_inv, const double *A) {
    double det=A[0]*(A[4]*A[8]-A[5]*A[7])-A[1]*(A[3]*A[8]-A[5]*A[6])+A[2]*(A[3]*A[7]-A[4]*A[6]);
    if(fabs(det)<1e-15)return 0;
    A_inv[0]=(A[4]*A[8]-A[5]*A[7])/det;A_inv[1]=(A[2]*A[7]-A[1]*A[8])/det;A_inv[2]=(A[1]*A[5]-A[2]*A[4])/det;
    A_inv[3]=(A[5]*A[6]-A[3]*A[8])/det;A_inv[4]=(A[0]*A[8]-A[2]*A[6])/det;A_inv[5]=(A[2]*A[3]-A[0]*A[5])/det;
    A_inv[6]=(A[3]*A[7]-A[4]*A[6])/det;A_inv[7]=(A[1]*A[6]-A[0]*A[7])/det;A_inv[8]=(A[0]*A[4]-A[1]*A[3])/det;
    return 1;
}
