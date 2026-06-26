/**
 * @file nav_gnss.c
 * @brief GNSS Measurement Models and PVT Solution
 *
 * L2: GNSS positioning, pseudorange model
 * L5: WLS positioning, satellite position, Klobuchar iono model
 * L6: Standalone GNSS PVT
 */

#include "nav_gnss.h"
#include "nav_kalman.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

#define GPS_GM 3.986005e14
#define GPS_OMEGA_DOT 7.2921151467e-5

void nav_gnss_sat_position(NAV_PRECISION pos[3], NAV_PRECISION vel[3],
                            NAV_PRECISION *clk,
                            const nav_gnss_ephemeris_t *eph, NAV_PRECISION t_tx) {
    if (!pos || !eph) return;
    NAV_PRECISION A = eph->sqrt_a * eph->sqrt_a;
    NAV_PRECISION n0 = sqrt(GPS_GM / (A*A*A));
    NAV_PRECISION n = n0 + eph->delta_n;
    NAV_PRECISION tk = t_tx - eph->toe;
    if (tk > 302400.0) tk -= 604800.0;
    if (tk < -302400.0) tk += 604800.0;
    NAV_PRECISION M = eph->M0 + n * tk;
    NAV_PRECISION E = M;
    for (int i = 0; i < 10; i++) {
        NAV_PRECISION dE = (M - E + eph->ecc*sin(E)) / (1.0 - eph->ecc*cos(E));
        E += dE;
        if (fabs(dE) < 1e-12) break;
    }
    NAV_PRECISION cosE = cos(E), sinE = sin(E);
    NAV_PRECISION v = atan2(sqrt(1.0-eph->ecc*eph->ecc)*sinE, cosE-eph->ecc);
    NAV_PRECISION phi = v + eph->w;
    NAV_PRECISION sin2p=sin(2.0*phi), cos2p=cos(2.0*phi);
    NAV_PRECISION du = eph->cus*sin2p + eph->cuc*cos2p;
    NAV_PRECISION dr = eph->crs*sin2p + eph->crc*cos2p;
    NAV_PRECISION di = eph->cis*sin2p + eph->cic*cos2p;
    NAV_PRECISION u = phi + du;
    NAV_PRECISION r = A*(1.0-eph->ecc*cosE) + dr;
    NAV_PRECISION incl = eph->i0 + di + eph->idot*tk;
    NAV_PRECISION Omega = eph->omega0 + (eph->omegadot-GPS_OMEGA_DOT)*tk
                          - GPS_OMEGA_DOT*eph->toe;
    NAV_PRECISION cosO=cos(Omega), sinO=sin(Omega);
    NAV_PRECISION cosu=cos(u), sinu=sin(u);
    NAV_PRECISION cosi=cos(incl), sini=sin(incl);
    pos[0] = r*cosu*cosO - r*sinu*cosi*sinO;
    pos[1] = r*cosu*sinO + r*sinu*cosi*cosO;
    pos[2] = r*sinu*sini;
    if (clk) {
        *clk = (eph->af0 + eph->af1*tk + eph->af2*tk*tk
                - 4.442807633e-10*eph->ecc*eph->sqrt_a*sinE) * NAV_C_LIGHT;
    }
    if (vel) memset(vel, 0, 3*sizeof(NAV_PRECISION));
}

void nav_gnss_azel(NAV_PRECISION *az, NAV_PRECISION *el,
                   const NAV_PRECISION rx[3], const NAV_PRECISION sat[3]) {
    if (!az||!el||!rx||!sat) return;
    nav_geodetic_t geo;
    nav_ecef_to_geodetic(&geo, rx);
    NAV_PRECISION R[9];
    nav_enu_to_ecef_matrix(R, geo.latitude, geo.longitude);
    NAV_PRECISION dx=sat[0]-rx[0], dy=sat[1]-rx[1], dz=sat[2]-rx[2];
    NAV_PRECISION e = R[0]*dx+R[3]*dy+R[6]*dz;
    NAV_PRECISION n = R[1]*dx+R[4]*dy+R[7]*dz;
    NAV_PRECISION u = R[2]*dx+R[5]*dy+R[8]*dz;
    NAV_PRECISION d = sqrt(e*e+n*n+u*u);
    if (d<1e-10){*az=0;*el=0;return;}
    *el=asin(u/d); *az=atan2(e,n);
    if(*az<0)*az+=2.0*NAV_PI;
}

void nav_iono_klobuchar(NAV_PRECISION *delay_m,
                         NAV_PRECISION lat, NAV_PRECISION lon,
                         NAV_PRECISION az, NAV_PRECISION el,
                         NAV_PRECISION tow,
                         const NAV_PRECISION a[4], const NAV_PRECISION b[4]) {
    if(!delay_m)return;
    *delay_m=0.0;
    if(!a||!b)return;
    NAV_PRECISION psi=0.0137/(el/NAV_PI+0.11)-0.022;
    NAV_PRECISION li=lat/NAV_PI+psi*cos(az);
    if(li>0.416)li=0.416;
    if(li<-0.416)li=-0.416;
    NAV_PRECISION loi=lon/NAV_PI+psi*sin(az)/cos(li*NAV_PI);
    NAV_PRECISION lm=li+0.064*cos((loi-1.617)*NAV_PI);
    NAV_PRECISION A=a[0]+a[1]*lm+a[2]*lm*lm+a[3]*lm*lm*lm;
    if(A<0)A=0;
    NAV_PRECISION P=b[0]+b[1]*lm+b[2]*lm*lm+b[3]*lm*lm*lm;
    if(P<72000)P=72000;
    NAV_PRECISION x=2.0*NAV_PI*(tow-50400.0)/P;
    NAV_PRECISION F=1.0+16.0*pow(0.53-el/NAV_PI,3.0);
    NAV_PRECISION Iz=5.0e-9;
    if(fabs(x)<1.57) Iz+=A*(1.0-x*x/2.0+x*x*x*x/24.0);
    *delay_m=Iz*NAV_C_LIGHT*F;
}

void nav_tropo_saastamoinen(NAV_PRECISION *d, NAV_PRECISION el, NAV_PRECISION h) {
    if(!d)return;
    NAV_PRECISION p0=1013.25*pow(1.0-2.2557e-5*h,5.2568);
    NAV_PRECISION T0=288.15-0.0065*h;
    NAV_PRECISION e0=11.66*exp(-4.0e-4*h);
    NAV_PRECISION se=sin(el); if(se<0.05)se=0.05;
    NAV_PRECISION zhd=0.0022768*p0/(1.0-0.00266*cos(2.0*el)-2.8e-7*h);
    NAV_PRECISION zwd=0.0022768*(1255.0/T0+0.05)*e0;
    NAV_PRECISION m=1.001/sqrt(0.002001+se*se);
    *d=zhd*m+zwd*m;
}

NAV_PRECISION nav_gnss_correct_pseudorange(NAV_PRECISION raw, NAV_PRECISION clk,
                                            NAV_PRECISION iono, NAV_PRECISION tropo) {
    return raw - clk + iono + tropo;
}

void nav_gnss_compute_dop(NAV_PRECISION *gd, NAV_PRECISION *pd,
                           NAV_PRECISION *hd, NAV_PRECISION *vd,
                           NAV_PRECISION *td, const NAV_PRECISION *G, int n) {
    if(!gd||!pd||!hd||!vd||!td||!G||n<4){
        if(gd)*gd=99;
        if(pd)*pd=99;
        if(hd)*hd=99;
        if(vd)*vd=99;
        if(td)*td=99;
        return;
    }
    NAV_PRECISION H[16]={0};
    for(int i=0;i<4;i++) for(int j=0;j<4;j++)
        for(int k=0;k<n;k++) H[i*4+j]+=G[k*4+i]*G[k*4+j];
    /* 4x4 inverse via adjugate */
    NAV_PRECISION det = H[0]*(H[5]*(H[10]*H[15]-H[11]*H[14])-H[6]*(H[9]*H[15]-H[11]*H[13])+H[7]*(H[9]*H[14]-H[10]*H[13]))
                       -H[1]*(H[4]*(H[10]*H[15]-H[11]*H[14])-H[6]*(H[8]*H[15]-H[11]*H[12])+H[7]*(H[8]*H[14]-H[10]*H[12]))
                       +H[2]*(H[4]*(H[9]*H[15]-H[11]*H[13])-H[5]*(H[8]*H[15]-H[11]*H[12])+H[7]*(H[8]*H[13]-H[9]*H[12]))
                       -H[3]*(H[4]*(H[9]*H[14]-H[10]*H[13])-H[5]*(H[8]*H[14]-H[10]*H[12])+H[6]*(H[8]*H[13]-H[9]*H[12]));
    if(fabs(det)<1e-15){*gd=99;*pd=99;*hd=99;*vd=99;*td=99;return;}
    NAV_PRECISION A[16];
    for(int r=0;r<4;r++) for(int c=0;c<4;c++){
        NAV_PRECISION s[9]; int si=0;
        for(int i=0;i<4;i++) if(i!=r) for(int j=0;j<4;j++) if(j!=c) s[si++]=H[i*4+j];
        NAV_PRECISION cf=s[0]*(s[4]*s[8]-s[5]*s[7])-s[1]*(s[3]*s[8]-s[5]*s[6])+s[2]*(s[3]*s[7]-s[4]*s[6]);
        A[c*4+r]=((r+c)%2? -1.0:1.0)*cf/det;
    }
    *gd=sqrt(A[0]+A[5]+A[10]+A[15]);
    *pd=sqrt(A[0]+A[5]+A[10]);
    *hd=sqrt(A[0]+A[5]);
    *vd=sqrt(A[10]);
    *td=sqrt(A[15]);
}

int nav_gnss_wls_pvt(nav_gnss_solution_t *sol, const nav_gnss_sv_t *svs,
                     int n_svs, const nav_gnss_config_t *cfg) {
    if (!sol || !svs || n_svs < 4) return -1;
    memset(sol, 0, sizeof(nav_gnss_solution_t));
    NAV_PRECISION xr[4] = {0,0,0,0};
    for (int iter = 0; iter < 10; iter++) {
        NAV_PRECISION G[64], dy[16];
        int nu = 0;
        for (int i = 0; i < n_svs && nu < 16; i++) {
            const nav_gnss_sv_t *sv = &svs[i];
            NAV_PRECISION s[3] = {sv->sat_x, sv->sat_y, sv->sat_z};
            NAV_PRECISION dx = s[0]-xr[0], dy2 = s[1]-xr[1], dz = s[2]-xr[2];
            NAV_PRECISION rng = sqrt(dx*dx+dy2*dy2+dz*dz);
            if (rng < 1.0) continue;
            G[nu*4+0] = -dx/rng;
            G[nu*4+1] = -dy2/rng;
            G[nu*4+2] = -dz/rng;
            G[nu*4+3] = 1.0;
            NAV_PRECISION pc = nav_gnss_correct_pseudorange(sv->pseudorange,
                                sv->sat_clk_bias, sv->iono_delay, sv->tropo_delay);
            dy[nu] = pc - rng - xr[3];
            nu++;
        }
        if (nu < 4) return -1;
        NAV_PRECISION GtG[16]={0}, Gtdy[4]={0};
        for (int i=0;i<4;i++) {
            for (int j=0;j<4;j++)
                for(int k=0;k<nu;k++) GtG[i*4+j] += G[k*4+i]*G[k*4+j];
            for(int k=0;k<nu;k++) Gtdy[i] += G[k*4+i]*dy[k];
        }
        NAV_PRECISION GtGi[16];
        if (nav_matrix_inverse_spd(GtGi, GtG, 4) != 0) return -1;
        NAV_PRECISION dxs[4]={0};
        for(int i=0;i<4;i++)
            for(int j=0;j<4;j++) dxs[i] += GtGi[i*4+j]*Gtdy[j];
        for(int i=0;i<4;i++) xr[i] += dxs[i];
        if (fabs(dxs[0])+fabs(dxs[1])+fabs(dxs[2]) < 0.001) break;
    }
    NAV_PRECISION ecef[3] = {xr[0], xr[1], xr[2]};
    nav_ecef_to_geodetic(&sol->pos, ecef);
    sol->clock_bias = xr[3];
    sol->num_svs = n_svs;
    (void)cfg;
    return 0;
}

void nav_ecef_to_geodetic(nav_geodetic_t *geo, const NAV_PRECISION e[3]) {
    if(!geo||!e)return;
    NAV_PRECISION x=e[0],y=e[1],z=e[2];
    NAV_PRECISION p=sqrt(x*x+y*y);
    NAV_PRECISION lon=atan2(y,x);
    NAV_PRECISION lat=atan2(z,p*(1.0-NAV_EARTH_ECCENTRICITY2));
    NAV_PRECISION alt=0.0;
    for(int i=0;i<5;i++){
        NAV_PRECISION sl=sin(lat);
        NAV_PRECISION Rn=NAV_EARTH_RADIUS_M/sqrt(1.0-NAV_EARTH_ECCENTRICITY2*sl*sl);
        alt=p/cos(lat)-Rn;
        lat=atan2(z,p*(1.0-NAV_EARTH_ECCENTRICITY2*Rn/(Rn+alt)));
    }
    geo->latitude=lat; geo->longitude=lon; geo->altitude=alt;
}

void nav_geodetic_to_ecef(NAV_PRECISION e[3], const nav_geodetic_t *g) {
    if(!e||!g)return;
    NAV_PRECISION sl=sin(g->latitude),cl=cos(g->latitude);
    NAV_PRECISION Rn=NAV_EARTH_RADIUS_M/sqrt(1.0-NAV_EARTH_ECCENTRICITY2*sl*sl);
    NAV_PRECISION Nph=(Rn+g->altitude)*cl;
    e[0]=Nph*cos(g->longitude); e[1]=Nph*sin(g->longitude);
    e[2]=(Rn*(1.0-NAV_EARTH_ECCENTRICITY2)+g->altitude)*sl;
}

void nav_enu_to_ecef_matrix(NAV_PRECISION R[9], NAV_PRECISION lat, NAV_PRECISION lon) {
    if(!R)return;
    NAV_PRECISION sl=sin(lat),cl=cos(lat),sLo=sin(lon),cLo=cos(lon);
    R[0]=-sLo;    R[1]=-sl*cLo; R[2]=cl*cLo;
    R[3]=cLo;     R[4]=-sl*sLo; R[5]=cl*sLo;
    R[6]=0.0;     R[7]=cl;      R[8]=sl;
}
