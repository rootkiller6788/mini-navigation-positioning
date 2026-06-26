/**
 * @file nav_ins_gnss_tight.c
 * @brief Tightly Coupled INS/GNSS Integration
 *
 * L6 Canonical Problem: Tight integration with pseudorange/Doppler in 17-state EKF.
 * Reference: Groves, Chapter 14
 */

#include "nav_integration.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

int nav_tight_init(nav_tight_integration_t *integ,
                    const nav_tight_config_t *config) {
    if (!integ || !config) return -1;
    memset(integ, 0, sizeof(nav_tight_integration_t));
    integ->config = *config;
    integ->ekf = nav_ekf_alloc(17, config->max_svs * 2, 0);
    if (!integ->ekf) return -1;
    integ->ins.pos.latitude = 0.0;
    integ->ins.pos.longitude = 0.0;
    integ->ins.pos.altitude = 0.0;
    nav_quat_identity(&integ->ins.quat);
    nav_quat_to_dcm(&integ->ins.dcm, &integ->ins.quat);
    NAV_PRECISION Qclk[289] = {0};
    for (int i = 0; i < 17; i++) Qclk[i*18] = 1e-8;
    Qclk[15*18] = 1e-2;
    Qclk[16*18] = 1e-4;
    nav_kf_set_Q(&integ->ekf->base, Qclk);
    return 0;
}

void nav_tight_predict(nav_tight_integration_t *integ,
                        const nav_imu_meas_t *imu) {
    if (!integ || !imu) return;
    NAV_PRECISION gyro_c[3], accel_c[3];
    NAV_PRECISION gyro_r[3] = {imu->gyro_x, imu->gyro_y, imu->gyro_z};
    NAV_PRECISION accel_r[3] = {imu->accel_x, imu->accel_y, imu->accel_z};
    nav_imu_correct_bias(gyro_c, accel_c,
                          gyro_r, accel_r,
                          integ->ins.gyro_bias, integ->ins.accel_bias);
    nav_ins_mechanize(&integ->ins, gyro_c, accel_c, imu->dt);
    if (integ->ekf) {
        NAV_PRECISION Phi[225], Qd[225];
        nav_ins_error_transition(Phi, &integ->ins, imu->dt);
        nav_ins_process_noise(Qd, &integ->ins, 1e-6, 1e-4, 1e-10, 1e-6, imu->dt);
        NAV_PRECISION F17[289];
        memset(F17, 0, 289*sizeof(NAV_PRECISION));
        for (int i = 0; i < 15; i++)
            for (int j = 0; j < 15; j++)
                F17[i*17+j] = Phi[i*15+j];
        F17[15*17+15] = 1.0; F17[15*17+16] = 1.0;
        F17[16*17+16] = 1.0;
        nav_kf_set_F(&integ->ekf->base, F17);
        nav_kf_predict(&integ->ekf->base);
    }
    integ->dt_accum += imu->dt;
}

void nav_los_vector(NAV_PRECISION los[3],
                    const NAV_PRECISION rx[3], const NAV_PRECISION sat[3]) {
    if (!los || !rx || !sat) return;
    NAV_PRECISION dx=sat[0]-rx[0], dy=sat[1]-rx[1], dz=sat[2]-rx[2];
    NAV_PRECISION r=sqrt(dx*dx+dy*dy+dz*dz);
    if (r<1e-10) { los[0]=los[1]=los[2]=0; return; }
    los[0]=dx/r; los[1]=dy/r; los[2]=dz/r;
}

void nav_tight_pr_innovation(NAV_PRECISION *innov,
                              const nav_gnss_sv_t *sv,
                              const NAV_PRECISION rx[3],
                              NAV_PRECISION clk_bias) {
    if (!innov || !sv || !rx) return;
    NAV_PRECISION dx=sv->sat_x-rx[0], dy=sv->sat_y-rx[1], dz=sv->sat_z-rx[2];
    NAV_PRECISION pr_pred=sqrt(dx*dx+dy*dy+dz*dz);
    NAV_PRECISION pr_corr=nav_gnss_correct_pseudorange(
        sv->pseudorange, sv->sat_clk_bias, sv->iono_delay, sv->tropo_delay);
    *innov = pr_corr - pr_pred - clk_bias;
}

void nav_tight_prr_innovation(NAV_PRECISION *innov,
                               const nav_gnss_sv_t *sv,
                               const NAV_PRECISION rx[3],
                               const NAV_PRECISION rv[3],
                               NAV_PRECISION clk_drift) {
    if(!innov||!sv||!rx||!rv)return;
    NAV_PRECISION dx=sv->sat_x-rx[0],dy=sv->sat_y-rx[1],dz=sv->sat_z-rx[2];
    NAV_PRECISION r=sqrt(dx*dx+dy*dy+dz*dz);
    if(r<1e-10){*innov=0;return;}
    NAV_PRECISION dvx=sv->sat_vx-rv[0],dvy=sv->sat_vy-rv[1],dvz=sv->sat_vz-rv[2];
    *innov=sv->pseudorange_rate - (dx*dvx+dy*dvy+dz*dvz)/r - clk_drift;
}

int nav_tight_update(nav_tight_integration_t *integ,
                      const nav_gnss_sv_t *svs, int n_svs,
                      const nav_gnss_ephemeris_t *eph, int n_eph) {
    if (!integ || !svs || n_svs <= 0) return -1;
    (void)eph; (void)n_eph;
    NAV_PRECISION rx[3];
    nav_geodetic_to_ecef(rx, &integ->ins.pos);
    NAV_PRECISION la[3];
    nav_lever_arm_compensate(la, rx, integ->ins.dcm.m,
                              integ->config.lever_arm,
                              integ->ins.pos.latitude,
                              integ->ins.pos.longitude);
    NAV_PRECISION rv[3];
    NAV_PRECISION Re[9];
    nav_enu_to_ecef_matrix(Re, integ->ins.pos.latitude, integ->ins.pos.longitude);
    NAV_PRECISION ve=integ->ins.vel_ned.y, vn=integ->ins.vel_ned.x, vu=-integ->ins.vel_ned.z;
    rv[0]=Re[0]*ve+Re[1]*vn+Re[2]*vu;
    rv[1]=Re[3]*ve+Re[4]*vn+Re[5]*vu;
    rv[2]=Re[6]*ve+Re[7]*vn+Re[8]*vu;
    int mp = integ->config.use_doppler ? 2 : 1;
    int m = n_svs * mp;
    NAV_PRECISION *z = (NAV_PRECISION*)calloc(m, sizeof(NAV_PRECISION));
    NAV_PRECISION *H = (NAV_PRECISION*)calloc(m*17, sizeof(NAV_PRECISION));
    NAV_PRECISION *Rm = (NAV_PRECISION*)calloc(m*m, sizeof(NAV_PRECISION));
    if (!z || !H || !Rm) { free(z); free(H); free(Rm); return -1; }
    for (int i = 0; i < n_svs; i++) {
        NAV_PRECISION los[3];
        NAV_PRECISION sat_p[3] = {svs[i].sat_x, svs[i].sat_y, svs[i].sat_z};
        nav_los_vector(los, la, sat_p);
        nav_tight_pr_innovation(&z[i*mp], &svs[i], la, integ->x_err[15]);
        for (int j = 0; j < 3; j++)
            H[i*mp*17 + j] = -los[j];
        H[i*mp*17 + 15] = 1.0;
        Rm[i*mp*m + i*mp] = integ->config.pr_noise;
        if (integ->config.use_doppler) {
            nav_tight_prr_innovation(&z[i*mp+1], &svs[i], la, rv, integ->x_err[16]);
            H[(i*mp+1)*17 + 16] = 1.0;
            for (int j = 0; j < 3; j++)
                H[(i*mp+1)*17 + (3+j)] = -los[j];
            Rm[(i*mp+1)*m + (i*mp+1)] = integ->config.prr_noise;
        }
    }
    integ->n_active = n_svs;
    nav_kf_set_H(&integ->ekf->base, H);
    nav_kf_set_R(&integ->ekf->base, Rm);
    int ret = nav_kf_update(&integ->ekf->base, z);
    free(z); free(H); free(Rm);
    if (ret == 0) {
        integ->x_err[15] = integ->ekf->base.x[15];
        integ->x_err[16] = integ->ekf->base.x[16];
        nav_tight_correct(integ);
    }
    return ret;
}

void nav_tight_correct(nav_tight_integration_t *integ) {
    if (!integ) return;
    const NAV_PRECISION *x = nav_kf_get_x(&integ->ekf->base);
    if (!x) return;
    NAV_PRECISION Rm = nav_meridian_radius(integ->ins.pos.latitude);
    NAV_PRECISION Rn = nav_transverse_radius(integ->ins.pos.latitude);
    NAV_PRECISION cl = cos(integ->ins.pos.latitude);
    integ->ins.pos.latitude  -= x[0] / (Rm + integ->ins.pos.altitude);
    integ->ins.pos.longitude -= x[1] / ((Rn + integ->ins.pos.altitude) * cl);
    integ->ins.pos.altitude  += x[2];
    integ->ins.vel_ned.x -= x[3]; integ->ins.vel_ned.y -= x[4];
    integ->ins.vel_ned.z -= x[5];
    nav_vector3_t psi; psi.x=x[6]; psi.y=x[7]; psi.z=x[8];
    nav_quat_apply_correction(&integ->ins.quat, &psi);
    nav_quat_to_dcm(&integ->ins.dcm, &integ->ins.quat);
    for (int i=0;i<3;i++) {
        integ->ins.gyro_bias[i] += x[9+i];
        integ->ins.accel_bias[i] += x[12+i];
    }
    memset(integ->ekf->base.x, 0, 17*sizeof(NAV_PRECISION));
}

void nav_tight_get_solution(nav_solution_t *sol,
                             const nav_tight_integration_t *integ) {
    if (!sol || !integ) return;
    memset(sol, 0, sizeof(nav_solution_t));
    sol->pos = integ->ins.pos;
    sol->vel_ned = integ->ins.vel_ned;
    { nav_euler_t eul; nav_quat_to_euler(&eul, &integ->ins.quat);
      sol->roll = eul.roll; sol->pitch = eul.pitch; sol->yaw = eul.yaw; }
    memcpy(sol->quat, &integ->ins.quat, 4*sizeof(NAV_PRECISION));
    memcpy(sol->gyro_bias, integ->ins.gyro_bias, 3*sizeof(NAV_PRECISION));
    memcpy(sol->accel_bias, integ->ins.accel_bias, 3*sizeof(NAV_PRECISION));
}

void nav_lever_arm_compensate(NAV_PRECISION ge[3],
                               const NAV_PRECISION ie[3],
                               const NAV_PRECISION Rnb[9],
                               const NAV_PRECISION la[3],
                               NAV_PRECISION lat, NAV_PRECISION lon) {
    if (!ge || !ie || !Rnb || !la) return;
    NAV_PRECISION ln[3];
    for (int i=0;i<3;i++)
        ln[i]=Rnb[i*3]*la[0]+Rnb[i*3+1]*la[1]+Rnb[i*3+2]*la[2];
    NAV_PRECISION Rne[9];
    nav_enu_to_ecef_matrix(Rne, lat, lon);
    ge[0]=ie[0]+Rne[0]*ln[1]+Rne[1]*ln[0]+Rne[2]*(-ln[2]);
    ge[1]=ie[1]+Rne[3]*ln[1]+Rne[4]*ln[0]+Rne[5]*(-ln[2]);
    ge[2]=ie[2]+Rne[6]*ln[1]+Rne[7]*ln[0]+Rne[8]*(-ln[2]);
}
