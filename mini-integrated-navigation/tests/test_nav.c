#include "nav_common.h"
#include "nav_rotation.h"
#include "nav_kalman.h"
#include "nav_imu.h"
#include "nav_gnss.h"
#include "nav_integration.h"
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <string.h>

#define TOL 1e-10
static int tests_run = 0, tests_passed = 0;

#define TEST(name) do { tests_run++; printf("  %s... ", name); } while(0)
#define PASS() do { tests_passed++; printf("PASS\n"); } while(0)
#define FAIL(m) do { printf("FAIL: %s\n", m); } while(0)
#define CHECK(c,m) do { if(!(c)){FAIL(m);return;} } while(0)
#define CF(a,b,t,m) do { if(fabs((a)-(b))>(t)){printf("FAIL: %s g%.6f e%.6f\n",m,a,b);return;} } while(0)

static void t_quat_id(void) {
    TEST("quat identity");
    nav_quat_t q; nav_quat_identity(&q);
    CHECK(q.w==1.0&&q.x==0.0&&q.y==0.0&&q.z==0.0,"bad id");
    PASS();
}

static void t_quat_norm(void) {
    TEST("quat normalize");
    nav_quat_t q; nav_quat_set(&q,2,2,2,2);
    CHECK(nav_quat_normalize(&q)==0,"norm fail");
    CF(nav_quat_norm(&q),1.0,TOL,"norm!=1");
    PASS();
}

static void t_quat_mul(void) {
    TEST("quat multiply");
    nav_quat_t a,b,r;
    nav_quat_set(&a,1,0,0,0); nav_quat_set(&b,0,1,0,0);
    nav_quat_multiply(&r,&a,&b);
    CF(r.x,1.0,TOL,"i*1");
    PASS();
}

static void t_quat_dcm(void) {
    TEST("quat-DCM roundtrip");
    nav_euler_t e={0.1,0.2,0.3},eb;
    nav_quat_t q,qb; nav_dcm_t d;
    nav_euler_to_quat(&q,&e);
    nav_quat_to_dcm(&d,&q);
    nav_dcm_to_quat(&qb,&d);
    nav_quat_to_euler(&eb,&qb);
    CF(eb.roll,e.roll,1e-9,"roll");
    CF(eb.pitch,e.pitch,1e-9,"pitch");
    CF(eb.yaw,e.yaw,1e-9,"yaw");
    PASS();
}

static void t_slerp(void) {
    TEST("SLERP");
    nav_quat_t q0,q1,r;
    nav_quat_identity(&q0);
    nav_euler_t e={0,NAV_PI/4,0};
    nav_euler_to_quat(&q1,&e);
    nav_quat_slerp(&r,&q0,&q1,0.5);
    CF(nav_quat_norm(&r),1.0,TOL,"slerp norm");
    PASS();
}

static void t_triad(void) {
    TEST("TRIAD");
    nav_vector3_t b1={1,0,0},b2={0,1,0},r1={1,0,0},r2={0,1,0};
    nav_dcm_t d;
    CHECK(nav_triad(&d,&b1,&b2,&r1,&r2)==0,"triad fail");
    for(int i=0;i<9;i++) CF(d.m[i],(i%4==0)?1.0:0.0,1e-9,"triad not I");
    PASS();
}

static void t_kf(void) {
    TEST("Kalman filter");
    nav_kf_t *kf=nav_kf_alloc(2,1);
    CHECK(kf!=NULL,"alloc");
    NAV_PRECISION x0[2]={0,0},P0[4]={1,0,0,1};
    nav_kf_set_x(kf,x0); nav_kf_set_P(kf,P0);
    NAV_PRECISION H[2]={1,0},R[1]={0.01};
    nav_kf_set_H(kf,H); nav_kf_set_R(kf,R);
    nav_kf_predict(kf);
    NAV_PRECISION z[1]={5.0};
    CHECK(nav_kf_update(kf,z)==0,"update fail");
    CF(nav_kf_get_x(kf)[0],5.0,0.1,"KF x[0]");
    nav_kf_free(kf);
    PASS();
}

static void t_align(void) {
    TEST("IMU alignment");
    nav_ins_state_t s;
    NAV_PRECISION a[3]={0,0,9.81},g[3]={0,0,7.29e-5};
    CHECK(nav_ins_coarse_alignment(&s,a,g,0.0)==0,"align fail");
    CF(s.pos.latitude,0.0,TOL,"lat");
    PASS();
}

static void t_gnss_conv(void) {
    TEST("GNSS coord conv");
    nav_geodetic_t g={0,0,0},g2;
    NAV_PRECISION e[3],e2[3];
    nav_geodetic_to_ecef(e,&g);
    nav_ecef_to_geodetic(&g2,e);
    nav_geodetic_to_ecef(e2,&g2);
    for(int i=0;i<3;i++) CF(e[i],e2[i],1e-6,"ecef r/t");
    PASS();
}

static void t_ci(void) {
    TEST("Cov intersection");
    NAV_PRECISION xa[2]={0,0},Pa[4]={1,0,0,1};
    NAV_PRECISION xb[2]={1,1},Pb[4]={1,0,0,1};
    NAV_PRECISION xf[2],Pf[4];
    nav_covariance_intersection(xf,Pf,xa,Pa,xb,Pb,2);
    CF(xf[0],0.5,0.1,"CI x[0]");
    CF(xf[1],0.5,0.1,"CI x[1]");
    PASS();
}

static void t_infof(void) {
    TEST("Info filter");
    nav_info_filter_t *f=nav_info_filter_alloc(2);
    CHECK(f!=NULL,"alloc");
    NAV_PRECISION H[4]={1,0,0,1},R[4]={0.1,0,0,0.1},z[2]={5,-3};
    nav_info_filter_update(f,H,R,z,2);
    NAV_PRECISION x[2];
    nav_info_filter_get_x(x,f);
    CF(x[0],5.0,0.5,"info x[0]");
    nav_info_filter_free(f);
    PASS();
}

static void t_ma(void) {
    TEST("Moving average");
    extern void nav_moving_average_3(NAV_PRECISION[3],const NAV_PRECISION[3],const NAV_PRECISION[3],const NAV_PRECISION[3]);
    NAV_PRECISION o[3],p[3]={1,2,3},c[3]={2,3,4},n[3]={3,4,5};
    nav_moving_average_3(o,p,c,n);
    CF(o[0],2.0,TOL,"MA0"); CF(o[1],3.0,TOL,"MA1"); CF(o[2],4.0,TOL,"MA2");
    PASS();
}

static void t_loose(void) {
    TEST("Loose init");
    nav_loose_integration_t integ;
    nav_loose_config_t cfg;
    memset(&cfg,0,sizeof(cfg));
    cfg.type=NAV_INTEG_LOOSE; cfg.use_gnss_velocity=1;
    cfg.pos_noise_n=cfg.pos_noise_e=1.0; cfg.pos_noise_d=4.0;
    cfg.vel_noise_n=cfg.vel_noise_e=0.01; cfg.vel_noise_d=0.04;
    CHECK(nav_loose_init(&integ,&cfg)==0,"loose init");
    CHECK(integ.ekf!=NULL,"no ekf");
    nav_ekf_free(integ.ekf);
    PASS();
}

static void t_tight(void) {
    TEST("Tight init");
    nav_tight_integration_t integ;
    nav_tight_config_t cfg;
    memset(&cfg,0,sizeof(cfg));
    cfg.type=NAV_INTEG_TIGHT; cfg.max_svs=12;
    cfg.pr_noise=100.0; cfg.prr_noise=0.01; cfg.use_doppler=1;
    CHECK(nav_tight_init(&integ,&cfg)==0,"tight init");
    CHECK(integ.ekf!=NULL,"no ekf");
    nav_ekf_free(integ.ekf);
    PASS();
}

static void t_outlier(void) {
    TEST("Outlier detect");
    extern int nav_outlier_detect(const NAV_PRECISION[3],const NAV_PRECISION[][3],int,NAV_PRECISION);
    NAV_PRECISION h[3][3]={{1,2,3},{1,2,3},{1,2,3}},d[3]={1,2,3};
    CHECK(nav_outlier_detect(d,h,3,3.0)==0,"false positive");
    PASS();
}

static void t_raim(void) {
    TEST("RAIM");
    extern int nav_raim_residual(int*,NAV_PRECISION*,NAV_PRECISION*,const NAV_PRECISION*,const NAV_PRECISION*,int,NAV_PRECISION);
    NAV_PRECISION G[24]={0,0,1,1, 0,1,0,1, 1,0,0,1, 0,-1,0,1, -1,0,0,1};
    NAV_PRECISION r[5]={0.1,-0.05,0.02,-0.03,0.01};
    int f; NAV_PRECISION s,t;
    CHECK(nav_raim_residual(&f,&s,&t,G,r,5,1e-5)==0,"RAIM fail");
    PASS();
}

static void t_hav(void) {
    TEST("Haversine");
    extern NAV_PRECISION nav_horizontal_distance(NAV_PRECISION,NAV_PRECISION,NAV_PRECISION,NAV_PRECISION);
    NAV_PRECISION d=nav_horizontal_distance(0,0,0,NAV_PI/180.0);
    CHECK(d>110000.0&&d<112000.0,"1deg != ~111km");
    PASS();
}

int main(void) {
    printf("=== mini-integrated-navigation Test Suite ===\n\n");
    t_quat_id(); t_quat_norm(); t_quat_mul(); t_quat_dcm();
    t_slerp(); t_triad(); t_kf(); t_align(); t_gnss_conv();
    t_ci(); t_infof(); t_ma(); t_loose(); t_tight();
    t_outlier(); t_raim(); t_hav();
    printf("\n=== %d/%d tests passed ===\n", tests_passed, tests_run);
    return (tests_passed==tests_run)?0:1;
}
