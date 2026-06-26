/**
 * @file nav_ahrs.c
 * @brief Attitude and Heading Reference System (AHRS)
 *
 * Implements Madgwick, Mahony, and complementary filter algorithms.
 * L5: AHRS algorithms   L6: IMU+magnetometer attitude determination
 *
 * Reference: Madgwick (2010), Mahony (2008)
 */

#include "nav_rotation.h"
#include "nav_imu.h"
#include <math.h>
#include <string.h>

void nav_ahrs_madgwick(nav_quat_t *q,
                        const NAV_PRECISION gyro[3],
                        const NAV_PRECISION accel[3],
                        const NAV_PRECISION mag[3],
                        NAV_PRECISION dt, NAV_PRECISION beta) {
    if (!q || !gyro || !accel) return;
    NAV_PRECISION qw=q->w, qx=q->x, qy=q->y, qz=q->z;
    NAV_PRECISION ax=accel[0], ay=accel[1], az=accel[2];
    NAV_PRECISION na=sqrt(ax*ax+ay*ay+az*az);
    if (na<1e-10) return;
    ax/=na; ay/=na; az/=na;
    NAV_PRECISION mx=0, my=0, mz=0;
    int um=0;
    if (mag) {
        mx=mag[0]; my=mag[1]; mz=mag[2];
        NAV_PRECISION nm=sqrt(mx*mx+my*my+mz*mz);
        if (nm>1e-10) { mx/=nm; my/=nm; mz/=nm; um=1; }
    }
    NAV_PRECISION f[6], J[24];
    f[0]=2.0*(qx*qz - qw*qy) - ax;
    f[1]=2.0*(qw*qx + qy*qz) - ay;
    f[2]=2.0*(0.5 - qx*qx - qy*qy) - az;
    J[0]=-2.0*qy; J[1]=2.0*qz; J[2]=-2.0*qw; J[3]=2.0*qx;
    J[4]=2.0*qx;  J[5]=2.0*qw; J[6]=2.0*qz;  J[7]=2.0*qy;
    J[8]=0.0;     J[9]=-4.0*qx;J[10]=-4.0*qy;J[11]=0.0;
    int ne=3;
    if (um) {
        ne=6;
        NAV_PRECISION hx=2.0*mx*(0.5-qy*qy-qz*qz)+2.0*my*(qx*qy-qw*qz)+2.0*mz*(qx*qz+qw*qy);
        NAV_PRECISION hy=2.0*mx*(qx*qy+qw*qz)+2.0*my*(0.5-qx*qx-qz*qz)+2.0*mz*(qy*qz-qw*qx);
        NAV_PRECISION bx=sqrt(hx*hx+hy*hy);
        NAV_PRECISION bz=2.0*mx*(qx*qz-qw*qy)+2.0*my*(qy*qz+qw*qx)+2.0*mz*(0.5-qx*qx-qy*qy);
        f[3]=2.0*bx*(0.5-qy*qy-qz*qz)+2.0*bz*(qx*qz-qw*qy)-mx;
        f[4]=2.0*bx*(qx*qy-qw*qz)+2.0*bz*(qw*qx+qy*qz)-my;
        f[5]=2.0*bx*(qw*qy+qx*qz)+2.0*bz*(0.5-qx*qx-qy*qy)-mz;
        J[12]=-bx*2.0*qy+bz*2.0*qz; J[13]=-bx*2.0*qz-bz*2.0*qy;
        J[14]=bx*2.0*qz+2.0*bz*qx;  J[15]=bx*2.0*qy-4.0*bx*qy;
        J[16]=-bx*2.0*qz+bz*2.0*qy; J[17]=bx*2.0*qy+bz*2.0*qz;
        J[18]=bx*2.0*qy+2.0*bz*qw;  J[19]=bx*2.0*qz-4.0*bx*qx;
        J[20]=bx*2.0*qy+bz*2.0*qx;  J[21]=bx*2.0*qz+bz*2.0*qw;
        J[22]=bx*2.0*qx-4.0*bx*qx;  J[23]=-bx*2.0*qw+bz*2.0*qz;
    }
    NAV_PRECISION grad[4]={0,0,0,0};
    for (int i=0;i<4;i++)
        for (int j=0;j<ne;j++)
            grad[i] += J[j*4+i] * f[j];
    NAV_PRECISION ng=sqrt(grad[0]*grad[0]+grad[1]*grad[1]+grad[2]*grad[2]+grad[3]*grad[3]);
    if (ng>1e-10) { grad[0]/=ng; grad[1]/=ng; grad[2]/=ng; grad[3]/=ng; }
    NAV_PRECISION gx=gyro[0], gy=gyro[1], gz=gyro[2];
    NAV_PRECISION dqw=0.5*(-qx*gx-qy*gy-qz*gz)-beta*grad[0];
    NAV_PRECISION dqx=0.5*(qw*gx+qy*gz-qz*gy)-beta*grad[1];
    NAV_PRECISION dqy=0.5*(qw*gy+qz*gx-qx*gz)-beta*grad[2];
    NAV_PRECISION dqz=0.5*(qw*gz+qx*gy-qy*gx)-beta*grad[3];
    q->w+=dqw*dt; q->x+=dqx*dt; q->y+=dqy*dt; q->z+=dqz*dt;
    nav_quat_normalize(q);
}

void nav_ahrs_mahony(nav_quat_t *q,
                      const NAV_PRECISION gyro[3],
                      const NAV_PRECISION accel[3],
                      const NAV_PRECISION mag[3],
                      NAV_PRECISION dt, NAV_PRECISION kp, NAV_PRECISION ki) {
    if (!q || !gyro || !accel) return;
    NAV_PRECISION ax=accel[0], ay=accel[1], az=accel[2];
    NAV_PRECISION na=sqrt(ax*ax+ay*ay+az*az);
    if (na<1e-10) return;
    ax/=na; ay/=na; az/=na;
    NAV_PRECISION qw=q->w, qx=q->x, qy=q->y, qz=q->z;
    NAV_PRECISION vx=2.0*(qx*qz - qw*qy);
    NAV_PRECISION vy=2.0*(qw*qx + qy*qz);
    NAV_PRECISION vz=2.0*(qw*qw + qz*qz - 0.5);
    NAV_PRECISION ex=ay*vz-az*vy, ey=az*vx-ax*vz, ez=ax*vy-ay*vx;
    if (mag) {
        NAV_PRECISION mx=mag[0], my=mag[1], mz=mag[2];
        NAV_PRECISION nm=sqrt(mx*mx+my*my+mz*mz);
        if (nm>1e-10) {
            mx/=nm; my/=nm; mz/=nm;
            NAV_PRECISION hx=2.0*mx*(0.5-qy*qy-qz*qz)+2.0*my*(qx*qy-qw*qz)+2.0*mz*(qx*qz+qw*qy);
            NAV_PRECISION hy=2.0*mx*(qx*qy+qw*qz)+2.0*my*(0.5-qx*qx-qz*qz)+2.0*mz*(qy*qz-qw*qx);
            NAV_PRECISION bx=sqrt(hx*hx+hy*hy);
            NAV_PRECISION bz=2.0*mx*(qx*qz-qw*qy)+2.0*my*(qy*qz+qw*qx)+2.0*mz*(0.5-qx*qx-qy*qy);
            NAV_PRECISION wx=2.0*bx*(0.5-qy*qy-qz*qz)+2.0*bz*(qx*qz-qw*qy);
            NAV_PRECISION wy=2.0*bx*(qx*qy-qw*qz)+2.0*bz*(qw*qx+qy*qz);
            NAV_PRECISION wz=2.0*bx*(qw*qy+qx*qz)+2.0*bz*(0.5-qx*qx-qy*qy);
            ex+=my*wz-mz*wy; ey+=mz*wx-mx*wz; ez+=mx*wy-my*wx;
        }
    }
    static NAV_PRECISION ei[3]={0,0,0};
    ei[0]+=ex*ki*dt; ei[1]+=ey*ki*dt; ei[2]+=ez*ki*dt;
    NAV_PRECISION gx=gyro[0]+kp*ex+ei[0];
    NAV_PRECISION gy=gyro[1]+kp*ey+ei[1];
    NAV_PRECISION gz=gyro[2]+kp*ez+ei[2];
    nav_vector3_t omega; omega.x=gx; omega.y=gy; omega.z=gz;
    nav_quat_kinematics(q, &omega, dt);
}

void nav_ahrs_complementary(nav_euler_t *euler,
                             const NAV_PRECISION gyro[3],
                             const NAV_PRECISION accel[3],
                             const NAV_PRECISION mag[3],
                             NAV_PRECISION dt, NAV_PRECISION alpha) {
    if (!euler || !gyro || !accel) return;
    NAV_PRECISION ax=accel[0], ay=accel[1], az=accel[2];
    NAV_PRECISION na=sqrt(ax*ax+ay*ay+az*az);
    if (na<1e-10) return;
    NAV_PRECISION ra=atan2(ay, az);
    NAV_PRECISION pa=atan2(-ax, sqrt(ay*ay+az*az));
    NAV_PRECISION cp=cos(euler->pitch), sp=sin(euler->pitch);
    NAV_PRECISION cr=cos(euler->roll), sr=sin(euler->roll);
    NAV_PRECISION rd=gyro[0]+sr*sp/cp*gyro[1]+cr*sp/cp*gyro[2];
    NAV_PRECISION pd=cr*gyro[1]-sr*gyro[2];
    NAV_PRECISION yd=sr/cp*gyro[1]+cr/cp*gyro[2];
    euler->roll = alpha*(euler->roll+rd*dt) + (1.0-alpha)*ra;
    euler->pitch= alpha*(euler->pitch+pd*dt) + (1.0-alpha)*pa;
    euler->yaw  = euler->yaw + yd*dt;
    if (mag) {
        NAV_PRECISION mx=mag[0], my=mag[1], mz=mag[2];
        NAV_PRECISION mX=mx*cos(euler->pitch)+my*sr*sp+mz*cr*sp;
        NAV_PRECISION mY=my*cr-mz*sr;
        NAV_PRECISION ym=atan2(-mY, mX);
        euler->yaw=alpha*euler->yaw+(1.0-alpha)*ym;
    }
    euler->roll=nav_wrap_pi(euler->roll);
    euler->pitch=nav_wrap_pi(euler->pitch);
    euler->yaw=nav_wrap_pi(euler->yaw);
}
