#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
/**
 * mini-uwb-localization: UWB Channel Models
 * IEEE 802.15.4a SV model, path loss, multipath generation
 * L3 Math Structures, L4 Friis for UWB, L6 Indoor channel modeling
 */
#include "uwb_channel.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

static const struct {
    double Lambda, lambda, Gamma, gamma, sigma_cluster, sigma_ray;
} sv_params[] = {
    [UWB_ENV_RESIDENTIAL_LOS]   = {0.047, 1.54, 22.61, 12.53, 2.75, 3.0},
    [UWB_ENV_RESIDENTIAL_NLOS]  = {0.118, 1.77, 26.27, 17.50, 2.93, 3.0},
    [UWB_ENV_OFFICE_LOS]        = {0.067, 1.01, 14.62,  7.81, 3.52, 3.0},
    [UWB_ENV_OFFICE_NLOS]       = {0.143, 1.55, 33.03, 16.18, 3.23, 3.0},
    [UWB_ENV_INDUSTRIAL_LOS]    = {0.019, 1.32,  7.76,  9.84, 3.12, 4.0},
    [UWB_ENV_INDUSTRIAL_NLOS]   = {0.274, 1.59, 18.58,  6.04, 3.97, 4.0},
    [UWB_ENV_OPEN_OUTDOOR_LOS]  = {0.033, 1.21, 17.27,  8.03, 1.90, 3.0},
    [UWB_ENV_OPEN_OUTDOOR_NLOS] = {0.099, 1.63, 32.61, 13.50, 2.61, 3.0},
};

void uwb_channel_sv_init(uwb_sv_channel_t *channel,
                         uwb_channel_environment_t env,
                         double distance_m, double fc_hz, double bw_hz) {
    if(!channel)return;
    memset(channel,0,sizeof(*channel));
    channel->environment=env;
    channel->distance_m=distance_m;
    channel->center_frequency_hz=fc_hz;
    channel->bandwidth_hz=bw_hz;
    if(env<=UWB_ENV_OPEN_OUTDOOR_NLOS){
        channel->cluster_decay_ns=sv_params[env].Gamma;
        channel->ray_decay_ns=sv_params[env].gamma;
        channel->cluster_arrival_rate_1_per_ns=sv_params[env].Lambda;
        channel->ray_arrival_rate_1_per_ns=sv_params[env].lambda;
    }
}

void uwb_channel_generate_paths(uwb_sv_channel_t *channel, uint32_t seed) {
    int ci,ri;
    double cluster_time=0.0;
    int path_idx=0;
    if(!channel)return;
    srand(seed);
    channel->num_clusters=0; channel->num_paths=0;
    for(ci=0;ci<UWB_CHANNEL_MAX_CLUSTERS&&channel->num_paths<UWB_CHANNEL_MAX_PATHS-10;ci++){
        double cluster_amplitude=exp(-cluster_time/channel->cluster_decay_ns)*
            exp(((double)rand()/RAND_MAX-0.5)*2.0*
                sv_params[channel->environment].sigma_cluster*0.23026);
        double ray_time=0.0;
        for(ri=0;ri<20&&path_idx<UWB_CHANNEL_MAX_PATHS;ri++){
            double ray_amplitude=cluster_amplitude*
                exp(-ray_time/channel->ray_decay_ns)*
                exp(((double)rand()/RAND_MAX-0.5)*2.0*
                    sv_params[channel->environment].sigma_ray*0.23026);
            channel->paths[path_idx].amplitude=ray_amplitude;
            channel->paths[path_idx].delay_ns=cluster_time+ray_time;
            channel->paths[path_idx].phase_rad=((double)rand()/RAND_MAX)*2.0*M_PI;
            channel->paths[path_idx].cluster_index=ci;
            path_idx++;
            ray_time+=(-log(1.0-(double)rand()/(RAND_MAX+1.0)))/
                channel->ray_arrival_rate_1_per_ns;
        }
        channel->num_clusters++;
        cluster_time+=(-log(1.0-(double)rand()/(RAND_MAX+1.0)))/
            channel->cluster_arrival_rate_1_per_ns;
    }
    channel->num_paths=path_idx;
}

double uwb_channel_path_loss_db(double distance_m, const uwb_path_loss_model_t *model) {
    if(!model||distance_m<=0.0)return 200.0;
    return model->pl0_db+10.0*model->path_loss_exponent*log10(distance_m/1.0);
}

void uwb_channel_get_path_loss_model(uwb_channel_environment_t env,
                                     uwb_path_loss_model_t *model) {
    if(!model)return;
    switch(env){
    case UWB_ENV_RESIDENTIAL_LOS: model->pl0_db=43.9;model->path_loss_exponent=1.79;model->shadowing_stddev_db=2.22;break;
    case UWB_ENV_RESIDENTIAL_NLOS:model->pl0_db=48.7;model->path_loss_exponent=4.58;model->shadowing_stddev_db=3.51;break;
    case UWB_ENV_OFFICE_LOS:      model->pl0_db=36.6;model->path_loss_exponent=1.63;model->shadowing_stddev_db=1.90;break;
    case UWB_ENV_OFFICE_NLOS:     model->pl0_db=51.4;model->path_loss_exponent=3.07;model->shadowing_stddev_db=3.90;break;
    case UWB_ENV_INDUSTRIAL_LOS:  model->pl0_db=56.7;model->path_loss_exponent=1.20;model->shadowing_stddev_db=6.00;break;
    case UWB_ENV_INDUSTRIAL_NLOS: model->pl0_db=56.7;model->path_loss_exponent=2.15;model->shadowing_stddev_db=6.00;break;
    default: model->pl0_db=44.0;model->path_loss_exponent=2.0;model->shadowing_stddev_db=3.0;break;
    }
}

void uwb_channel_compute_statistics(const uwb_sv_channel_t *channel,
                                    uwb_channel_statistics_t *stats) {
    int i;
    double sum_e=0.0,sum_t=0.0,sum_t2=0.0;
    if(!channel||!stats)return;
    memset(stats,0,sizeof(*stats));
    for(i=0;i<channel->num_paths;i++){
        double e=channel->paths[i].amplitude*channel->paths[i].amplitude;
        double t=channel->paths[i].delay_ns;
        sum_e+=e; sum_t+=e*t; sum_t2+=e*t*t;
    }
    if(sum_e<1e-20)return;
    stats->total_energy=sum_e;
    stats->mean_excess_delay_ns=sum_t/sum_e;
    double var=sum_t2/sum_e-stats->mean_excess_delay_ns*stats->mean_excess_delay_ns;
    stats->rms_delay_spread_ns=(var>0.0)?sqrt(var):0.0;
    stats->coherence_bandwidth_hz=(stats->rms_delay_spread_ns>0.0)?
        1.0/(5.0*stats->rms_delay_spread_ns*1e-9):1e12;
    stats->nakagami_m=uwb_channel_nakagami_m(channel->environment);
}

void uwb_channel_apply_shadowing(uwb_sv_channel_t *channel,
                                 double shadowing_std_db, uint32_t seed) {
    int i; double shadow,scale;
    if(!channel)return;
    srand(seed);
    shadow=((double)rand()/RAND_MAX-0.5)*2.0*shadowing_std_db;
    scale=pow(10.0,shadow/20.0);
    for(i=0;i<channel->num_paths;i++) channel->paths[i].amplitude*=scale;
}

double uwb_channel_nakagami_m(uwb_channel_environment_t env) {
    switch(env){case UWB_ENV_RESIDENTIAL_LOS:return 3.5;case UWB_ENV_OFFICE_LOS:return 2.5;default:return 1.5;}
}

void uwb_channel_to_cir(const uwb_sv_channel_t *sv_channel,
                        uwb_cir_t *cir, double sampling_period_ps, int num_samples) {
    int i;
    if(!sv_channel||!cir)return;
    memset(cir,0,sizeof(*cir));
    cir->sampling_period_ps=sampling_period_ps;
    cir->num_samples=(num_samples<UWB_CIR_MAX_SAMPLES)?num_samples:UWB_CIR_MAX_SAMPLES;
    for(i=0;i<sv_channel->num_paths;i++){
        int idx=(int)(sv_channel->paths[i].delay_ns*1000.0/sampling_period_ps);
        if(idx>=0&&idx<cir->num_samples){
            double amp=sv_channel->paths[i].amplitude;
            cir->samples[idx].i+=(int32_t)(amp*cos(sv_channel->paths[i].phase_rad)*1000);
            cir->samples[idx].q+=(int32_t)(amp*sin(sv_channel->paths[i].phase_rad)*1000);
        }
    }
}

double uwb_friis_rx_power_dbm(double tx_power_dbm, double tx_gain_dbi,
                               double rx_gain_dbi, double distance_m,
                               const uwb_path_loss_model_t *model) {
    return tx_power_dbm+tx_gain_dbi+rx_gain_dbi-uwb_channel_path_loss_db(distance_m,model);
}

void uwb_channel_apply_doppler(uwb_sv_channel_t *channel,
                               double velocity_ms, double time_s) {
    int i; double fd_max;
    if(!channel||channel->center_frequency_hz<=0.0)return;
    fd_max=velocity_ms*channel->center_frequency_hz/UWB_C;
    for(i=0;i<channel->num_paths;i++)
        channel->paths[i].phase_rad+=2.0*M_PI*fd_max*time_s;
}

double uwb_bandwidth_for_accuracy(double required_accuracy_m, double snr_linear) {
    if(required_accuracy_m<=0.0||snr_linear<=0.0)return 1e12;
    return UWB_C/(2.0*M_PI*sqrt(2.0*snr_linear)*required_accuracy_m);
}
