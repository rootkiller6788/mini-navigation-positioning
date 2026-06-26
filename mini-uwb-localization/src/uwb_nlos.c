#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
/**
 * mini-uwb-localization: NLOS Detection and Mitigation
 * L5: NLOS detection, kurtosis/skewness, leading edge detection
 * L8: ML-based NLOS classification (logistic regression, decision tree)
 */
#include "uwb_nlos.h"
#include <math.h>
#include <string.h>

void nlos_extract_features(const uwb_cir_t *cir, nlos_features_t *features) {
    double sum_energy=0.0,sum_tau=0.0,sum_tau2=0.0;
    double fp_energy=0.0,peak_energy=0.0,total_energy=0.0;
    int i;
    if(!cir||!features)return;
    memset(features,0,sizeof(*features));
    for(i=0;i<cir->num_samples;i++){
        double mag2=(double)cir->samples[i].i*cir->samples[i].i+
                    (double)cir->samples[i].q*cir->samples[i].q;
        double tau=cir->sampling_period_ps*i*1e-3;
        total_energy+=mag2;
        sum_energy+=mag2;
        sum_tau+=mag2*tau;
        sum_tau2+=mag2*tau*tau;
        if(i==cir->first_path_index) fp_energy=mag2;
        if(i==cir->peak_path_index) peak_energy=mag2;
    }
    if(total_energy<1e-20)return;
    double mean_tau=sum_tau/total_energy;
    double mean_tau2=sum_tau2/total_energy;
    double var_tau=mean_tau2-mean_tau*mean_tau;
    features->rms_delay_spread_ns=(var_tau>0)?sqrt(var_tau):0.0;
    features->mean_excess_delay_ns=mean_tau;
    features->max_excess_delay_ns=cir->sampling_period_ps*cir->num_samples*1e-3;
    features->total_energy=total_energy;
    features->first_path_energy_ratio=(total_energy>0)?fp_energy/total_energy:0.0;
    features->peak_to_energy_ratio=(total_energy>0)?peak_energy/total_energy:0.0;
    features->fp_to_total_energy_ratio=features->first_path_energy_ratio;
    features->cir_growth_rate=0.0; /* from ranging meas, not CIR */
    features->rise_time_ps=cir->sampling_period_ps*5.0;
    features->fp_amplitude_decay=0.0;
    features->peak_before_fp_db=0.0;
    features->noise_std_energy=cir->noise_floor*cir->noise_floor*cir->num_samples;
    features->kurtosis=cir_compute_kurtosis(cir);
    features->skewness=cir_compute_skewness(cir);
    features->path_loss_db=(total_energy>1e-20)?-10.0*log10(total_energy):100.0;
}

void nlos_lrt_detect(const nlos_features_t *features,
                     nlos_detection_result_t *result) {
    double score;
    if(!features||!result)return;
    memset(result,0,sizeof(*result));
    result->classifier_used=NLOS_CLASSIFIER_JOINT_LRT;
    score=0.0;
    if(features->kurtosis>5.0) score+=0.4;
    else if(features->kurtosis>4.0) score+=0.2;
    else if(features->kurtosis>3.5) score+=0.1;
    if(features->mean_excess_delay_ns>25.0) score+=0.4;
    else if(features->mean_excess_delay_ns>15.0) score+=0.2;
    else if(features->mean_excess_delay_ns>10.0) score+=0.1;
    if(features->rms_delay_spread_ns>20.0) score+=0.2;
    result->is_nlos=(score>0.5)?1:0;
    result->nlos_probability=score;
    result->confidence=(score>0.8||score<0.2)?0.9:0.5;
    memcpy(&result->features,features,sizeof(nlos_features_t));
}

void nlos_decision_tree_detect(const nlos_features_t *features,
                               nlos_detection_result_t *result) {
    int is_nlos=0;
    if(!features||!result)return;
    memset(result,0,sizeof(*result));
    result->classifier_used=NLOS_CLASSIFIER_DECISION_TREE;
    if(features->kurtosis>4.5&&features->rms_delay_spread_ns>15.0) is_nlos=1;
    else if(features->skewness>0.5&&features->mean_excess_delay_ns>20.0) is_nlos=1;
    else if(features->first_path_energy_ratio<0.1&&
            features->total_energy>features->noise_std_energy*10) is_nlos=1;
    else if(features->cir_growth_rate>0.8) is_nlos=1;
    result->is_nlos=is_nlos;
    result->nlos_probability=is_nlos?0.85:0.15;
    result->confidence=0.75;
    memcpy(&result->features,features,sizeof(nlos_features_t));
}

void nlos_logistic_detect(const nlos_features_t *features,
                          const double *weights, double bias,
                          nlos_detection_result_t *result) {
    double z;
    if(!features||!weights||!result)return;
    memset(result,0,sizeof(*result));
    result->classifier_used=NLOS_CLASSIFIER_LOGISTIC;
    z=weights[0]*features->kurtosis+weights[1]*features->skewness+
      weights[2]*features->rms_delay_spread_ns+weights[3]*features->mean_excess_delay_ns+
      weights[4]*features->first_path_energy_ratio+weights[5]*features->cir_growth_rate+bias;
    double prob=1.0/(1.0+exp(-z));
    result->nlos_probability=prob;
    result->is_nlos=(prob>0.5)?1:0;
    result->confidence=fabs(prob-0.5)*2.0;
    memcpy(&result->features,features,sizeof(nlos_features_t));
}

void nlos_mitigate_range(uwb_ranging_meas_t *meas,
                         const nlos_detection_result_t *nlos_result,
                         const nlos_mitigation_config_t *config) {
    if(!meas||!nlos_result||!config)return;
    if(!nlos_result->is_nlos||nlos_result->nlos_probability<config->nlos_prob_threshold)
        return;
    switch(config->strategy){
    case NLOS_MITIGATE_REJECT:
        meas->quality=UWB_RANGE_QUALITY_REJECT; break;
    case NLOS_MITIGATE_WEIGHT:
        meas->distance_variance*=config->nlos_variance_scale;
        if(meas->quality<UWB_RANGE_QUALITY_POOR) meas->quality=UWB_RANGE_QUALITY_POOR;
        break;
    case NLOS_MITIGATE_CORRECT:
        meas->distance_m-=config->nlos_range_bias_m;
        if(meas->distance_m<0.0) meas->distance_m=0.0;
        meas->distance_variance*=config->nlos_variance_scale;
        break;
    case NLOS_MITIGATE_SMOOTH:
        meas->distance_variance*=config->smoothing_factor; break;
    }
}

double nlos_estimate_range_bias(const nlos_features_t *features) {
    if(!features)return 0.0;
    double log_rms=(features->rms_delay_spread_ns>1.0)?
                    log10(features->rms_delay_spread_ns):0.0;
    double bias=0.15+0.08*features->kurtosis+0.12*features->skewness+0.10*log_rms;
    return (bias>0.0)?bias:0.0;
}

double cir_compute_rms_delay_spread(const uwb_cir_t *cir) {
    double sum_e=0.0,sum_t=0.0,sum_t2=0.0;
    int i;
    if(!cir||cir->num_samples==0)return 0.0;
    for(i=0;i<cir->num_samples;i++){
        double mag2=(double)cir->samples[i].i*cir->samples[i].i+
                    (double)cir->samples[i].q*cir->samples[i].q;
        double tau=cir->sampling_period_ps*i*1e-3;
        sum_e+=mag2; sum_t+=mag2*tau; sum_t2+=mag2*tau*tau;
    }
    if(sum_e<1e-20)return 0.0;
    double mean=sum_t/sum_e, var=sum_t2/sum_e-mean*mean;
    return (var>0.0)?sqrt(var):0.0;
}

double cir_compute_kurtosis(const uwb_cir_t *cir) {
    double sum=0.0,sum2=0.0,sum4=0.0;
    int i;
    if(!cir||cir->num_samples<4)return 3.0;
    for(i=0;i<cir->num_samples;i++){
        double mag=sqrt((double)cir->samples[i].i*cir->samples[i].i+
                         (double)cir->samples[i].q*cir->samples[i].q);
        sum+=mag; sum2+=mag*mag;
    }
    double n=cir->num_samples;
    double mean=sum/n, var=sum2/n-mean*mean;
    if(var<1e-20)return 3.0;
    for(i=0;i<cir->num_samples;i++){
        double mag=sqrt((double)cir->samples[i].i*cir->samples[i].i+
                         (double)cir->samples[i].q*cir->samples[i].q);
        double diff=mag-mean;
        sum4+=diff*diff*diff*diff;
    }
    return (sum4/n)/(var*var);
}

double cir_compute_skewness(const uwb_cir_t *cir) {
    double sum=0.0,sum2=0.0,sum3=0.0;
    int i;
    if(!cir||cir->num_samples<4)return 0.0;
    for(i=0;i<cir->num_samples;i++){
        double mag=sqrt((double)cir->samples[i].i*cir->samples[i].i+
                         (double)cir->samples[i].q*cir->samples[i].q);
        sum+=mag; sum2+=mag*mag;
    }
    double n=cir->num_samples;
    double mean=sum/n, var=sum2/n-mean*mean;
    if(var<1e-20)return 0.0;
    double sigma=sqrt(var);
    for(i=0;i<cir->num_samples;i++){
        double mag=sqrt((double)cir->samples[i].i*cir->samples[i].i+
                         (double)cir->samples[i].q*cir->samples[i].q);
        double diff=mag-mean;
        sum3+=diff*diff*diff;
    }
    return (sum3/n)/(sigma*sigma*sigma);
}

int cir_detect_leading_edge(const uwb_cir_t *cir, double threshold_factor,
                            uint16_t *leading_edge_index) {
    int i;
    double threshold;
    if(!cir||!leading_edge_index||cir->num_samples==0)return 0;
    threshold=cir->noise_floor*threshold_factor;
    for(i=0;i<cir->num_samples;i++){
        double mag2=(double)cir->samples[i].i*cir->samples[i].i+
                    (double)cir->samples[i].q*cir->samples[i].q;
        if(sqrt(mag2)>threshold){
            *leading_edge_index=i;
            return 1;
        }
    }
    return 0;
}
