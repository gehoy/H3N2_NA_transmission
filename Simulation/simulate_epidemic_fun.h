#ifndef SIMULATE_EPIDEMIC_FUN_PERF_PKG__H
#define SIMULATE_EPIDEMIC_FUN_PERF_PKG__H

#include <Rcpp.h>
using namespace Rcpp;

extern double mIncub;
extern double sdIncub;
extern double maxPCRDetectability;
extern double m_incub_g;
extern double sd_incub_g;
extern double shape_incub_g;
extern double scale_incub_g;


extern double mGamma;
extern double vGamma;
extern double shift;
extern double shape_infectivity;
extern double scale_infectivity;

extern double shape_gt;
extern double scale_gt;

extern double mu_f;
extern double sigma_f; 
extern double alpha_f;
extern double tau;

extern double cf;


RcppExport void test();
RcppExport double incubPeriod_g();
RcppExport double incubPeriod();
RcppExport double cumuldensityPeriod_g();
/*RcppExport double detectionPeriod();
RcppExport double detectionPeriod2(double infectionDate, double stFU); 
RcppExport double rIncub(double d);
RcppExport double rInfectionAsymptomatic(double d);*/
RcppExport NumericVector foi(
    double t, 
    double dt, 
    double lastDate,  

    NumericVector incub_periods,
    NumericVector symptomOnset, 
    NumericVector infection_onset,
    IntegerVector age, 
    IntegerVector haihk14_titer_q,
    IntegerVector na_titer_q,
    IntegerVector infectionStatus,
    IntegerVector stalk_titer_q, 
    double alpha,
    double beta,
    double delta,
    double muchild,
    double musHIUpper,
    double muNAUpper,
    double muStalkUpper,
    double muiHIUpper,
    double muiNAUpper,
    double muiStalkUpper,
    double hhsize,
    double mainHHSize,
    int display
);

#endif
