// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppGSL.h>
#include <Rcpp.h>

using namespace Rcpp;

// aevacc2_by_condition
IntegerVector aevacc2_by_condition(float sd, float theta, float drift, int non_decision_time, int maxdur, NumericVector update, Function fixation_model, int nr_reps, int timestep);
RcppExport SEXP addmtoolbox_aevacc2_by_condition(SEXP sdSEXP, SEXP thetaSEXP, SEXP driftSEXP, SEXP non_decision_timeSEXP, SEXP maxdurSEXP, SEXP updateSEXP, SEXP fixation_modelSEXP, SEXP nr_repsSEXP, SEXP timestepSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< float >::type sd(sdSEXP );
        Rcpp::traits::input_parameter< float >::type theta(thetaSEXP );
        Rcpp::traits::input_parameter< float >::type drift(driftSEXP );
        Rcpp::traits::input_parameter< int >::type non_decision_time(non_decision_timeSEXP );
        Rcpp::traits::input_parameter< int >::type maxdur(maxdurSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type update(updateSEXP );
        Rcpp::traits::input_parameter< Function >::type fixation_model(fixation_modelSEXP );
        Rcpp::traits::input_parameter< int >::type nr_reps(nr_repsSEXP );
        Rcpp::traits::input_parameter< int >::type timestep(timestepSEXP );
        IntegerVector __result = aevacc2_by_condition(sd, theta, drift, non_decision_time, maxdur, update, fixation_model, nr_reps, timestep);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// aevacc2_by_trial
int aevacc2_by_trial(float sd, float theta, float drift, int non_decision_time, int maxdur, int mindur, int cur_decision, NumericVector update, IntegerVector fixpos, IntegerVector fixdur, int nr_reps, int timestep);
RcppExport SEXP addmtoolbox_aevacc2_by_trial(SEXP sdSEXP, SEXP thetaSEXP, SEXP driftSEXP, SEXP non_decision_timeSEXP, SEXP maxdurSEXP, SEXP mindurSEXP, SEXP cur_decisionSEXP, SEXP updateSEXP, SEXP fixposSEXP, SEXP fixdurSEXP, SEXP nr_repsSEXP, SEXP timestepSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< float >::type sd(sdSEXP );
        Rcpp::traits::input_parameter< float >::type theta(thetaSEXP );
        Rcpp::traits::input_parameter< float >::type drift(driftSEXP );
        Rcpp::traits::input_parameter< int >::type non_decision_time(non_decision_timeSEXP );
        Rcpp::traits::input_parameter< int >::type maxdur(maxdurSEXP );
        Rcpp::traits::input_parameter< int >::type mindur(mindurSEXP );
        Rcpp::traits::input_parameter< int >::type cur_decision(cur_decisionSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type update(updateSEXP );
        Rcpp::traits::input_parameter< IntegerVector >::type fixpos(fixposSEXP );
        Rcpp::traits::input_parameter< IntegerVector >::type fixdur(fixdurSEXP );
        Rcpp::traits::input_parameter< int >::type nr_reps(nr_repsSEXP );
        Rcpp::traits::input_parameter< int >::type timestep(timestepSEXP );
        int __result = aevacc2_by_trial(sd, theta, drift, non_decision_time, maxdur, mindur, cur_decision, update, fixpos, fixdur, nr_reps, timestep);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// aevacc2_full_output
NumericVector aevacc2_full_output(float sd, float theta, float drift, int non_decision_time, int maxdur, NumericVector update, Function fixation_model, int nr_reps, int timestep);
RcppExport SEXP addmtoolbox_aevacc2_full_output(SEXP sdSEXP, SEXP thetaSEXP, SEXP driftSEXP, SEXP non_decision_timeSEXP, SEXP maxdurSEXP, SEXP updateSEXP, SEXP fixation_modelSEXP, SEXP nr_repsSEXP, SEXP timestepSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< float >::type sd(sdSEXP );
        Rcpp::traits::input_parameter< float >::type theta(thetaSEXP );
        Rcpp::traits::input_parameter< float >::type drift(driftSEXP );
        Rcpp::traits::input_parameter< int >::type non_decision_time(non_decision_timeSEXP );
        Rcpp::traits::input_parameter< int >::type maxdur(maxdurSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type update(updateSEXP );
        Rcpp::traits::input_parameter< Function >::type fixation_model(fixation_modelSEXP );
        Rcpp::traits::input_parameter< int >::type nr_reps(nr_repsSEXP );
        Rcpp::traits::input_parameter< int >::type timestep(timestepSEXP );
        NumericVector __result = aevacc2_full_output(sd, theta, drift, non_decision_time, maxdur, update, fixation_model, nr_reps, timestep);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// aevacc_by_condition
IntegerVector aevacc_by_condition(float sd, float theta, float drift, int non_decision_time, int maxdur, NumericVector update, Function fixation_model, int nr_reps, int timestep);
RcppExport SEXP addmtoolbox_aevacc_by_condition(SEXP sdSEXP, SEXP thetaSEXP, SEXP driftSEXP, SEXP non_decision_timeSEXP, SEXP maxdurSEXP, SEXP updateSEXP, SEXP fixation_modelSEXP, SEXP nr_repsSEXP, SEXP timestepSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< float >::type sd(sdSEXP );
        Rcpp::traits::input_parameter< float >::type theta(thetaSEXP );
        Rcpp::traits::input_parameter< float >::type drift(driftSEXP );
        Rcpp::traits::input_parameter< int >::type non_decision_time(non_decision_timeSEXP );
        Rcpp::traits::input_parameter< int >::type maxdur(maxdurSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type update(updateSEXP );
        Rcpp::traits::input_parameter< Function >::type fixation_model(fixation_modelSEXP );
        Rcpp::traits::input_parameter< int >::type nr_reps(nr_repsSEXP );
        Rcpp::traits::input_parameter< int >::type timestep(timestepSEXP );
        IntegerVector __result = aevacc_by_condition(sd, theta, drift, non_decision_time, maxdur, update, fixation_model, nr_reps, timestep);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// aevacc_by_condition_memnoise
IntegerVector aevacc_by_condition_memnoise(float sd, float theta, float drift, int non_decision_time, float items_seen_bias, float items_seen_noise_bias, int maxdur, NumericVector update, Function fixation_model, int nr_reps, int timestep);
RcppExport SEXP addmtoolbox_aevacc_by_condition_memnoise(SEXP sdSEXP, SEXP thetaSEXP, SEXP driftSEXP, SEXP non_decision_timeSEXP, SEXP items_seen_biasSEXP, SEXP items_seen_noise_biasSEXP, SEXP maxdurSEXP, SEXP updateSEXP, SEXP fixation_modelSEXP, SEXP nr_repsSEXP, SEXP timestepSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< float >::type sd(sdSEXP );
        Rcpp::traits::input_parameter< float >::type theta(thetaSEXP );
        Rcpp::traits::input_parameter< float >::type drift(driftSEXP );
        Rcpp::traits::input_parameter< int >::type non_decision_time(non_decision_timeSEXP );
        Rcpp::traits::input_parameter< float >::type items_seen_bias(items_seen_biasSEXP );
        Rcpp::traits::input_parameter< float >::type items_seen_noise_bias(items_seen_noise_biasSEXP );
        Rcpp::traits::input_parameter< int >::type maxdur(maxdurSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type update(updateSEXP );
        Rcpp::traits::input_parameter< Function >::type fixation_model(fixation_modelSEXP );
        Rcpp::traits::input_parameter< int >::type nr_reps(nr_repsSEXP );
        Rcpp::traits::input_parameter< int >::type timestep(timestepSEXP );
        IntegerVector __result = aevacc_by_condition_memnoise(sd, theta, drift, non_decision_time, items_seen_bias, items_seen_noise_bias, maxdur, update, fixation_model, nr_reps, timestep);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// aevacc_by_trial
int aevacc_by_trial(int nr_reps, int maxdur, int mindur, int cur_decision, float sd, float theta, float drift, int non_decision_time, int timestep, NumericVector update, IntegerVector fixpos, IntegerVector fixdur);
RcppExport SEXP addmtoolbox_aevacc_by_trial(SEXP nr_repsSEXP, SEXP maxdurSEXP, SEXP mindurSEXP, SEXP cur_decisionSEXP, SEXP sdSEXP, SEXP thetaSEXP, SEXP driftSEXP, SEXP non_decision_timeSEXP, SEXP timestepSEXP, SEXP updateSEXP, SEXP fixposSEXP, SEXP fixdurSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< int >::type nr_reps(nr_repsSEXP );
        Rcpp::traits::input_parameter< int >::type maxdur(maxdurSEXP );
        Rcpp::traits::input_parameter< int >::type mindur(mindurSEXP );
        Rcpp::traits::input_parameter< int >::type cur_decision(cur_decisionSEXP );
        Rcpp::traits::input_parameter< float >::type sd(sdSEXP );
        Rcpp::traits::input_parameter< float >::type theta(thetaSEXP );
        Rcpp::traits::input_parameter< float >::type drift(driftSEXP );
        Rcpp::traits::input_parameter< int >::type non_decision_time(non_decision_timeSEXP );
        Rcpp::traits::input_parameter< int >::type timestep(timestepSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type update(updateSEXP );
        Rcpp::traits::input_parameter< IntegerVector >::type fixpos(fixposSEXP );
        Rcpp::traits::input_parameter< IntegerVector >::type fixdur(fixdurSEXP );
        int __result = aevacc_by_trial(nr_reps, maxdur, mindur, cur_decision, sd, theta, drift, non_decision_time, timestep, update, fixpos, fixdur);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// aevacc_by_trial_memnoise
int aevacc_by_trial_memnoise(float sd, float theta, float drift, int non_decision_time, float items_seen_bias, float items_seen_noise_bias, int maxdur, int mindur, int cur_decision, NumericVector update, IntegerVector fixpos, IntegerVector fixdur, int nr_reps, int timestep);
RcppExport SEXP addmtoolbox_aevacc_by_trial_memnoise(SEXP sdSEXP, SEXP thetaSEXP, SEXP driftSEXP, SEXP non_decision_timeSEXP, SEXP items_seen_biasSEXP, SEXP items_seen_noise_biasSEXP, SEXP maxdurSEXP, SEXP mindurSEXP, SEXP cur_decisionSEXP, SEXP updateSEXP, SEXP fixposSEXP, SEXP fixdurSEXP, SEXP nr_repsSEXP, SEXP timestepSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< float >::type sd(sdSEXP );
        Rcpp::traits::input_parameter< float >::type theta(thetaSEXP );
        Rcpp::traits::input_parameter< float >::type drift(driftSEXP );
        Rcpp::traits::input_parameter< int >::type non_decision_time(non_decision_timeSEXP );
        Rcpp::traits::input_parameter< float >::type items_seen_bias(items_seen_biasSEXP );
        Rcpp::traits::input_parameter< float >::type items_seen_noise_bias(items_seen_noise_biasSEXP );
        Rcpp::traits::input_parameter< int >::type maxdur(maxdurSEXP );
        Rcpp::traits::input_parameter< int >::type mindur(mindurSEXP );
        Rcpp::traits::input_parameter< int >::type cur_decision(cur_decisionSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type update(updateSEXP );
        Rcpp::traits::input_parameter< IntegerVector >::type fixpos(fixposSEXP );
        Rcpp::traits::input_parameter< IntegerVector >::type fixdur(fixdurSEXP );
        Rcpp::traits::input_parameter< int >::type nr_reps(nr_repsSEXP );
        Rcpp::traits::input_parameter< int >::type timestep(timestepSEXP );
        int __result = aevacc_by_trial_memnoise(sd, theta, drift, non_decision_time, items_seen_bias, items_seen_noise_bias, maxdur, mindur, cur_decision, update, fixpos, fixdur, nr_reps, timestep);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// aevacc_full_output_memnoise
NumericVector aevacc_full_output_memnoise(float sd, float theta, float drift, int non_decision_time, float items_seen_bias, float items_seen_noise_bias, int maxdur, NumericVector update, Function fixation_model, int nr_reps, int timestep);
RcppExport SEXP addmtoolbox_aevacc_full_output_memnoise(SEXP sdSEXP, SEXP thetaSEXP, SEXP driftSEXP, SEXP non_decision_timeSEXP, SEXP items_seen_biasSEXP, SEXP items_seen_noise_biasSEXP, SEXP maxdurSEXP, SEXP updateSEXP, SEXP fixation_modelSEXP, SEXP nr_repsSEXP, SEXP timestepSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< float >::type sd(sdSEXP );
        Rcpp::traits::input_parameter< float >::type theta(thetaSEXP );
        Rcpp::traits::input_parameter< float >::type drift(driftSEXP );
        Rcpp::traits::input_parameter< int >::type non_decision_time(non_decision_timeSEXP );
        Rcpp::traits::input_parameter< float >::type items_seen_bias(items_seen_biasSEXP );
        Rcpp::traits::input_parameter< float >::type items_seen_noise_bias(items_seen_noise_biasSEXP );
        Rcpp::traits::input_parameter< int >::type maxdur(maxdurSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type update(updateSEXP );
        Rcpp::traits::input_parameter< Function >::type fixation_model(fixation_modelSEXP );
        Rcpp::traits::input_parameter< int >::type nr_reps(nr_repsSEXP );
        Rcpp::traits::input_parameter< int >::type timestep(timestepSEXP );
        NumericVector __result = aevacc_full_output_memnoise(sd, theta, drift, non_decision_time, items_seen_bias, items_seen_noise_bias, maxdur, update, fixation_model, nr_reps, timestep);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// dynamicaddm
double dynamicaddm(float sd, float theta, float drift, int non_decision_time, int decision, NumericVector valuations, NumericVector fixpos, NumericVector fixdur, int rt, float stateStep);
RcppExport SEXP addmtoolbox_dynamicaddm(SEXP sdSEXP, SEXP thetaSEXP, SEXP driftSEXP, SEXP non_decision_timeSEXP, SEXP decisionSEXP, SEXP valuationsSEXP, SEXP fixposSEXP, SEXP fixdurSEXP, SEXP rtSEXP, SEXP stateStepSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< float >::type sd(sdSEXP );
        Rcpp::traits::input_parameter< float >::type theta(thetaSEXP );
        Rcpp::traits::input_parameter< float >::type drift(driftSEXP );
        Rcpp::traits::input_parameter< int >::type non_decision_time(non_decision_timeSEXP );
        Rcpp::traits::input_parameter< int >::type decision(decisionSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type valuations(valuationsSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type fixpos(fixposSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type fixdur(fixdurSEXP );
        Rcpp::traits::input_parameter< int >::type rt(rtSEXP );
        Rcpp::traits::input_parameter< float >::type stateStep(stateStepSEXP );
        double __result = dynamicaddm(sd, theta, drift, non_decision_time, decision, valuations, fixpos, fixdur, rt, stateStep);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
