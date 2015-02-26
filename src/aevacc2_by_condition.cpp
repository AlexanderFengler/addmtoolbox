// [[Rcpp::depends(RcppZiggurat)]]
#include <Rcpp.h>
#include <Ziggurat.h>
#include <algorithm>
#include <stdlib.h>
using namespace Rcpp;
static Ziggurat::Ziggurat::Ziggurat zigg;

//' Runs evidence accumulation function (2 item case) for one trial condition
//' @author Alexander Fengler, \email{alexanderfengler@@gmx.de}
//' @title evidence accumulation (2 item case)
//' \code{aevacc2_by_condition} Returns data.table with log likelihoods and corresponding parameter combinations
//' @return Returns a vector that stores decisions and rts for each simulation run
//' @param sd standard deviation used for drift diffusion process
//' @param theta theta used for drift diffusion process
//' @param drift drift-rate used for drift diffusion process
//' @param non_decision_time non decision time used for drift diffusion process
//' @param timestep timestep in ms associated with each step in the drift diffusion process
//' @param nr_reps number of repitions (simulation runs)
//' @param maxdur maximum duration in ms that the process is allowed to simulate
//' @param update Vector that stores the item valuations for the trial conditon simulated
//' @param fixation_model a user supplied fixation model that will be utilized to supply fixation locations and potentially fixation durations
//' @export
// [[Rcpp::export]]
IntegerVector aevacc2_by_condition(float sd,
                                   float theta,
                                   float drift,
                                   int non_decision_time,
                                   int timestep,
                                   int nr_reps,
                                   int maxdur,
                                   NumericVector update,
                                   Function fixation_model){

  // Set seed for random sampler ------------------------------------------------------------------
  NumericVector seed(1);
  seed = floor(runif(1,-100000,100000));
  zigg.setSeed(seed[0]);
  // ----------------------------------------------------------------------------------------------

  // Initialize Variable that collects output -----------------------------------------------------
  IntegerVector out(2*nr_reps);
  // ----------------------------------------------------------------------------------------------

  // Initialize Variables needed to propagate model -----------------------------------------------
  int nr_items = update.size();
    float rdv;
  bool decision = 0;
  int cur_rt = 0;
  int out_cnt = -2; // index for output vector
  int out_plus = 0;
  int cur_fix_cnt = 0;
  int cur_fixpos_indice = 0;

  NumericVector cur_update(nr_items);
  for (int i = 0; i < nr_items; ++i){
    update[i] = update[i]*drift;
    cur_update[i] = theta*update[i];
  }

  NumericMatrix fixdat = fixation_model();
  // ------------------------------------------------------------------------------------------------

  // Outer loop cycles through simulation numbers --------------------------------------------------
  for (int rep_cnt = 0; rep_cnt < nr_reps;++rep_cnt){

    // Reset Variables before entering simulation ------------------------------------------------
    cur_rt = 0;
    cur_fix_cnt = 0;
    out_cnt += 2;
    out_plus = out_cnt + 1;
    decision = 0;
    rdv = 0;

    // Compute fixation path
    fixdat = fixation_model();
    // -------------------------------------------------------------------------------------------

    // Propagate model through simulation run --------------------------------------------------
    for (int fix_cnt = 0; fix_cnt < 1000; ++fix_cnt){  //

      // Make updates according to new fixation location ---------------------------------------
      cur_fixpos_indice = fixdat(fix_cnt,0) - 1;
      cur_update[cur_fixpos_indice] = update[cur_fixpos_indice];
      // ---------------------------------------------------------------------------------------

      // Propagate Model through current fixation ----------------------------------------------
      for (int rt_cur_fix = 0; rt_cur_fix < fixdat(fix_cnt,1);rt_cur_fix += timestep){

        // Accumulate rdvence for current timestep
        rdv += (cur_update[0] - cur_update[1]) + sd*zigg.norm();

        //update RT
        cur_rt += timestep;

        // check whether decision made
        if (abs(rdv) >= 1){
          decision = 1;
        }

        // if decision taken or maximum duration reached we break out of loop
        if ((decision == 1) || (cur_rt > maxdur)){
          break;
        }
      }

      // revert back cur_update to neutral state (times theta for update everywhere)
      cur_update[cur_fixpos_indice] = cur_update[cur_fixpos_indice]*theta;

      cur_fix_cnt +=  1;

      if ((decision == 1) || (cur_rt >= maxdur)){
        break;
      }
    }
    // ------------------------------------------------------------------------------------------

    // store decision
    if (rdv < 0){
      out[out_cnt] = 2;
    } else{
      out[out_cnt] = 1;
    }
    // store reaction time
    if(decision == 0){
      out[out_plus] = maxdur + non_decision_time;
    } else {
      out[out_plus] = cur_rt + non_decision_time;
    }
    // --------------------------------------------------------------------------------------------
  }
  // ----------------------------------------------------------------------------------------------
  return out;
}

/// [[Rcpp::export]]
//void myzsetseed(unsigned long int x) {
//  zigg.setSeed(x);
//return; }
