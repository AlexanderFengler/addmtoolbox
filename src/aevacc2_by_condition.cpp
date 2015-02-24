//' Runs evidence accumulation function (2 item case) for one trial condition
//' \code{aevacc2_by_condition} Returns data.table with log likelihoods and corresponding parameter combinations
//' @author Alexander Fengler, \email{alexanderfengler@@gmx.de}
//' @title evidence accumulation (2 item case)
//' @return Returns a vector that stores decisions and rts for each simulation run
//' @export
//' @param sd standard deviation used for drift diffusion process
//' @param theta theta used for drift diffusion process
//' @param non_decision_time non decision time used for drift diffusion process
//' @param timestep timestep in ms associated with each step in the drift diffusion process
//' @param nr_reps number of repitions (simulation runs)
//' @param maxdur maximum duration in ms that the process is allowed to simulate
//' @param update Vector that stores the item valuations for the trial conditon simulated
//' @param fixpos Vector that stores the locations for a supplied fixed fixation pathway
//' @param fixdur Vector that stores the fixation durations for a supplied fixed fixation pathway
//' @param fixdursamples Vector from which fixation duration can be sampled once supplied fixations run out
//' @param fixation_model a user supplied fixation model that will be utilized to supply fixation locations and potentially fixation durations

#include <Rcpp.h>
#include <Ziggurat.h>
#include <algorithm>
#include <stdlib.h>

using namespace Rcpp;
static Ziggurat::Ziggurat::Ziggurat zigg;

// [[Rcpp::depends(RcppZiggurat)]]
// [[Rcpp::export]]

IntegerVector aevacc2_by_condition(float sd,
                                   float theta,
                                   float drift,
                                   int non_decision_time,
                                   int timestep,
                                   int nr_reps,
                                   int maxdur,
                                   NumericVector update,
                                   IntegerVector fixpos,
                                   IntegerVector fixdur,
                                   IntegerVector fixdursamples,
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
  int num_fixpos = fixpos.size();
  IntegerVector cur_fixpos(1);
  IntegerVector cur_fixdur(1);
  IntegerVector temp_fixpos(1);
  int cur_fix_cnt = 0;
  int cur_fixpos_indice = 0;

  IntegerVector potential_fixpos(2);
  potential_fixpos[0] = 1;
  potential_fixpos[1] = 2;

  NumericVector cur_update(nr_items);
  for (int i = 0; i < nr_items; ++i){
    update[i] = update[i]*drift;
    cur_update[i] = theta*update[i];
  }
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

    // -------------------------------------------------------------------------------------------

    // Propagate model through simulation run --------------------------------------------------
    for (int fix_cnt = 0; fix_cnt < 1000; ++fix_cnt){  //

      // Handle next fixation ------------------------------------------------------------------
      // identifying current fixation position and duration
      // If empirical fixations still to be supplied fine,
      // otherwise sample from eligible fixations positions
      // and sample duration from fixdursamples vector.
      if (fix_cnt <= (num_fixpos - 1)){
        cur_fixpos[0] = fixpos[fix_cnt];
        cur_fixdur[0] = fixdur[fix_cnt];

        // get fixation duration
        if (fix_cnt == num_fixpos - 1){
          cur_fixdur = fixation_model(fixdursamples,1);
        }
      } else if (fix_cnt > (num_fixpos - 1)){

        temp_fixpos[0] = cur_fixpos[0];
        // sample fixation and make sure the new fixation is not the same as the old one
        while (cur_fixpos[0] == temp_fixpos[0]){
          cur_fixpos = fixation_model(potential_fixpos,1);
        }
        cur_fixdur = fixation_model(fixdursamples,1);
      }
      // ---------------------------------------------------------------------------------------

      // Make updates according to new fixation location ---------------------------------------
      cur_fixpos_indice = cur_fixpos[0] - 1;
      cur_update[cur_fixpos_indice] = update[cur_fixpos_indice];
      // ---------------------------------------------------------------------------------------

      // Propagate Model through current fixation ----------------------------------------------
      for (int rt_cur_fix = 0; rt_cur_fix < fixdur[fix_cnt];rt_cur_fix += timestep){

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
