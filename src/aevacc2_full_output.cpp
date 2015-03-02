// [[Rcpp::depends(RcppZiggurat)]]
#include <Rcpp.h>
#include <algorithm>
#include <stdlib.h>
#include <Ziggurat.h>
using namespace Rcpp;
static Ziggurat::Ziggurat::Ziggurat zigg;

//' Evidence accumulation function for one trial/condition, and returns detailed model output (2 items)
//' @author Alexander Fengler, \email{alexanderfengler@@gmx.de}
//' @title Evidence accumulation with detailed output (2 items)
//' \code{aevacc2_full_output}
//' @return Returns a vector that stores detailed output by simulation run
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
NumericVector aevacc2_full_output(float sd,
                                  float theta,
                                  float drift,
                                  int non_decision_time,
                                  int maxdur,
                                  NumericVector update,
                                  Function fixation_model,
                                  int nr_reps,
                                  int timestep){

  // Set seed for random sampler ------------------------------------------------------------------
  NumericVector seed(1);
  seed = floor(runif(1,-100000,100000));
  zigg.setSeed(seed[0]);
  // ----------------------------------------------------------------------------------------------

  // Output Collection Variables ------------------------------------------------------------------
  int nr_items = update.size();
  NumericVector out((6+2*nr_items)*nr_reps);
  int maxpos = 0;
  bool decision = 1;
  int cur_rt = 0;
  int item_last_attended = 0;
  float value_last_attended = 0;
  int last_attended_chosen = 0;
  NumericVector Durations(nr_items);
  NumericVector Fixations(nr_items);

  NumericMatrix fixdat = fixation_model();
  // ----------------------------------------------------------------------------------------------

  // Initialization of Varibales needed in loop ---------------------------------------------------
  float rdv = 0;
  int out_cnt = -6 - nr_items - nr_items;   // output counter
  int cur_fixpos_indice = 0;             // storing indice of current fixation (-1 to make it vector indice)
  int cur_fix_cnt = 0;                   // counts number of fixations used

  NumericVector cur_update(nr_items);    // storing current updates
  // ------------------------------------------------------------------------------------------------

  for (int i = 0; i < nr_items; ++i){
    update[i] = update[i]*drift;
    cur_update[i] = theta*update[i];
  }
  // ------------------------------------------------------------------------------------------------

  // Loop cycling through simulations ---------------------------------------------------------------
  for (int rep_cnt = 0; rep_cnt < nr_reps;++rep_cnt){
    // reset ----------------------------------------------------------------------------------------
    cur_rt = 0;
    cur_fix_cnt = 0;
    out_cnt += 6 + 2*nr_items;
    rdv = 0;
    decision = 0;

    for (int i = 0; i < nr_items; ++i){
      Durations[i] = 0;
      Fixations[i] = 0;
    }
    // -----------------------------------------------------------------------------------------------

    // Generate fixations for current simulation run -------------------------------------------------
    fixdat = fixation_model();
    // -----------------------------------------------------------------------------------------------

    // Enter simulation run --------------------------------------------------------------------------
    for (int fix_cnt = 0; fix_cnt < 1000; ++fix_cnt){  // the < 1000 condition is arbitrary, "maxdur" is used to break out of loop
      // Make updates according to new fixation location ---------------------------------------------

      // adjust cur_update vector to reflect current fixation position
      cur_fixpos_indice = (fixdat(fix_cnt,0) - 1);
      cur_update[cur_fixpos_indice] = update[cur_fixpos_indice];

      // Update Fixations vector
      Fixations[fixdat(fix_cnt,0)-1]++;
      // -----------------------------------------------------------------------------------------------

      // Propagate Model -------------------------------------------------------------------------------
      for (int rt_cur_fix = 0; rt_cur_fix < fixdat(fix_cnt,1);rt_cur_fix += timestep){
        // update Duration on Item
        Durations[fixdat(fix_cnt,0)-1] += timestep;

        // Accumulate rdvence for current timestep
        rdv += (cur_update[0] - cur_update[1]) + sd*zigg.norm();

        //update RT
        cur_rt += timestep;

        // check whether decision made
        if (rdv >= 1 || rdv <= -1){
          decision = 1;
        }

        // if decision taken or maximum duration reached we break out of loop
        if ((decision == 1) || (cur_rt > maxdur)){
          break;
        }
      }
      // -------------------------------------------------------------------------------------------------

      // revert back cur_update to neutral state (times theta for update everywhere)
      cur_update[cur_fixpos_indice] = cur_update[cur_fixpos_indice]*theta;

      // update total fixation count
      cur_fix_cnt +=  1;

      // Check whether decision has been taken or maximum duration has been reached at this point --------
      if ((decision == 1) || (cur_rt >= maxdur)){
        // Check which side has been chosen --------------------------------------------------------------
        if (rdv < 0){
          maxpos = 1;
        } else {
          maxpos = 0;
        }
        // -----------------------------------------------------------------------------------------------
        item_last_attended = fixdat(fix_cnt,0);
        value_last_attended = update[fixdat(fix_cnt,0) - 1];
        last_attended_chosen = 0;
        if (item_last_attended == (maxpos + 1)){
          last_attended_chosen = 1;
        }
        break;
      }
    }
    // --------------------------------------------------------------------------------------------------

    // Storing output -----------------------------------------------------------------------------------
    out[out_cnt] = maxpos + 1;
    out[out_cnt + 1] = cur_fix_cnt;

    if(decision == 0){
      out[out_cnt + 2] = maxdur;
    }
    else {
      out[out_cnt + 2] = cur_rt + non_decision_time;
    }

    out[out_cnt + 3] = item_last_attended;
    out[out_cnt + 4] = value_last_attended/drift;
    out[out_cnt + 5] = last_attended_chosen;

    for(int i = 0; i < nr_items;i++){
      out[out_cnt + 6 + i] = Durations[i];
      out[out_cnt + 6 + nr_items + i] = Fixations[i];
    }
    // ---------------------------------------------------------------------------------------------------
  }
  // -----------------------------------------------------------------------------------------------------
  return out;
}

//// [[Rcpp::export]]
//void myzsetseed(unsigned long int x) {
//    zigg.setSeed(x);
//return; }
