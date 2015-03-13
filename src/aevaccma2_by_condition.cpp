// [[Rcpp::depends(RcppZiggurat)]]
#include <Rcpp.h>
#include <Ziggurat.h>
#include <algorithm>
#include <stdlib.h>
using namespace Rcpp;
static Ziggurat::Ziggurat::Ziggurat zigg;

//' Simulate aDDM process by unique trial condition (2 items, multiattribute)
//' @author Alexander Fengler, \email{alexanderfengler@@gmx.de}
//' @title Simulate aDDM process (by condition, 2 items)
//' \code{aevaccma2_by_condition()}
//' @return vector that stores decisions and rts for each simulation run
//' @param parameters vector that stores the parameters used for the simulations (Order: [non.decision.time, drift, sd, theta, gamma])
//' @param timestep timestep in ms associated with each step in the drift diffusion process
//' @param nr_reps number of repitions (simulation runs)
//' @param maxdur maximum duration in ms that the process is allowed to simulate
//' @param update Vector that stores the item valuations for the trial conditon simulated
//' @param fixation_model a user supplied fixation model that will be utilized to supply fixation locations and potentially fixation durations
//' @export
// [[Rcpp::export]]
IntegerVector aevaccma2_by_condition(NumericVector parameters,
                                     int maxdur = 0,
                                     NumericVector update = 0,
                                     int nr_attributes = 1,
                                     Function fixation_model = 0,
                                     int nr_reps = 0,
                                     int timestep = 0){

  // Set seed for random sampler ------------------------------------------------------------------
  NumericVector seed(1);
  seed = floor(runif(1,-100000,100000));
  zigg.setSeed(seed[0]);
  // ----------------------------------------------------------------------------------------------

  // Initialize parameters ------------------------------------------------------------------------
  int non_decision_time = parameters[0];
  float drift = parameters[1];
  float sd = parameters[2];
  float theta = parameters[3];
  float gamma = parameters[4];
  // ----------------------------------------------------------------------------------------------

  // Initialize Variable that collects output -----------------------------------------------------
  IntegerVector out(2*nr_reps);
  // ----------------------------------------------------------------------------------------------

  // Initialize Variables needed to propagate model -----------------------------------------------
  int nr_screen_positions = update.size();
  int nr_items = nr_screen_positions / nr_attributes;
  float rdv;
  bool decision = 0;
  int cur_rt = 0;
  int out_cnt = -2; // index for output vector
  int out_plus = 0;
  int cur_fix_cnt = 0;
  int cur_fixpos_indice = 0;

  NumericMatrix fixdat = fixation_model(update);
  // ------------------------------------------------------------------------------------------------

  // Precompute drifts for each fixation position ----------------------------------------------------
  NumericMatrix updatemat(nr_attributes,2);
  NumericVector drifts(2*nr_attributes);

  // generate matrix that stores attribute wise values (rows) by item (columns)
  for(int i = 0; i < nr_items; i++){
    for (int j = 0; j < nr_attributes; j++){
      updatemat(j,i) = update[((nr_attributes*i)+j)];
    }
  }

  // use matrix to fill in drifts vector
  for (int i = 0; i < nr_items; i++){
    for (int j = 0; j < nr_attributes; j++){
      for (int k = 0; k <nr_attributes; k++){
        if (j == k){
          if (i == 0){
            drifts[((nr_attributes*i)+j)] += updatemat(j,i) - theta*updatemat(j,1);
          } else {
            drifts[((nr_attributes*i)+j)] += theta*updatemat(j,0) - updatemat(j,i);
          }
        } else {
          if (i == 0){
            drifts[((nr_attributes*i)+j)] += gamma*(updatemat(j,i) - theta*updatemat(j,1));
          } else {
            drifts[((nr_attributes*i)+j)] += gamma*(theta*updatemat(j,0) - updatemat(j,i));
          }
        }
      }
    }
  }

  // scale values in drifts vector according to current drift rate
  for (int i = 0; i < drifts.size(); i++){
    drifts[i] = drifts[i]*drift;
  }
  // -------------------------------------------------------------------------------------------------

  // Outer loop cycles through simulation numbers ---------------------------------------------------
  for (int rep_cnt = 0; rep_cnt < nr_reps;++rep_cnt){

    // Reset Variables before entering simulation ---------------------------------------------------
    cur_rt = 0;
    cur_fix_cnt = 0;
    out_cnt += 2;
    out_plus = out_cnt + 1;
    decision = 0;
    rdv = 0;

    // Compute fixation path
    fixdat = fixation_model(update);
    // ---------------------------------------------------------------------------------------------

    // Propagate model through simulation run ------------------------------------------------------
    for (int fix_cnt = 0; fix_cnt < 1000; ++fix_cnt){

      // Make updates according to new fixation location ---------------------------------------
      cur_fixpos_indice = fixdat(fix_cnt,0) - 1;
      // ---------------------------------------------------------------------------------------

      // Propagate Model through current fixation ----------------------------------------------
      for (int rt_cur_fix = 0; rt_cur_fix < fixdat(fix_cnt,1);rt_cur_fix += timestep){

        // Accumulate rdvence for current timestep
        rdv += drifts[cur_fixpos_indice] + sd*zigg.norm();

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

