// [[Rcpp::depends(RcppZiggurat)]]
#include <Rcpp.h>
#include <Ziggurat.h>
#include <algorithm>
#include <stdlib.h>
using namespace Rcpp;
static Ziggurat::Ziggurat::Ziggurat zigg;

//' Simulate DDM process by unique trial condition (2 items)
//' @author Alexander Fengler, \email{alexanderfengler@@gmx.de}
//' @title Simulate DDM process (by condition, 2 items)
//' \code{evacc2_by_condition()}
//' @return vector that stores decisions and rts for each simulation run
//' @param parameters vector that stores the parameters used for the simulations (Order: [non.decision.time, drift, sd])
//' @param timestep timestep in ms associated with each step in the drift diffusion process
//' @param nr_reps number of repitions (simulation runs)
//' @param maxdur maximum duration in ms that the process is allowed to simulate
//' @param update Vector that stores the item valuations for the trial conditon simulated
//' @param fixation_model placeholder for interface consistency / see attentional versions for specification
//' @param nr_attributes placeholder for interface consistency / see multiattribute versions for specification
//' @export
// [[Rcpp::export]]
IntegerVector evacc2_by_condition(NumericVector parameters,
                                  int maxdur,
                                  NumericVector update,
                                  int nr_attributes,
                                  Function fixation_model,
                                  int nr_reps,
                                  int timestep){

  // Set seed for random sampler ------------------------------------------------------------------
  NumericVector seed(1);
  seed = floor(runif(1,-100000,100000));
  zigg.setSeed(seed[0]);
  // ----------------------------------------------------------------------------------------------

  // Initialize Variable that collects output -----------------------------------------------------
  IntegerVector out(2*nr_reps);
  // ----------------------------------------------------------------------------------------------

  // Initialize parameters ------------------------------------------------------------------------
  int non_decision_time = parameters[0];
  float drift = parameters[1];
  float sd = parameters[2];
  // ----------------------------------------------------------------------------------------------

  // Initialize Variables needed to propagate model -----------------------------------------------
  int nr_items = update.size();
  float rdv;
  bool decision = 0;
  int cur_rt = 0;
  int out_cnt = -2; // index for output vector
  int out_plus = 0;
  float cur_drift = 0;

  for (int i = 0; i < nr_items; ++i){
    update[i] = update[i]*drift;
  }

  cur_drift = update[0] - update[1];
  // ----------------------------------------------------------------------------------------------

  // Outer loop cycles through simulation numbers -------------------------------------------------
  for (int rep_cnt = 0; rep_cnt < nr_reps;++rep_cnt){

    // Reset Variables before entering simulation -------------------------------------------------
    cur_rt = 0;
    out_cnt += 2;
    out_plus = out_cnt + 1;
    decision = 0;
    rdv = 0;
    // --------------------------------------------------------------------------------------------

    // Propagate model through simulation run -----------------------------------------------------
    while(cur_rt < maxdur){  //
        // Accumulate rdvence for current timestep
        rdv += cur_drift + sd*zigg.norm();

        //update RT
        cur_rt += timestep;

        // check whether decision made
        if (rdv >= 1 || rdv <= -1){
          decision = 1;
          break;
        }
      }
    // --------------------------------------------------------------------------------------------

    // Collect outputs ----------------------------------------------------------------------------
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
