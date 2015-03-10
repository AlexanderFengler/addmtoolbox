// [[Rcpp::depends(RcppZiggurat)]]

#include <Rcpp.h>
#include <Ziggurat.h>
#include <algorithm>
#include <stdlib.h>
using namespace Rcpp;
static Ziggurat::Ziggurat::Ziggurat zigg;

//' Simulate DDM process by unique trial (2 items)
//' \code{evacc2_by_trial()}
//' @author Alexander Fengler, \email{alexanderfengler@@gmx.de}
//' @title Simulate DDM process (by trial, 2 items)
//' @return numeric variable that provides a success count (runs that predicted a reaction time in the correct rt-bin and simultaneously the correct decision)
//' @param sd standard deviation used for drift diffusion process
//' @param theta placeholder for interface consistency / see attentional versions for specification
//' @param drift drift-rate used for drift diffusion process
//' @param non_decision_time non decision time used for drift diffusion process
//' @param timestep timestep in ms associated with each step in the drift diffusion process
//' @param nr_reps number of repitions (simulation runs)
//' @param maxdur numeric variable that supplies the maximum reaction time considered a success in simulations
//' @param mindur numeric variable that supplies the minimum reaction time considered a succes in simulations
//' @param cur_decision numeric variable that provides the empirical decision taken in trial
//' @param update Vector that stores the item valuations for the trial conditon simulated
//' @param fixpos Vector placeholder for interface consistency / see attentional versions for specification
//' @param fixdur Vector placeholder for interface consistency / see attentional versions for specification
//' @param cur_maxfix integer that provides number of fixation in trial
//' @param gamma placeholder for interface consistency / see multiattribute versions for specification
//' @param nr_attributes placeholder for interface consistency / see multiattribute versions for specification
//' @export
// [[Rcpp::export]]
int evacc2_by_trial(float sd = 0,
                     float theta = 0,
                     float gamma = 0,
                     float drift = 0,
                     int non_decision_time = 0,
                     int maxdur = 0,
                     int mindur = 0,
                     int cur_decision = 0,
                     NumericVector update = 0,
                     int nr_attributes = 0,
                     IntegerVector fixpos = 0,
                     IntegerVector fixdur = 0,
                     int nr_reps = 0,
                     int timestep = 0){

  // Set seed for random sampler --------------------------------------------------------------------
  NumericVector seed(1);
  seed = floor(runif(1,-100000,100000));
  zigg.setSeed(seed[0]);
  // ------------------------------------------------------------------------------------------------

  // Variable collecting success counts -------------------------------------------------------------
  int out = 0;
  // ------------------------------------------------------------------------------------------------

  // Initialize Variables need for model propagation ------------------------------------------------
  int nr_items = update.size();
  float rdv = 0;
  bool decision = 1;
  int maxpos = 0;
  int cur_rt = 0;
  float cur_drift = 0;

  for (int i = 0; i < nr_items; ++i){
    update[i] = update[i]*drift;
  }

  cur_drift = update[0] - update[1];
  // ------------------------------------------------------------------------------------------------

  // Cycle through simulations ----------------------------------------------------------------------
  for (int rep_cnt = 0; rep_cnt < nr_reps;++rep_cnt){

    // reset variables before entering next simulation ----------------------------------------------
    cur_rt = 0;
    decision = 0;
    rdv = 0;
    // ----------------------------------------------------------------------------------------------

    // Propagate model through simulation run -------------------------------------------------------
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
    // ----------------------------------------------------------------------------------------------

    // determine whether outcome is success (correct rt and decision) -------------------------------
    if (decision == 1){
      if (rdv < 0){
        maxpos = 1;
      } else {
        maxpos = 0;
      }
      if (maxpos == cur_decision){
        if ((cur_rt + non_decision_time <= maxdur) && (cur_rt + non_decision_time >= mindur)){
          ++out;
        }
      }
    }
    // ----------------------------------------------------------------------------------------------
  }
  // ------------------------------------------------------------------------------------------------
  return out;
}

//// [[Rcpp::export]]
//void myzfastsetseed(unsigned long int x){
//    zigg.setSeed(x);
//return;}

