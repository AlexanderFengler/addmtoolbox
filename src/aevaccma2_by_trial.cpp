// [[Rcpp::depends(RcppZiggurat)]]

#include <Rcpp.h>
#include <Ziggurat.h>
#include <algorithm>
#include <stdlib.h>
using namespace Rcpp;
static Ziggurat::Ziggurat::Ziggurat zigg;

//' Simulate aDDM process by unique trial (2 items // multiattribute)
//' \code{aevaccma2_by_trial()}
//' @author Alexander Fengler, \email{alexanderfengler@@gmx.de}
//' @title Simulate aDDM process (by trial, 2 items, multiattribute)
//' @return numeric variable that provides a success count (runs that predicted a reaction time in the correct rt-bin and simultaneously the correct decision)
//' @param parameters vector that stores the parameters used for the simulations (Order: [non.decision.time, drift, sd, theta, gamma])
//' @param timestep timestep in ms associated with each step in the drift diffusion process
//' @param nr_reps number of repitions (simulation runs)
//' @param maxdur numeric variable that supplies the maximum reaction time considered a success in simulations
//' @param mindur numeric variable that supplies the minimum reaction time considered a succes in simulations
//' @param cur_decision numeric variable that provides the empirical decision taken in trial
//' @param update Vector that stores the item valuations for the trial conditon simulated
//' @param fixpos Vector that stores the locations for a supplied fixed fixation pathway
//' @param fixdur Vector that stores the fixation durations for a supplied fixed fixation pathway
//' @param cur_maxfix integer that provides number of fixation in trial
//' @export
// [[Rcpp::export]]
int aevaccma2_by_trial(NumericVector parameters,
                       int maxdur,
                       int mindur,
                       int cur_decision,
                       NumericVector update,
                       int nr_attributes,
                       IntegerVector fixpos,
                       IntegerVector fixdur,
                       int nr_reps,
                       int timestep){

  // Set seed for random sampler --------------------------------------------------------------------
  NumericVector seed(1);
  seed = floor(runif(1,-100000,100000));
  zigg.setSeed(seed[0]);
  // ------------------------------------------------------------------------------------------------

  // Initialize parameters ------------------------------------------------------------------------
  int non_decision_time = parameters[0];
  float drift = parameters[1];
  float sd = parameters[2];
  float theta = parameters[3];
  float gamma = parameters[4];
  // ----------------------------------------------------------------------------------------------

  // Variable collecting success counts -------------------------------------------------------------
  int out = 0;
  // ------------------------------------------------------------------------------------------------

  // Initialize Variables need for model propagation ------------------------------------------------
  int nr_screen_positions = update.size();
  int nr_items = nr_screen_positions / nr_attributes;
  float rdv = 0;
  bool decision = 1;
  int maxpos = 0;
  int cur_rt = 0;
  int cur_fixpos_indice = 0;
  int num_fixpos = fixpos.size();
  int cur_fix_cnt = 0;
  // -------------------------------------------------------------------------------------------------

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

  // Cycle through simulations -----------------------------------------------------------------------
  for (int rep_cnt = 0; rep_cnt < nr_reps;++rep_cnt){

    // reset variables before entering next simulation -----------------------------------------------
    cur_rt = 0;
    cur_fix_cnt = 0;
    decision = 0;
    rdv = 0;
    // ----------------------------------------------------------------------------------------------

    // run a simulation round -----------------------------------------------------------------------
    for (int fix_cnt = 0; fix_cnt < num_fixpos; ++fix_cnt){  // maybe start fix_cnt at one?

      // get current fixation position
      cur_fixpos_indice = fixpos[fix_cnt] - 1;


      // propagate model ----------------------------------------------------------------------------
      for (int rt_cur_fix = 0; rt_cur_fix < fixdur[fix_cnt];rt_cur_fix += timestep){
        //update rt
        cur_rt += timestep;

        // accumulate evidence
        rdv += drifts[cur_fixpos_indice] + sd*zigg.norm();

        // check whether decision made
        if (rdv >= 1 || rdv <= -1){
          decision = 1;
        }

        // if decision taken or maximum duration reached we break out of loop
        if ((decision == 1) || (cur_rt > maxdur)){
          break;
        }
      }

      // count fixation upwards
      cur_fix_cnt +=  1;

      // if decision taken or maximum duration reached we break out of loop
      if ((decision == 1) || (cur_rt > maxdur)){
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
