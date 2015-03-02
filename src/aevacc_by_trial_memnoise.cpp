// [[Rcpp::depends(RcppZiggurat)]]
#include <Rcpp.h>
#include <Ziggurat.h>
#include <algorithm>
#include <stdlib.h>
using namespace Rcpp;
static Ziggurat::Ziggurat::Ziggurat zigg;

//' Runs evidence accumulation for one trial (item general case). Flexible version allowing for memory effects.
//' \code{aevacc_by_trial_memnoise} Returns success counts for model simulation
//' @author Alexander Fengler, \email{alexanderfengler@@gmx.de}
//' @title Flexible evidence accumulation (item general) allowing for memory effects
//' @return Returns a numeric variable that provides a success count (runs that predicted a reaction time in the correct rt-bin and simultaneously the correct decision)
//' @param sd standard deviation used for drift diffusion process
//' @param theta theta used for drift diffusion process
//' @param drift drift-rate used for drift diffusion process
//' @param non_decision_time non decision time used for drift diffusion process
//' @param timestep timestep in ms associated with each step in the drift diffusion process
//' @param nr_reps number of repitions (simulation runs)
//' @param maxdur numeric variable that supplies the maximum reaction time considered a success in simulations
//' @param mindur numeric variable that supplies the minimum reaction time considered a succes in simulations
//' @param cur_decision numeric variable that provides the empirical decision taken in trial
//' @param update Vector that stores the item valuations for the trial conditon simulated
//' @param fixpos Vector that stores the locations for a supplied fixed fixation pathway
//' @param fixdur Vector that stores the fixation durations for a supplied fixed fixation pathway
//' @param items_seen_bias Numeric Variable storing the relative amount of drift that unseen items receive
//' @param items_seen_noise_bias Numeric Variable storing the relative noise sd that unseen items receive
//' @export
// [[Rcpp::export]]
int aevacc_by_trial_memnoise(float sd,
                             float theta,
                             float drift,
                             int non_decision_time,
                             float items_seen_bias,
                             float items_seen_noise_bias,
                             int maxdur,
                             int mindur,
                             int cur_decision,
                             NumericVector update,
                             IntegerVector fixpos,
                             IntegerVector fixdur,
                             int nr_reps,
                             int timestep){

  // Set seed for random sampler ------------------------------------------------------------------
  NumericVector seed(1);
  seed = floor(runif(1,-100000,100000));
  zigg.setSeed(seed[0]);
  // ----------------------------------------------------------------------------------------------

  // Variable collecting success counts -----------------------------------------------------------
  int out = 0;
  // ----------------------------------------------------------------------------------------------


  // Initialize Variables need for model propagation ----------------------------------------------
  int nr_items = update.size();
  NumericVector Evid(nr_items);
  int maxpos = 0;
  float temp = 0;
  bool decision = 1;
  int cur_rt = 0;
  int cur_fixpos_indice = 0;
  int num_fixpos = fixpos.size();
  int cur_fixpos = 0;
  int cur_fix_cnt = 0;
  NumericVector cur_update(nr_items);
  NumericVector temp_update(nr_items);
  IntegerVector items_seen(nr_items);
  IntegerVector items_seen_noise(nr_items);

  for (int i = 0; i < nr_items; ++i){
    update[i] = update[i]*drift;
    cur_update[i] = theta*update[i];
    items_seen_noise[i] = items_seen_noise_bias;
    items_seen[i] = items_seen_bias;
  }
  // ----------------------------------------------------------------------------------------------

  // cycle through simulations --------------------------------------------------------------------
  for (int rep_cnt = 0; rep_cnt < nr_reps;++rep_cnt){

    // reset variables before entering next simulation --------------------------------------------
    cur_rt = 0;
    cur_fix_cnt = 0;

    for (int i = 0; i < nr_items; ++i){
      Evid[i] = 0;
    }
    // --------------------------------------------------------------------------------------------

    // run simulation round -----------------------------------------------------------------------
    for (int fix_cnt = 0; fix_cnt < num_fixpos; ++fix_cnt){

      // identifying current fixation position
      cur_fixpos = fixpos[fix_cnt];

      cur_fixpos_indice = cur_fixpos - 1;

      // adjust cur_update to reflect current fixation bias (we take the multiplication times theta away for currently fixated item)
      cur_update[cur_fixpos_indice] = update[cur_fixpos_indice];

      // update items_seen vector to account for current fixation
      items_seen[cur_fixpos_indice] = 1;
      items_seen_noise[cur_fixpos_indice] = 1;

      // update temp_update vector to reflect current items seen
      for(int i = 0; i < nr_items; i++){
        temp_update[i] = cur_update[i]*items_seen[i];
      }

      for (int rt_cur_fix = 0; rt_cur_fix < fixdur[fix_cnt];rt_cur_fix += timestep){

        // update decision_term --> performance isues lead me to solve the decision_taken issue this way
        decision = 1;

        // accumulate Evidence --> In aDDM this happens through switch statement because fixation position matters

        for (int i = 0; i < nr_items; ++i){
          Evid[i] += temp_update[i] + sd*zigg.norm()*items_seen_noise[i];
        }

        //update RT
        cur_rt += timestep;

        // check whether decision taken -------------------------------------------------------------

        // get current max signal
        for (int j = 0; j < nr_items; ++j){
          temp = Evid[maxpos] - Evid[j];

          if (temp < 0){
            maxpos = j;
          }
        }

        // compute rdv
        for (int k = 0; k < nr_items; ++k){
          if (k != maxpos){
            temp = Evid[maxpos] - Evid[k];
            if (temp < 1){
              decision = 0;
              break;
            }
          }
        }
        // ------------------------------------------------------------------------------------------

        // if decision taken or maximum duration reacher, break out of loop
        if ((decision == 1) || (cur_rt > maxdur)){
          break;
        }
      }
      // --------------------------------------------------------------------------------------------

      // revert back cur_update to neutral state (times theta for update everywhere)
      cur_update[cur_fixpos_indice] = cur_update[cur_fixpos_indice]*theta;

      // count fixation upwards
      cur_fix_cnt +=  1;

      // if decision taken or maximum duration reached we break out of loop
      if ((decision == 1) || (cur_rt > maxdur)){
        break;
      }
    }
    // ----------------------------------------------------------------------------------------------

    // check whether outcome is success (right rt bin and decision)
    if (decision == 1){
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
//void myzsetseed(unsigned long int x) {
//    zigg.setSeed(x);
//return; }
