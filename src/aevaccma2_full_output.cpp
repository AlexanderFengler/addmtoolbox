// [[Rcpp::depends(RcppZiggurat)]]
#include <Rcpp.h>
#include <algorithm>
#include <stdlib.h>
#include <Ziggurat.h>
using namespace Rcpp;
static Ziggurat::Ziggurat::Ziggurat zigg;

//' Simulate aDDM process (2 items, multiattribute) with detailed output
//' @author Alexander Fengler, \email{alexanderfengler@@gmx.de}
//' @title Simulate aDDM process (2 items) with detailed output
//' \code{aevaccma2_full_output()}
//' @return vector that stores detailed output by simulation run
//' @param sd standard deviation used for drift diffusion process
//' @param theta theta (attentional bias) used for drift diffusion process
//' @param drift drift-rate used for drift diffusion process
//' @param non_decision_time non decision time used for drift diffusion process
//' @param timestep timestep in ms associated with each step in the drift diffusion process
//' @param nr_reps number of repitions (simulation runs)
//' @param maxdur maximum duration in ms that the process is allowed to simulate
//' @param update Vector that stores the item valuations for the trial conditon simulated
//' @param fixation_model a user supplied fixation model that will be utilized to supply fixation locations and potentially fixation durations
//' @param gamma placeholder for interface consistency / see multiattribute versions for specification
//' @param nr_attributes placeholder for interface consistency / see multiattribute versions for specification
//' @export
// [[Rcpp::export]]
NumericVector aevaccma2_full_output(float sd,
                                    float theta,
                                    float gamma,
                                    float drift,
                                    int non_decision_time,
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

  // Output Collection Variables ------------------------------------------------------------------
  int nr_screen_positions = update.size();
  int nr_items = nr_screen_positions / nr_attributes;
  NumericVector out((6+2*nr_screen_positions)*nr_reps);
  int maxpos = 0;
  bool decision = 1;
  int cur_rt = 0;
  int item_last_attended = 0;
  float value_last_attended = 0;
  int last_attended_chosen = 0;
  NumericVector Durations(nr_screen_positions);
  NumericVector Fixations(nr_screen_positions);

  NumericMatrix fixdat = fixation_model(update, nr_screen_positions);
  // ----------------------------------------------------------------------------------------------

  // Initialization of Varibales needed in loop ---------------------------------------------------
  float rdv = 0;
  int out_cnt = -6 - nr_screen_positions - nr_screen_positions;   // output counter
  int cur_fixpos_indice = 0;             // storing indice of current fixation (-1 to make it vector indice)
  int cur_fix_cnt = 0;                   // counts number of fixations used

  NumericVector cur_update(nr_screen_positions);    // storing current updates
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

  // Loop cycling through simulations ---------------------------------------------------------------
  for (int rep_cnt = 0; rep_cnt < nr_reps;++rep_cnt){
    // reset ----------------------------------------------------------------------------------------
    cur_rt = 0;
    cur_fix_cnt = 0;
    out_cnt += 6 + 2*nr_screen_positions;
    rdv = 0;
    decision = 0;

    for (int i = 0; i < nr_screen_positions; ++i){
      Durations[i] = 0;
      Fixations[i] = 0;
    }
    // -----------------------------------------------------------------------------------------------

    // Generate fixations for current simulation run -------------------------------------------------
    fixdat = fixation_model(update, nr_screen_positions);
    // -----------------------------------------------------------------------------------------------

    // Enter simulation run --------------------------------------------------------------------------
    for (int fix_cnt = 0; fix_cnt < 1000; ++fix_cnt){  // the < 1000 condition is arbitrary, "maxdur" is used to break out of loop
      // Make updates according to new fixation location ---------------------------------------------

      // adjust cur_update vector to reflect current fixation position
      cur_fixpos_indice = (fixdat(fix_cnt,0) - 1);

      // Update Fixations vector
      Fixations[fixdat(fix_cnt,0)-1]++;
      // -----------------------------------------------------------------------------------------------

      // Propagate Model -------------------------------------------------------------------------------
      for (int rt_cur_fix = 0; rt_cur_fix < fixdat(fix_cnt,1);rt_cur_fix += timestep){
        // update Duration on Item
        Durations[fixdat(fix_cnt,0)-1] += timestep;

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
      // -------------------------------------------------------------------------------------------------

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
    out[out_cnt + 4] = value_last_attended;
    out[out_cnt + 5] = last_attended_chosen;

    for(int i = 0; i < nr_screen_positions;i++){
      out[out_cnt + 6 + i] = Durations[i];
      out[out_cnt + 6 + nr_screen_positions + i] = Fixations[i];
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
