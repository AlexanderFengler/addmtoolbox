// [[Rcpp::depends(RcppZiggurat)]]
#include <Rcpp.h>
#include <algorithm>
#include <stdlib.h>
#include <Ziggurat.h>
using namespace Rcpp;
static Ziggurat::Ziggurat::Ziggurat zigg;

//' Runs evidence accumulation function (item general case) for one trial, and returns detailed model output
//' @author Alexander Fengler, \email{alexanderfengler@@gmx.de}
//' @title Evidence accumulation by condition (item general) with detailed output
//' \code{aevacc_by_condition}
//' @return Returns a vector that stores detailed output by simulation run
//' @param sd standard deviation used for drift diffusion process
//' @param theta theta used for drift diffusion process
//' @param drift drift-rate used for drift diffusion process
//' @param non_decision_time non decision time used for drift diffusion process
//' @param timestep timestep in ms associated with each step in the drift diffusion process
//' @param nr_reps number of repitions (simulation runs)
//' @param maxdur maximum duration in ms that the process is allowed to simulate
//' @param update Vector that stores the item valuations for the trial conditon simulated
//' @param fixpos Vector that stores the locations for a supplied fixed fixation pathway
//' @param fixdur Vector that stores the fixation durations for a supplied fixed fixation pathway
//' @param fixdursamples Vector from which fixation duration can be sampled once supplied fixations run out
//' @param fixation_model a user supplied fixation model that will be utilized to supply fixation locations and potentially fixation durations
//' @export
// [[Rcpp::export]]
NumericVector aevacc_full_output_memnoise(float sd,
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
  // ----------------------------------------------------------------------------------------------

  // Initialization of Varibales needed in loop ---------------------------------------------------
  NumericVector Evid(nr_items);
  int out_cnt = -7+2*nr_items;           // output counter
  float temp = 0;                        // loop internal (used at point of current-rdv calculation)
  int num_fixpos = fixpos.size();        // number of real fixations supplied
  IntegerVector temp_fixpos(1);
  IntegerVector cur_fixpos(1);           // storing current fixation position
  int cur_fixpos_indice;                 // storing indice of current fixation (-1 to make it vector indice)
  IntegerVector cur_fixdur(1);           // storing current fixation duration
  int cur_fix_cnt = 0;                   // counts number of fixations used
  IntegerVector eligible(nr_items);      // binary indication of which fixation position is valid for sampling

  for (int i = 0; i < nr_items; i++){
    eligible[i] = i + 1;
  }

  NumericVector cur_update(nr_items);    // storing current updates
  NumericVector temp_update(nr_items);   // storing current updates adjusted for items_seen and items_seen noise
  NumericVector items_seen(nr_items);    // vector that scales drift for unseen items
  NumericVector items_seen_noise(nr_items); // vector than scales noise for unseen items
  float items_seen_noise_scalar = 0;


  for (int i = 0; i < nr_items; ++i){
    update[i] = update[i]*drift;
    cur_update[i] = theta*update[i];
    items_seen_noise[i] = items_seen_noise_scalar;
  }
  // ----------------------------------------------------------------------------------------------

  // Loop cycling through simulations -------------------------------------------------------------
  for (int rep_cnt = 0; rep_cnt < nr_reps;++rep_cnt){
    // reset ----------------------------------------------------------------------------------
    cur_rt = 0;
    cur_fix_cnt = 0;
    out_cnt += 7 + 2*nr_items;

    for (int i = 0; i < nr_items; ++i){
      Evid[i] = 0;
      Durations[i] = 0;
      Fixations[i] = 0;
    }
    // -----------------------------------------------------------------------------------------

    // Enter simulation run --------------------------------------------------------------------
    for (int fix_cnt = 0; fix_cnt < 1000; ++fix_cnt){  // the < 1000 condition is arbitrary, "maxdur" is used to break out of loop

      // Get next fixation positions -----------------------------------------------------------
      // identifying current fixation position and duration
      // If empirical fixations still to be supplied, fine otherwise sample from eligible fixation positions
      // and sample duration from fixdursamples vector
      if (fix_cnt <= (num_fixpos - 1)){
        cur_fixpos[0] = fixpos[fix_cnt];
        cur_fixdur[0] = fixdur[fix_cnt];
        // If we reached the last fixation, here we actually sample the duration from normal refixations!
        // Rationale: In case there is any "shorter-final-fixation" bias in the data, the algorithm should be given the chance to reproduce such an effect
        // NOT be supplied with such an effect !
        if (fix_cnt == num_fixpos - 1){
          cur_fixdur = fixation_model(fixdursamples,1);
        }
      } else if (fix_cnt > (num_fixpos - 1)){
        temp_fixpos[0] = cur_fixpos[0];
        // use sample function
        while (cur_fixpos[0] == temp_fixpos[0]){
          cur_fixpos = fixation_model(eligible,1);
        }
        cur_fixdur = fixation_model(fixdursamples,1);
      }
      // -----------------------------------------------------------------------------------------

      // Make updates according to new fixation location -----------------------------------------

      // adjust cur_update vector to reflect current fixation position
      cur_fixpos_indice = (cur_fixpos[0] - 1);
      cur_update[cur_fixpos_indice] = update[cur_fixpos_indice];


      // update items_seen vector to account for current fixation
      items_seen[cur_fixpos_indice] = 1;
      items_seen_noise[cur_fixpos_indice] = 1;

      // update temp_update vector to reflect current items seen
      for(int i = 0; i < nr_items; i++){
        temp_update[i] = cur_update[i]*items_seen[i];
      }

      // Update Fixations vector
      Fixations[(cur_fixpos[0]-1)]++;
      // -----------------------------------------------------------------------------------------


      // Propagate Model -------------------------------------------------------------------------
      for (int rt_cur_fix = 0; rt_cur_fix < cur_fixdur[0];++rt_cur_fix){
        // update Duration on Item
        Durations[(cur_fixpos[0]-1)]++;

        // update decision_term --> performance isues lead me to solve the decision_taken issue this way
        decision = 1;

        // accumulate Evidence for timestep
        for (int i = 0; i < nr_items; ++i){
          Evid[i] += temp_update[i] + sd*zigg.norm()*items_seen_noise[i];
        }

        //update RT
        cur_rt++;

        // Find signal with maximum value
        for (int j = 0; j < nr_items; ++j){
          temp = Evid[maxpos] - Evid[j];

          if (temp < 0){
            maxpos = j;
          }
        }

        // compute RDV
        for (int k=0; k < nr_items; ++k){
          if (k != maxpos)
            temp = Evid[maxpos] - Evid[k];
          if (temp < 1){
            decision = 0;
            break;
          }
        }

        if (decision == 1){
          break;
        }
      }
      // ---------------------------------------------------------------------------------------

      // revert back cur_update to neutral state (times theta for update everywhere)
      cur_update[cur_fixpos_indice] = cur_update[cur_fixpos_indice]*theta;

      // update total fixation count
      cur_fix_cnt +=  1;

      // Check whether decision has been taken or maximum duration has been reached at this point
      if ((decision == 1) || (cur_rt >= maxdur)){
        item_last_attended = cur_fixpos[0];
        value_last_attended = update[(cur_fixpos[0] - 1)];
        last_attended_chosen = 0;
        if (item_last_attended == (maxpos + 1)){
          last_attended_chosen = 1;
        }
        break;
      }
    }
    // ------------------------------------------------------------------------------------------

    // Storing output -----------------------------------------------------------------------------
    out[out_cnt] = maxpos + 1;
    out[out_cnt+1] = cur_fix_cnt;

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
      out[out_cnt + 6 + i] = Durations[i]*timestep;
      out[out_cnt + 6 + nr_items + i] = Fixations[i];
    }
    // --------------------------------------------------------------------------------------------
  }
  // ----------------------------------------------------------------------------------------------
  return out;
}

//// [[Rcpp::export]]
//void myzsetseed(unsigned long int x) {
//    zigg.setSeed(x);
//return; }
