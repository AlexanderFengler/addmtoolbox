// [[Rcpp::depends(RcppZiggurat)]]

#include <Rcpp.h>
#include <Ziggurat.h>
#include <algorithm>
#include <stdlib.h>
using namespace Rcpp;
static Ziggurat::Ziggurat::Ziggurat zigg;

//' Simulate aDDM process by unique trial condition (>2 items)
//' @author Alexander Fengler, \email{alexanderfengler@@gmx.de}
//' @title Simulate aDDM process by unique trial condition (>2 items)
//' \code{aevacc_by_condition()}
//' @return vector that stores decisions and rts for each simulation run
//' @param sd standard deviation used for drift diffusion process
//' @param theta theta (attentional bias) used for drift diffusion process
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
IntegerVector aevacc_by_condition(float sd,
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

  // Initialize Variable that collects output -----------------------------------------------------
  IntegerVector out(2*nr_reps);
  // ----------------------------------------------------------------------------------------------

  // Initialize Variables needed to propagate model -----------------------------------------------
  int nr_items = update.size();
  NumericVector Evid(nr_items);
  int maxpos = 0;
  float temp = 0;
  bool decision = 1;
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

       for (int i = 0; i < nr_items; ++i){
         Evid[i] = 0;
       }

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

             decision = 1;

             // Accumulate Evidence for current timestep
             for (int i = 0; i < nr_items; ++i){
                         Evid[i] += cur_update[i] + sd*zigg.norm();
                         }

             //update RT
             cur_rt += timestep;

             // find signal with maximum value
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

             // revert back cur_update to neutral state (times theta for update everywhere)
             cur_update[cur_fixpos_indice] = cur_update[cur_fixpos_indice]*theta;

             cur_fix_cnt +=  1;

             if ((decision == 1) || (cur_rt >= maxdur)){
               break;
             }
         }
         // ------------------------------------------------------------------------------------------

       // store decision
       out[out_cnt] = maxpos+1; // maxpos plus one to convert from indice to item

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
