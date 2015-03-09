// [[Rcpp::depends(RcppZiggurat)]]
#include <Rcpp.h>
#include <Ziggurat.h>
#include <algorithm>
#include <stdlib.h>
using namespace Rcpp;
static Ziggurat::Ziggurat::Ziggurat zigg;

//' Simulate aDDM process by unique trial (>2 items)
//' \code{aevacc_by_trial()}
//' @author Alexander Fengler, \email{alexanderfengler@@gmx.de}
//' @title Simulate aDDM process by unique trial (>2 items)
//' @return Returns a numeric variable that provides a success count (runs that predicted a reaction time in the correct rt-bin and simultaneously the correct decision)
//' @param sd standard deviation used for drift diffusion process
//' @param theta theta (attentional bias) used for drift diffusion process
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
//' @param gamma placeholder for interface consistency / see multiattribute versions for specification
//' @param nr_attributes placeholder for interface consistency / see multiattribute versions for specification
//' @export
// [[Rcpp::export]]
int aevacc_by_trial(float sd,
                    float theta,
                    float gamma,
                    float drift,
                    int non_decision_time,
                    int maxdur,
                    int mindur,
                    int cur_decision,
                    NumericVector update,
                    int nr_attributes,
                    IntegerVector fixpos,
                    IntegerVector fixdur,
                    int nr_reps,
                    int timestep){

  // Set seed for random sampler ------------------------------------------------------------------
  NumericVector seed(1);
  seed = floor(runif(1,-100000,100000));
  zigg.setSeed(seed[0]);
  // ----------------------------------------------------------------------------------------------


   // Variable collecting success counts ----------------------------------------------------------
   int out = 0;
   // ---------------------------------------------------------------------------------------------

   // Initialize Variables need for model propagation ------------------------------------------------
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
   for (int i = 0; i < nr_items; ++i){
     update[i] = update[i]*drift;
     cur_update[i] = theta*update[i];
   }
   // ---------------------------------------------------------------------------------------------

   // Cycle through simulations -------------------------------------------------------------------
   for (int rep_cnt = 0; rep_cnt < nr_reps;++rep_cnt){

     // reset variables before entering next simulation -------------------------------------------
     cur_rt = 0;
     cur_fix_cnt = 0;

     for (int i = 0; i < nr_items; ++i){
       Evid[i] = 0;
     }
     // -------------------------------------------------------------------------------------------

     // run a simulation round --------------------------------------------------------------------
     for (int fix_cnt = 0; fix_cnt < num_fixpos; ++fix_cnt){  // maybe start fix_cnt at one?

       // get current fixation position
       cur_fixpos = fixpos[fix_cnt];
       cur_fixpos_indice = cur_fixpos - 1;

       // adjust cur_update to reflect effect of current fixation
       cur_update[cur_fixpos_indice] = update[cur_fixpos_indice];

       // propagate model -------------------------------------------------------------------------
       for (int rt_cur_fix = 0; rt_cur_fix < fixdur[fix_cnt];rt_cur_fix += timestep){

          decision = 1;

         // accumulate evidence
         for (int i = 0; i < nr_items; ++i){
           Evid[i] += cur_update[i] + sd*zigg.norm();
         }

         //update rt
         cur_rt += timestep;

         // check whether decision has been taken ---------------------------------------------------

         // find current max signal
         for (int j = 0; j < nr_items; ++j){
           temp = Evid[maxpos] - Evid[j];

           if (temp < 0){
             maxpos = j;
           }
         }

         // compute current rdv
         for (int k = 0; k < nr_items; ++k){
           if (k != maxpos){
             temp = Evid[maxpos] - Evid[k];
             if (temp < 1){
               decision = 0;
               break;
             }
           }
         }
         // -------------------------------------------------------------------------------------------

         // if decision taken or maximum duration reached we break out of loop
         if ((decision == 1) || (cur_rt> maxdur)){
           break;
         }
       }

       // revert back cur_update to neutral state (times theta for update everywhere)
       cur_update[cur_fixpos_indice] = cur_update[cur_fixpos_indice]*theta;

       // count fixation upwards
       cur_fix_cnt +=  1;

       // if decision taken or maximum duration reached we break out of loop
       if ((decision == 1) || (cur_rt > maxdur)){
         break;
       }
     }
     // -------------------------------------------------------------------------------------------

     // determine whether outcome is success (correct rt and decision) ----------------------------
     if (decision == 1){
       if (maxpos == cur_decision){
         if ((cur_rt + non_decision_time <= maxdur) && (cur_rt + non_decision_time >= mindur)){
           ++out;
         }
       }
     }
     // ---------------------------------------------------------------------------------------------
   }
   // ---------------------------------------------------------------------------------------------
   return out;
}

//// [[Rcpp::export]]
//void myzfastsetseed(unsigned long int x){
//    zigg.setSeed(x);
//return;}
