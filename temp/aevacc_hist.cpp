// SOME IMPORTANT THINGS TO NOTE:
// In the aDDM the maximally possible RT as outcome of the simulation is --> max.RT in the choice set plus the longest fixation in the choice set
// This is the result of two effects
// First: In this model we check whether current RTn is higher than maximal duration only after each fixation ends --> performance considerations
// Second  --> follows from First: Potentially this can lead to the combined effects of: The end of the second last fixation simulated is exactly one ms short of "maxdur" and
// therefore the model allows for one more fixation. The last fixation duration that is simulated is the maximum fixation duration of the current subject
// and moreover, the simulation is far from reaching the choice barrier --> given that the decision is not reached during the last fixation the maximum time that can be simulated
// will end up to be the Maximum RT of the current subject plus the Longest fixation of this subject


#include <Rcpp.h>
#include <Ziggurat.h>
using namespace Rcpp;
#include <algorithm>

#include <stdlib.h>

static Ziggurat::Ziggurat::Ziggurat zigg;

// [[Rcpp::depends(RcppZiggurat)]]
// [[Rcpp::export]]

IntegerVector aevacc_hist(int nr_reps,
                            int maxdur,
                            float cur_sd,
                            float theta,
                            NumericVector update,
                            IntegerVector fixpos,
                            IntegerVector fixdur,
                            IntegerVector fixdursamples,
                            int nr_items,
                            Function s){

  // Set seed for random sampler
   NumericVector seed(1);
   seed = floor(runif(1,-100000,100000));
   zigg.setSeed(seed[0]);

  // IntegerVector "out" collects all relevant information per repetition (Choice,RT,Nr. Fixations)
   IntegerVector out(2*nr_reps);

  //int RT = 0;
  // Evid will be updated per ms ... it collects the current evidence for each item
   NumericVector Evid(nr_items);

   int maxpos = 0;
   float temp = 0;
   bool decision = 1;
   int cur_rt = 0;
   int out_cnt = -2;
   int out_plus = 0;

   // Define the number of real fixations and respective durations for current condition
   int num_fixpos = fixpos.size();

   // moreover we have to add a variable storing the current fixation position
   IntegerVector cur_fixpos(1);
   IntegerVector cur_fixdur(1);
   IntegerVector temp_fixpos(1);

   // moreover we need a variable that counts the amount of fixations we use
   int cur_fix_cnt = 0;

   // we need one variable that basically is a converted fixation position (-1) to use the current position for accessing a vector (which starts at zero)
   int cur_fixpos_indice = 0;

   // Initialize vector that stores all fixation locations that can be sampled from later (see usage of vector below)
   IntegerVector eligible(nr_items);
   for (int i = 0; i < nr_items; i++){
     eligible[i] = i + 1;
   }

   // We generate the attention-dependent evidence updates here to calculate only once and afterwards access
    // At the start of every fixation we will update the respective position in the cur_update vector to account for the fixation bias
      NumericVector cur_update(nr_items);

      for (int i = 0; i < nr_items; ++i){
        cur_update[i] = theta*update[i];
      }

   // outer loop counts the simluations
     for (int rep_cnt = 0; rep_cnt < nr_reps;++rep_cnt){
         cur_rt = 0;
         cur_fix_cnt = 0;
         out_cnt += 2;
         out_plus = out_cnt + 1;

      for (int i = 0; i < nr_items; ++i){
          Evid[i] = 0;
         }

       // inner for loops run until decision taken
         // loop that counts until we reached last fixations
         for (int fix_cnt = 0; fix_cnt < 1000; ++fix_cnt){  // maybe start fix_cnt at one?

            // identifying current fixation position and duration
           // If empirical fixations still to be supplied fine otherwise sample from eligible fixations positions
           // and sample duration from fixdursamples vector
           if (fix_cnt <= (num_fixpos - 1)){

             cur_fixpos[0] = fixpos[fix_cnt];
             cur_fixdur[0] = fixdur[fix_cnt];
             // If we reached the last fixation, here we actually sample the duration from normal refixations!
             // Rationale: In case there is any "shorter-final-fixation" bias in the data, the algorithm should be given the chance to reproduce such an effect
             // NOT be supplied with such an effect !
             if (fix_cnt == num_fixpos - 1){
               cur_fixdur = s(fixdursamples,1);
             }
           }
           else if (fix_cnt > (num_fixpos - 1)){
             temp_fixpos[0] = cur_fixpos[0];
             // use sample function
             while (cur_fixpos[0] == temp_fixpos[0]){
               cur_fixpos = s(eligible,1);
             }
             cur_fixdur = s(fixdursamples,1);
           }


         cur_fixpos_indice = cur_fixpos[0] - 1;
         cur_update[cur_fixpos_indice] = update[cur_fixpos_indice];

           for (int rt_cur_fix = 0; rt_cur_fix < fixdur[fix_cnt];++rt_cur_fix){

             // update decision_term --> performance isues lead me to solve the decision_taken issue this way
             decision = 1;

             // accumulate Evidence --> In aDDM this happens through switch statement because fixation position matters

             for (int i = 0; i < nr_items; ++i){
                         Evid[i] += cur_update[i] + cur_sd*zigg.norm();
                         }

               //update RT
               cur_rt++;

            // now checking for whther decision has been taken

               for (int j = 0; j < nr_items; ++j){
                  temp = Evid[maxpos] - Evid[j];
                  if (temp < 0){
                    maxpos = j;
                  }
                }

            // computing RDV

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


       // it is max_pos + 1 because of difference in vector indexing in R and C++
       // stores decision
       out[out_cnt] = maxpos+1;

       //storing RT
       if(decision == 0){
        out[out_plus] = maxdur;
       }
       else {
        out[out_plus] = cur_rt;
       }

     }

   return out;
}

/// [[Rcpp::export]]
//void myzsetseed(unsigned long int x) {
  //  zigg.setSeed(x);
//return; }
