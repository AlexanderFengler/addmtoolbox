// Author: Alexander Fengler
// Date: Jan 25h 2015
// Purpose: Full output version of the addm evidence accumulation function


// This function now also allows for passing an R function that predicts next fixation locations and corresponding fixation durations
// If basically no fixation pathways are supplied the function will completely simulate the fixation path !

#include <Rcpp.h>
using namespace Rcpp;
#include <algorithm>
#include <stdlib.h>
#include <Ziggurat.h>

static Ziggurat::Ziggurat::Ziggurat zigg;

// [[Rcpp::depends(RcppZiggurat)]]
// [[Rcpp::export]]
NumericVector aevacc_full_output_withmem(int nr_reps,
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
   NumericVector out(22*nr_reps);

  //int RT = 0;
  // Evid will be updated per ms ... it collects the current evidence for each item
   NumericVector Evid(nr_items);

   int maxpos = 0;
   float temp = 0;
   bool decision = 1;
   int cur_rt = 0;
   int out_cnt = -22;
   IntegerVector temp_fixpos(1);
   IntegerVector cur_fixdur(1);


   // Other relevant output variables
   int item_last_attended = 0;
   float value_last_attended = 0;
   int last_attended_chosen = 0;

   NumericVector Durations(8);

         Durations[0] = 0;
         Durations[1] = 0;
         Durations[2] = 0;
         Durations[3] = 0;
         Durations[4] = 0;
         Durations[5] = 0;
         Durations[6] = 0;
         Durations[7] = 0;

   NumericVector Fixations(8);

         Fixations[0] = 0;
         Fixations[1] = 0;
         Fixations[2] = 0;
         Fixations[3] = 0;
         Fixations[4] = 0;
         Fixations[5] = 0;
         Fixations[6] = 0;
         Fixations[7] = 0;

   // Define the number of real fixations and respective durations for current condition
   int num_fixpos = fixpos.size();
   //int num_fixpos_sim = fixpos_sim.size()/nr_reps;

   // moreover we have to add a variable storing the current fixation position
   IntegerVector cur_fixpos(1);

   // moreover we need a variable that counts the amount of fixations we use
   int cur_fix_cnt = 0;

   // cur_fixpos_indice is a variable that tracks the current fixation position -- but deducts one to translate it into vector indice

   int cur_fixpos_indice;

  // Initialize vector that stores all fixation locations that can be sampled from later (see usage of vector below)
   IntegerVector eligible(nr_items);
   for (int i = 0; i < nr_items; i++){
     eligible[i] = i + 1;
   }

   // We generate the attention-dependent evidence updates here to calculate only once and afterwards access
    // At the start of every fixation we will update the respective position in the cur_update vector to account for the fixation bias
      NumericVector cur_update(nr_items);
      NumericVector temp_update(nr_items); // temp update is a new variable that is needed when we want to adjust the updates according to whether items have been seen at all. Will be used below.
      IntegerVector items_seen(nr_items); // vector tracking which items have been seen. will be used later

      for (int i = 0; i < nr_items; ++i){
        cur_update[i] = theta*update[i];
      }


   // outer loop counts the simluations
     for (int rep_cnt = 0; rep_cnt < nr_reps;++rep_cnt){
         cur_rt = 0;
         cur_fix_cnt = 0;
         out_cnt += 22;


         for (int i = 0; i < nr_items; ++i){
            Evid[i] = 0;
            Durations[i] = 0;
            Fixations[i] = 0;
         }

       // inner for loops run until decision taken

         // loop that counts until we reached last fixations
         for (int fix_cnt = 0; fix_cnt < 1000; ++fix_cnt){  // the < 1000 condition is just a placeholde, the maxdur condition is the real condition on breaking out of the for loop

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

           cur_fixpos_indice = (cur_fixpos[0] - 1);

           cur_update[cur_fixpos_indice] = update[cur_fixpos_indice];


           // update items_seen vector to account for current fixation
           items_seen[cur_fixpos_indice] = 1;

           // update temp_update vector to reflect current items seen
           for(int i = 0; i < nr_items; i++){
             temp_update[i] = cur_update[i]*items_seen[i];
           }

           // Update Fixations vector
           Fixations[(cur_fixpos[0]-1)]++;

           for (int rt_cur_fix = 0; rt_cur_fix < cur_fixdur[0];++rt_cur_fix){
             // update Duration on Item
             Durations[(cur_fixpos[0]-1)]++;

             // update decision_term --> performance isues lead me to solve the decision_taken issue this way
             decision = 1;

            // accumulate Evidence --> In aDDM this happens through switch statement because fixation position matters
            for (int i = 0; i < nr_items; ++i){
                Evid[i] += temp_update[i] + cur_sd*zigg.norm();
            }

            //update RT
            cur_rt++;

            // now checking for whether decision has been taken
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
               item_last_attended = cur_fixpos[0];
               value_last_attended = update[(cur_fixpos[0] - 1)];
               last_attended_chosen = 0;
               if (item_last_attended == (maxpos + 1)){
                 last_attended_chosen = 1;
               }
               break;
             }
         }

       // it is max_pos + 1 because of difference in vector indexing in R and C++
       // stores decision
       out[out_cnt] = maxpos+1;

       //storing Nr. Fixations
       out[out_cnt+1] = cur_fix_cnt;

       //storing RT
       if(decision == 0){
        out[out_cnt + 2] = maxdur;
       }
       else {
        out[out_cnt + 2] = cur_rt;
       }

       out[out_cnt + 3] = item_last_attended;
       out[out_cnt + 4] = value_last_attended;
       out[out_cnt + 5] = last_attended_chosen;

       out[out_cnt + 6] = Durations[0];
       out[out_cnt + 7] = Durations[1];
       out[out_cnt + 8] = Durations[2];
       out[out_cnt + 9] = Durations[3];
       out[out_cnt + 10] = Durations[4];
       out[out_cnt + 11] = Durations[5];
       out[out_cnt + 12] = Durations[6];
       out[out_cnt + 13] = Durations[7];


       out[out_cnt + 14] = Fixations[0];
       out[out_cnt + 15] = Fixations[1];
       out[out_cnt + 16] = Fixations[2];
       out[out_cnt + 17] = Fixations[3];
       out[out_cnt + 18] = Fixations[4];
       out[out_cnt + 19] = Fixations[5];
       out[out_cnt + 20] = Fixations[6];
       out[out_cnt + 21] = Fixations[7];
     }
   return out;
}

//// [[Rcpp::export]]
//void myzsetseed(unsigned long int x) {
//    zigg.setSeed(x);
//return; }
