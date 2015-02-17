// aDDM algorithm for model fitting
// Output are simple "success counts" meaning whether the model spit out the right decision within the right RT-bin according to the original data as input
// No further, more detailed, output provided

// MODEL VARIANT: Theta 0 until item seen, then theta according to what is provided as input parameter

#include <Rcpp.h>
#include <Ziggurat.h>
using namespace Rcpp;
#include <algorithm>

#include <stdlib.h>

static Ziggurat::Ziggurat::Ziggurat zigg;

// [[Rcpp::depends(RcppZiggurat)]]
// [[Rcpp::export]]

int aevacc_fast_withmem_zeronoise(int nr_reps,
                            int maxdur,
                            int mindur,
                            int cur_decision,
                            float cur_sd, 
                            float theta,
                            NumericVector update, 
                            IntegerVector fixpos, 
                            IntegerVector fixdur,
                            int nr_items){
   
   // Set seed for random sampler
   NumericVector seed(1);
   seed = floor(runif(1,-100000,100000));
   zigg.setSeed(seed[0]);
   
   // IntegerVector "out" collects all relevant information per repetition (Choice,RT,Nr. Fixations)
   int out = 0;
  
   // Evid will be updated per timestep ... it collects the current evidence for each item
   NumericVector Evid(nr_items);
   
   int maxpos = 0;
   float temp = 0;
   bool decision = 1;
   int cur_rt = 0;
   int cur_fixpos_indice = 0;
   
   // Define the number of real fixations and respective durations for current condition
   int num_fixpos = fixpos.size();
   
   // moreover we have to add a variable storing the current fixation position
   int cur_fixpos = 0;
   
   // moreover we need a variable that counts the amount of fixations we use
   int cur_fix_cnt = 0;
   
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
         
         // reset RT, Fixation count and Evidence for new simulation round
         cur_rt = 0;
         cur_fix_cnt = 0;
         
         for (int i = 0; i < nr_items; ++i){
          Evid[i] = 0;
         }
         
       // inner for loops run until decision taken
      
         // count by fixation
         for (int fix_cnt = 0; fix_cnt < num_fixpos; ++fix_cnt){  // maybe start fix_cnt at one?
         
         // identifying current fixation position
         cur_fixpos = fixpos[fix_cnt];
         
         cur_fixpos_indice = cur_fixpos - 1;
             
         // adjust cur_update to reflect current fixation bias (we take the multiplication times theta away for currently fixated item)
         cur_update[cur_fixpos_indice] = update[cur_fixpos_indice];
         
         // update items_seen vector to account for current fixation
         items_seen[cur_fixpos_indice] = 1;
         
         // update temp_update vector to reflect current items seen
         for(int i = 0; i < nr_items; i++){
           temp_update[i] = cur_update[i]*items_seen[i];
         }
         
           for (int rt_cur_fix = 0; rt_cur_fix < fixdur[fix_cnt];++rt_cur_fix){
           
             // update decision_term --> performance isues lead me to solve the decision_taken issue this way
             decision = 1;
             
             // accumulate Evidence --> In aDDM this happens through switch statement because fixation position matters   
             
             for (int i = 0; i < nr_items; ++i){
                         Evid[i] += temp_update[i] + cur_sd*zigg.norm()*items_seen[i];
                         }
             
               //update RT
               cur_rt++;
               
            // now checking for whether decision has been taken
            
            // First figure out element that has max Evidence
               for (int j = 0; j < nr_items; ++j){
                  temp = Evid[maxpos] - Evid[j];
                  
                  if (temp < 0){
                    maxpos = j;
                  }
                }
                
            // Second compute current RDV
               for (int k = 0; k < nr_items; ++k){
                 if (k != maxpos){
                  temp = Evid[maxpos] - Evid[k];
                   if (temp < 1){
                     decision = 0;
                     break;
                   }
                 }
               }
                 // if decision taken or maximum duration reached we break out of loop
                 if ((decision == 1) || (cur_rt > maxdur)){
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
         
         
         // if decision is the same as the real one and falls within a specified RT.bin circumscribing the real RT we add one to success counts (out) 
        if (decision == 1){
         if (maxpos == cur_decision){
           if ((cur_rt <= maxdur) && (cur_rt >= mindur)){
              ++out;    
           } 
         } 
        }
     }
   
   return out;
}

//// [[Rcpp::export]]
//void myzsetseed(unsigned long int x) {
//    zigg.setSeed(x);
//return; }

