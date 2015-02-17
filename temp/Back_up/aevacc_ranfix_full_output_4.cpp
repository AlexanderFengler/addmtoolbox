// SOME IMPORTANT THINGS TO NOTE:
  // In the aDDM the maximally possible RT as outcome of the simulation is --> max.RT in the choice set plus the longest fixation in the choice set
// This is the result of two effects 
// First: In this model we check whether current RTn is higher than maximal duration only after each fixation ends --> performance considerations
// Second  --> follows from First: Potentially this can lead to the combined effects of: The end of the second last fixation simulated is exactly one ms short of "maxdur" and 
// therefore the model allows for one more fixation. The last fixation duration that is simulated is the maximum fixation duration of the current subject 
// and moreover, the simulation is far from reaching the choice barrier --> given that the decision is not reached during the last fixation the maximum time that can be simulated
// will end up to be the Maximum RT of the current subject plus the Longest fixation of this subject


#include <Rcpp.h>
using namespace Rcpp;
#include <algorithm>

#include <stdlib.h>

// [[Rcpp::export]]
NumericVector aevacc_ranfix_full_output_4(int nr_reps,
                                   int maxdur, 
                                   float cur_sd, 
                                   float theta,
                                   NumericVector update, 
                                   NumericVector fixpos_sim, 
                                   NumericVector fixdur_sim){
  
  
  
  // IntegerVector "out" collects all relevant information per repetition (Choice,RT,Nr. Fixations)
  NumericVector out(22*nr_reps);
  // Noise collects the "noise-terms" per ms
  NumericVector noise(4); //= rnorm(4,0,cur_sd);
  
  //int RT = 0;
  // Evid will be updated per ms ... it collects the current evidence for each item
  NumericVector Evid(4);
  
  int maxpos = 0;
  int vec_len = 4;
  float temp = 0;
  bool decision = 1;
  int cur_rt = 0;
  int out_cnt = -22;
  
  
  // Other relevant output variables
  
  int item_last_attended = 0;
  float value_last_attended = 0;
  int last_attended_chosen = 0;
  
  NumericVector Durations(8);
  
  Durations[0] = 0;
  Durations[1] = 0;
  Durations[2] = 0;
  Durations[3] = 0;
  Durations[4] = -1;
  Durations[5] = -1;
  Durations[6] = -1;
  Durations[7] = -1;
  
  NumericVector Fixations(8);
  
  Fixations[0] = 0;
  Fixations[1] = 0;
  Fixations[2] = 0;
  Fixations[3] = 0;
  Fixations[4] = -1;
  Fixations[5] = -1;
  Fixations[6] = -1;
  Fixations[7] = -1;
  
  
  int num_fixpos_sim = fixpos_sim.size()/nr_reps;

  // moreover we have to add a variable storing the current fixation position
  
  int cur_fixpos = 0;
  
  // moreover we need a variable that counts the amount of fixations we use
  
  int cur_fix_cnt = 0;
  
  // We generate the attention-dependent evidence updates here to calculate only once and afterwards access
  
  NumericVector update_1(4);
  NumericVector update_2(4);
  NumericVector update_3(4);
  NumericVector update_4(4);
  
  // Fixating 1
  update_1[0] = update[0];
  update_1[1] = theta*update[1];
  update_1[2] = theta*update[2];
  update_1[3] = theta*update[3];
  
  // Fixating 2
  update_2[0] = theta*update[0];
  update_2[1] = update[1];
  update_2[2] = theta*update[2];
  update_2[3] = theta*update[3];
  
  // Fixating 3
  update_3[0] = theta*update[0];
  update_3[1] = theta*update[1];
  update_3[2] = update[2];
  update_3[3] = theta*update[3];
  
  // Fixating 4
  update_4[0] = theta*update[0];
  update_4[1] = theta*update[1];
  update_4[2] = theta*update[2];
  update_4[3] = update[3];
  
  // The following numeric_vector "cur_update" will be filled for at the start of every fixation 
  NumericVector cur_update(4);
  
  
  // lastly we need counts for the simulated fixations and simulated durations
  
  int total_fixpos_sim_cnt = 0;
  
  int total_fixdur_sim_cnt = 0;
  
  // outer loop counts the simluations 
  for (int rep_cnt = 0; rep_cnt < nr_reps;++rep_cnt){
   
    cur_rt = 0;
    cur_fix_cnt = 0;
    out_cnt += 22;
    cur_fixpos = 0;
    
    Evid[0] = 0;
    Evid[1] = 0;
    Evid[2] = 0;
    Evid[3] = 0;
    
    Durations[0] = 0;
    Durations[1] = 0;
    Durations[2] = 0;
    Durations[3] = 0;
    
    
    Fixations[0] = 0;
    Fixations[1] = 0;
    Fixations[2] = 0;
    Fixations[3] = 0;
    
    // inner for loops run until decision taken
    
    
    for (int fix_cnt_sim = 0; fix_cnt_sim < num_fixpos_sim; ++fix_cnt_sim){  // maybe start fix_cnt at one?
                                                                                 decision = 0;
                                                                                 if (fixpos_sim[total_fixpos_sim_cnt] != cur_fixpos){
                                                                                 // identifying current fixation position
                                                                                 cur_fixpos = fixpos_sim[total_fixpos_sim_cnt];
                                                                                 
                                                                                 // Update Fixations vector
                                                                                 Fixations[(cur_fixpos-1)]++;
                                                                                 
                                                                                 switch(cur_fixpos){
                                                                                   
                                                                                   case 1:   cur_update = update_1; 
                                                                                   break;
                                                                                   
                                                                                   case 2:   cur_update = update_2; 
                                                                                   break;
                                                                                   
                                                                                   case 3:   cur_update = update_3; 
                                                                                   break;
                                                                                   
                                                                                   case 4:   cur_update = update_4; 
                                                                                   break;         
                                                                                 }
                                                                                 
                                                                                 for (int rt_cur_fix = 0; rt_cur_fix < fixdur_sim[total_fixdur_sim_cnt];++rt_cur_fix){
                                                                                   
                                                                                   // update Duration on Item
                                                                                   Durations[(cur_fixpos-1)]++;
                                                                                   
                                                                                   // update decision_term --> performance isues lead me to solve the decision_taken issue this way
                                                                                   decision = 1;
                                                                                   
                                                                                   //load noise-terms
                                                                                   noise = rnorm(4,0,cur_sd);
                                                                                   
                                                                                   // accumulate Evidence --> In aDDM this happens through switch statement because fixation position matters   
                                                                                   
                                                                                   for (int i = 0; i < vec_len; ++i){
                                                                                     Evid[i] += cur_update[i] + noise[i];
                                                                                   }
                                                                                   //update RT
                                                                                   cur_rt += 1;
                                                                                   
                                                                                   // now checking for whther decision has been taken
                                                                                   
                                                                                   for (int j = 0; j < vec_len; ++j){
                                                                                     temp = Evid[maxpos] - Evid[j];
                                                                                     
                                                                                     if (temp < 0){
                                                                                       maxpos = j;
                                                                                     }
                                                                                     
                                                                                    }
                                                                                   
                                                                                   // alternatively the above code-chunk (for loop) can be avoided by using which_max() ... TEST
                                                                                   
                                                                                   // computing RDV
                                                                                   
                                                                                   for (int k=0; k < vec_len; ++k){
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
                                                                                 
                                                                                 total_fixdur_sim_cnt += 1;
                                                                                 cur_fix_cnt += 1;
                                                                             }
                                                                              
                                                                               if ((decision == 1) || (cur_rt >= maxdur)){
                                                                                 item_last_attended = cur_fixpos;
                                                                                 value_last_attended = cur_update[cur_fixpos - 1];
                                                                                   last_attended_chosen = 0;
                                                                                 if (item_last_attended == (maxpos + 1)){
                                                                                   last_attended_chosen = 1;
                                                                                 }
                                                                                 break;
                                                                               }
                                                                               
                                                                                total_fixpos_sim_cnt++; 
                                                                                
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