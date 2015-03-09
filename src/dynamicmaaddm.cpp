#include <Rcpp.h>
using namespace Rcpp;

//' Simulate aDDM process by unique trial condition (2 items)
//' @author Alexander Fengler, \email{alexanderfengler@@gmx.de}
//' @title Simulate aDDM process (by condition, 2 items)
//' \code{dynamicaddm()}
//' @return numeric variable storing likelihood value
//' @param sd standard deviation used for drift diffusion process
//' @param theta theta used for drift diffusion process [0,1]
//' @param gamma secondary attentional bias which refers to attributes no inspected [0,1]
//' @param drift drift-rate used for drift diffusion process
//' @param non_decision_time non decision time used for drift diffusion process
//' @param decision integer which gives choice in current trial (1 or 2)
//' @param rt reaction time of provided trial
//' @param valuations vector that stores the item valuations for the provided trial
//' @param nr_attributes integer which gives the number of attributes by item
//' @param fixpos vector with all empirical fixations positions encountered in provided trial
//' @param fixdur vector with all empirical fixation durations encounrtered in provided trial
//' @param stateStep numeric variable between [0,1] that indicates how finegrained the vertical grid of the model space shall be computed
//' @export
// [[Rcpp::export]]
double dynamicmaaddm(float sd = 0,
                     float theta = 0,
                     float gamma = 0,
                     float drift = 0,
                     int non_decision_time = 0,
                     int decision = 0,
                     NumericVector valuations = 0,
                     int nr_attributes = 0,
                     NumericVector fixpos = 0,
                     NumericVector fixdur = 0,
                     int rt = 0,
                     float stateStep = 0){

  // INITIALIZATIONS -----------------------------------------------------------------------------------

  // intialize loglik output variable
  double loglik = 0;

  // get number of fixations to consider
  int fixnum = fixpos.size();
  int nr_screen_positions = valuations.size();
  int nr_items = nr_screen_positions / nr_attributes;

  // Precompute drifts for each fixation position ----------------------------------------------------
  NumericMatrix valuationsmat(nr_attributes,2);
  NumericVector drifts(2*nr_attributes);

  // generate matrix that stores attribute wise values (rows) by item (columns)
  for(int i = 0; i < nr_items; i++){
    for (int j = 0; j < nr_attributes; j++){
      valuationsmat(j,i) = valuations[((nr_attributes*i)+j)];
    }
  }

  // use matrix to fill in drifts vector
  for (int i = 0; i < nr_items; i++){
    for (int j = 0; j < nr_attributes; j++){
      for (int k = 0; k <nr_attributes; k++){
        if (j == k){
          if (i == 0){
            drifts[((nr_attributes*i)+j)] += valuationsmat(j,i) - theta*valuationsmat(j,1);
          } else {
            drifts[((nr_attributes*i)+j)] += theta*valuationsmat(j,0) - valuationsmat(j,i);
          }
        } else {
          if (i == 0){
            drifts[((nr_attributes*i)+j)] += gamma*(valuationsmat(j,i) - theta*valuationsmat(j,1));
          } else {
            drifts[((nr_attributes*i)+j)] += gamma*(theta*valuationsmat(j,0) - valuationsmat(j,i));
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


  // define barriers
  int barrierUp = 1;
  int barrierDown = -1;

  // define barriers by timestep
  NumericVector barrierTimeUp(rt);
  NumericVector barrierTimeDown(rt);

  for (int i = 0; i < rt; i++){
    barrierTimeUp[i] = barrierUp;
    barrierTimeDown[i] = barrierDown;
  }

  // Initialize potential decay
  //float decay = 0;
  //for (int i = 1; i < rt; i++){
  //   barrierTimeUp[i] = barrierUp /(1+decay*(i));
  //  barrierTimeDown[i] = barrierDown /(1+decay*(i));
  //}

  // Define positions in vertical grid
  int numStates = ((barrierUp - barrierDown)/stateStep) + 1;
  NumericVector states(numStates);

  for (int i = 0; i < numStates; i++){
    states[i] = barrierUp - i*stateStep;
  }

  // Initialize probability states
  NumericVector prStates(numStates);

  // Initialize starting point with probability 1
  int zeroPos = ((numStates-1)/2); // note the '-1' is to account for indexing starting at 0
  prStates[zeroPos] = 1;

  // Initialize the vectors that collect the barrier crossing probabilities by timestep
  NumericVector upCrossing(rt);
  NumericVector downCrossing(rt);

  // Now main loop that propagates the model through the grid
  // Initialize a few variables that will be used in the loop
  NumericMatrix pCrossBarrierUp(numStates,nr_items*nr_attributes);
  NumericMatrix pCrossBarrierDown(numStates,nr_items*nr_attributes);

  for (int j = 0; j < drifts.size(); j++){
    for (int i = 0; i < numStates;i++){
      pCrossBarrierUp(i,j) = R::pnorm(stateStep*i,drifts[j],sd,0,0);
      pCrossBarrierDown(i,j) = R::pnorm((numStates-i-1)*stateStep,-drifts[j],sd,0,0);
    }
  }

  // fresh probability states vector
  NumericVector PrStatesNew(numStates);

  // storage variables used inside the for loop
  double tempUpCross = 0;
  double tempDownCross = 0;
  float sumIn = 0;
  float sumCurrent = 0;
  float ratioNewOld = 0;

  float temp = 0; // state variable used inside the for loop

  // we precompute pnorm probabilities used for barrier crossings
  NumericMatrix pChange(numStates*2,nr_items*nr_attributes);

  int curmax = numStates*2;
  for(int j = 0; j < drifts.size(); j++){
    for(int i = 0; i < curmax;i++){
      pChange(i,j) = R::dnorm(stateStep*(i-numStates),drifts[j],sd,0);
    }
  }
  //---------------------------------------------------------------------------------------------------------

  // Run the grid -------------------------------------------------------------------------------------------
  int t = -1;
  int cur_fixpos = 0;

  for (int fixcnt = 0; fixcnt < fixnum; fixcnt++){
    cur_fixpos = fixpos[fixcnt];
    for (int cur_fixdur = 0; cur_fixdur < fixdur[fixcnt]; cur_fixdur++){
      t += 1;
      PrStatesNew[0] = 0;
      PrStatesNew[numStates-1] = 0;
      for (int s = 0;s < numStates;s++){
        PrStatesNew[s] = 0;
        if (s > 0 & s < (numStates - 1)){
          // store new probabilities
          for(int i=0; i < numStates;i++){
            // temp is the current indice that we acces in the pChange vector which stores the transition probabilities
            // for all possible state transitions
            // example: current state ("to" in Antonios code) is 3 -- then we access the pChange vector in order 3,2,1,0,1,2,3,4,5....
            // P from three to three is a 0 vertical step changes, we start form the top vertical position and go down to bottom
            // abs() because we have to have the indices as psotive numbers
            temp = i-s + numStates;
            PrStatesNew[s] += stateStep*prStates[i]*pChange(temp,cur_fixpos-1);
          }
        }
      }

      // store the barrier crossings
      tempUpCross = 0;
      tempDownCross = 0;
      for(int i = 0; i < numStates;i++){
        tempUpCross += prStates[i]*pCrossBarrierUp(i,cur_fixpos-1);
        tempDownCross += prStates[i]*pCrossBarrierDown(i,cur_fixpos-1);
      }

      // renormalization
      sumIn = 0;
      sumCurrent = 0;

      for (int i=0;i < numStates;i++){
        sumIn+=prStates[i];
        sumCurrent+=PrStatesNew[i];
      }

      sumCurrent += tempUpCross + tempDownCross;
      ratioNewOld = sumIn/sumCurrent;

      tempUpCross = tempUpCross * ratioNewOld;
      tempDownCross = tempDownCross * ratioNewOld;

      for (int i=0;i<numStates;i++){
        //PrStatesNew[i] = PrStatesNew[i] * ratioNewOld;
        // make new states the old ones for next iteration
        prStates[i] = PrStatesNew[i] * ratioNewOld;
      }

      upCrossing[t] = tempUpCross;
      downCrossing[t] = tempDownCross;
    }
  }
  //---------------------------------------------------------------------------------------------------------------

  // Return correct likelihood depending on supplied decision -----------------------------------------------------
  if (decision == 1){
    return(loglik = upCrossing[t]);
  }  else {
    return(loglik = downCrossing[t]);
  }
  //---------------------------------------------------------------------------------------------------------------
}

//    return List::create(_["upCrossing"] = upCrossing,
//                        _["downCrossing"] = downCrossing,
//                        _["pCrossBarrierUp"] = pCrossBarrierUp,
//                        _["pCrossBarrierDown"] = pCrossBarrierDown);

