#include <Rcpp.h>
using namespace Rcpp;

//' Simulate aDDM process by unique trial condition (2 items)
//' @author Alexander Fengler, \email{alexanderfengler@@gmx.de}
//' @title Simulate aDDM process (by condition, 2 items)
//' \code{dynamicaddm()}
//' @return numeric variable storing likelihood value
//' @param sd standard deviation used for drift diffusion process
//' @param theta theta used for drift diffusion process
//' @param drift drift-rate used for drift diffusion process
//' @param non_decision_time non decision time used for drift diffusion process
//' @param rt reaction time of provided trial
//' @param valuations vector that stores the item valuations for the provided trial
//' @param fixpos vector with all empirical fixations positions encountered in provided trial
//' @param fixdur vector with all empirical fixation durations encounrtered in provided trial
//' @param stateStep numeric variable between [0,1] that indicates how finegrained the vertical grid of the model space shall be computed
//' @export
// [[Rcpp::export]]
int dynamicaddm(float sd,
                float theta,
                float drift,
                int non_decision_time,
                NumericVector valuations,
                NumericVector fixpos,
                NumericVector fixdur,
                int rt,
                float stateStep){

  // get number of fixations to consider
  int fixnum = fixpos.size();

  // generate the two drifts corresponding to potential fixation locations
  NumericVector drifts(2);
  drifts[0] = drift*(valuations[0] - theta*valuations[1]);
  drifts[1] = drift*(theta*valuations[0] - valuations[1]);

  // define barriers
  int barrierUp = 1;
  int barrierDown = -1;

  // Barriers by timesteps
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

  // Define grid vertically
  int numStates = ((barrierUp - barrierDown)/stateStep) + 1;
  NumericVector states(numStates);

  for (int i = 0; i < numStates; i++){
    states[i] = barrierUp - i*stateStep;
  }

  // Initialize probability states (state in the middle zero for the moment)
  NumericVector prStates(numStates);
  int zeroPos = ((numStates-1)/2); // assuming odd number of states this provides the indice of the middle state (-1 in the end to account for vectors-indices starting at 0 in c++)
  prStates[zeroPos] = 1;

  // Now initialize the vectors that collect the barrier crossing probabilities by timestep
  NumericVector upCrossing(rt);
  NumericVector downCrossing(rt);

  // Now main loop that propagates the model through the grid
  // initialize a few variables that will be used in the loop
  NumericMatrix pCrossBarrierUp(numStates,2);
  NumericMatrix pCrossBarrierDown(numStates,2);

  for (int i = 0; i < numStates;i++){
    pCrossBarrierUp(i,0) = R::pnorm(stateStep*i,drifts[0],sd,0,0);
    pCrossBarrierDown(i,0) = R::pnorm((numStates-i-1)*stateStep,-drifts[0],sd,0,0);
    pCrossBarrierUp(i,1) = R::pnorm(stateStep*i,drifts[1],sd,0,0);
    pCrossBarrierDown(i,1) = R::pnorm((numStates-i-1)*stateStep,-drifts[1],sd,0,0);
  }

  // fresh probability states vector
  NumericVector PrStatesNew(numStates);
  // a vector that stores current states we test for in inner for loop
  //float to = 0;
  // vector that stores current distance from all all states to the "to" states
  NumericVector change(numStates);
  // a temp variable may or may not be used for multiple purposes
  float temp = 0;

  // storage variables used inside the for loop
  //float changeUp = 0;
  //float changeDown = 0;
  double tempUpCross = 0;
  double tempDownCross = 0;
  float sumIn = 0;
  float sumCurrent = 0;
  float ratioNewOld = 0;

  // we precompute pnorm probabilities used for barrier crossings
  NumericMatrix pChange(numStates*2,2);

  int curmax = numStates*2;
  for(int i = 0; i < curmax;i++){
    pChange(i,0) = R::dnorm(stateStep*(i-numStates),drifts[0],sd,0);
    pChange(i,1) = R::dnorm(stateStep*(i-numStates),drifts[1],sd,0);
  }

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

      // now we store the barrier crossings and do some normalizations
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
        PrStatesNew[i] = PrStatesNew[i] * ratioNewOld;
        // make new states the old ones for next iteration
        prStates[i] = PrStatesNew[i];
      }

      upCrossing[t] = tempUpCross;
      downCrossing[t] = tempDownCross;
    }
  }

//   return List::create(_["upCrossing"] = upCrossing,
//                       _["downCrossing"] = downCrossing,
//                       _["pCrossBarrierUp"] = pCrossBarrierUp,
//                       _["pCrossBarrierDown"] = pCrossBarrierDown);
return(1);
}
