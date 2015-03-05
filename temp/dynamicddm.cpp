#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List dynamicddm(){

  // initialize mean and sd
  float sd = 0.3;
  float mean = 0.03;

  // how coars do you want the steps?
  float stateStep = 0.05;

  // how much time do you want to give the process?
  int timeMax = 200;

  // define barriers
  int barrierUp = 1;
  int barrierDown = -1;

  // Barriers by timesteps
  NumericVector barrierTimeUp(timeMax);
  NumericVector barrierTimeDown(timeMax);

  for (int i = 0; i < timeMax; i++){
    barrierTimeUp[i] = barrierUp;
    barrierTimeDown[i] = barrierDown;
  }

  // Initialize potential decay
  //float decay = 0;

  //for (int i = 1; i < timeMax; i++){
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
  NumericVector upCrossing(timeMax);
  NumericVector downCrossing(timeMax);

  // Now main loop that propagates the model through the grid
  // initialize a few variables that will be used in the loop
  NumericVector pCrossBarrierUp(numStates);
  NumericVector pCrossBarrierDown(numStates);

  for (int i = 0; i < numStates;i++){
    pCrossBarrierUp[i] = R::pnorm(stateStep*i,mean,sd,0,0);
    pCrossBarrierDown[i] = R::pnorm((numStates-i-1)*stateStep,-mean,sd,0,0);
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
  float tempUpCross = 0;
  float tempDownCross = 0;
  float sumIn = 0;
  float sumCurrent = 0;
  float ratioNewOld = 0;

  // we precompute pnorm probabilities used for barrier crossings
  NumericVector pChange(numStates*2);

  int curmax = numStates*2;
  for(int i = 0; i < curmax;i++){
    pChange[i] = R::dnorm(stateStep*(i-numStates),mean,sd,0);
  }

  for (int t = 0; t < timeMax;t++){
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
          PrStatesNew[s] += stateStep*prStates[i]*pChange[temp];
        }
      }
    }

    // now we store the barrier crossings and do some normalizations
    tempUpCross = 0;
    tempDownCross = 0;
    for(int i = 0; i < numStates;i++){
      tempUpCross += prStates[i]*pCrossBarrierUp[i];
      tempDownCross += prStates[i]*pCrossBarrierDown[i];
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

  return List::create(_["upCrossing"] = upCrossing,
                      _["downCrossing"] = downCrossing,
                      _["pCrossBarrierUp"] = pCrossBarrierUp,
                      _["pCrossBarrierDown"] = pCrossBarrierDown);
}
