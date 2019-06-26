/* Zero Range process with periodic boundary conditions
1 dimensional system on the interval [0,1]
g(k)=k^power
Function recieves as input: pointer to the configuration of the system, dT and power
*/

#include <iostream>
#include <vector>
#include <random>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

using namespace std;

void ZRP_KMC_fixedBoundary(vector<int>* eta_adr, double dT)
{
  //Parameters of ZRP process
  double power=2.0;                        //Power pf the ZRP process: g(k)=k^power
 
  //Evolution of the ZRP via a KMC algorithm
  double t = 0.0;                      //Time variable
  int Nbins = (*eta_adr).size();       //Number of bins
  vector<double> rate_cumsum(Nbins);   //Cumulated jump rate
  //random_device rd;
  //srand(rd());
  random_device rd;
  mt19937_64 seed(rd());
  uniform_real_distribution<double> distribution(0.0,1.0);
  
  double rand_no_time =0;                   //Random number for time update
  double rand_no_rate =0;                   //Random number for finding bin
  double rand_no_direction =0;              //Random number for direction
  int index;                           //Index of the bin to carry out the jump
  int direction;                       //Direction of jump: 1(right),-1(left)
  int KMC_steps = 0;
  while(t <=dT)
  {
      //Compute cumulates sum for the rates
      rate_cumsum.at(0)=pow((double)(*eta_adr).at(0),power);
      for(int j=1; j<Nbins; j++){
          rate_cumsum.at(j) = rate_cumsum.at(j-1) + pow((double)(*eta_adr).at(j),power);
      }
      //Update time
      rand_no_time = distribution(seed);
      t -= log(rand_no_time)/rate_cumsum.at(Nbins-1);
      if(t <= dT)
      {
          KMC_steps++;
	  //Find index of the bin to carry out the jump
          rand_no_rate = distribution(seed);
          index=0;
          while(rate_cumsum.at(index) < rate_cumsum.at(Nbins-1)*rand_no_rate){
              index++;
          }
	  //Find direction to move
          rand_no_direction = distribution(seed);
          if(rand_no_direction < 0.5){
              direction = -1;
          }
          else{
              direction = 1;
          }
//	  //Update configuration with periodic boundary conditions
//          if(index==0 && direction==-1){
//              (*eta_adr).at(0)--;
//              (*eta_adr).at(Nbins-1)++;
//          }
//          else if(index==Nbins-1 && direction==1){
//              (*eta_adr).at(Nbins-1)--;
//              (*eta_adr).at(0)++;
//          }
//          else{
//              (*eta_adr).at(index)--;
//              (*eta_adr).at(index+direction)++;
//          }
      //Update configuration with fixed boundary conditions
          if(index == 1 && direction == -1){
              (*eta_adr).at(1)--;
          }
          else if(index == Nbins-2 && direction == 1){
              (*eta_adr).at(Nbins-2)--;
          }
          else if(index == 0 && direction == -1){ //start again
              t += log(rand_no_time)/rate_cumsum.at(Nbins-1);
          }
          else if(index == 0 && direction == 1){
              (*eta_adr).at(index+direction)++;
          }
          else if(index == Nbins-1 && direction == 1){ //start again
              t += log(rand_no_time)/rate_cumsum.at(Nbins-1);
          }
          else if(index == Nbins-1 && direction == -1){
              (*eta_adr).at(index+direction)++;
          }
          else{
              (*eta_adr).at(index)--;
              (*eta_adr).at(index+direction)++;
          }
      }
  }

  cout << "KMC_steps = " << KMC_steps << endl;
  return;
}
