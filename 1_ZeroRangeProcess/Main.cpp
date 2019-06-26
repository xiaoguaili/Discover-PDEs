#include <iostream>
#include <fstream>
#include <ctime>
#include <string>
#include <stdio.h>
#include <vector>
#include <math.h>
#include <cmath>

using namespace std;

#include "InitialProfile.h"
#include "ZRP_KMC_fixedBoundary.h"
#include "Postprocessing.h"
#include "SavingParameters.h"

int main()
{
  //Configuration parameters and variables
  static int Nbin = 5000;                //Number of bins
  static int Npart = 20000;              //Number of particles
  vector<int> rho(Nbin);                 //Configuration at current time
  string profile = "linear";             //Initial profile
  string process = "ZRP";                //Zero range process
  static double slope = 0.0;              //Slope of linear profile
  static double center = 4.00;           //Value at center of linear profile
  double x[Nbin];                        //Coordinate at the center of each bin assuming total length 1
  double dx = 1.0/(double)Nbin;
  x[0] = 0.0;
  for(int j=1; j<Nbin; j++)
      x[j] = x[j-1]+dx;
    
  //Time parameters
  static double t_equilibration = 100/pow(5000.0,2.0);    //Time for equilibration
  static double t_randomize = 0.1/pow(5000.0,2.0);        //Time to create different initial conditions, after t_equilibration
  static double h = 0.001/pow(5000.0,2.0);               //Time to calculations actually used
  static double t_factor = pow((double)Nbin,2.0);               //t_micro = t_macro*t_factor
  
  //Sampling parameters
  static int R2 = 200;                  //Number of realizations; R2 in paper
  static int R1 = 5;                   //R1 in paper
    
  //Filter function parameters: two groups of basis functions are used
  //First group of basis functions
  static int Ngamma1 = 20;               //Number of basis functions
  double x_basis1[Ngamma1];              //Concentration points of the basis functions
  for(int i=0; i<Ngamma1; i++){
      x_basis1[i] = i*1.0/(double)Ngamma1;
  }
  //Second group of basis functions
  static int Ngamma2 = 40;               //Number of basis functions
  double x_basis2[Ngamma2];              //Concentration points of the basis functions
  for(int i=0; i<Ngamma2; i++){
      x_basis2[i] = i*1.0/(double)Ngamma2;
  }

  //Output variables
  ofstream myfile;
  char filename[250];
  string full_filename1;
  string full_filename2;
  //If nessasary, change the size of following vectors
  static double rho_gamma_initial1[100000];
  static double rho_gamma_final1[100000];
  static double rho_gamma_initial2[100000];
  static double rho_gamma_final2[100000];
  string file_path1="../Ngamma_20/";     //Output of different groups of basis functions are saved in different folders.
  string file_path2="../Ngamma_40/";
    
  //Timer
  clock_t t=clock();
  
  //Saving parameters on summary file
  SavingParameters(file_path1, Nbin, Npart, profile, slope, center, process, t_equilibration, t_randomize, h, t_factor, R2, Ngamma1, x_basis1, x);
  SavingParameters(file_path2, Nbin, Npart, profile, slope, center, process, t_equilibration, t_randomize, h, t_factor, R2, Ngamma2, x_basis2, x);

  for(int req=0; req<R1; req++)
  {
      //Initialization
      InitialProfile(x, &rho, Npart, profile, slope, center);
       	      
      //Equilibration configuration ZRP
      if(process=="ZRP")
          ZRP_KMC_fixedBoundary(&rho, t_equilibration*t_factor);
      else
          cout<< "Invalid process" << endl;
      vector<int> rho_eq(rho);
      t=clock()-t;
      printf("It took me %f seconds to equilibrate.\n", ((float)t)/CLOCKS_PER_SEC);
      
      //Saving parameters on data file
      sprintf(filename, "Data_%d.m",req);
      full_filename1 = file_path1 + filename;
      full_filename2 = file_path2 + filename;
      
      myfile.open(full_filename1);
      if(myfile.is_open()){
	  myfile << "\%It took " << ((float)t)/CLOCKS_PER_SEC << " seconds to equilibrate\n";
      myfile << "\%<rho_initial,gamma1>, <rho_final,gamma1>, <rho_initial,gamma2>, <rho_final,gamma2>, ...\n";
      myfile << "rho_gamma_initial_final = [";
	  myfile.close();
      }
        
      myfile.open(full_filename2);
      if(myfile.is_open()){
          myfile << "\%It took " << ((float)t)/CLOCKS_PER_SEC << " seconds to equilibrate\n";
          myfile << "\%<rho_initial,gamma1>, <rho_final,gamma1>, <rho_initial,gamma2>, <rho_final,gamma2>, ...\n";
          myfile << "rho_gamma_initial_final = [";
          myfile.close();
      }
        
      //Various realizations
      vector<int> rho_initial(rho);
      if(process=="ZRP")
	  {
          for(int r=0; r<R2; r++)
	      {
	          rho=rho_eq;
	          ZRP_KMC_fixedBoundary(&rho, t_randomize*t_factor);
	          rho_initial=rho;
	          ZRP_KMC_fixedBoundary(&rho, h*t_factor);
	          Postprocessing(x, &rho_initial, &rho, rho_gamma_initial1, rho_gamma_final1, x_basis1, r, R2, Ngamma1, full_filename1);
              Postprocessing(x, &rho_initial, &rho, rho_gamma_initial2, rho_gamma_final2, x_basis2, r, R2, Ngamma2, full_filename2);
	      }
	  }
      else
	      cout<< "Invalid process" << endl;
      
      t=clock()-t;
      printf("It took me %f seconds to do %d samples.\n", ((float)t)/CLOCKS_PER_SEC, R2);
      myfile.open(full_filename1, ios::app);
      if(myfile.is_open()){
	      myfile << "\%It took " << ((float)t)/CLOCKS_PER_SEC << " seconds to do " << R2 << " samples";
	      myfile.close();
      }
      myfile.open(full_filename2, ios::app);
      if(myfile.is_open()){
          myfile << "\%It took " << ((float)t)/CLOCKS_PER_SEC << " seconds to do " << R2 << " samples";
          myfile.close();
      }

  }
  
  //Terminate the program
  return 0;
}
