#include <iostream>
#include <fstream>
#include <ctime>
#include <vector>
#include <string>
#include <math.h>
#include <iomanip>
using namespace std;

void SavingParameters(string file_path, int Nbin, int Npart, string profile, double slope, double center, string process, double t_equilibration, double t_randomize, double h, double t_factor, int R2, int Ngamma, double x_basis[], double x[])
{
  ofstream myfile;
  char filename[200];
  double rhoLowBound = center-slope/2.0;
  double rhoUpBound = center+slope/2.0;
  sprintf(filename,"Summary.m");
  myfile.open(file_path+filename);
  if(myfile.is_open())
  {
      myfile << "\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\n";
      myfile << "\%Parameters of the simulation, " << profile << " profile, " << process << " process \n\n";
      myfile << "\%Slope is " << slope << " \n";
      myfile << "\%Density range is [" << rhoLowBound << ","<< rhoUpBound << "] \n";
      myfile << "\%Geometry\n";
      myfile << "Nbin = " << Nbin << ";\t\t\t \%Number of bins \n";
      myfile << "Npart = " << Npart << ";\t\t\t \%Number of particles \n\n";
      myfile << "\%Times\n";
      myfile << "t_prep - t_ini = " << t_equilibration << ";\t\t \%Time for equilibration \n";
      myfile << "t_0 - t_prep = " << t_randomize << ";\t\t \%Time to create different initial conditions\n";
      myfile << "h = " << h << ";\t\t\t \% Short interval of actual calculations\n";
      myfile << "t_factor = " << t_factor << ";\t\t \%T=t_factor*t (ratio between micro and macro scale)\n\n";
      myfile << "\%Sampling parameters\n";
      myfile << "R2 = " << R2 << ";\t\t\t\t \%Number of realizations\n\n";
      myfile << "\%Positions of bins are\n";
      myfile << "x = [" ;
      myfile << std::fixed;
      myfile << setprecision(6) << x[0];
      for(int j=1; j<Nbin; j++){
          myfile << ", ";
          myfile << setprecision(6) << x[j];
      }
      myfile << "];\n";
      myfile << "\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\n\n";
      myfile << "\%Filter function paratmeters\n";
      myfile << "Ngamma = " << Ngamma << ";\t\t\t \%Number of basis functions\n";
      myfile << "x_basis = [" << x_basis[0];
      for(int j=1; j<Ngamma; j++)
          myfile << ", " << x_basis[j];
      myfile << "];\n";
      myfile << "\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\n\n";
      myfile.close();
  }

  return;
}
