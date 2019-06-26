#include <iostream>
#include <fstream>
#include <ctime>
#include <vector>
#include <string>
#include <math.h>
using namespace std;

void Postprocessing(double x[], vector<int>* rho_initial_adr, vector<int>* rho_adr, double* rho_gamma_initial, double* rho_gamma_final, double* x_basis, int r, int R, int Ngamma, string filename1)
{
//  Parameter declaration
    ofstream myfile;
    int Nbin = (*rho_adr).size();      //Number of bins
    double dx = 1/(double)Nbin;          //size of bin
    double gamma_x;
//  Postprocessing and saving
    myfile.open(filename1, ios::app);
    myfile.precision(11);
    for(int j=0; j<Ngamma; j++)
    {
        int k=0;
        rho_gamma_initial[r*Ngamma+j]=0;
        rho_gamma_final[r*Ngamma+j]=0;
        while(x[k] < (x_basis[j] - 1.0/(double)Ngamma)){
            k++;
        }
	
        while(x[k] < (x_basis[j] + 1.0/(double)Ngamma)){
            gamma_x = max(0., 1.0-(double)Ngamma*fabs(x[k]-x_basis[j]));
            rho_gamma_initial[r*Ngamma+j] = rho_gamma_initial[r*Ngamma+j] + gamma_x*((*rho_initial_adr)[k])*dx;
            rho_gamma_final[r*Ngamma+j] = rho_gamma_final[r*Ngamma+j] + gamma_x*((*rho_adr)[k])*dx;
            k++;
        }
        
        myfile << rho_gamma_initial[r*Ngamma+j] << ", " << rho_gamma_final[r*Ngamma+j];
        if(j == Ngamma-1){
            if(r == R-1)
                myfile << "];";
            else
                myfile << ";\n";
        }
        else
            myfile << ",";
    }
    myfile.close();

}