#define _USE_MATH_DEFINES  //Needed to get pi: M_PI

#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <random>

using namespace std;

void InitialProfile(double x[], vector<int>* rho, int Npart, string profile, double slope, double center)
{
    int Nbins = (*rho).size();
    double average = (double)Npart/(double)Nbins;  //Average number of particles per bin
    random_device rd;
    mt19937_64 seed(rd());
    uniform_real_distribution<double> distribution (0.0,1.0);
    double rand_no;
    vector<int>::iterator it;
    int j;
    int cumsum;
    if(profile == "flat"){ //Flat profile
        fill((*rho).begin(),(*rho).end(),(int)average);
    }
    else if(profile == "linear") //Linear profile
    {
        *((*rho).begin()) = slope*x[0] - 1.0/2.0*slope + center;
        *((*rho).end()-1) = slope*x[Nbins-1] - 1.0/2.0*slope + center;
        for(it=(*rho).begin()+1,j=1; it!=(*rho).end()-1; it++, j++)
        {
            cumsum = 0.0;
            for(int k=0; k<j; k++){
                cumsum += (*rho)[k];
            }
            (*it) = (int)((double)Nbins*(1.0/2.0*slope*pow(x[j],2) + (center-slope/2.0)*x[j])) - cumsum;
            rand_no = distribution(seed);
            if (rand_no < 0.5){
                (*it) -= 1;
            }
            else{
                (*it) += 1;
            }
            if((*it) < 0){
                (*it) = 0;
            }
        }
    }
    else if(profile == "triangle")  //Triangle profile
    {
        *((*rho).begin()) = slope*x[0] - 1.0/4.0*slope + center;
        *((*rho).end()-1) = -slope*x[Nbins-1] + 3.0/4.0*slope + center;
        for(it=(*rho).begin()+1,j=1; it!=(*rho).end()-1; it++, j++)
        {
            cumsum = 0.0;
            for(int k=0; k<j; k++){
                cumsum += (*rho)[k];
            }
            if(x[j] <= 0.5){
                (*it) = (int)((double)Nbins*(1.0/2.0*slope*pow(x[j],2) + (center-slope/4.0)*x[j])) - cumsum;
            }
            else{
                (*it) = (int)((double)Nbins*(center/2.0-1.0/2.0*slope*(pow(x[j],2)-1.0/4.0) + (center+slope*3.0/4.0)*(x[j]-1.0/2.0)) - cumsum);
            }
            rand_no = distribution(seed);
            if (rand_no < 0.5){
                (*it) -= 1;
            }
            else{
                (*it) += 1;
            }
            if((*it) < 0){
                (*it) = 0;
            }
        }
    }
    else
    {
        cout << "Warning - unknown profile " << profile << ".\n";
    }
    return;
}
