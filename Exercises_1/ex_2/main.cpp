#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>

#include "random.h"

using namespace std;
 
int main (int argc, char *argv[]){

    Random rnd;      //random numbers generator

    int M = 100000;      //number of values to plot
   
    //setting the generator
    rnd.Initialize();

    //to store values to plot
    vector<double> S_unif;
    vector<double> S_exp;
    vector<double> S_lorentz;

    
    //N=1

    for(int i=0; i<M; i++){
        S_unif.push_back(rnd.Rannyu());
        S_exp.push_back(rnd.Exp(1.));
        S_lorentz.push_back(rnd.Lorentz(0.,1.));
    }

    ofstream out("plot_1.out");
    for(int i=0; i<M; i++)
        out << S_unif[i] << "\t" << S_exp[i] << "\t" << S_lorentz[i] << endl;
    out.close();

    
    //N=2

    for(int i=0; i<M; i++){
        S_unif[i] = (rnd.Rannyu()+rnd.Rannyu())/2.;
        S_exp[i] = (rnd.Exp(1.)+rnd.Exp(1.))/2.;
        S_lorentz[i] = (rnd.Lorentz(0.,1.)+rnd.Lorentz(0.,1.))/2.;
    }

    out.open("plot_2.out");
    for(int i=0; i<M; i++)
        out << S_unif[i] << "\t" << S_exp[i] << "\t" << S_lorentz[i] << endl;
    out.close();


    double sum_unif;
    double sum_exp;
    double sum_lorentz;


    //N=10

    for(int i=0; i<M; i++){

        sum_unif = 0.;
        sum_exp = 0.;
        sum_lorentz = 0.;

        for(int j=0; j<10; j++){
            sum_unif += rnd.Rannyu();
            sum_exp += rnd.Exp(1.);
            sum_lorentz += rnd.Lorentz(0.,1.);
        }
        S_unif[i] = sum_unif/10.;
        S_exp[i] = sum_exp/10.;
        S_lorentz[i] = sum_lorentz/10.;

    }

    out.open("plot_10.out");
    for(int i=0; i<M; i++)
        out << S_unif[i] << "\t" << S_exp[i] << "\t" << S_lorentz[i] << endl;
    out.close();


    //N=100

    for(int i=0; i<M; i++){

        sum_unif = 0.;
        sum_exp = 0.;
        sum_lorentz = 0.;

        for(int j=0; j<100; j++){
            sum_unif += rnd.Rannyu();
            sum_exp += rnd.Exp(1.);
            sum_lorentz += rnd.Lorentz(0.,1.);
        }
        S_unif[i] = sum_unif/100.;
        S_exp[i] = sum_exp/100.;
        S_lorentz[i] = sum_lorentz/100.;

    }

    out.open("plot_100.out");
    for(int i=0; i<M; i++)
        out << S_unif[i] << "\t" << S_exp[i] << "\t" << S_lorentz[i] << endl;
    out.close();

  
    rnd.SaveSeed();


    return 0;

}
