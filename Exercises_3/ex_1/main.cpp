#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>

#include "random.h"

using namespace std;

 
int main (int argc, char *argv[]){

    Random rnd;          //generator of uniform random numbers in [0,1)

    int M = 100000;      //number of throws
    int N = 100;         //maximum number of blocks
    int L = M/N;         //number of steps per block

   
    //setting the generator
    rnd.Initialize();

    //parameters
    double S_0 = 100.;  //initial asset price
    double T = 1.;      //delivery time
    double K = 100.;    //prescribed price
    double r = 0.1;     //risk-free interest rate
    double v = 0.25;    //volatility

    vector<double> blk_ave_C;           //for call option
    vector<double> blk_ave2_C;
    vector<double> sum_prog_C;
    vector<double> sum2_prog_C;
    vector<double> err_prog_C;
    double sum_C;

    vector<double> blk_ave_P;           //for put option
    vector<double> blk_ave2_P;
    vector<double> sum_prog_P;
    vector<double> sum2_prog_P;
    vector<double> err_prog_P;
    double sum_P;

    for(int i=0;i<N;i++){

            blk_ave_C.push_back(0.);
            blk_ave2_C.push_back(0.);
            sum_prog_C.push_back(0.);
            sum2_prog_C.push_back(0.);
            err_prog_C.push_back(0.);

            blk_ave_P.push_back(0.);
            blk_ave2_P.push_back(0.);
            sum_prog_P.push_back(0.);
            sum2_prog_P.push_back(0.);
            err_prog_P.push_back(0.);

   }

    double S;   //price at time t


    //DIRECT SAMPLING

    for(int i=0; i<N; i++){                         //storing average and squared average per block

        sum_C = 0., sum_P = 0.;
        for(int j=0; j<L; j++){
            S = S_0 * exp( (r-pow(v,2)/2.)*T + v*rnd.Gauss(0,T) );      //computing price directly at time T
            sum_C += exp(-r*T) * max(0.,S-K);
            sum_P += exp(-r*T) * max(0.,K-S);
        }
        blk_ave_C[i] = sum_C / L;
        blk_ave2_C[i] = blk_ave_C[i]*blk_ave_C[i];
        blk_ave_P[i] = sum_P / L;
        blk_ave2_P[i] = blk_ave_P[i]*blk_ave_P[i];

    }

    for(int i=0; i<N; i++){

        for(int j=0; j<i+1; j++){
            sum_prog_C[i] += blk_ave_C[j];      //sum from j=0 to j=i (for i+1 blocks)
            sum2_prog_C[i] += blk_ave2_C[j];     //sum2 from j=0 to j=i (for i+1 blocks)
            sum_prog_P[i] += blk_ave_P[j];
            sum2_prog_P[i] += blk_ave2_P[j];
        }

        sum_prog_C[i]/=(i+1);                 //progressive average for i blocks
        sum2_prog_C[i]/=(i+1);                //progressive squared average for i blocks
        sum_prog_P[i]/=(i+1);
        sum2_prog_P[i]/=(i+1);
        if(i==0)                             //progressive error for i blocks (if i=0 err=0)
            err_prog_C[i] = 0., err_prog_P[i] = 0.;
        else{
            err_prog_C[i] = sqrt((sum2_prog_C[i] - pow(sum_prog_C[i],2))/i);
            err_prog_P[i] = sqrt((sum2_prog_P[i] - pow(sum_prog_P[i],2))/i);
        }
    }

    ofstream out;
    out.open("direct_sampling.out");
    for(int i=0; i<N; i++)
        out << i+1 << "\t" << sum_prog_C[i] << "\t" << err_prog_C[i] << "\t" << sum_prog_P[i] << "\t" << err_prog_P[i] << endl;
        //printing call option-error-put option-error per increasing number of blocks i+1
    out.close();


    for(int i=0;i<N;i++){

            blk_ave_C[i] = 0.;
            blk_ave2_C[i] = 0.;
            sum_prog_C[i] = 0.;
            sum2_prog_C[i] = 0.;
            err_prog_C[i] = 0.;

            blk_ave_P[i] = 0.;
            blk_ave2_P[i] = 0.;
            sum_prog_P[i] = 0.;
            sum2_prog_P[i] = 0.;
            err_prog_P[i] = 0.;

   }


    //DISCRETIZED SAMPLING

    double t, t_0;

    for(int i=0; i<N; i++){                         //storing average and squared average per block

        sum_C = 0., sum_P = 0.;
        for(int j=0; j<L; j++){
            S = S_0;
            t_0 = 0.;
            for(int l=1; l<=100; l++){                //computing price at time T after 100 steps in time
                t = t_0 + T/100.;               
                S = S * exp( (r-pow(v,2)/2.)*(t-t_0) + v*rnd.Gauss(0,1)*sqrt(t-t_0) );
                t_0 = t;
            }
            sum_C += exp(-r*T) * max(0.,S-K);
            sum_P += exp(-r*T) * max(0.,K-S);
        }
        blk_ave_C[i] = sum_C / L;
        blk_ave2_C[i] = blk_ave_C[i]*blk_ave_C[i];
        blk_ave_P[i] = sum_P / L;
        blk_ave2_P[i] = blk_ave_P[i]*blk_ave_P[i];

    }

    for(int i=0; i<N; i++){

        for(int j=0; j<i+1; j++){
            sum_prog_C[i] += blk_ave_C[j];      //sum from j=0 to j=i (for i+1 blocks)
            sum2_prog_C[i] += blk_ave2_C[j];     //sum2 from j=0 to j=i (for i+1 blocks)
            sum_prog_P[i] += blk_ave_P[j];
            sum2_prog_P[i] += blk_ave2_P[j];
        }

        sum_prog_C[i]/=(i+1);                 //progressive average for i blocks
        sum2_prog_C[i]/=(i+1);                //progressive squared average for i blocks
        sum_prog_P[i]/=(i+1);
        sum2_prog_P[i]/=(i+1);
        if(i==0)                             //progressive error for i blocks (if i=0 err=0)
            err_prog_C[i] = 0., err_prog_P[i] = 0.;
        else{
            err_prog_C[i] = sqrt((sum2_prog_C[i] - pow(sum_prog_C[i],2))/i);
            err_prog_P[i] = sqrt((sum2_prog_P[i] - pow(sum_prog_P[i],2))/i);
        }
    }    

    out.open("discretized_sampling.out");
    for(int i=0; i<N; i++)
        out << i+1 << "\t" << sum_prog_C[i] << "\t" << err_prog_C[i] << "\t" << sum_prog_P[i] << "\t" << err_prog_P[i] << endl;
        //printing call option-error-put option-error per increasing number of blocks i+1
    out.close();


  
    rnd.SaveSeed();      

    return 0;

}












