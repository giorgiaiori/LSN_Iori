#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

#include "random.h"

using namespace std;
 
int main (int argc, char *argv[]){

   Random rnd;      //generator of uniform random numbers in [0,1)

   int M = 100000;      //number of throws
   int N = 100;         //maximum number of blocks
   int L = M/N;         //number of steps per block

   
   //setting the generator
   rnd.Initialize();

   //temporary variables
   double x, y;


   //part 1 - mean method

   vector<double> blk_ave;
   vector<double> blk_ave2;
   vector<double> sum_prog;
   vector<double> sum2_prog;
   vector<double> err_prog;

   for(int i=0;i<N;i++){
        blk_ave.push_back(0.);
        blk_ave2.push_back(0.);
        sum_prog.push_back(0.);
        sum2_prog.push_back(0.);
        err_prog.push_back(0.);
   }
   double sum;

    for(int i=0; i<N; i++){                         //storing average and squared average per block

        sum = 0.;
        for(int j=0; j<L; j++){
            x = rnd.Rannyu();
            sum += M_PI/2. * cos( M_PI/2.*x );
        }
        
        blk_ave[i] = sum/L;
        blk_ave2[i] = blk_ave[i]*blk_ave[i];
        
    }

    for(int i=0; i<N; i++){

        for(int j=0; j<i+1; j++){
            sum_prog[i] += blk_ave[j];      //sum from j=0 to j=i (for i+1 blocks)
            sum2_prog[i] += blk_ave2[j];     //sum2 from j=0 to j=i (for i+1 blocks)
        }

        sum_prog[i]/=(i+1);                 //progressive average for i blocks
        sum2_prog[i]/=(i+1);                //progressive squared average for i blocks
        if(i==0)                             //progressive error for i blocks (if i=0 err=0)
            err_prog[i] = 0.;
        else
            err_prog[i] = sqrt((sum2_prog[i] - pow(sum_prog[i],2))/i);
    }

    ofstream out;
    out.open("mean_method.out");
    for(int i=0; i<N; i++)
        out << i+1 << "\t" << sum_prog[i] << "\t" << err_prog[i] << endl;     //printing number of blocks - global average - global error
        
    out.close();

    
    //part 2 - importance sampling

    
    // 1st attempt: p(x) = 2*(1-x) with inversion of the cumulative function

    for(int i=0;i<N;i++){
        blk_ave[i] = 0.;
        blk_ave2[i] = 0.;
        sum_prog[i] = 0.;
        sum2_prog[i] = 0.;
        err_prog[i] = 0.;
   }

    for(int i=0; i<N; i++){                         //storing average and squared average per block

        sum = 0.;
        for(int j=0; j<L; j++){
            x = 1. - sqrt(1.-rnd.Rannyu());
            sum += M_PI * cos(M_PI*x/2.) / (4.*(1.-x));
        }
        
        blk_ave[i] = sum/L;
        blk_ave2[i] = blk_ave[i]*blk_ave[i];
        
    }

    for(int i=0; i<N; i++){

        for(int j=0; j<i+1; j++){
            sum_prog[i] += blk_ave[j];      //sum from j=0 to j=i (for i+1 blocks)
            sum2_prog[i] += blk_ave2[j];     //sum2 from j=0 to j=i (for i+1 blocks)
        }

        sum_prog[i]/=(i+1);                 //progressive average for i blocks
        sum2_prog[i]/=(i+1);                //progressive squared average for i blocks
        if(i==0)                             //progressive error for i blocks (if i=0 err=0)
            err_prog[i] = 0.;
        else
            err_prog[i] = sqrt((sum2_prog[i] - pow(sum_prog[i],2))/i);
    }

    out.open("importance_sampling1.out");
    for(int i=0; i<N; i++)
        out << i+1 << "\t" << sum_prog[i] << "\t" << err_prog[i] << endl;     //printing number of blocks - global average - global error

    out.close();



    // 2nd attempt: p(x) = 3/2 * (1-x^2) with rejection technique

     for(int i=0;i<N;i++){
        blk_ave[i] = 0.;
        blk_ave2[i] = 0.;
        sum_prog[i] = 0.;
        sum2_prog[i] = 0.;
        err_prog[i] = 0.;
   }

    double p;       //to store the chosen probability density

    for(int i=0; i<N; i++){                         //storing average and squared average per block

        sum = 0.;
        for(int i=0; i<L; i++){

            x = rnd.Rannyu();
            y = rnd.Rannyu();

            p = 3./2. * ( 1.-pow(x,2) );

            if( y < 2./3.* p )                                          //accept
                sum += M_PI * cos(M_PI*x/2.) / (3.*(1.-pow(x,2)));
            else                                                        //reject
                i--;

        }
        
        blk_ave[i] = sum/L;
        blk_ave2[i] = blk_ave[i]*blk_ave[i];
        
    }

    for(int i=0; i<N; i++){

        for(int j=0; j<i+1; j++){
            sum_prog[i] += blk_ave[j];      //sum from j=0 to j=i (for i+1 blocks)
            sum2_prog[i] += blk_ave2[j];     //sum2 from j=0 to j=i (for i+1 blocks)
        }

        sum_prog[i]/=(i+1);                 //progressive average for i blocks
        sum2_prog[i]/=(i+1);                //progressive squared average for i blocks
        if(i==0)                             //progressive error for i blocks (if i=0 err=0)
            err_prog[i] = 0.;
        else
            err_prog[i] = sqrt((sum2_prog[i] - pow(sum_prog[i],2))/i);
    }

    out.open("importance_sampling2.out");
    for(int i=0; i<N; i++)
        out << i+1 << "\t" << sum_prog[i] << "\t" << err_prog[i] << endl;     //printing number of blocks - global average - global error

    out.close();




   rnd.SaveSeed();      


   return 0;

}
