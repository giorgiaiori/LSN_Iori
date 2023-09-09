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


   //ex. 1.1

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
        for(int j=0; j<L; j++)
            sum += rnd.Rannyu();
        
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
    out.open("part1.out");
    for(int i=0; i<N; i++)
        out << i+1 << "\t" << sum_prog[i] << "\t" << err_prog[i] << endl;     //printing number of blocks - global average - global error
        
    out.close();

    //ex. 1.2

    //same procedure as above, but changing the measured values

    for(int i=0;i<N;i++){
        blk_ave[i] = 0.;
        blk_ave2[i] = 0.;
        sum_prog[i] = 0.;
        sum2_prog[i] = 0.;
        err_prog[i] = 0.;
   }

    for(int i=0; i<N; i++){                         //storing average and squared average per block

        sum = 0.;
        for(int j=0; j<L; j++)
            sum += pow( rnd.Rannyu()-0.5,2 );
        
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

    out.open("part2.out");
    for(int i=0; i<N; i++)
        out << i+1 << "\t" << sum_prog[i] << "\t" << err_prog[i] << endl;     //printing number of blocks - global average - global error

    out.close();



    //ex. 1.3

   int m = 100;                 //number of sub-intervals in [0,1)
   int n = 10000;               //number of throws for each chi-square
   double expected = n/m;       //expected number of events in each sub-interval

   int n_chi = 100;             //number of chi-square computations 

   double chi_square;     
   vector<int> n_events;        //n_events[k] contains the number of events occured in the sub-interval k of [0,1), which is [k/m,(k+1)/m)      
   double x;                    //to store casual numbers

    for(int k=0; k<m; k++)
            n_events.push_back(0);

   out.open("part3.out");

   for(int j=0; j<n_chi; j++){

        chi_square = 0.;
        for(int k=0; k<m; k++)
            n_events[k] = 0;

        for(int i=0; i<n; i++){

            x = rnd.Rannyu();

            for(int k=0; k<m; k++)
                if( x >= (double) k/m and x < (double) (k+1)/m ){
                    n_events[k]++;
                    break;
                }

        }

        for(int k=0; k<m; k++)
            chi_square += pow( n_events[k]-expected,2 );

        chi_square = chi_square/expected;
        out << j+1 << "\t" << chi_square << endl;       //printing number of computation - value of chi-square  

   }

   out.close(); 

   rnd.SaveSeed();      


   return 0;

}
