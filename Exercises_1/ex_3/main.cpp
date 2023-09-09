#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

#include "random.h"

using namespace std;
 
int main (int argc, char *argv[]){

   Random rnd;      //generator of uniform random numbers in [0,1)

   int M = 100000;      //number of throws of the needle
   int N = 100;         //maximum number of blocks
   int L = M/N;         //number of throws per block

   double d = 1.0;  //distance between the lines
   double l = 0.8;  //length of the needle

   
   //setting the generator
   rnd.Initialize();

   //parameters I will need
   double y_1, y_2;
   int counter = 0; //all throws
   int inter_counter = 0;   //intersections
   double prob;

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

    for(int i=0; i<N; i++){                         //storing average and squared average per block

        counter = 0;
        inter_counter = 0;
        for(int j=0; j<L; j++){
            y_1 = rnd.Rannyu(0,d);
            if(y_1==0 or y_1==d)
                inter_counter++;
            else{
            double a = rnd.Rannyu(-1,1);
            y_2 = a * l;
            if(y_2<=0 or y_2>=d)
                inter_counter++;
            }
            counter++;
        }
        prob = (double) (inter_counter) / (double) (counter);
        
        blk_ave[i] = 2.0 * l / (prob * d);
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
    out.open("pi.out");
    for(int i=0; i<N; i++)
        out << i+1 << "\t" << sum_prog[i] << "\t" << err_prog[i] << endl;     //printing number of blocks - global average - global error
        
    out.close();


   rnd.SaveSeed();      


   return 0;

}
