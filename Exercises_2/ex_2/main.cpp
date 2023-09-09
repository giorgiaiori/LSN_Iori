#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

#include "random.h"

using namespace std;
 
int main (int argc, char *argv[]){

    Random rnd;      //generator of uniform random numbers in [0,1)

    int M = 10000;      //number of throws
    int N = 100;         //maximum number of blocks
    int L = M/N;         //number of steps per block

    int T = 100;        //number of time steps in each RW

   
    //setting the generator
    rnd.Initialize();

    
    vector<double> sum;     //to accumulate the distances of the j=1,...,100 RW at fixed step i; its dimension must then be 100
    vector<double> sum2;    //to accumulate the squared distances of the j=1,...,100 RW at fixed step i; dim=100
    vector<double> err;     //to store the uncertainty associated to steps j=1,...,100; dim=100

    //data blocking
    vector<double> blk_ave;
    for(int j=0; j<T; j++)
        blk_ave.push_back(0.);
    
    for(int j=0; j<T; j++){
        sum.push_back(0.);
        sum2.push_back(0.);
        err.push_back(0.);
    }

    //random walk in the continuum

    double x0,y0,z0;                    //initial position 
    double x,y,z;                       //position at the next step                       

    double theta, phi;                  //angles to define the direction
    double a = 1.;                      //length of each step

    for(int iblk=1; iblk<=N; iblk++){       //in each block I have L RWs of T=100 steps each

        for(int ithrow=0; ithrow<L; ithrow++){              //data blocking
    
            x0 = 0.;        //starting from the origin
            y0 = 0.;
            z0 = 0.;

            for(int j=0; j<T; j++){     //j is the step in each block

                theta = rnd.Rannyu(0.,M_PI);
                phi = rnd.Rannyu(0.,2*M_PI);

                x = x0 + a * sin(theta) * cos(phi);
                y = y0 + a * sin(theta) * sin(phi);
                z = z0 + a * cos(theta);
    
                blk_ave[j] += sqrt(x*x + y*y + z*z);      //keeping track of the position in step j for each throw

                x0 = x;     //updating the current position which the next step will start from
                y0 = y;
                z0 = z;

            }

        }
        
        for(int j=0; j<T; j++){             //for each block I accumulate the mean of the L RWs at fixed time j in sum[j]
            sum[j] += blk_ave[j]/L;              
            sum2[j] += sum2[j]*sum2[j];
            blk_ave[j] = 0.;                //cleaning up blk_ave for the next block
        }

    }

    for(int j=0; j<T; j++){         //for each fixed step j I have N values of the position: 
        sum[j] = sum[j]/N;          //I take their mean
        sum2[j] = sum2[j]/N;
        err[j] = sqrt(fabs(sum2[j] - pow(sum[j],2))/N);  //and their uncertainty
    }

    ofstream out;
    out.open("RW_continuum.out");
    out << 0 << "\t" << 0 << "\t" << 0 << endl; //printing the starting point (position=0, error=0)   
    for(int j=0; j<L; j++)
        out << j+1 << "\t" << sum[j] << "\t" << err[j] << endl; //printing step in time-mean-error

    out.close();


    for(int j=0; j<T; j++){
        sum[j] = 0.;
        sum2[j] = 0.;
        err[j] = 0.;
    }

    //random walk in a cubic lattice

    double num;      //to choose the direction

    for(int iblk=1; iblk<=N; iblk++){       //in each block I have L RWs of T=100 steps each
    
        for(int ithrow=0; ithrow<L; ithrow++){              //data blocking

            x0 = 0.;        //starting from the origin
            y0 = 0.;
            z0 = 0.;

            for(int j=0; j<T; j++){     //j is the step in each block

                num = rnd.Rannyu();

                if(num>=0 and num<1./6.){                   //moving in direction x with probability 1/3 (+a prob. 1/6, -a prob. 1/6)
                    x = x0 + a;
                    y = y0;
                    z = z0;
                }

                if(num>=1./6. and num<2./6.){
                    x = x0 - a;
                    y = y0;
                    z = z0;
                }

                if(num>=2./6. and num<3./6.){               //moving in direction y with probability 1/3 (+a prob. 1/6, -a prob. 1/6)
                    x = x0;
                    y = y0 + a;
                    z = z0;
                }

                if(num>=3./6. and num<4./6.){
                    x = x0;
                    y = y0 - a;
                    z = z0;
                }

                if(num>=4./6. and num<5./6.){               //moving in direction z with probability 1/3 (+a prob. 1/6, -a prob. 1/6)
                    x = x0;
                    y = y0;
                    z = z0 + a;
                }

                if(num>=5./6. and num<1){
                    x = x0;
                    y = y0;
                    z = z0 - a;
                }
    
                blk_ave[j] += sqrt(x*x + y*y + z*z);      //keeping track of the position in step j for each throw

                x0 = x;     //updating the current position which the next step will start from
                y0 = y;
                z0 = z;

            }

        }

        for(int j=0; j<T; j++){             //for each block I accumulate the mean of the L RWs at fixed time j in sum[j]
            sum[j] += blk_ave[j]/L;              
            sum2[j] += sum2[j]*sum2[j];
            blk_ave[j] = 0.;                //cleaning up blk_ave for the next block
        }

    }

    for(int j=0; j<T; j++){         //for each fixed step j I have N values of the position: 
        sum[j] = sum[j]/N;          //I take their mean
        sum2[j] = sum2[j]/N;
        err[j] = sqrt(fabs(sum2[j] - pow(sum[j],2))/N);  //and their uncertainty
    }    
    
    out.open("RW_cubiclattice.out");
    out << 0 << "\t" << 0 << "\t" << 0 << endl; //printing the starting point (position=0, error=0)   
    for(int j=0; j<L; j++)
        out << j+1 << "\t" << sum[j] << "\t" << err[j] << endl; //printing step in time-mean-error

    out.close();
    


   rnd.SaveSeed();      


   return 0;

}
