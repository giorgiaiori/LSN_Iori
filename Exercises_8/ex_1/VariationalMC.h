/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#ifndef __varMC__
#define __varMC__

//Random numbers
#include "random.h"
int seed[4];
Random rnd;

//parameters, observables
double mu, sigma;
double walker;
double psisq[1000];

// averages
double blk_av, blk_norm, accepted, attempted;
double glob_av, glob_av2;
double stima_ene;
double err_ene;

//configuration
double x;
double x_old;

// simulation
int nstep, nblk, restart, nstart;
double delta;
int nbins;
double bin_size;
int count_histo;
double sample_range = 3.0;

//functions
void Input(void);
void Reset(int);
void Accumulate(void);
void Averages(int);
void Move(bool);
void ConfFinal(void);
void Measure(void);
double Error(double,double,int);
double local_energy(double);
double potential(double);
double psi_squared(double);
void PrintConf(void);
int ComputeBin(double);
void HistoPsiSq(void);

#endif

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
