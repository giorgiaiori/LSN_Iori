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
double walker;
double temp, temp_start, temp_min, alpha;
double beta;
double ene, error_ene, ene_old, error_ene_old;
double psisq[1000];

// averages
double blk_av, blk_norm, accepted, attempted, accepted_par, attempted_par;
double glob_av, glob_av2;
double stima_ene;
double err_ene;

//configuration
double x;
double x_old;
double mu, sigma;
double mu0, sigma0;
double mu_old, sigma_old;

// simulation
int nstep, nblk, restart, nstart, ntemp, nanneal;
double delta, delta_mu, delta_sigma;
int nbins;
double bin_size;
int count_histo;
double sample_range = 6.0;
bool do_i_compute_error;

//functions
void Input(void);
void Reset(int);
void Accumulate(void);
//void Averages(int);
void Move(void);
void ConfFinal(void);
void Measure(void);
double Error(double,double,int);
double local_energy(double);
double potential(double);
double psi_squared(double);
int ComputeBin(double);
void HistoPsiSq(void);
void SimulatedAnnealing(void);
void PrintCurrent(int);
void DecreaseTemp(void);
void Move_Params(bool);
void ComputeEnergy(bool,bool);

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
