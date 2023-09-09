#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include "VariationalMC.h"

using namespace std;

int main(){

    Input();    //Initialization

    for(int itemp=1; itemp<=ntemp; itemp++){    //Cicle over temperature steps

        SimulatedAnnealing();
        PrintCurrent(itemp);
        DecreaseTemp();

    }

    ConfFinal();

    ComputeEnergy(true,true);   //for the final configuration, compute energy with error and psi_squared
    HistoPsiSq();

    return 0;

}









void Input(void){

    ifstream ReadInput;

    cout << "Variational Monte Carlo for 1 particle in 1D with simulated annealing" << endl;
    cout << "Potential: V(x) = x^4 - 5/2 * x^2" << endl;
    cout << "Trial ansatz: sum of two gaussians centered in +-mu with erf sigma" << endl;
    cout << "The program uses h_bar = mass = K_Boltz = 1 units " << endl;

    //Read input informations

    ReadInput.open("input.in");

    //Read seed for random numbers
    int p1, p2;
    ifstream Primes("Primes");
    Primes >> p1 >> p2 ;
    Primes.close();

    ReadInput >> restart;

    ifstream Seed;
    if(restart) 
        Seed.open("seed.out");
    else 
        Seed.open("seed.in");
    Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
    rnd.SetRandom(seed,p1,p2);
    Seed.close();
   
    ReadInput >> delta;      //delta Metropolis that needs to be regulated in order to have an acceptance of 50%
    ReadInput >> nblk;       //number of blocks
    ReadInput >> nstep;      //number of Monte Carlo steps in a block
    ReadInput >> nstart;     //number of equilibration steps

    cout << endl;
    cout << "Information about computing the energy" << endl;
    cout << "The program performs Metropolis moves with uniform translations" << endl;
    cout << "Moves parameter = " << delta << endl;
    cout << "Number of blocks = " << nblk << endl;
    cout << "Number of steps in one block = " << nstep << endl << endl;

    ReadInput >> nbins;
    cout << "Number of bins for squared wave function = " << nbins << endl;
    bin_size = sample_range/double(nbins);
    cout << "Size of bins = " << bin_size << endl;

    ReadInput >> delta_mu;
    ReadInput >> delta_sigma;
    ReadInput >> nanneal;

    cout << endl;
    cout << "Information about simulated annealing" << endl;
    cout << "The program performs Metropolis moves with uniform translations" << endl;
    cout << "Moves parameter for mu = " << delta_mu << endl;
    cout << "Moves parameter for sigma = " << delta_sigma << endl;
    cout << "Number of steps for each temperature = " << nanneal << endl;
    
    ReadInput >> temp_start;
    ReadInput >> temp_min;
    ReadInput >> alpha;
    temp = temp_start;
    beta = 1.0/temp;
    ntemp = int( log(temp_min/temp_start)/log(alpha) );
    cout << "Initial temperature = " << temp_start << endl;
    cout << "Minimum final temperature = " << temp_min << endl;
    cout << "Number of temperatures explored = " << ntemp << endl;
    cout << "Decreasing rate = " << alpha << endl;

    if(restart==0){
        ReadInput >> mu0;
        ReadInput >> sigma0;
    }

    ReadInput.close();
    
    //Initial configuration
    cout << endl;
    cout << "Reading initial configuration... ";

    ifstream ReadConf;

    if(restart){
        ReadConf.open("config.out");
        ReadConf >> mu;
        ReadConf >> sigma;
    }
    else{
        mu = mu0;
        sigma = sigma0;
    }

    ReadConf.close();

    mu_old = mu;
    sigma_old = sigma;

    cout << "Done!" << endl;

    cout << "Starting parameters for trial ansatz:" << endl;
    cout << "mu = " << mu << endl;
    cout << "sigma = " << sigma << endl;

    //Initializing the bins of the squared wave function
    count_histo = 0;
    for(int k=0; k<nbins; k++)
        psisq[k] = 0.0;

    //initial configuration for x
    x = mu0;
    x_old = x;

    do_i_compute_error = false;

    ComputeEnergy(do_i_compute_error,false);
    ene_old = ene;
    error_ene = 0.; //now I don't care about the error
    error_ene_old = error_ene;

    return;

}


void Move(bool psi_sample){

    double p_old, p_new, p_acc;

    p_old = psi_squared(x_old);

    //Metropolis move
    x = x_old + delta*(rnd.Rannyu() - 0.5);
    
    p_new = psi_squared(x);

    p_acc = p_new/p_old;

    if(p_acc>=rnd.Rannyu()){
        x_old = x;
        accepted = accepted + 1.;
    }
    else
        x = x_old;

    attempted = attempted + 1.;

    if(psi_sample){
        int bin = ComputeBin(x);
        if(bin>=0 and bin<nbins){   // histogram only in sample_range
            psisq[bin] += 1;
            count_histo++;
        }
    }

}


void Reset(int iblk){ //Reset block averages
   
    if(iblk == 1){
        glob_av = 0.;
        glob_av2 = 0.;
    }

    blk_av = 0.;

    blk_norm = 0;
    attempted = 0;
    accepted = 0;

}


void Accumulate(void){ //Update block averages

    blk_av = blk_av + walker;
    blk_norm = blk_norm + 1.0;

}


double Error(double sum, double sum2, int iblk){
    return sqrt(fabs(sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)iblk);
}


void ConfFinal(void){

    ofstream WriteConf, WriteSeed;

    cout << "Print final configuration to file config.out" << endl << endl;
    WriteConf.open("config.out");
    WriteConf << mu << endl << sigma << endl;
    WriteConf.close();

    rnd.SaveSeed();

}


void Measure(){ //Properties measurement

    walker = local_energy(x);

    return;

}


double local_energy(double x){

    double factor = 1./pow(sigma,2);
    double A = exp(- factor/2.*pow(x-mu,2));
    double B = exp(- factor/2.*pow(x+mu,2));

    double kin_term = factor/2. * ( 1. - factor * ( pow(x,2) + pow(mu,2) - 2.*mu*x*(A-B)/(A+B) ) );
    double pot_term = potential(x);

    return kin_term + pot_term;

}


double potential(double x){
    return pow(x,4)-5./2.*pow(x,2);
}


double psi_squared(double x){

    double factor = 1./pow(sigma,2);
    double A = exp(- factor/2.*pow(x-mu,2));
    double B = exp(- factor/2.*pow(x+mu,2));
    double psi = A + B;

    return psi*psi;

}


int ComputeBin(double x){

    double x_trasl = x + sample_range/2.;  // shifting x in order to get only positive values
    int bin = (int) (x_trasl/bin_size);

    return bin;

}


void HistoPsiSq(void){

    for(int k=0; k<nbins; k++)      // normalizing the histogram
        psisq[k] /= (count_histo*bin_size);

    ofstream Histo;
    const int wd=12;
    Histo.open("psi_squared_final.dat",ios::app);

    double x_bin;
    for(int k=0; k<nbins; k++){ 
        x_bin = (k+0.5)*bin_size - sample_range/2.; // shifting x to its right place
        Histo << setw(wd) << x_bin << setw(wd) << psisq[k] << endl;
    }

    Histo.close();

    return;

}


// FOR SIMULATED ANNEALING

void DecreaseTemp(void){
    temp = temp*alpha;
    beta = 1.0/temp;
    return;
}


void SimulatedAnnealing(){

    do_i_compute_error = false;

    for(int iequi=1; iequi<=nstart; iequi++)    // Equilibration
        Move_Params(do_i_compute_error);

    for(int n=1; n<=nanneal; n++){
        if(n==nanneal)
            do_i_compute_error = true;
        Move_Params(do_i_compute_error);
    }

    return;

}


void Move_Params(bool compute_error){

    double p_old, p_new, p_acc;

    //Metropolis move
    mu = mu_old + delta_mu*(rnd.Rannyu() - 0.5);
    sigma = sigma_old + delta_sigma*(rnd.Rannyu() - 0.5);
    ComputeEnergy(compute_error,false);   //computes new energy and registers it in ene

    p_old = exp(-beta*ene_old);
    p_new = exp(-beta*ene);

    p_acc = p_new/p_old;

    if(p_acc>=rnd.Rannyu()){
        mu_old = mu;
        sigma_old = sigma;
        ene_old = ene;
        error_ene_old = error_ene;
        accepted_par = accepted_par + 1.;
    }
    else{
        mu = mu_old;
        sigma = sigma_old;
        ene = ene_old;
        error_ene = error_ene_old;
    }
    attempted_par = attempted_par + 1.;

}


void PrintCurrent(int itemp){

    ofstream Energy, Params;
    const int wd=12;
    
    cout << endl;
    cout << "Temperature-step number " << itemp << endl;
    cout << "Acceptance rate " << accepted_par/attempted_par << endl << endl;

    accepted_par = 0.0;
    attempted_par = 0.0;
    
    Energy.open("energy.dat",ios::app);
    Params.open("parameters.dat",ios::app);

    Energy << setw(wd) << itemp << setw(wd) << temp << setw(wd) << ene << setw(wd) << error_ene << endl;
    Params << setw(wd) << itemp << setw(wd) << temp << setw(wd) << mu << setw(wd) << sigma << endl;

    cout << "----------------------------" << endl << endl;

    Energy.close();
    Params.close();

    return;

}


void ComputeEnergy(bool compute_error, bool psi_sample){

    int nstart_x = 2000;

    for(int i=0; i<nstart_x; i++)   //Equilibration
        Move(false);

    if (psi_sample){    //in the final computation of energy do it well even if before it had not so much statystics
        nblk = 100;
        nstep = 1000;
    }

    for(int iblk=1; iblk<=nblk; iblk++){    //Cicle over blocks

        Reset(iblk);    //Reset block averages

        for(int istep=1; istep<=nstep; istep++){    //Cicle of steps in each block
            Move(psi_sample); //Metropolis move
            Measure();  //Measure the average of the hamiltonian
            Accumulate(); //Update block averages
        }

        ofstream Acc;   // monitoring the acceptance for x moves
        Acc.open("acceptance.dat",ios::app);
        Acc << "Block number " << iblk << endl;
        Acc << "Acceptance rate " << accepted/attempted << endl << endl; 

        stima_ene = blk_av/blk_norm;
        glob_av += stima_ene;
        glob_av2 += stima_ene*stima_ene;
        if(compute_error)
            err_ene = Error(glob_av,glob_av2,iblk);

        if(psi_sample){
            ofstream Energy;
            Energy.open("final_energy.dat",ios::app);
            const int wd=12;
            Energy << setw(wd) << iblk << setw(wd) << stima_ene << setw(wd) << glob_av/(double)iblk << setw(wd) << err_ene << endl;
            Energy.close();
        }

    }

    ene = glob_av/(double)nblk;
    error_ene = err_ene;

}