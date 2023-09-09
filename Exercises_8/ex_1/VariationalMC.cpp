#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include "VariationalMC.h"

using namespace std;

int main(){

    Input();    //Initialization: defining the trial ansatz

    if(restart==0){
        for(int istep=1; istep <= nstart; istep++)  //Equilibration
            Move(false);
    }

    for(int iblk=1; iblk<=nblk; iblk++){    //Cicle over blocks

        Reset(iblk);    //Reset block averages

        //ChangeStart(iblk,2);  // Change starting point if iblk is a multiple of n

        for(int istep=1; istep<=nstep; istep++){    //Cicle of steps in each block
            Move(true); //Metropolis move
            Measure();  //Measure the average of the hamiltonian
            Accumulate(); //Update block averages
            PrintConf();
        }

        Averages(iblk); //Print results of block iblk

    }

    ConfFinal();
    HistoPsiSq();

    return 0;

}









void Input(void){

    ifstream ReadInput;

    cout << "Variational Monte Carlo for 1 particle in 1D" << endl;
    cout << "Potential: V(x) = x^4 - 5/2 * x^2" << endl;
    cout << "Trial ansatz: sum of two gaussians centered in +-mu with erf sigma" << endl;
    cout << "The program uses h_bar = mass = 1 units " << endl;

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
    ReadInput >> nstart;

    cout << endl;
    cout << "The program performs Metropolis moves with uniform translations" << endl;
    cout << "Moves parameter = " << delta << endl;
    cout << "Number of blocks = " << nblk << endl;
    cout << "Number of steps in one block = " << nstep << endl << endl;

    ReadInput >> mu;
    ReadInput >> sigma;

    cout << endl;
    cout << "Starting parameters for trial ansatz:" << endl;
    cout << "mu = " << mu << endl;
    cout << "sigma = " << sigma << endl;

    ReadInput >> nbins;
    cout << "Number of bins for squared wave function = " << nbins << endl;
    bin_size = sample_range/double(nbins);
    cout << "Size of bins = " << bin_size << endl;

    ReadInput.close();

    //Initializing the bins of the squared wave function
    count_histo = 0;
    for(int k=0; k<nbins; k++)
        psisq[k] = 0.0;
    
    //Initial configuration
    cout << endl;
    cout << "Reading initial configuration... ";

    ifstream ReadConf;

    if(restart){
        ReadConf.open("config.out");
        ReadConf >> x;
    }
    else{
        x = mu;
    }

    ReadConf.close();

    x_old = x;
    
    cout << "Done!" << endl;

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


void Averages(int iblk){ //Print results for current block
    
    ofstream Energy;
    const int wd=12;
    
    cout << endl;
    cout << "Block number " << iblk << endl;
    cout << "Acceptance rate " << accepted/attempted << endl << endl;
    
    Energy.open("output_energy.dat",ios::app);

    stima_ene = blk_av/blk_norm;
    glob_av += stima_ene;
    glob_av2 += stima_ene*stima_ene;
    err_ene = Error(glob_av,glob_av2,iblk);

    Energy << setw(wd) << iblk << setw(wd) << stima_ene << setw(wd) << glob_av/(double)iblk << setw(wd) << err_ene << endl;

    cout << "----------------------------" << endl << endl;

    Energy.close();

}


double Error(double sum, double sum2, int iblk){
    return sqrt(fabs(sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)iblk);
}


void ConfFinal(void){

    ofstream WriteConf, WriteSeed;

    cout << "Print final configuration to file config.out" << endl << endl;
    WriteConf.open("config.out");
    WriteConf << x << endl;
    WriteConf.close();

    rnd.SaveSeed();

}


void Measure(){ //Properties measurement

    walker = local_energy(x);

    return;

}


double local_energy(double x){

    double a = x-mu;
    double b = x+mu;
    double A = exp( - pow(a,2)/(2.*pow(sigma,2)) );
    double B = exp( - pow(b,2)/(2.*pow(sigma,2)) );

    double second_derivative = - 1./(2.*pow(sigma,4)) * ( (pow(a,2)-pow(sigma,2))*A + (pow(b,2)-pow(sigma,2))*B );
    double psi = A + B;

    double kin_term = second_derivative/psi;

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


void PrintConf(){

    ofstream ConfNow;
    ConfNow.open("now_config.dat",ios::app);
    ConfNow << x << endl;
    ConfNow.close();

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
    Histo.open("psi_squared.dat",ios::app);

    double x_bin;
    for(int k=0; k<nbins; k++){ 
        x_bin = (k+0.5)*bin_size - sample_range/2.; // shifting x to its right place
        Histo << setw(wd) << x_bin << setw(wd) << psisq[k] << endl;
    }

    Histo.close();

    return;

}







