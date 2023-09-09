#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <iomanip>

#include "random.h"
#include "classes.h"
#include "globals.h"

using namespace std;

int main(){ 


  // reading input information
  Input();
  

  // creating generation zero from scratch (redo it until it's okay)
  do{
    Create_GenZero();
    if(_ancestor.Check()==false) cout << "ancestor not valid... making a new one" << endl;
  }while(_ancestor.Check()==false);
  cout << "ancestor created!" << endl << endl;

  // computing lengths of the paths and sorting (from the minimum length to the maximum)
  _ancestor.ComputeAllLengths();

  _ancestor.PrintPopulation("ancestor.out");  // print to check

  _ancestor.Sort();

  // printing the average of the best half and the minimum length for Gen0
  PrintLength(0,_ancestor);


  // creating _n_gen new generations
  for(int num_gen=1; num_gen<=_n_gen; num_gen++){

    // creating the new generation (redo it until it's okay)
    do{
      Create_NewGen();
      if(_descendant.Check()==false) cout << "descendant not valid" << endl;
    }while(_descendant.Check()==false);
    cout << "Descendant " << num_gen << " created!" << endl;

    //_descendant.PrintPopulation("descendant.out");  // print to check

    // putting the current generation in the old one: it is the new ancestor for the next generation
    UpdateAncestor();

    // printing the average of the best half and the minimum length for the current generation
    PrintLength(num_gen,_ancestor);

  }
  

  return 0;
  

}













// Functions


void Input(){

    //reading input information

    ifstream ReadInput;
    ReadInput.open("input.in");
    ReadInput >> _ncities;  // # cities in a path
    ReadInput >> _build;    // how to build the cities (0=circle,1=square)
    ReadInput >> _npaths;   // # paths in a population
    ReadInput >> _n_gen;    // # generations (excluding Gen0)
    ReadInput >> _p;        // exponent for the selection operator
    ReadInput >> _p_swap;   // probability for swap
    ReadInput >> _p_perm;   // probability for permute
    ReadInput >> _p_shift;  // probability for shift
    ReadInput >> _p_invert; // probability for invert
    ReadInput >> _p_crossover;  // probability for crossover
    ReadInput.close();

    cout << "Travelling salesman problem:" << endl;
    cout << _ncities << " cities built ";
    if(_build) cout << "in a 1x1 square" << endl;
    else cout << "on a circumference of radius 1" << endl;
    cout << "Performing genetic algorithm with " << _n_gen << " generations" << endl;
    cout << "Each generation is made of " << _npaths << " paths" << endl << endl;
  
    //initializing the random generator
    _rnd.Initialize();


    //generating _ncities cities on a circumference if build=0 or in a square if build=1 printing the configuration of generated cities

    ofstream Conf;
    Conf.open("cities_conf.out");
    const int wd = 12;
    for(int i=0; i<_ncities; i++){
      city pluto;
      if(_build)
        pluto.Build_InSquare(i+1,_rnd);
      else
        pluto.Build_OnCircle(i+1,_rnd);
      Conf << pluto.GetIndex() << setw(wd) << pluto.GetX() << setw(wd) << pluto.GetY() << endl; 
    }

    Conf.close(); 


    return;

}


void Create_GenZero(){

  _ancestor.Initialize(_npaths,_ncities,"cities_conf.out",true,_rnd);   // random initialization: read cities from the file and then make a random number of swaps

  int num_cities, num_shifts, num_start;
  num_cities = _ncities-1;
  num_start = 1;
  double prob = 0.5;

  for(int i=0; i<_npaths; i++)  // now with probability 50% make some shifts of all the cities after the first one
    if(_ancestor.CallMutation(prob,_rnd)){
      num_shifts = int(_rnd.Rannyu(1,_ncities-1));
      _ancestor.Shift(i,num_cities,num_shifts,num_start);
    }

}


// print average length of best half and minimum length after sorting!

void PrintLength(int num_gen, population generation){    //num_gen is the number of generation
  ofstream LenOut;
  LenOut.open("length.out",ios::app);
  const int wd = 12;
  LenOut << setw(wd) << num_gen << setw(wd) << generation.GetBestHalf() << setw(wd) << generation.GetPath(0).GetLength() << endl;
  LenOut.close();   // printing for every row # generation - average length of best half - minimum length
}


void Create_NewGen(){

  _descendant.Initialize(_npaths,_ncities,"cities_conf.out",false,_rnd);    // just creating the vector of paths, not giving paths to it

  int i_select;   //selecting individuals for next generation
  for(int i=0; i<_npaths; i++){
    i_select = _ancestor.Select(_p,_rnd);
    _descendant.SetPath(i,_ancestor.GetPath(i_select));   // initializing the generation from the old one with the selected paths 
  }

  CrossoverRoutine();   // crossover
  MutationRoutine();    // mitations

  return;

}


void MutationRoutine(){

  for(int i=0; i<_npaths; i++){   // for every path
    // swap cities
    if(_descendant.CallMutation(_p_swap,_rnd)){   // calling the mutations according to their probability
      int a = int(_rnd.Rannyu(1,_ncities));
      int b = int(_rnd.Rannyu(1,_ncities));
      _descendant.Swap(i,a,b);
    }

    // shift cities
    if(_descendant.CallMutation(_p_shift,_rnd)){
      int start = int(_rnd.Rannyu(1,_ncities-3));
      int m = int(_rnd.Rannyu(3,_ncities-start));
      int n = int(_rnd.Rannyu(1,_ncities-2));
      _descendant.Shift(i,m,n,start);
    }

    // permute cities
    if(_descendant.CallMutation(_p_perm,_rnd)){
      int m = int(_rnd.Rannyu(2,(_ncities-1)/2));
      int a = int(_rnd.Rannyu(1,_ncities-2*m));
      int b = int(_rnd.Rannyu(a+m,_ncities-m));
      _descendant.Permute(i,m,a,b);
    }

    // invert cities
    if(_descendant.CallMutation(_p_invert,_rnd)){
      int m = int(_rnd.Rannyu(1,_ncities-1));
      int start = int(_rnd.Rannyu(1,_ncities-m));
      _descendant.Invert(i,m,start);
    }

  }

}


void CrossoverRoutine(){

  int num_cross = int(_rnd.Rannyu(1,_npaths/2));  // number of possible crossovers

  int mum_select, dad_select;
  for(int i=0; i<num_cross; i+=2){
    if(_descendant.CallMutation(_p_crossover,_rnd)){    // for every possible crossover call it according to its probability
      mum_select = _ancestor.Select(_p,_rnd); // selecting individuals from the ancestor
      dad_select = _ancestor.Select(_p,_rnd);
      _descendant.Crossover(_ancestor.GetPath(mum_select),_ancestor.GetPath(dad_select),i,i+1,_rnd);  // crossover creating two elemets of the descendant
    }
  }

}


void UpdateAncestor(){

  for(int i=0; i<_npaths; i++)    // reinizialize paths of the ancestor copying those of the descendant
    _ancestor.SetPath(i,_descendant.GetPath(i));

  _ancestor.ComputeAllLengths();    // compute length of all paths

  _ancestor.Sort();   // sort descendant

}