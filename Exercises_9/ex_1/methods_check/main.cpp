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


  //reading input informations

    ifstream ReadInput;
    ReadInput.open("input.in");
    ReadInput >> _ncities;
    ReadInput >> _build;
    ReadInput >> _npaths;
    ReadInput >> _n_gen;
    ReadInput >> _p;
    ReadInput >> _p_swap;
    ReadInput >> _p_perm;
    ReadInput >> _p_shift;
    ReadInput >> _p_invert;
    ReadInput >> _p_crossover;
    ReadInput.close();

    cout << "Travelling salesman problem:" << endl;
    cout << _ncities << " cities built ";
    if(_build) cout << "in a 1x1 square" << endl;
    else cout << "on a circumference of radius 1" << endl;
    cout << "Performing genetic algorithm with " << _n_gen << " generations" << endl;
  
    //initializing the random generator
    _rnd.Initialize();


    //generating N cities on a circumference if build=0 or in a square if build=1 printing the configuration of generated cities

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

    
  // check if the random initialization of the path and the check method works

  path prova;
  prova.SetCities(_ncities,"cities_conf.out"); // initializing the cities (indexes are ordered in the file)
  prova.Initialize(_rnd);
  prova.PrintPath("path.out");
  if(prova.Check())
    cout << "first path is fine" << endl;
  else
    cout << "path not valid" << endl;

  // check if the second path is different from the first one
  path prova2;
  prova2.SetCities(_ncities,"cities_conf.out"); // initializing the cities (indexes are ordered in the file)
  prova2.Initialize(_rnd);
  prova2.PrintPath("path2.out");
  if(prova2.Check())
    cout << "second path is fine" << endl;
  else
    cout << "path not valid" << endl;


  // check my methods for the mutations

  prova.InvertGroup(10,5);  // invert 10 cities starting from index 5 in the vector
  prova.ShiftCities(3,2,1); // shift 3 cities of 2 positions starting from index 1 in the vector
  prova.PermuteGroups(4,3,12);  // permute two groups of 4 cities that start from index 3 and 12 in the vector
  prova.PrintPath("path_mutation.out");


  // check SetCity and the check method for wrong paths

  prova.SetCity(4,prova.GetCity(3));
  prova.PrintPath("path4.out");
  if(prova.Check())
    cout << "third path is fine" << endl;
  else
    cout << "path not valid" << endl;


  // build a path of three cities and check the computation of the length

  city place_A;
  place_A.SetIndex(35);
  place_A.SetX(2.0);
  place_A.SetY(1.0);

  city place_B;
  place_B.SetIndex(36);
  place_B.SetX(0.0);
  place_B.SetY(-2.0);

  city place_C;
  place_C.SetIndex(37);
  place_C.SetX(-1.0);
  place_C.SetY(1.0);

  ofstream NewPath;
  NewPath.open("new.out",ios::app);
  NewPath << place_A.GetIndex() << setw(wd) << place_A.GetX() << setw(wd) << place_A.GetY() << endl; 
  NewPath << place_B.GetIndex() << setw(wd) << place_B.GetX() << setw(wd) << place_B.GetY() << endl; 
  NewPath << place_C.GetIndex() << setw(wd) << place_C.GetX() << setw(wd) << place_C.GetY() << endl; 
  NewPath.close(); 

  path prova3;
  prova3.SetCities(3,"new.out");
  prova3.ComputeLength();
  NewPath.open("new.out",ios::app);
  NewPath << "length of the path " << prova3.GetLength() << endl;
  NewPath.close();



  // check population methods: random initialization, check, computation of lengths, minimum length, sorting, average of the best half

  population ancestor;
  ancestor.Initialize(_npaths,_ncities,"cities_conf.out",true,_rnd);

  if(ancestor.Check())
    cout << "first population is fine" << endl;
  else
    cout << "population not valid" << endl;
  
  ancestor.ComputeAllLengths();

  //printing the starting population 
  ofstream pop;
  pop.open("check_ancestor.out");
  pop << fixed << setprecision(6);
  for(int i=0; i<_npaths; i++){
    pop << "Path " << i+1 << ":\t";
    for(int j=0; j<_ncities; j++)
      pop << ancestor.GetPath(i).GetIndex(j) << " " << ancestor.GetPath(i).GetCity(j).GetX() << " " << ancestor.GetPath(i).GetCity(j).GetY() << "\t";
    //ancestor.GetPath(i).ComputeLength();
    pop << ancestor.GetPath(i).GetLength() << endl;
  }
  pop << endl << "minimum length for path " << ancestor.MinLength(0,_npaths-1)+1 << endl;
  pop << "minimum length between path 5 and 9 for path " << ancestor.MinLength(4,8)+1 << endl;

  pop << endl << "sorting population..." << endl;
  ancestor.Sort();
  for(int i=0; i<_npaths; i++){
    pop << "Path " << i+1 << ":\t";
    for(int j=0; j<_ncities; j++)
      pop << ancestor.GetPath(i).GetIndex(j) << " " << ancestor.GetPath(i).GetCity(j).GetX() << " " << ancestor.GetPath(i).GetCity(j).GetY() << "\t";
    //ancestor.GetPath(i).ComputeLength();
    pop << ancestor.GetPath(i).GetLength() << endl;
  }
  pop << endl << "minimum length for path " << ancestor.MinLength(0,_npaths-1)+1 << endl;
  pop << "minimum length between path 5 and 9 for path " << ancestor.MinLength(4,8)+1 << endl;
  pop << "average length of best half " << ancestor.GetBestHalf() << endl;

  // check selection operator
  pop << "selecting some individuals... " << endl;
  pop << "p = " << _p << " : ";
  for(int i=0; i<20; i++)
    pop << ancestor.Select(_p,_rnd) << " ";
  pop << endl;
  pop << "p = " << 0.2 << " : ";
  for(int i=0; i<20; i++)
    pop << ancestor.Select(0.2,_rnd) << " ";
  pop << endl;
  pop << "p = " << 0.8 << " : ";
  for(int i=0; i<20; i++)
    pop << ancestor.Select(0.8,_rnd) << " ";
  pop << endl;
  pop << "p = " << 1.4 << " : ";
  for(int i=0; i<20; i++)
    pop << ancestor.Select(1.4,_rnd) << " ";
  pop << endl;

  // check mutation operators
  pop << endl << "inverting a group of cities in a path..." << endl;
  ancestor.Invert(9,3,1);
  for(int i=0; i<_npaths; i++){
    pop << "Path " << i+1 << ":\t";
    for(int j=0; j<_ncities; j++)
      pop << ancestor.GetPath(i).GetIndex(j) << " " << ancestor.GetPath(i).GetCity(j).GetX() << " " << ancestor.GetPath(i).GetCity(j).GetY() << "\t";
    //ancestor.GetPath(i).ComputeLength();
    //pop << ancestor.GetPath(i).GetLength() << endl;
    pop << endl;
  }

  pop << endl << "permuting two groups of cities in a path..." << endl;
  ancestor.Permute(8,3,4,9);
  for(int i=0; i<_npaths; i++){
    pop << "Path " << i+1 << ":\t";
    for(int j=0; j<_ncities; j++)
      pop << ancestor.GetPath(i).GetIndex(j) << " " << ancestor.GetPath(i).GetCity(j).GetX() << " " << ancestor.GetPath(i).GetCity(j).GetY() << "\t";
    //ancestor.GetPath(i).ComputeLength();
    //pop << ancestor.GetPath(i).GetLength() << endl;
    pop << endl;
  }

  pop << endl << "shifting cities in a path..." << endl;
  ancestor.Shift(7,4,3,2);
  for(int i=0; i<_npaths; i++){
    pop << "Path " << i+1 << ":\t";
    for(int j=0; j<_ncities; j++)
      pop << ancestor.GetPath(i).GetIndex(j) << " " << ancestor.GetPath(i).GetCity(j).GetX() << " " << ancestor.GetPath(i).GetCity(j).GetY() << "\t";
    //ancestor.GetPath(i).ComputeLength();
    //pop << ancestor.GetPath(i).GetLength() << endl;
    pop << endl;
  }

  pop << endl << "swapping two cities in a path..." << endl;
  ancestor.Swap(1,0,33);
  for(int i=0; i<_npaths; i++){
    pop << "Path " << i+1 << ":\t";
    for(int j=0; j<_ncities; j++)
      pop << ancestor.GetPath(i).GetIndex(j) << " " << ancestor.GetPath(i).GetCity(j).GetX() << " " << ancestor.GetPath(i).GetCity(j).GetY() << "\t";
    //ancestor.GetPath(i).ComputeLength();
    //pop << ancestor.GetPath(i).GetLength() << endl;
    pop << endl;
  }
  if(ancestor.Check())
    cout << "population is fine" << endl;
  else
    cout << "population not valid" << endl;

  pop << endl << "making the population valid again..." << endl;
  path percorso;
  percorso.SetCities(_ncities,"cities_conf.out");
  percorso = ancestor.GetPath(0);
  ancestor.SetPath(1,percorso);
  for(int i=0; i<_npaths; i++){
    pop << "Path " << i+1 << ":\t";
    for(int j=0; j<_ncities; j++)
      pop << ancestor.GetPath(i).GetIndex(j) << " " << ancestor.GetPath(i).GetCity(j).GetX() << " " << ancestor.GetPath(i).GetCity(j).GetY() << "\t";
    //ancestor.GetPath(i).ComputeLength();
    //pop << ancestor.GetPath(i).GetLength() << endl;
    pop << endl;
  }
  if(ancestor.Check())
    cout << "population is fine" << endl;
  else
    cout << "population not valid" << endl;

  // check calling mutation depending on the probability
  pop << endl << "calling a mutation on a path..." << endl;
  double accepted = 0.;
  double all = 0.;
  for(int i=0; i<100; i++){
    if(ancestor.CallMutation(_p_invert,_rnd))
      accepted += 1.0;
    all += 1.0;
  }
  pop << "should be " << _p_invert << " and I have " << accepted/all << endl;

  pop.close();

  /*path mamma = path(_ncities);
  path papa = path(_ncities);
  mamma = _ancestor.GetPath(8);
  papa = _ancestor.GetPath(9);
  _ancestor.Crossover(mamma,papa,0,1,_rnd);
  _ancestor.PrintPopulation("crossover.out");*/     // crossover operator doesn't change anything in the population

  return 0;

}

