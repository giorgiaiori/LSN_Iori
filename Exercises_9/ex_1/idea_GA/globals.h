#ifndef __globals_h__
#define __globals_h__

#include "random.h"
#include "classes.h"

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <iomanip>

using namespace std;

// global variables

Random _rnd; // random generator
int _ncities;  //number of cities
int _npaths, _n_gen; //number of paths in a population, maximum number of generations
int _build;  //0 if cities are on a circumference, 1 if they are in a square
double _p;  //parameter for the selection operator
double _p_swap, _p_perm, _p_shift, _p_invert, _p_crossover;  //probabilities for calling genetic mutations

// population: first generation and derived generation
population _ancestor;
population _descendant;


// global functions

void Input(void);
void Create_GenZero(void);
void PrintLength(int,population);
void Create_NewGen(void);
void MutationRoutine(void);
void CrossoverRoutine(void);
void UpdateAncestor(void);





#endif
