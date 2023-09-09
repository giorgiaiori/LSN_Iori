#ifndef __classes_h__
#define __classes_h__

#include "random.h"

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <iomanip>

using namespace std;


// classes

class city{

public:

    //constructor
    city();

    //destructor
    ~city();

    //methods
    void Build_OnCircle(int index, Random& rnd);
    void Build_InSquare(int index, Random& rnd);
    void SetAll(int index, double x, double y);
    void SetIndex(int index);
    void SetX(double x);
    void SetY(double y);
    double GetX() const;
    double GetY() const;
    int GetIndex() const;
        

private:

    int m_index;        //index defining the city
    double m_x, m_y;    //cartesian coordinates of the city

};

    

class path{

public:

	//constructor
    path();
	path(int N);
	
	//destructor
	~path();

	//methods
    void Initialize(Random& rnd);
    void SetCities(int num_cities, string filename);
    void SwapCities(int i, int j);
    void ShiftCities(int m, int n, int i_start);
    void PermuteGroups(int m, int a, int b);
    void InvertGroup(int m, int i_start);
	void ComputeLength();
    bool Check();
    double GetLength() const;
    city GetCity(int i) const;
    void SetCity(int i, city place);
    vector<int> GetIndexes() const;
    int GetIndex(int i) const;
    void PrintPath(string filename);
	

private:
    int m_ncities;                    //number of cities in the path
    double m_L;                 //length of the path (considering that I come back to the starting city)
    vector<int> m_indexes;      //array of the indexes defining the cities (the order defines the path)
    vector<city> m_cities;               //array of the cities (the order defines the path)

};




class population{

public:

	//constructor
    population();
	
	//destructor
	~population();

	//methods
    void Initialize(int num_paths, int num_cities, string filename, bool random_initialize, Random& rnd);
    bool Check();
    void SwapPaths(int i, int j);
    int MinLength(int beg, int end);
    void Sort();
    double GetBestHalf();
    int Select(double p, Random& rnd);
    void Crossover(path mum, path dad, int ison, int idaughter, Random& rnd);
    path GetPath(int i) const;
    void SetPath(int i, path a);
    void ComputeAllLengths();
    void Swap(int path_index, int i, int j);
    void Permute(int path_index, int m, int a, int b);
    void Invert(int path_index, int m, int i_start);
    void Shift(int path_index, int m, int n, int i_start);
    bool CallMutation(double prob, Random& rnd);
    void PrintPopulation(string filename);
	

private:
    int m_npaths;            //number of paths in the population
    int m_ncities;            //number of cities in each path
    vector<path> m_paths;     //array of paths composing the population

};







// functions in main

void Input(void);





#endif
