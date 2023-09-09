#include "classes.h"
#include "random.h"

#include <iostream>
#include <cmath>
#include <vector>

using namespace std;

//CITY

        city::city(){
        }

        city::~city(){
        }


        //methods
        void city::Build_OnCircle(int index, Random& rnd){
            m_index = index;
            double theta = rnd.Rannyu(0.,2*M_PI);
            m_x = cos(theta);
            m_y = sin(theta);
        }

        void city::Build_InSquare(int index, Random& rnd){  
            m_index = index;
            m_x = rnd.Rannyu();
            m_y = rnd.Rannyu();
        }

        void city::SetAll(int index, double x, double y){
            m_index = index;
            m_x = x;
            m_y = y;
        }

        void city::SetIndex(int index){
            m_index = index;
        }

        void city::SetX(double x){
            m_x = x;
        }

        void city::SetY(double y){
            m_y = y;
        }

        double city::GetX() const{    
            return m_x;
        }

        double city::GetY() const{
            return m_y;
        }

        int city::GetIndex() const{
            return m_index;
        }


//PATH

        path::path(){
        }

        path::~path(){
            m_cities.clear();
            m_indexes.clear();
        }


        //methods
        //initialization
        void path::Initialize(Random& rnd){
            m_indexes[0] = 1;
            int num = int(rnd.Rannyu(1.,10.));   //the number of swaps to be done (numero massimo 10 è messo a caso!!!)
            int n, m;
            for(int i=0; i<num; i++){
                n = int(rnd.Rannyu(1,m_ncities-1));
                m = int(rnd.Rannyu(1,m_ncities-1));
                SwapCities(n,m);
            }    
        }       //to generate the starting population I take ordered paths (1,2,3,...,N) and make a random number of swaps between two random cities, except for the first one

        void path::SetCities(int num_cities, string filename){
            m_ncities = num_cities;
            ifstream ReadCities;
            ReadCities.open(filename);
            int read_index;
            double read_x, read_y;
            for(int i=0; i<m_ncities; i++){
                ReadCities >> read_index >> read_x >> read_y;
                city pluto;
                m_cities.push_back(pluto);
                m_cities[i].SetAll(read_index,read_x,read_y);
                m_indexes.push_back(read_index);
            }
            ReadCities.close();
        }

        //swap two cities
        void path::SwapCities(int i, int j){
            //swapping indexes 
            int temp = m_indexes[i];
            m_indexes[i] = m_indexes[j];
            m_indexes[j] = temp;
            //swapping cities
            city a = m_cities[i];
            m_cities[i] = m_cities[j];
            m_cities[j] = a;
        }

        //shift m contiguous cities for n positions
        void path::ShiftCities(int m, int n, int i_start){
            int i_end = i_start + m - 1;
            vector<city> shifted;
            for(int num=1; num<=n; num++){
                shifted.push_back(m_cities[i_end]);
                for(int i=i_start, j=1; i<=i_end-1 and j<=m-1; i++, j++){    //shifting cities for 1 position in shifted
                    shifted.push_back(m_cities[i]);
                }
                for(int i=i_start, j=0; i<=i_end and j<=m-1; i++, j++){      //copying the new order in m_seq
                    m_cities[i].SetAll(shifted[j].GetIndex(),shifted[j].GetX(),shifted[j].GetY());
                    m_indexes[i] = m_cities[i].GetIndex();
                }
                shifted.clear();
            }   //doing the shift of one position n times

        }

        //permute two different groups of m contiguous cities
        void path::PermuteGroups(int m, int a, int b){
            //assuming that b>a+m and b+m<m_N
            for(int i=a, j=b; i<a+m and j<b+m; i++, j++)
                SwapCities(i,j);
        }

        //invert the order of a sequence of m contiguous cities except for the first one
        void path::InvertGroup(int m, int i_start){
            int i_end = i_start + m - 1; //the -1 is because i_end has to be the last one of the sequence, it is not out of it
            for(int i=0; i<m/2; i++)            // m/2 integer division!
                SwapCities(i_start+i,i_end-i);
        }

        //length
        void path::ComputeLength(){
            double L = 0.;
            double sqrd_x, sqrd_y;
            // distance (first-last city)
            for(int i=0; i<m_ncities-1; i++){
                sqrd_x = pow( m_cities[i].GetX() - m_cities[i+1].GetX() , 2 );
                sqrd_y = pow( m_cities[i].GetY() - m_cities[i+1].GetY() , 2 );
                L += sqrt( sqrd_x + sqrd_y );
            }    
            // distance between last and first (always index 0 in the matrix) city
            sqrd_x = pow( m_cities[m_ncities-1].GetX() - m_cities[0].GetX() , 2 );
            sqrd_y = pow( m_cities[m_ncities-1].GetY() - m_cities[0].GetY() , 2 );
            L += sqrt( sqrd_x + sqrd_y );
    
            m_L = L;

        }

        //check operator
        bool path::Check(){
            bool fine = true;
            vector<int> indexes = m_indexes;

            // sorting indexes from the minimum to the maximum, from position 1 to m_N-1 (0 is always the first city with index 1)
            for (int i=1; i<m_ncities-1; i++) {
                int minIndex = i;
                for (int j=i+1; j<m_ncities; j++) {
                    if (indexes[j] < indexes[minIndex])
                        minIndex = j;
                }
            int temp = indexes[i];
            indexes[i] = indexes[minIndex];
            indexes[minIndex] = temp;
            }

            for(int i=0; i<m_ncities; i++){   //I also check the position 0 to make sure it contains the first city
                if(indexes[i]!=i+1){
                    fine = false;
                    break;
                }
            }
            return fine;
        }       //returns true if the path is okay, false otherwise

        double path::GetLength() const{
            return m_L;
        }

        city path::GetCity(int i) const{
            return m_cities[i];
        }

        void path::SetCity(int i, city place){
            m_indexes[i] = place.GetIndex();
            m_cities[i] = place;
        }

        vector<int> path::GetIndexes() const{
            return m_indexes;
        }

        int path::GetIndex(int i) const{
            return m_indexes[i];
        }

        void path::PrintPath(string filename){
            ofstream PathOut;
            PathOut.open(filename);
            const int wd = 12;
            for(int i=0; i<m_ncities; i++)
                PathOut << m_cities[i].GetIndex() << setw(wd) << m_cities[i].GetX() << setw(wd) << m_cities[i].GetY() << endl;
            PathOut.close();
        }




//POPULATION

        population::population(){
        }

        population::~population(){
            m_paths.clear();
        }


        //methods
        void population::Initialize(int num_paths, int num_cities, string filename, bool random_initialize, Random& rnd){
            m_npaths = num_paths;
            m_ncities = num_cities;
            for(int i=0; i<m_npaths; i++){
                path pippo;
                m_paths.push_back(pippo);
                if(random_initialize){
                    m_paths[i].SetCities(num_cities,filename);
                    m_paths[i].Initialize(rnd);
                }
            }
        }

        //check the population
        bool population::Check(){
            bool fine = true;
            for(int i=0; i<m_npaths; i++){
                if(m_paths[i].Check()==false){
                    fine = false;
                    break;
                }
            }
            return fine;
        }       //returns true if every path is okay, false otherwise

        //swap two paths i and j
        void population::SwapPaths(int i, int j){
            path a = m_paths[i];
            m_paths[i] = m_paths[j];
            m_paths[j] = a;
        }

        //find path with minimum length within a range of indexes beg-end
        int population::MinLength(int beg, int end){
            vector<double> len;
            for(int i=0; i<m_npaths; i++){
                m_paths[i].ComputeLength();
                len.push_back( m_paths[i].GetLength() );
            }    
            double min = len[beg];
	        int min_index = beg; 
	        for(int i=beg; i<=end; i++){
		        if(len[i]<min){
			        min = len[i];
			        min_index = i;
		        }
	        }
	        return min_index;
        }

        //sort (from the shortest to the longest)
        void population::Sort(){
            int min_index;
	        for(int i=0; i<m_npaths-1; i++){
		        min_index = MinLength(i,m_npaths-1);
		        SwapPaths(i,min_index);
	        }
        }

        //compute the average length of the best half of the population (after sorting!)
        double population::GetBestHalf(){
            int half = m_npaths/2;
            double sum = 0.;
            for(int i=0; i<half; i++)
                sum += m_paths[i].GetLength();
            sum = sum/half;
            return sum;
        }

        //selection operator
        int population::Select(double p, Random& rnd){
            double r = rnd.Rannyu();
            int j = int( m_npaths * pow(r,p) ) + 1;
            return j;
        }       //returning the index of the selected path in the array m_paths

        //crossover operator: takes as input two paths (parents) and creates sons in this generation at index ison and idaughter
        void population::Crossover(path mum, path dad, int ison, int idaughter, Random& rnd){
            
            int cut_index = int( rnd.Rannyu(1,m_ncities-2) );   //I will cut right after the position given by the index (ex. if cut_index=2, 12345 becomes 123|45)
            for(int i=0; i<=cut_index; i++){
                m_paths[ison].SetCity(i,mum.GetCity(i));
                m_paths[idaughter].SetCity(i,dad.GetCity(i));
            }

            bool not_taken = true;
            int j = 0;
            for(int i=cut_index+1; i<m_ncities; i++){
                for(;j<m_ncities and not_taken; j++){
                    for(int k=0; k<=cut_index; k++){
                        if(m_paths[ison].GetIndex(k)==mum.GetIndex(j)){
                            not_taken = false;
                            break;      //dovrebbe essere inutile
                        }
                    }  
                }
                m_paths[ison].SetCity(i,mum.GetCity(j));
                not_taken = true;
            }

            j = 0;
            for(int i=cut_index+1; i<m_ncities; i++){
                for(;j<m_ncities and not_taken; j++){
                    for(int k=0; k<=cut_index; k++){
                        if(m_paths[idaughter].GetIndex(k)==dad.GetIndex(j)){
                            not_taken = false;
                            break;
                        }
                    }  
                }
                m_paths[idaughter].SetCity(i,dad.GetCity(j));
                not_taken = true;
            }

        }

        path population::GetPath(int i) const{
            return m_paths[i];
        }

        void population::SetPath(int i, path a){
            m_paths[i] = a;
        }

        void population::ComputeAllLengths(){
            for(int i=0; i<m_npaths; i++)
                m_paths[i].ComputeLength();
        }

        // methods for genetic mutations

        void population::Swap(int path_index, int i, int j){
            m_paths[path_index].SwapCities(i,j);
        }

        void population::Permute(int path_index, int m, int a, int b){
            m_paths[path_index].PermuteGroups(m,a,b);
        }

        void population::Invert(int path_index, int m, int i_start){
            m_paths[path_index].InvertGroup(m,i_start);
        }

        void population::Shift(int path_index, int m, int n, int i_start){
            m_paths[path_index].ShiftCities(m,n,i_start);
        }

        bool population::CallMutation(double prob, Random& rnd){
            double r = rnd.Rannyu();
            if(r <= prob)
                return true;
            else
                return false;
        }

        void population::PrintPopulation(string filename){

            ofstream pop;
            pop.open(filename);
            //pop << fixed << setprecision(6);
            for(int i=0; i<m_npaths; i++){
                pop << "Path " << i+1 << ":\t";
                for(int j=0; j<m_ncities; j++)
                    pop << m_paths[i].GetIndex(j) << " " << m_paths[i].GetCity(j).GetX() << " " << m_paths[i].GetCity(j).GetY() << "\t";
                pop << m_paths[i].GetLength() << endl;
            }

            pop.close();

        }

