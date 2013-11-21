// 
//  SparseMat.h
//  
//
//  Created by Ben Yee on 11/5/13.
//
//

#ifndef ____SparseMat__
#define ____SparseMat__

#include <iostream>
#include <fstream>
#include <cstddef>
#include <string>
#include <vector>
#include <stdlib.h>
#include "Utilities.h"
#include "GMRES_sol.h"


using namespace std;

class Utilities;

class SparseMat{
public:
    static const unsigned int MAX_IT = 100;
    
    SparseMat();
    ~SparseMat();
    
    vector<double> getA(){return A;}
    vector<unsigned int> getIA(){return IA;}
    vector<unsigned int> getJA(){return JA;}
    
    void readFullMatrix(string inputfile = "example1.txt",char delim = '\t');
    
    
    vector< vector<double> > convertToFullMatrix(bool rowFirst = true){return convertToFullMatrix(IA.size()-1,IA.size()-1,rowFirst);} //Default numRows and numCols to the implied num. of rows given by the size of IA
    vector< vector<double> > convertToFullMatrix(unsigned int numRows, unsigned int numCols,bool rowFirst = true);
    //rowFirst determines whether the output vector has rows or columns as its outer structure
    //rowFirst = true means that [i][j] refers to the i-th row and j-th column
    
    vector<double> smvp(const vector<double> &vec);
    
    //Compute (A*v_0, A^2*v_0, ... A^s*v_0)
    vector< vector<double> > matrixPowers(const vector<double> &v_0, /*const vector<double> &th, */unsigned int s);
    /*void matrixPowers_fixed(const vector<double> &v_0, double (&V)[Utilities::A_SIZE][Utilities::s]);
    void matrixPowers_fixednorm(const vector<double> &v_0, double V[Utilities::A_SIZE][Utilities::s]);*/
    
    //Classical GMRES algorithm:
    struct GMRES_sol classicalGMRES(const vector<double> &b, double tol = 1.0E-6, unsigned int max_it = MAX_IT);
    struct GMRES_sol classicalGMRES(const vector<double> &b, vector<double> x, double tol = 1.0E-6, unsigned int max_it = MAX_IT);
    
    void print_matrix();
    
private:
    vector<double> A;
    vector<unsigned int> IA;
    vector<unsigned int> JA;
    
};
#endif /* defined(____SparseMat__) */
