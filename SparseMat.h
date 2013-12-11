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
    
    SparseMat();
    ~SparseMat();
    
    vector<double> getA(){return A;}
    vector<unsigned int> getIA(){return IA;}
    vector<unsigned int> getJA(){return JA;}
    
    void readFullMatrix(string inputfile = "example1.txt",char delim = '\t');
    void readSparseMatrix(string inputfile = "example1.txt", bool matflag = false,char delim = '\t');
    
    vector< vector<double> > convertToFullMatrix(bool rowFirst = true){return convertToFullMatrix(IA.size()-1,IA.size()-1,rowFirst);} //Default numRows and numCols to the implied num. of rows given by the size of IA
    vector< vector<double> > convertToFullMatrix(unsigned int numRows, unsigned int numCols,bool rowFirst = true);
    //rowFirst determines whether the output vector has rows or columns as its outer structure
    //rowFirst = true means that [i][j] refers to the i-th row and j-th column
    
    vector<double> smvp(const vector<double> &vec);
    void smvp(const vector<double> &vec, vector<double> &outb);
    
    //Compute (A*v_0, A^2*v_0, ... A^s*v_0)
    void regMatrixPowers(vector<vector<double> > &V, const unsigned int ind[2]);
    
    //CA matrix powers attempt:
    void matrixPowersMapper();
    void matrixPowers(const vector<double> &v_0, /*const vector<double> &th, */vector<vector<double> > &V);
    //Takes in two inputs.  v_0 is the starting vector.  V is the matrix to store the matrix powers in.
    
    
    
    /*void matrixPowers_fixed(const vector<double> &v_0, double (&V)[Utilities::A_SIZE][Utilities::s]);
    void matrixPowers_fixednorm(const vector<double> &v_0, double V[Utilities::A_SIZE][Utilities::s]);*/
    
    //Classical GMRES algorithm:
    struct GMRES_sol classicalGMRES(const vector<double> &b, double tol = 1.0E-6, unsigned int max_it = 100);
    struct GMRES_sol classicalGMRES(const vector<double> &b, vector<double> x, double tol = 1.0E-6,unsigned int max_it = 100);
    
    //Communication avoiding GMRES:
    struct GMRES_sol ca_GMRES(const vector<double> &b, double tol = 1.0E-6, unsigned int max_it = Utilities::RESTART);
    struct GMRES_sol ca_GMRES(const vector<double> &b, vector<double> x, double tol = 1.0E-6,unsigned int max_it = Utilities::RESTART);
    
    void print_matrix();
    
private:
    vector<double> A;
    vector<unsigned int> IA;
    vector<unsigned int> JA;
    
    unsigned int map[Utilities::A_SIZE*Utilities::s][2];
    
};
#endif /* defined(____SparseMat__) */
