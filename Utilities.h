//
//  Utilities.h
//  
//
//  Created by Ben Yee on 11/7/13.
//
//

#ifndef ____Utilities__
#define ____Utilities__

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <stdlib.h>
#include "GMRES_sol.h"
#include "SparseMat.h"

using namespace std;

class Utilities{
public:
    static vector<double> EMPTY_DOUBLE_VECTOR;
    Utilities();
    ~Utilities();
    
    static const unsigned int MAX_IT = 100;
    
    static void printDVector(vector<double> vec); //Prints out a vector of doubles
    static void printFullMatrix(vector< vector<double> > mat, bool rowFirst = true); //Prints out a vector of double vectors
    
    //Make a vector of zeros
    static vector<double> zeros(unsigned int rows);
    //Make a matrix of zeros
    static vector< vector<double> > zeros(unsigned int rows, unsigned int cols);
    //Expand a vector vec until mat is of size rows
    static vector<double>  expandVec(vector<double> vec, unsigned int rows);
    //Expand a RECTANGULAR matrix mat until mat is of size rows,cols
    static vector< vector<double> > expandMat(vector< vector<double> > mat, unsigned int rows, unsigned int cols);
    
    static vector<double> readVectorFile(string inputfile = "example2.txt",char delim = '\t');
    static vector<double> axpy(vector<double> x, vector<double> y);
    static vector<double> axpy(vector<double> x, double a);
    static vector<double> axpy(vector<double> x, double a, vector<double> y);
    static vector<double> matvec(vector< vector<double> > A, vector<double> x, bool rowFirst = true);
    static vector< vector<double> > transpose(vector< vector<double> > A);
    static double dotProd(vector<double> x, vector<double> y);
    
    //Get a piece of a matrix:
    static vector< vector<double> > subMatrix(vector< vector<double> > A, pair<unsigned int,unsigned int> ind1, pair<unsigned int,unsigned int> ind2);
    //Stack matrices vertically:
    static vector< vector<double> > stackMat(vector<vector<double> > A, vector<vector<double> > B, bool rowFirst = true);
    
    //Back-substitution for solving upper triangular matrix:
    static vector<double> backSub(vector< vector<double> > A, vector<double> b, bool rowFirst = true);
    
    //mat should be given as mat[row][col]
    //NOTE THAT THIS RETURNS Q-transpose, NOT Q
    static pair<vector< vector<double> 	>, vector< vector<double> > > mgs(vector< vector<double> > mat);
    static pair<vector< vector<double> 	>, vector< vector<double> > > mgs(vector< vector<double> > mat, unsigned int numcols, vector<vector<double> > R, vector<vector<double> > Q);
    
    //Use the normal equations to solve Least Squares:  (Large condition #, may result in inaccuracy)
    static vector<double> leastSquaresNormal(vector< vector<double> > A, vector<double> y);
    //Use a QR (mgs) factorization to solve least squares: (not optimal, but more stable)
    static vector<double> leastSquaresQR(vector< vector<double> > A, vector<double> y);
    //Use tall-skinny QR to factorize a matrix
    static vector < vector <double> > tsQR(vector< vector<double> > A,unsigned int blksiz);
    
    
    //******************Fixed block size stuff***************************************
    static const unsigned int NUMCOLS = 15;
    static const unsigned int BLOCK_SIZE2 = 7200/NUMCOLS;
    static const unsigned int BLOCK_SIZE = BLOCK_SIZE2 - NUMCOLS;
    static void mgs(double At[NUMCOLS][BLOCK_SIZE2], unsigned int numcols, double R[NUMCOLS][NUMCOLS], double Q[NUMCOLS][BLOCK_SIZE2]);
    static double twoNorm(double x[BLOCK_SIZE2]);
    static double dotProd(double x[BLOCK_SIZE2], double y[BLOCK_SIZE2]);
    static void axpy(double x[BLOCK_SIZE2], double y[BLOCK_SIZE2],double out[BLOCK_SIZE2]);
    static void axpy(double x[BLOCK_SIZE2], double a, double out[BLOCK_SIZE2]);
    static void axpy(double x[BLOCK_SIZE2], double a, double y[BLOCK_SIZE2],double out[BLOCK_SIZE2]);
    static void RAtoAi(vector<vector<double> >  A, vector<vector<double> >  R, double Ai[NUMCOLS][BLOCK_SIZE2], double ind1);
    static void RAtoAi(vector<vector<double> >  A, double R[NUMCOLS][NUMCOLS], double Ai[NUMCOLS][BLOCK_SIZE2], double ind1);
    //Use TSQR where the array size is already determined
    static vector < vector <double> > tsQR_fixed(vector< vector<double> > A);
    //************************************************************************************************
    
    
    //Classical GMRES algorithm:
    static struct GMRES_sol classicalGMRES(SparseMat* A, vector<double> b, double tol = 1.0E-6, unsigned int max_it = MAX_IT);
    static struct GMRES_sol classicalGMRES(SparseMat* A, vector<double> b, vector<double> x, double tol = 1.0E-6, unsigned int max_it = MAX_IT);
    
    static double twoNorm(vector<double> x);
    static double infNorm(vector<double> x);
};

#endif /* defined(____Utilities__) */
