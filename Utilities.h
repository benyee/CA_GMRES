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

using namespace std;

class Utilities{
public:
    static vector<double> EMPTY_DOUBLE_VECTOR;
    Utilities();
    ~Utilities();
    
    static const unsigned int s = 15; // How many vectors to compute at a time for matrix powers.  Also determines the width of your TSQR matrix (I think...)
    static const unsigned int RESTART = 60; //This doesn't do anything for now.
    static const unsigned int A_SIZE = 2500;  //This is the size of your matrix
    static const unsigned int BLOCK_SIZE2 = 7200/s; //Block size for TSQR + s.  This value should always be equal to BLOCK_SIZE + s
    static const unsigned int BLOCK_SIZE = BLOCK_SIZE2 - s; //Block size for TSQR
    
    static void printDVector(const vector<double> &vec); //Prints out a vector of doubles
    static void printFullMatrix(const vector< vector<double> > &mat, bool rowFirst = true); //Prints out a vector of double vectors
    
    //Make a vector of zeros
    static vector<double> zeros(unsigned int rows);
    static vector<unsigned int> unsignedint_zeros(unsigned int rows);
    //Make a matrix of zeros
    static vector< vector<double> > zeros(unsigned int rows, unsigned int cols);
    //Expand a vector vec until mat is of size rows
    static void expandVec(vector<double> &vec, unsigned int rows);
    //Expand a RECTANGULAR matrix mat until mat is of size rows,cols
    static void expandMat(vector< vector<double> > &mat, unsigned int rows, unsigned int cols);
    
    static vector<double> readVectorFile(string inputfile = "example2.txt",char delim = '\t');
    static vector<double> axpy(const vector<double> &x, const vector<double> &y);
    static vector<double> axpy(const vector<double> &x, double a);
    static vector<double> axpy(const vector<double> &x, double a, const vector<double> &y);
    static vector<double> matvec(const vector< vector<double> > &A, const vector<double> &x, bool rowFirst = true);
    static vector< vector<double> > transpose(const vector< vector<double> > &A);
    static double dotProd(vector<double> x, vector<double> y);
    
    //Get a piece of a matrix:
    static vector< vector<double> > subMatrix(const vector< vector<double> > &A, pair<unsigned int,unsigned int> ind1, pair<unsigned int,unsigned int> ind2);
    //Stack matrices vertically:
    static vector< vector<double> > stackMat(vector<vector<double> > A, vector<vector<double> > B, bool rowFirst = true);
    
    //Back-substitution for solving upper triangular matrix:
    static vector<double> backSub(const vector< vector<double> > &A,const vector<double> &b, bool rowFirst = true);
    
    //mat should be given as mat[row][col]
    //NOTE THAT THIS RETURNS Q-transpose, NOT Q
    static pair<vector< vector<double> 	>, vector< vector<double> > > mgs(const vector< vector<double> > &mat);
    static pair<vector< vector<double> 	>, vector< vector<double> > > mgs(const vector< vector<double> > &mat, unsigned int numcols, vector<vector<double> > R, vector<vector<double> > Q);
    
    
    static pair<vector< vector<double> 	>, vector< vector<double> > > mgs_col(const vector< vector<double> > &mat);
    static pair<vector< vector<double> 	>, vector< vector<double> > > mgs_col(const vector< vector<double> > &mat, unsigned int numcols, vector<vector<double> > R, vector<vector<double> > Q);
    
    static void mgs(const vector< vector<double> > &mat, vector<vector<double> > &R, vector<vector<double> > &Q);
    
    //Use the normal equations to solve Least Squares:  (Large condition #, may result in inaccuracy)
    static vector<double> leastSquaresNormal(vector< vector<double> > A, vector<double> y);
    //Use a QR (mgs) factorization to solve least squares: (not optimal, but more stable)
    static vector<double> leastSquaresQR(const vector< vector<double> > &A,const vector<double> &y);
    //Use tall-skinny QR to factorize a matrix
    static vector< vector<double> > tsQR(const vector< vector<double> > &A,unsigned int blksiz);
    
    //******************Fixed block size stuff***************************************
    static void mgs(double At[s][BLOCK_SIZE2], double R[s][s], double Q[s][BLOCK_SIZE2]);
    //static void mgs(vector<vector<double> > &At, unsigned int numcols, double Q[NUMCOLS][BLOCK_SIZE2], unsigned int ind);
//    static void mgs(const vector< vector<double> > &mat,  double R[s][s], double Q[s][A_SIZE]);
    static double twoNorm(double x[],unsigned int len = BLOCK_SIZE2);
    static double dotProd(double x[], double y[],unsigned int len = BLOCK_SIZE2);
    static void axpy(double x[], double y[],double out[],unsigned int len = BLOCK_SIZE2);
    static void axpy(double x[], double a, double out[],unsigned int len = BLOCK_SIZE2);
    static void axpy(double x[], double a, double y[],double out[],unsigned int len = BLOCK_SIZE2);
    static void RAtoAi(const vector<vector<double> >  &A, vector<vector<double> >  R, double Ai[s][BLOCK_SIZE2], double ind1);
    static void RAtoAi(const vector<vector<double> >  &A, double R[s][s], double Ai[s][BLOCK_SIZE2], double ind1);
    //Use TSQR where the array size is already determined
    static vector < vector <double> > tsQR_fixed(const vector< vector<double> > &A);
    //****************Column based fixed block size stuff************************
    static void RAtoAi_col(const vector<vector<double> >  &A, vector<vector<double> >  R, double Ai[s][BLOCK_SIZE2], double ind1);
    static void RAtoAi_col(const vector<vector<double> >  &A, double R[s][s], double Ai[s][BLOCK_SIZE2], double ind1);
    static vector < vector <double> > tsQR_col(const vector< vector<double> > &A);
    //************************************************************************************************

    
    
    static double twoNorm(const vector<double> &x);
    static double infNorm(const vector<double> &x);
};

#endif /* defined(____Utilities__) */
