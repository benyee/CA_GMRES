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

using namespace std;

class Utilities{
public:
    static vector<double> EMPTY_DOUBLE_VECTOR;
    Utilities();
    ~Utilities();
    
    static void printDVector(vector<double> vec); //Prints out a vector of doubles
    static void printFullMatrix(vector< vector<double> > mat); //Prints out a vector of double vectors
    
    static vector<double> readVectorFile(string inputfile = "example2.txt",char delim = '\t');
    static vector<double> axpy(vector<double> x, vector<double> y);
    static vector<double> axpy(vector<double> x, double a);
    static vector<double> axpy(vector<double> x, double a, vector<double> y);
    static vector<double> matvec(vector< vector<double> > A, vector<double> x, bool rowFirst = true);
    static double dotProd(vector<double> x, vector<double> y);
    
    //mat should be given as mat[row][col]
    static pair<vector< vector<double> 	>, vector< vector<double> > > mgs(vector< vector<double> > mat);
    
    //Use the normal equations to solve Least Squares:  (Large condition #, may result in inaccuracy)
    static vector<double> leastSquaresNormal(vector< vector<double> > A, vector<double> y);
    
    static double twoNorm(vector<double> x);
    static double infNorm(vector<double> x);
};

#endif /* defined(____Utilities__) */
