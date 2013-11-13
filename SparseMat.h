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

using namespace std;

class SparseMat{
public:
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
    
    vector<double> smvp(vector<double> vec);
    
    void print_matrix();
    
private:
    vector<double> A;
    vector<unsigned int> IA;
    vector<unsigned int> JA;
    
};
#endif /* defined(____SparseMat__) */
