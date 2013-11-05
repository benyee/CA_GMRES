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

#include "SparseVec.h"

using namespace std;

class SparseMat{
public:
    SparseMat();
    ~SparseMat();
    
    vector<double> getA(){return A;}
    vector<unsigned int> getIA(){return IA;}
    vector<unsigned int> getJA(){return JA;}
    
    void readFullMatrix(string inputfile = "example1.txt",char delim = '\t');
    
    SparseVec* smvp(SparseVec *vec);
    
    void print_matrix();
    
private:
    vector<double> A;
    vector<unsigned int> IA;
    vector<unsigned int> JA;
    
};
#endif /* defined(____SparseMat__) */
