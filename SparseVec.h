//
//  SparseVec.h
//  
//
//  Created by Ben Yee on 11/5/13.
//
//

#ifndef ____SparseVec__
#define ____SparseVec__

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using namespace std;

#endif /* defined(____SparseVec__) */

class SparseVec{
public:
    SparseVec();
    ~SparseVec();
    
    void readFullVector(string inputfile = "example2.txt",char delim = '\t');
    
    void print_vector();
    
private:
    vector <double> b;
    vector <unsigned int> Ib;
};