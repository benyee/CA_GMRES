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
#include <cmath>
#include <stdlib.h>

using namespace std;


class SparseVec{
public:
    SparseVec();
    SparseVec(SparseVec* input);
    SparseVec(vector<double> bIn, vector<unsigned int> IbIn){setb(bIn); setIb(IbIn);}
    ~SparseVec();
    
    vector<double> getb(){return b;}
    vector<unsigned int> getIb(){return Ib;}
    void setb(vector<double> bIn){b.clear(); b.insert(b.end(),bIn.begin(),bIn.end());}
    void setIb(vector<unsigned int> IbIn){Ib.clear(); Ib.insert(Ib.end(),IbIn.begin(),IbIn.end());}
    
    void readFullVector(string inputfile = "example2.txt",char delim = '\t');
    
    void print_vector();
    
    double dotProduct(SparseVec *vec2);
    SparseVec* axpy(double a = 1,SparseVec *y = NULL);
    
    double infNorm();
    double twoNorm();
    
private:
    vector <double> b;
    vector <unsigned int> Ib;
};

#endif /* defined(____SparseVec__) */