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
#include <vector>

using namespace std;

#endif /* defined(____SparseMat__) */
class SparseMat{
public:
    SparseMat();
    ~SparseMat();
    
    vector<double> getA(){return A;}
    vector<int> getIA(){return IA;}
    vector<unsigned int> getJA(){return JA;}
    
    void print_matrix();
    
    
private:
    vector<double> A;
    vector<int> IA;
    vector<unsigned int> JA;
    
};