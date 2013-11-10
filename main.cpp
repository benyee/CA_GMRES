//
//  main.cpp
//
//
//  Created by Ben Yee on 9/20/13.
//
//

#include <iostream>
#include <string>
#include <vector>

#include "SparseMat.h"

using namespace std;

int main ()
{
    cout << "Hello world!"<<endl;
    
    SparseMat *sample = new SparseMat();
    sample->readFullMatrix();
    sample->print_matrix();
    
    SparseVec *samplevec = new SparseVec();
    samplevec->readFullVector();
    samplevec->print_vector();
    
    SparseVec *samplevec2 = new SparseVec();
    samplevec2->readFullVector("example3.txt",'\t');
    samplevec2->print_vector();
    cout<< "Dot product is: "<<samplevec->dotProduct(samplevec2)<<endl;
    
    cout<< "A times vector 1 gives..."<<endl;
    (sample->smvp(samplevec))->print_vector();
    cout<< "A times vector 2 gives..."<<endl;
    (sample->smvp(samplevec2 ))->print_vector();
    
    cout<<"5*v1 is..."<<endl;
    (samplevec->axpy(5.0))->print_vector();
    cout<<"v1+v2 is..."<<endl;
    (samplevec->axpy(1,samplevec2))->print_vector();
    cout<<"5*v1+v2 is..."<<endl;
    (samplevec->axpy(5.0,samplevec2))->print_vector();
    
    cout<<"Inf Norm: "<<samplevec->infNorm()<<endl;
    cout<<"Two Norm: "<<samplevec->twoNorm()<<endl;
    
    cout << "Goodbye world!"<<endl;
    return 0;
}
