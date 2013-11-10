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
#include "Utilities.h"

using namespace std;

int main ()
{
    cout << "Hello world!"<<endl;
    
    SparseMat *sample = new SparseMat();
    sample->readFullMatrix();
    sample->print_matrix();
    
    vector<double> samplevec = Utilities::readVectorFile();
    Utilities::printDVector(samplevec);
    
    vector<double> samplevec2 = Utilities::readVectorFile("example3.txt");
    Utilities::printDVector(samplevec2);
    cout<< "Dot product is: "<<Utilities::dotProd(samplevec, samplevec2)<<endl;
    
    cout<< "A times vector 1 gives..."<<endl;
    Utilities::printDVector(sample->smvp(samplevec));
    cout<< "A times vector 2 gives..."<<endl;
    Utilities::printDVector(sample->smvp(samplevec2));
    
    cout<<"5*v1 is..."<<endl;
    Utilities::printDVector(Utilities::axpy(samplevec,5.0));
    cout<<"v1+v2 is..."<<endl;
    Utilities::printDVector(Utilities::axpy(samplevec,samplevec2));
    cout<<"5*v1+v2 is..."<<endl;
    Utilities::printDVector(Utilities::axpy(samplevec,5,samplevec2));
    
    cout<<"Inf Norm: "<<Utilities::infNorm(samplevec)<<endl;
    cout<<"Two Norm: "<<Utilities::twoNorm(samplevec)<<endl;
    
    cout << "Goodbye world!"<<endl;
    return 0;
}
