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
    
    SparseMat *example = new SparseMat();
    vector<double> samplevec;
    if(Utilities::A_SIZE==16){
        example->readFullMatrix("thbyth.txt");
        samplevec= Utilities::readVectorFile();
        Utilities::expandVec(samplevec,Utilities::A_SIZE);
    }else if(Utilities::A_SIZE<=2500){
        example->readFullMatrix("example_matpow.txt");
        samplevec= Utilities::readVectorFile("example_matpowvec.txt");
    }
    
    GMRES_sol sol = example->ca_GMRES(samplevec);
    
    if(Utilities::A_SIZE<20){
        Utilities::printDVector(sol.x);
    }
    
    cout << "Goodbye world!"<<endl;
    return 0;
}
