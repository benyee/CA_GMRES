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
    example->readFullMatrix("thbyth.txt");
    
    vector<double> samplevec = Utilities::readVectorFile();
    Utilities::expandVec(samplevec,Utilities::A_SIZE);
    GMRES_sol sol = example->ca_GMRES(samplevec);
    
    Utilities::printDVector(sol.x);
    
    cout << "Goodbye world!"<<endl;
    return 0;
}
