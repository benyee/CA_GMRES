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

#include "SparseVec.h"
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
    
    cout << "Goodbye world!"<<endl;
    return 0;
}
