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
    sample->print_matrix();
    
    cout << "Goodbye world!"<<endl;
    return 0;
}
