//
//  Utilities.h
//  
//
//  Created by Ben Yee on 11/7/13.
//
//

#ifndef ____Utilities__
#define ____Utilities__

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>


class Utilities{
public:
    Utilities();
    ~Utilities();
    
    vector<double> axpy(vector<double> x, vector<double> y = NULL, double a = 1);
    double dotProd(vector<double> x, vector<double> y);
};

#endif /* defined(____Utilities__) */
