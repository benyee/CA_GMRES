//
//  Utilities.cpp
//  
//
//  Created by Ben Yee on 11/7/13.
//
//

#include "Utilities.h"

Utilities::Utilities(){
}
    
vector<double> Utilities::axpy(vector<double> x, vector<double> y, double a){
    vector<double> out(x);
    if (a != 1){
        for(unsigned int i = 0; i<x.size();i++){
            out[i] = a*out[i];
        }
    }
    if (y != NULL){
        for(unsigned int i = 0; i<x.size();i++){
            out[i] = out[i] + y[i];
        }
    }
    return out;
}

double Utilities::dotProd(vector<double> x, vector<double> y){
    double sum = 0;
    for(unsigned int i = 0; i<x.size();i++){
        sum += x[i]*y[i];
    }
    return sum;
}

double Utilities::twoNorm(vector<double> x){
    return sqrt(dotProd(x,x));
}