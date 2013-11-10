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

void Utilities::printDVector(vector<double> vec){
    cout<<"[ ";
    for(unsigned int i = 0; i < vec.size(); i++){
        cout<<vec[i]<<" ";
    }
    cout<<"]"<<endl;
}

void Utilities::printFullMatrix(vector< vector<double> > mat){
    cout<<"---Matrix output begin--"<<endl;
    for(unsigned int i = 0; i < mat.size();i++){
        printDVector(mat[i]);
    }
    cout<<"---Matrix output end--"<<endl;
}

vector<double> Utilities::readVectorFile(string inputfile,char delim){
    ifstream inputFile;
    inputFile.open(inputfile.c_str());
    vector<double> out;
    
    string line;
    
    while(getline(inputFile,line)){ //while not at the end of the file
        double value = atof(line.c_str());
        out.push_back(value);
    }
    
    return out;
}

vector<double> Utilities::axpy(vector<double> x, vector<double> y){
    return axpy(x,1,y);
}
vector<double> Utilities::axpy(vector<double> x, double a){
    vector<double> y;
    return axpy(x,a,y);
}
vector<double> Utilities::axpy(vector<double> x, double a, vector<double> y){
    vector<double> out(x);
    if (a != 1){
        for(unsigned int i = 0; i<x.size();i++){
            out[i] = a*out[i];
        }
    }
    unsigned int minsize = min(x.size(),y.size());
    for(unsigned int i = 0; i<minsize;i++){
        out[i] = out[i] + y[i];
    }
    return out;
}

vector<double> Utilities::matvec(vector< vector<double> > A, vector<double> x, bool rowFirst){
    vector<double> out;
    if(rowFirst){
        unsigned int minSize = A.size();
        for(unsigned int i = 0;i<minSize;i++){
            out.push_back(dotProd(A[i],x));
        }
        return out;
    }
    unsigned int minSize = A[0].size();
    for(unsigned int i = 0;i<minSize;i++){
        out.push_back(0);
    }
    minSize = min(A.size(),x.size());
    for(unsigned int i = 0; i<minSize;i++){
        out = axpy(A[i],x[i],out);
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

pair<vector< vector<double> >, vector< vector<double> > > Utilities::mgs(vector< vector<double> > mat){
    vector< vector<double> > Q;
    vector< vector<double> > R;
    
    //Stuff goes here
    
    pair<vector< vector<double> >, vector< vector<double> > > output (Q,R);
    return output;
}

vector<double> Utilities::leastSquaresNormal(vector< vector<double> > A, vector<double> y){
    vector<double> x;
    
    //stuff goes here
    
    return x;
}

double Utilities::infNorm(vector<double> x){
    double max = abs(x[0]);
    
    for(unsigned int i=1; i<x.size(); i++){
        double val=abs(x[i]);
        if(val>max){
            max = val;
        }
    }
    return max;
}

double Utilities::twoNorm(vector<double> x){
    return sqrt(dotProd(x,x));
}