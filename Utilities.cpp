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

void Utilities::printFullMatrix(vector< vector<double> > mat, bool rowFirst){
    if(rowFirst){
        cout<<"---Matrix output begin--"<<endl;
        for(unsigned int i = 0; i < mat.size();i++){
            printDVector(mat[i]);
        }
        cout<<"---Matrix output end--"<<endl;
    }else{
        cout<<"---Matrix output begin--"<<endl;
        for(unsigned int i = 0; i < mat.back().size();i++){
            cout<<"[ ";
            for(unsigned int j = 0; j < mat.size();j++){
                cout<<mat[j][i]<<" ";
            }
            cout<<"] "<<endl;
        }
        cout<<"---Matrix output end--"<<endl;
    }
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
    unsigned int minSize = A.back().size();
    for(unsigned int i = 0;i<minSize;i++){
        out.push_back(0);
    }
    minSize = min(A.size(),x.size());
    for(unsigned int i = 0; i<minSize;i++){
        out = axpy(A[i],x[i],out);
    }
    return out;
}

vector< vector<double> > Utilities::transpose(vector< vector<double> > A){
    vector< vector<double> > At(A);
    
    for(unsigned int i = 0; i < A.size(); i++){
        for(unsigned int j = 0; j < A[i].size(); j++){
            At[i][j] = A[j][i];
        }
    }
    
    return At;
}

double Utilities::dotProd(vector<double> x, vector<double> y){
    double sum = 0;
    for(unsigned int i = 0; i<x.size();i++){
        sum += x[i]*y[i];
    }
    return sum;
}

vector<double> Utilities::backSub(vector< vector<double> > A, vector<double> b,bool rowFirst){
    if(!rowFirst){
        A = transpose(A);
    }
    
    vector<double> x(b);
    for(int i = A.size()-1; i>=0;i--){
        for(unsigned int j = i+1; j < A.size(); j++){
            x[i] = x[i] - A[i][j]*x[j];
        }
        x[i] = x[i]/A[i][i];
    }
    
    if(!rowFirst){
        A = transpose(A);
    }
    
    return x;
}

//NOTE THAT THIS RETURNS Q-transpose, NOT Q
pair<vector< vector<double> >, vector< vector<double> > > Utilities::mgs(vector< vector<double> > mat){
    vector< vector<double> > At = transpose(mat);
    
    vector< vector<double> > Q(At);
    vector< vector<double> > R;
    
    //Initialize R:
    unsigned int numrows = mat.size();
    unsigned int numcols = mat.back().size();
    for(unsigned int i = 0; i<numcols;i++){
        vector<double> temp;
        for(unsigned int j = 0; j<numcols; j++){
            temp.push_back(0);
        }
        R.push_back(temp);
    }
    
    for(unsigned int i = 0; i<numcols;i++){
        R[i][i] = twoNorm(At[i]);
        Q[i] = axpy(At[i],1.0/R[i][i]);
        for(unsigned int j = i+1; j<numcols;j++){
            R[i][j] = dotProd(Q[i],At[j]);
            At[j] = axpy(axpy(Q[i],R[i][j]),-1,At[j]);
        }
    }
    //Stuff goes here
    
    //NOTE THAT THIS RETURNS Q-transpose, NOT Q
    pair<vector< vector<double> >, vector< vector<double> > > output (Q,R);
    return output;
}


vector< vector<double> > Utilities::subMatrix(vector< vector<double> > A, pair<unsigned int,unsigned int> ind1, pair<unsigned int,unsigned int> ind2){
    vector< vector<double> > Aout;
    
    unsigned int i = ind1.first;
    while(i<A.size() && i < ind1.second){
        unsigned int j = ind2.first;
        vector<double> temp;
        while(j<A.size() && j < ind2.second){
            temp.push_back(A[i][j]);
            j++;
        }
        while(j<ind2.second){
            temp.push_back(0);
            j++;
        }
        Aout.push_back(temp);
        i++;
    }
    while(i<ind1.second){
        vector<double> temp;
        for(unsigned int j = ind1.first; j<ind2.second; j++){
            temp.push_back(0);
        }
        Aout.push_back(temp);
        i++;
    }
    return Aout;
}

vector<double> Utilities::leastSquaresNormal(vector< vector<double> > A, vector<double> y){
    vector<double> x;
    
    //Stuff goes here
    
    return x;
}

vector<double> Utilities::leastSquaresQR(vector< vector<double> > A, vector<double> y){
    pair<vector< vector<double> >, vector< vector<double> > >  QR = mgs(A);
    vector<double> x = matvec(QR.first,y);
    x = backSub(QR.second,x);
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