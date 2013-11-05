//
//  SparseVec.cpp
//  
//
//  Created by Ben Yee on 11/5/13.
//
//

#include "SparseVec.h"

SparseVec::SparseVec(){
}

SparseVec::SparseVec(SparseVec *input){
    setb(input->getb());
    setIb(input->getIb());
}

void SparseVec::readFullVector(string inputfile, char delim){
    ifstream inputFile;
    inputFile.open(inputfile.c_str());
    
    string line;
    
    unsigned int row = 0;
    while(getline(inputFile,line)){ //while not at the end of the file
        double value = atof(line.c_str());
        if (value != 0){
            b.push_back(value);
            Ib.push_back(row);
        }
        row++;
    }
    Ib.push_back(b.size());
}

void SparseVec::print_vector(){
    cout<<"b = [";
    for(unsigned int i = 0; i<b.size(); i++){
        cout<<b[i]<<" ";
    }
    cout<<"]"<<endl;
    
    cout<<"Ib = [";
    for(unsigned int i = 0; i<Ib.size(); i++){
        cout<<Ib[i]<<" ";
    }
    cout<<"]"<<endl;
}

double SparseVec::dotProduct(SparseVec *vec2){
    vector<double> b2 = vec2->getb();
    vector<unsigned int> Ib2 = vec2->getIb();
    
    unsigned int i = 0;
    unsigned int j = 0;
    
    double sum = 0;
    while( i < Ib.back() && j < Ib2.back() ){
        if(Ib[i] == Ib2[j]){
            sum += b[i]*b2[j];
            i++;
        }else if (Ib[i]< Ib2[j]){
            i++;
        }else{
            j++;
        }
    }
    
    return sum;
}

SparseVec* SparseVec::axpy(double a, SparseVec *y){
    vector<double> outb = getb();
    vector<unsigned int> outIb;
    SparseVec *out;
    if(a!=1){ //If not default a value, multiply x by a
        for(int i = 0; i<Ib.back(); i++){
            outb[i] = outb[i]*a;
        }
    }
    
    vector<double> outby;
    if(y != NULL){ //If we have a y argument
        vector<double> b2 = y->getb();
        vector<unsigned int> Ib2 = y->getIb();
        
        unsigned int i = 0;
        unsigned int j = 0;
        
        while( i < Ib.back() && j < Ib2.back() ){
            if(Ib[i] == Ib2[j]){
                outby.push_back(outb[i]+b2[j]);
                outIb.push_back(Ib[i]);
                i++; j++;
            }else if (Ib[i]< Ib2[j]){
                outby.push_back(outb[i]);
                outIb.push_back(Ib[i]);
                i++;
            }else{
                outby.push_back(b2[j]);
                outIb.push_back(Ib2[j]);
                j++;
            }
        }
        
        out = new SparseVec(outby, outIb);
        
    }else{ //Otherwise just return what we have now
        outIb = getIb();
        out = new SparseVec(outb, outIb);
    }
    
    return out;
}

double SparseVec::infNorm(){
    double max = abs(b[0]);
    
    for(unsigned int i=1; i<Ib.back(); i++){
        double val=abs(b[i]);
        if(val>max){
            max = val;
        }
    }
    return max;
}

double SparseVec::twoNorm(){
    return sqrt(dotProduct(this));
}
