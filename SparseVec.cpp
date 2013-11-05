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
        }else if (Ib[i]< Ib[j]){
            i++;
        }else{
            j++;
        }
    }
    
    return sum;
}
