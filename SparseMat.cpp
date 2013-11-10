//
//  SparseMat.cpp
//  
//
//  Created by Ben Yee on 11/5/13.
//
//

#include "SparseMat.h"

SparseMat::SparseMat(){
}

void SparseMat::readFullMatrix(string inputfile, char delim){
    ifstream inputFile;
    inputFile.open(inputfile.c_str());
    
    string line;
    
    unsigned int counter = 0;
    unsigned int col = 0;
    while(getline(inputFile,line)){ //while not at the end of the file
        IA.push_back(counter);
        
        size_t found_old = 0;
        size_t found = line.find_first_of(delim); //Find index of first tab chararacter
        
        while (true)
        {
            double value = atof(line.substr(found_old,found).c_str());
            if(value != 0){
                A.push_back(value);
                JA.push_back(col);
                counter++;
            }
            col++;
            found_old = found; //Remember last index
            found=line.find_first_of(delim,found+1); //Find the index of the next tab character
            
            if(found==std::string::npos){ //If can't find next tab character, break
                break;
            }
            //cout<<"counter = "<<counter<<"; col = "<<col<<endl;
        }
        
        //Need this f
        double value = atof(line.substr(found_old,found).c_str());
        if(value != 0){
            A.push_back(value);
            JA.push_back(col);
            counter++;
        }
        col = 0;
    }
    
    IA.push_back(counter); //Last entry of IA is the number of nonzeros
}

void SparseMat::print_matrix(){
    cout<<"A = [";
    for(unsigned int i = 0; i<A.size(); i++){
        cout<<A[i]<<" ";
    }
    cout<<"]"<<endl;
    
    cout<<"IA = [";
    for(unsigned int i = 0; i<IA.size(); i++){
        cout<<IA[i]<<" ";
    }
    cout<<"]"<<endl;
    
    cout<<"JA = [";
    for(unsigned int i = 0; i<JA.size(); i++){
        cout<<JA[i]<<" ";
    }
    cout<<"]"<<endl;
}

/*vector<double> SparseMat::smvp(SparseVec *vec){
    vector<double> outb;
    vector<unsigned int> outIb;
    
    vector<double> b = vec->getb();
    vector<unsigned int> Ib = vec->getIb();
    
    unsigned int j; //Tracks the index of A vector
    for(unsigned int row = 0; row<IA.size()-1; row++){
        unsigned int i = 0; //Tracks the index of b vector
        j = IA[row];
        
        //Dot product:
        double sum = 0;
        while( i < Ib.back() && j < IA[row+1]){
            if(Ib[i] == JA[j]){
                sum += b[i]*A[j];
                i++;
            }else if (Ib[i]< JA[j]){
                i++;
            }else{
                j++;
            }
        }
        
        //Store if nonzero:
        if(sum !=0){
            outb.push_back(sum);
            outIb.push_back(row);
        }
    }
    
    outIb.push_back(outb.size()); //Add in the number of nonzero elements
    
    SparseVec *outvec = new SparseVec(outb,outIb);
    return outvec;
}*/

vector<double> SparseMat::smvp(vector<double> vec){
    vector<double> outb;
    
    unsigned int j; //Tracks the index of A vector
    for(unsigned int row = 0; row<IA.size()-1; row++){
        j = IA[row];
        double sum = 0;
        while(j < IA[row+1]){
            sum += A[j]*vec[JA[j]];
            j++;
        }
        outb.push_back(sum);
    }
    
    return outb;
}


