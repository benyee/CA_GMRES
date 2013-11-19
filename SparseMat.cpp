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

vector< vector<double> > SparseMat::convertToFullMatrix(unsigned int numRows, unsigned int numCols, bool rowFirst){
    
    vector< vector<double> > outvec;
    
    if(rowFirst){
        unsigned int A_counter = 0;
        for(unsigned int i = 0; i < numRows; i++){
            vector<double> row;
            for(unsigned int j = 0; j< numCols;j++){
                if(j == JA[A_counter] && A_counter < IA[i+1]){
                    row.push_back(A[A_counter]);
                    A_counter++;
                }else{
                    row.push_back(0);
                }
            }
            outvec.push_back(row);
        }
        return outvec;
    }
    
    for(unsigned int j = 0; j< numCols;j++){
        vector<double> emptyrow;
        outvec.push_back(emptyrow);
        for(unsigned int i = 0; i < numRows;i++){
            outvec[i].push_back(0);
        }
    }
    unsigned int A_counter;
    for(unsigned int i = 0; i < numRows; i++){
        A_counter = IA[i];
        while(A_counter < IA[i+1]){
            outvec[JA[i+1]][i] = A[A_counter];
            A_counter++;
        }
    }
    return outvec;
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

vector<double> SparseMat::smvp(const vector<double> &vec){
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


//Matrix powers:
vector< vector<double> > SparseMat::matrixPowers(const vector<double> &v_0, /*const vector<double> &th, */unsigned int s){
    unsigned int length = v_0.size();
    
    unsigned int start = clock();
    vector<vector<double> > V = Utilities::zeros(length,s);
    cout<<"This took "<<clock()-start<<endl;
    
    //which vector am I working on right now?
    unsigned int whichvector = 0;
    
    vector<unsigned int> index = Utilities::unsignedint_zeros(s);
    //How far along each vector has been computed
    //i.e., vector i has been computed completely up to (not including) entry index[i]
    //This is also the row of A that I shoudl be working on for that vector
    
    while(index[s-1]<length){
        //cout<<"Checking first condition..."<<endl;
        if(index[whichvector] >= length && whichvector < s){
            whichvector++;
            //cout<<"Increasing whichvector because I'm done with this vector"<<endl;
            continue;
        }
        //cout<<"Checking second condition..."<<endl;
        if(whichvector < s-1 && index[whichvector]>JA[ IA[ index[whichvector+1]+1 ]-1 ]){
            whichvector++;
            //cout<<"Increasing whichvector because I can move on!"<<endl;
            continue;
        }
        //cout<<"Checking third condition..."<<endl;
        if(whichvector > 0 && index[whichvector-1] <= JA[ IA[ index[whichvector]+1 ]-1 ]){
            whichvector--;
            //cout<<"Decreasing whichvector because I need to"<<endl;
            continue;
        }
        //cout<<"Computing for whichvector = "<<whichvector<<" with row = "<<index[whichvector]<<endl;
        
        
        if(whichvector){
            double local = 0;
            for(unsigned int j = IA[index[whichvector]]; j<IA[index[whichvector]+1];j++){
                local += A[j]*V[JA[j]][whichvector-1];
            }
            V[index[whichvector]][whichvector] = local;
        }else{
            double local = 0;
            for(unsigned int j = IA[index[whichvector]]; j<IA[index[whichvector]+1];j++){
                //cout<<"Doing Aval = "<<Aval[j]<<" times "<<v_0[JA[j]]<<endl;
                local += A[j]*v_0[JA[j]];
            }
            V[index[0]][0] = local;
        }
        index[whichvector]++;
    }
    
    cout<<"Checking took "<<checkTime<<endl;
    cout<<"Computation took "<<compTime<<endl;
    return V;
}

