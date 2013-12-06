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
#include "Utilities.h"

using namespace std;

int main ()
{
    cout << "Hello world!"<<endl;
    
    //To change the matrix being run, do three things:
    //  (1) Change the input to example->readFullMatrix(...);
    //  (2) Adjust s, BLOCK_SIZE2, and A_SIZE appropriately in Utilities.h
    //          s is the # of Krylov vectors (goes up to A^{s-1}*v)
    //          BLOCK_SIZE2 = BLOCK_SIZE + s = block size for TSQR + s
    //          A_SIZE = size of A (which should be square)
    //  (3) Make sure you read in a vector of appropriate size. Use Utilities::expandVec if necessary to pad the vector with zeros.
    SparseMat *example = new SparseMat();
    vector<double> samplevec;
    cout<<"Reading in A"<<endl;
    example->readSparseMatrix("A.txt",true);
    if(Utilities::A_SIZE<20){
        example->print_matrix();
    }
    cout<<"Reading in b"<<endl;
    samplevec= Utilities::readVectorFile("b.txt");
    
    cout<<"Mapping A.."<<endl;
    example->matrixPowersMapper();
    vector<vector<double> > V_CA(Utilities::s,vector<double>(Utilities::A_SIZE,0));
    //CA matrix powers:
    cout<<"Starting CA matrix powers..."<<endl;
    int start = clock();
    example->matrixPowers(samplevec,V_CA);
    cout<<"CA matrix powers took"<<clock()-start<<endl;
    
    vector<vector<double> > V_reg(Utilities::sp1,vector<double>(Utilities::A_SIZE,0));
    V_reg[0] = samplevec;
    unsigned int ind[2] = {1,Utilities::sp1};
    //Regular matrix powers:
    cout<<"Starting regular matrix powers..."<<endl;
    start = clock();
    example->regMatrixPowers(V_reg,ind);
    cout<<"Regular matrix powers took"<<clock()-start<<endl;
    
    cout<<endl;
    cout<<endl;
    cout<<"Running tsqr..."<<endl;
    start = clock();
    Utilities::tsQR_col(V_reg);
    cout<<"tsqr took "<<clock()-start<<endl;
    cout<<"Running mgs..."<<endl;
    start = clock();
    Utilities::mgs_col(V_reg);
    cout<<"mgs took "<<clock()-start<<endl;
    
    
    
    /* CA GMRES stuff
    SparseMat *example = new SparseMat();
    vector<double> samplevec;
    if(Utilities::A_SIZE==16){
        example->readFullMatrix("thbyth.txt");
        samplevec= Utilities::readVectorFile();
        Utilities::expandVec(samplevec,Utilities::A_SIZE);
    }else if(Utilities::A_SIZE<=2500){
        example->readFullMatrix("example_matpow.txt");
        samplevec= Utilities::readVectorFile("example_matpowvec.txt");
    }
    
    GMRES_sol sol = example->ca_GMRES(samplevec);
    
    if(Utilities::A_SIZE<20){
        Utilities::printDVector(sol.x);
    }*/
    
    cout << "Goodbye world!"<<endl;
    return 0;
}
