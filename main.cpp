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
    
    bool fullDiagnosis = false;
    
    cout << "Testing reading in a sparse matrix..."<<endl;
    SparseMat *sample = new SparseMat();
    sample->readFullMatrix();
    sample->print_matrix();
    
    cout << "Testing converting a sparse matrix to a full matrix..."<<endl;
    vector< vector<double> > fullSample = sample->convertToFullMatrix();
    Utilities::printFullMatrix(fullSample);
    
    if(fullDiagnosis){
        cout << "Testing transpose on a full matrix..."<<endl;
        cout<<"The transpose of this is.."<<endl;
        Utilities::printFullMatrix(Utilities::transpose(fullSample));
        cout << "Testing submatrix of a full matrix..."<<endl;
        cout<<"A[1:2,2:3] ="<<endl;
        Utilities::printFullMatrix(Utilities::subMatrix(fullSample,make_pair(1,3),make_pair(2,4)));
    }
    
    cout << "Testing reading in two vectors..."<<endl;
    vector<double> samplevec = Utilities::readVectorFile();
    Utilities::printDVector(samplevec);
    vector<double> samplevec2 = Utilities::readVectorFile("example3.txt");
    Utilities::printDVector(samplevec2);
    cout<< "Dot product is: "<<Utilities::dotProd(samplevec, samplevec2)<<endl;
    
    if(fullDiagnosis){
        cout<< "A times vector 1 gives..."<<endl;
        Utilities::printDVector(sample->smvp(samplevec));
        cout<<"A times vector 1 with a full matvec gives..."<<endl;
        Utilities::printDVector(Utilities::matvec(fullSample,samplevec));
        cout<< "A times vector 2 gives..."<<endl;
        Utilities::printDVector(sample->smvp(samplevec2));
    
        cout<<"Ax=v1, x =  "<<endl;
        Utilities::printDVector(Utilities::leastSquaresQR(fullSample,samplevec));
        
        cout<<"5*v1 is..."<<endl;
        Utilities::printDVector(Utilities::axpy(samplevec,5.0));
        cout<<"v1+v2 is..."<<endl;
        Utilities::printDVector(Utilities::axpy(samplevec,samplevec2));
        cout<<"5*v1+v2 is..."<<endl;
        Utilities::printDVector(Utilities::axpy(samplevec,5,samplevec2));
        
        cout<<"Inf Norm: "<<Utilities::infNorm(samplevec)<<endl;
        cout<<"Two Norm: "<<Utilities::twoNorm(samplevec)<<endl;
    }
    
    cout<<endl<<"Testing mgs..."<<endl;
    pair<vector< vector<double> >, vector< vector<double> > >  QR = Utilities::mgs(fullSample);
    vector < vector < double > > R1 = QR.second;
    Utilities::printFullMatrix(R1);
    
    cout << endl <<"Testing TSQR..."<<endl;
    vector < vector <double> > R2 = Utilities::tsQR(fullSample, 2);
    Utilities::printFullMatrix(R2);
    
    if(fullDiagnosis){
        cout<<endl<<"Testing expand matrix function..."<<endl;
        vector<vector<double> > testfullmat(fullSample);
        Utilities::expandMat(testfullmat,10,10);
        Utilities::printFullMatrix(testfullmat);
        
        cout<<endl<<"Testing expand vector function..."<<endl;
        vector<double> testsamplevec(samplevec);
        Utilities::expandVec(testsamplevec,11);
        Utilities::printDVector(testsamplevec);
    }
    
    cout<<endl<<"Testing GMRES_sol struct and its print function..."<<endl;
    GMRES_sol sol = sample->classicalGMRES(samplevec);
    sol.print();
    
    
    SparseMat *bigmat = new SparseMat;
    bigmat->readFullMatrix("tallskinny.txt");
    vector<vector<double> > bigmatfull = bigmat->convertToFullMatrix(Utilities::A_SIZE,Utilities::s);
    cout<<"Running mgs on a tall skinny matrix..."<<endl;
    vector<vector<double> > Q(Utilities::s,vector<double>(Utilities::A_SIZE));
    vector<vector<double> > R(Utilities::s,vector<double>(Utilities::s));
    unsigned int start = clock();
    Utilities::mgs(bigmatfull);
    if(Utilities::A_SIZE ==4){
        Utilities::printFullMatrix(bigmatfull);
        cout<<"----R----"<<endl;
        for(unsigned int i = 0; i< Utilities::s; i++){
            for(unsigned int j =0; j<Utilities::s;j++){
                //cout<<R[i][j]<<'\t';
            }
            cout<<endl;
        }
        cout<<"----------"<<endl;
        return 0;
    }
    cout<<"Running tsQR..."<<clock()-start<<endl;
    unsigned int mgs = clock();
    Utilities::tsQR_fixed(bigmatfull);
    cout<<"Done!"<<clock()-mgs<<endl;
    
    if(fullDiagnosis){
        cout<<endl<<"Testing matvec for non-square matrices"<<endl;
        vector<vector<double> > nonsqmat;
        vector<double> temp;
        temp.push_back(1);
        temp.push_back(2);
        temp.push_back(3);
        vector<double> temp2;
        temp2.push_back(0);
        temp2.push_back(5);
        temp2.push_back(6);
        nonsqmat.push_back(temp);
        nonsqmat.push_back(temp2);
        vector<double> sample3vec;
        sample3vec.push_back(1);
        sample3vec.push_back(1);
        sample3vec.push_back(1);
        Utilities::printDVector(Utilities::matvec(nonsqmat,sample3vec));
    }
    

    
    if(Utilities::A_SIZE == 2500){
        cout<<endl<<"Testing matrix powers kernel"<<endl;
        SparseMat *secondsample = new SparseMat;
        secondsample->readFullMatrix("example_matpow.txt");
        vector<double> smvptest = Utilities::readVectorFile("example_matpowvec.txt");
        vector<vector<double> > V = Utilities::zeros(Utilities::A_SIZE,Utilities::s);
        secondsample->matrixPowersMapper();
        start = clock();
        secondsample->matrixPowers(smvptest,V);
        cout<<"matrix powers took "<<clock()-start<<endl;
        vector<vector<double> > matpowtest2 = Utilities::zeros(Utilities::s,smvptest.size());
        start = clock();
        for(unsigned int i = 0; i<Utilities::s;i++){
            matpowtest2[i] =secondsample->smvp(smvptest);
        }
        cout<<"regular mat pow took "<<clock()-start<<endl;
    }else if(Utilities::A_SIZE==4){/*
        cout<<endl<<"Testing matrix powers kernel"<<endl;
        double V[Utilities::A_SIZE][Utilities::s];
        //    vector<vector<double> > matpowtest =secondsample->matrixPowers_fixed(smvptest,15,2500,V);
        sample->matrixPowers_fixed(samplevec,V);
        cout<<endl;
        for(unsigned int i = 0; i< Utilities::A_SIZE;i++){
            for(unsigned int j = 0; j<Utilities::s;j++){
                cout<<V[i][j]<<'\t';
            }
            cout<<endl;
        }
        cout<<"Normal version of matrix powers.."<<endl;
        sample->matrixPowers_fixednorm(samplevec,V);
        cout<<endl;
        for(unsigned int i = 0; i< Utilities::A_SIZE;i++){
            for(unsigned int j = 0; j<Utilities::s;j++){
                cout<<V[i][j]<<' ';
            }
            cout<<endl;
        }*/
    }
    
    
    
    cout << "Goodbye world!"<<endl;
    return 0;
}
