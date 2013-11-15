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
    
    cout << "Testing reading in a sparse matrix..."<<endl;
    SparseMat *sample = new SparseMat();
    sample->readFullMatrix();
    sample->print_matrix();
    
    cout << "Testing converting a sparse matrix to a full matrix..."<<endl;
    vector< vector<double> > fullSample = sample->convertToFullMatrix();
    Utilities::printFullMatrix(fullSample);
    cout << "Testing transpose on a full matrix..."<<endl;
    cout<<"The transpose of this is.."<<endl;
    Utilities::printFullMatrix(Utilities::transpose(fullSample));
    cout << "Testing submatrix of a full matrix..."<<endl;
    cout<<"A[1:2,2:3] ="<<endl;
    Utilities::printFullMatrix(Utilities::subMatrix(fullSample,make_pair(1,3),make_pair(2,4)));
    
    cout << "Testing reading in two vectors..."<<endl;
    vector<double> samplevec = Utilities::readVectorFile();
    Utilities::printDVector(samplevec);
    vector<double> samplevec2 = Utilities::readVectorFile("example3.txt");
    Utilities::printDVector(samplevec2);
    cout<< "Dot product is: "<<Utilities::dotProd(samplevec, samplevec2)<<endl;
    
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
    
    cout<<endl<<"Testing mgs..."<<endl;
    pair<vector< vector<double> >, vector< vector<double> > >  QR = Utilities::mgs(fullSample);
    vector < vector < double > > R1 = QR.second;
    Utilities::printFullMatrix(R1);
    
    cout << endl <<"Testing TSQR..."<<endl;
    vector < vector <double> > R2 = Utilities::tsQR(fullSample, 2);
    Utilities::printFullMatrix(R2);
    
    cout<<endl<<"Testing expand matrix function..."<<endl;
    Utilities::printFullMatrix(Utilities::expandMat(fullSample,10,10));
    
    cout<<endl<<"Testing expand vector function..."<<endl;
    Utilities::printDVector(Utilities::expandVec(samplevec,11));
    
    cout<<endl<<"Testing GMRES_sol struct and its print function..."<<endl;
    GMRES_sol sol = Utilities::classicalGMRES(sample,samplevec);
    sol.print();
    
    
    SparseMat *bigmat = new SparseMat;
    bigmat->readFullMatrix("tallskinny.txt");
    vector<vector<double> > bigmatfull = bigmat->convertToFullMatrix(100000,Utilities::NUMCOLS);
    unsigned int start = clock();
    cout<<"Running mgs..."<<clock()-start<<endl;
    Utilities::mgs(bigmatfull);
    cout<<"Running tsQR..."<<clock()-start<<endl;
    unsigned int mgs = clock();
    Utilities::tsQR_fixed(bigmatfull);
    cout<<"Done!"<<clock()-mgs<<endl;
    
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
    
    
    
    cout << "Goodbye world!"<<endl;
    return 0;
}
