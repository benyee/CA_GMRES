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
    
    SparseMat *sample = new SparseMat();
    sample->readFullMatrix();
    sample->print_matrix();
    
    vector< vector<double> > fullSample = sample->convertToFullMatrix();
    Utilities::printFullMatrix(fullSample);
    cout<<"The transpose of this is.."<<endl;
    Utilities::printFullMatrix(Utilities::transpose(fullSample));
    cout<<"A[1:2,2:3] ="<<endl;
    Utilities::printFullMatrix(Utilities::subMatrix(fullSample,make_pair(1,3),make_pair(2,4)));
    
    
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
    
    pair<vector< vector<double> >, vector< vector<double> > >  QR = Utilities::mgs(fullSample);
    vector < vector < double > > R1 = QR.second;
    
    cout << endl <<"MGS" << endl;
    Utilities::printFullMatrix(R1);
    
    vector < vector <double> > R2 = Utilities::tsQR(fullSample, 2);
    
    cout << endl <<"TSQR"<<endl;
    Utilities::printFullMatrix(R2);

    GMRES_sol sol;
    sol.converged = true;
    sol.num_its = 5;
    sol.res = samplevec;
    sol.x = samplevec2;
    
    sol.print();
    
    cout << "Goodbye world!"<<endl;
    return 0;
}
