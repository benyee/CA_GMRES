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

#include "papi.h"
#include "SparseMat.h"
#include "Utilities.h"

using namespace std;

int main (int argc, char **argv)
{
    int retval;
  int Events[2];
  long long event_values[2];
  long long overhead_values[2];
  
  //Initialize the papi library
  retval = PAPI_library_init( PAPI_VER_CURRENT ); 

  if ( retval != PAPI_VER_CURRENT )
  {
     cout << "Failed to initialize PAPI!" << endl;
     return retval;
  }
  
  Events[0] = PAPI_FP_OPS;   //Floating Point Operations
  Events[1] = PAPI_TOT_CYC;  //Total number of clock cycles
  
  //Start the PAPI hardware counters 
  retval = PAPI_start_counters( (int *)Events, 2);
  if (retval != PAPI_OK )
  {
    cout << "Oops! Trouble starting counters!" << endl;
    return retval;
  }

  retval = PAPI_read_counters(overhead_values, 2);
  if (retval != PAPI_OK )
  {
    cout << "Oops! Trouble reading counter values!" << endl;
    return retval;
  }

  
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
    
    cout << "Goodbye world!"<<endl;
	
	
	
	//Read the current values from the counters
  retval = PAPI_read_counters(event_values, 2);
  if (retval != PAPI_OK )
  {
    cout << "Oops! Trouble reading counter values!" << endl;
    return retval;
  }
  cout << "PAPI counted " << overhead_values[0] << " floating point operations during initialization." << endl;
  cout << "PAPI counted for " << overhead_values[1] << " clock cycles during initialization." << endl << endl;
  cout << "PAPI counted " << event_values[0] - overhead_values[0] << " floating point operations." << endl;
  cout << "PAPI counted for " << event_values[1] - overhead_values[1] << " clock cycles." << endl;

  //Done using PAPI
  PAPI_shutdown();
  
    return 0;
}
