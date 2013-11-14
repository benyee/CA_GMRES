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
  long long load_values[2];
  long long mgs_values[2];
  long long tsqr_values[2];
  
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

  // Get overhead value for counters
  retval = PAPI_read_counters(overhead_values, 2);
  if (retval != PAPI_OK )
  {
    cout << "Oops! Trouble reading counter values!" << endl;
    return retval;
  }
  
  //Load Data
  cout << "Loading Data... " << endl;
    SparseMat *tmpA = new SparseMat();
    tmpA->readFullMatrix("tall_skinny_10000x100.txt");
	cout << "Data loaded!" << endl << "Converting to full matrix..." << endl;
    vector< vector<double> > tsA = tmpA->convertToFullMatrix(10000,100);
	cout << "Matrix converted!" << endl;

  //Get load value for counters
  retval = PAPI_read_counters(load_values, 2);
  if (retval != PAPI_OK )
  {
    cout << "Oops! Trouble reading counter values!" << endl;
    return retval;
  }
  
  //Compute MGS Solution
  cout << "Performing MGS... " << endl;
  pair<vector< vector<double> >, vector< vector<double> > >  QR = Utilities::mgs(tsA);
  vector < vector < double > > R1 = QR.second;
  cout << "MGS complete!" << endl;
  
  //Get MGS value for counters
  retval = PAPI_read_counters(mgs_values, 2);
  if (retval != PAPI_OK )
  {
    cout << "Oops! Trouble reading counter values!" << endl;
    return retval;
  }
  
 //Compute tsQR Solution
	cout << "Performing TSQR..." << endl;
   vector < vector <double> > R2 = Utilities::tsQR(tsA, 100);
   cout << "TSQR complete!" << endl;
  
  //Get tsQR value for counters
  retval = PAPI_read_counters(tsqr_values, 2);
  if (retval != PAPI_OK )
  {
    cout << "Oops! Trouble reading counter values!" << endl;
    return retval;
  }
  
  cout << endl;
  cout << "PAPI counted " << mgs_values[0] - load_values[0] - overhead_values[0] << " floating point operations during MGS." << endl;
  cout << "PAPI counted " << mgs_values[1] - load_values[1] - overhead_values[1] << " clock cycles during MGS." << endl;
  cout << endl;
  cout << "PAPI counted " << tsqr_values[0] - mgs_values[0] - overhead_values[0] << " floating point operations during TSQR." << endl;
  cout << "PAPI counted " << tsqr_values[1] - mgs_values[1] - overhead_values[1] << " clock cycles during TSQR." << endl;

    //Done using PAPI
  PAPI_shutdown();
  
    return 0;
}