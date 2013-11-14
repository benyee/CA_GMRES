//
//  GMRES_sol.h
//  
//
//  Created by Ben Yee on 11/13/13.
//
//

#ifndef _GMRES_sol_h
#define _GMRES_sol_h

#include <vector>
#include <iostream>

using namespace std;

//Structure for GMRES output
struct GMRES_sol{
    bool converged;
    unsigned int num_its;
    vector<double> res;
    vector<double> x;
    
    void print(){
        cout<<"---------------"<<endl;
        cout<<"GMRES output:"<<endl;
        if(converged){
            cout<<"GMRES converged in "<<num_its<<" iterations."<<endl;
        }else{
            cout<<"GMRES did not converge in "<<num_its<<" iterations."<<endl;
        }
        cout<<"Solution computed = "<<endl;
        cout<<"[ ";
        for(unsigned int i = 0; i < x.size(); i++){
            cout<<x[i]<<" ";
        }
        cout<<"]"<<endl;
        cout<<"Vector of residuals = "<<endl;
        cout<<"[ ";
        for(unsigned int i = 0; i < res.size(); i++){
            cout<<res[i]<<" ";
        }
        cout<<"]"<<endl;
        cout<<"--------------"<<endl;
    }
};

#endif
