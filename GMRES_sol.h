//
//  GMRES_sol.h
//  
//
//  Created by Ben Yee on 11/13/13.
//
//

#ifndef _GMRES_sol_h
#define _GMRES_sol_h

//Structure for GMRES output
struct GMRES_sol{
    unsigned int num_its;
    vector<double> res;
    vector<double> x;
};

#endif
