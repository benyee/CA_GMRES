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
            A_counter = IA[i+1];
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
        double sum = 0;
        for(unsigned int j = IA[row]; j<IA[row+1];j++){
            sum += A[j]*vec[JA[j]];
        }
        outb.push_back(sum);
    }
    
    return outb;
}

void SparseMat::regMatrixPowers(vector<vector<double> > &V, const unsigned int ind[2]){
    if(ind[0] == 0){
        return;
    }
    for(unsigned int i = ind[0]; i<ind[1];i++){
        V[i] = smvp(V[i-1]);
    }
}

void SparseMat::matrixPowersMapper(){
    unsigned int counter = 0;
    //which vector am I working on right now?
    unsigned int whichvector = 0;
    vector<unsigned int> index = Utilities::unsignedint_zeros(Utilities::s);
    while(index[Utilities::s-1]<Utilities::A_SIZE){
        //cout<<"Checking first condition..."<<endl;
        if(index[whichvector] >= Utilities::A_SIZE && whichvector < Utilities::s){
            whichvector++;
//            cout<<"Increasing whichvector because I'm done with vector "<<whichvector<<endl;
            continue;
        }
        //cout<<"Checking second condition..."<<endl;
        if(whichvector < Utilities::s-1 && index[whichvector]>JA[ IA[ index[whichvector+1]+1 ]-1 ]){
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
        
        map[counter][1] = index[whichvector];
        map[counter][2] = whichvector;
        index[whichvector]++;
        counter++;
    }
}


//Matrix powers:
void SparseMat::matrixPowers(const vector<double> &v_0, vector<vector<double> > &V){
    //which vector am I working on right now?
    for(unsigned int i = 0; i<Utilities::A_SIZE*Utilities::s; i++){
        unsigned int index = map[i][1];
        unsigned int whichvector = map[i][2];
        if(whichvector){
            double local = 0;
            for(unsigned int j = IA[index]; j<IA[index+1];j++){
                local += A[j]*V[JA[j]][whichvector-1];
            }
            V[index][whichvector] = local;
        }else{
            double local = 0;
            for(unsigned int j = IA[index]; j<IA[index+1];j++){
                //cout<<"Doing Aval = "<<Aval[j]<<" times "<<v_0[JA[j]]<<endl;
                local += A[j]*v_0[JA[j]];
            }
            V[index][0] = local;
        }
    }
}
/*

void SparseMat::matrixPowers_fixed(const vector<double> &v_0, double (&V)[Utilities::A_SIZE][Utilities::s]){

    const unsigned int length = Utilities::A_SIZE;
    const unsigned int s = Utilities::s;
    //which vector am I working on right now?
    unsigned int whichvector = 0;
    
    unsigned int index[s] = {0};
//    vector<unsigned int> index = Utilities::unsignedint_zeros(s);
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
    
}



void SparseMat::matrixPowers_fixednorm(const vector<double> &v_0,  double V[Utilities::A_SIZE][Utilities::s]){
    for(unsigned int i = 0; i < Utilities::s; i++){
        unsigned int j; //Tracks the index of A vector
        if(i == 0){
            for(unsigned int row = 0; row<IA.size()-1; row++){
                j = IA[row];
                double sum = 0;
                while(j < IA[row+1]){
                    sum += A[j]*v_0[JA[j]];
                    j++;
                }
                V[row][i] = sum;
            }
            continue;
        }
        for(unsigned int row = 0; row<IA.size()-1; row++){
            j = IA[row];
            double sum = 0;
            while(j < IA[row+1]){
                sum += A[j]*V[JA[j]][i-1];
                j++;
            }
            V[row][i] = sum;
        }
    }
    
}*/



//Classical GMRES algorithm:
struct GMRES_sol SparseMat::classicalGMRES(const vector<double> &b, double tol, unsigned int max_it){
    return classicalGMRES(b,Utilities::zeros(b.size()),tol,max_it);
}
struct GMRES_sol SparseMat::classicalGMRES(const vector<double> &b, vector<double> x, double tol, unsigned int max_it){
    GMRES_sol sol;
    sol.converged = false;
    vector<double> res;
    
    unsigned int j;
    unsigned int restart = 0;
    
    
    while(restart < max_it && !sol.converged){
        restart++;
        vector< vector<double> > h;
        vector< vector<double> > v;
        vector<double> y;
        
        vector<double> r = Utilities::axpy(smvp(x),-1,b);
        double beta = Utilities::twoNorm(r);
        res.push_back(beta);
        if(beta == 0){
            sol.converged = true;
            break;
        }
        v.push_back(Utilities::axpy(r,1./beta));
        
        vector<double> x_0(x);
        vector<double> e_1;
        e_1.push_back(beta);
        j = 0;
        //One thing to consider: Check every RESTART/s loops instead of every loop.
        while(j<Utilities::RESTART){
            cout<<"Running iteration "<<(restart-1)*Utilities::RESTART+j+1<<endl;
            v.push_back(smvp(v[j]));
            Utilities::expandMat(h,j+2,j+1);
            for(unsigned int i = 0; i<=j;i++){
                h[i][j] = Utilities::dotProd(v[j+1],v[i]);
                v[j+1] = Utilities::axpy(v[i],-h[i][j],v[j+1]);
            }
            h[j+1][j] = Utilities::twoNorm(v[j+1]);
            if(h[j+1][j] == 0){
                sol.converged = true;
                res.push_back(0);
                break;
            }
            v[j+1] = Utilities::axpy(v[j+1],1.0/h[j+1][j]);
            
            Utilities::expandVec(e_1,j+2);
            y = Utilities::leastSquaresQR(h,e_1);
            x = Utilities::axpy(x_0,Utilities::matvec(v,y,false));
            res.push_back(Utilities::twoNorm(Utilities::axpy(smvp(x),-1,b)));
            j++;
            if(res.back()<=tol){
                sol.converged = true;
                break;
            }
        }
    }
    
    sol.num_its = (restart-1)*Utilities::RESTART+j;
    sol.x = x;
    sol.res =res;
    return sol;
}

//Communication avoiding GMRES:
struct GMRES_sol SparseMat::ca_GMRES(const vector<double> &b, double tol, unsigned int max_it){
    return ca_GMRES(b,Utilities::zeros(b.size()),tol,max_it);
}
struct GMRES_sol SparseMat::ca_GMRES(const vector<double> &b, vector<double> x, double tol, unsigned int max_it){
    GMRES_sol sol;
    sol.converged = false;
    vector<double> res;
    
    const double s = Utilities::s;
    const double A_SIZE= Utilities::A_SIZE;
    
    unsigned int restart = 0;
    
    vector<double> r0 = Utilities::axpy(b,-1,smvp(x));
    double beta = Utilities::twoNorm(r0);
    vector<double> q1 = Utilities::axpy(r0,1.0/beta);
    
    //Form basis matrix B:
    vector<vector<double> > B(s+1,vector<double>(s,0));
    for(unsigned int i = 0;i<=s;i++){
        B[i+1][i] = 1;
    }
    
    //Initialize V:
    vector<vector<double> > V(s,vector<double>(A_SIZE,0));
    V[0] = b;
    
    //Matrix powers:
    unsigned int ind[2] = {1,s};
    regMatrixPowers(V, ind);
    
    //Initialize more matrices.  blackQ, Q and V are stored as transposes.  All other matrices stored as regular matrices.
    vector<vector<double> > R;
    vector<vector<double> > Rinv(s,vector<double>(s,0));
    vector<vector<double> > Q(V);
    vector<vector<double> > Qtemp(Q);
    vector<vector<double> > blackQ(Utilities::RESTART,vector<double>(A_SIZE,0));
    
    //Might need to resize this later.  This is blackH_k
    vector<vector<double> > blackh_k(Utilities::RESTART+1,vector<double>(Utilities::RESTART,0));
    
    vector<vector<double> > temp(s+1,vector<double>(s,0));
    
    //Step 8:
    //beta*e1:
    vector<double> Be1(Utilities::RESTART+1,0);
    Be1[0] = beta;
    
    double hk;
    
    for(unsigned int k = 0; k<Utilities::RESTART/s;k++){
        if(k==0){
            //Step 6:
            R = Utilities::tsQR_col(V,blackQ,Qtemp);
            
            //Calculate Rinv:  (Step7)
            Utilities::invertUpperT(R,Rinv);
            //Calculate blackh_k  (Step 7)
            Utilities::matmat(B,Rinv,temp);
            Utilities::matmat(R,temp,blackh_k);
            
//            double hk = blackh_k[]
        }else{
            
        }
    }
    
    sol.converged = true;
    sol.num_its = (restart+1)*Utilities::RESTART;
    sol.x = x;
    sol.res =res;
    return sol;
}




