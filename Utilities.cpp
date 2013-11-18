//
//  Utilities.cpp
//  
//
//  Created by Ben Yee on 11/7/13.
//
//

#include "Utilities.h"

Utilities::Utilities(){
}

void Utilities::printDVector(vector<double> vec){
    cout<<"[ ";
    for(unsigned int i = 0; i < vec.size(); i++){
        cout<<vec[i]<<" ";
    }
    cout<<"]"<<endl;
}

void Utilities::printFullMatrix(vector< vector<double> > mat, bool rowFirst){
    if(rowFirst){
        cout<<"---Matrix output begin--"<<endl;
        for(unsigned int i = 0; i < mat.size();i++){
            printDVector(mat[i]);
        }
        cout<<"---Matrix output end--"<<endl;
    }else{
        cout<<"---Matrix output begin--"<<endl;
        for(unsigned int i = 0; i < mat.back().size();i++){
            cout<<"[ ";
            for(unsigned int j = 0; j < mat.size();j++){
                cout<<mat[j][i]<<" ";
            }
            cout<<"] "<<endl;
        }
        cout<<"---Matrix output end--"<<endl;
    }
}

vector<double> Utilities::zeros(unsigned int rows){
    vector<double> out;
    for(unsigned int i = 0; i<rows; i++){
        out.push_back(0);
    }
    return out;
}
vector < vector<double> > Utilities::zeros(unsigned int rows, unsigned int cols){
    vector< vector<double> > out;
    for(unsigned int i = 0; i<rows; i++){
        out.push_back(zeros(cols));
    }
    return out;
}

vector<double> Utilities::expandVec(vector<double> vec, unsigned int rows){
    unsigned int oldNumRows = vec.size();
    if(oldNumRows == 0){
        return zeros(rows);
    }
    for(unsigned int i = oldNumRows; i<rows; i++){
        vec.push_back(0);
    }
    return vec;
}
vector< vector<double> > Utilities::expandMat(vector< vector<double> > mat, unsigned int rows, unsigned int cols){
    unsigned int oldNumRows = mat.size();
    if(oldNumRows == 0){
        return zeros(rows,cols);
    }
    unsigned int oldNumCols = mat[0].size();
    if(oldNumCols == 0){
        return zeros(rows,cols);
    }
    for(unsigned int i = oldNumRows; i<rows; i++){
        mat.push_back(zeros(cols));
    }
    for(unsigned int i = 0; i<oldNumRows; i++){
        for(unsigned int j = oldNumCols;j<cols;j++){
            mat[i].push_back(0);
        }
    }
    return mat;
}

vector<double> Utilities::readVectorFile(string inputfile,char delim){
    ifstream inputFile;
    inputFile.open(inputfile.c_str());
    vector<double> out;
    
    string line;
    
    while(getline(inputFile,line)){ //while not at the end of the file
        double value = atof(line.c_str());
        out.push_back(value);
    }
    
    return out;
}

vector<double> Utilities::axpy(vector<double> x, vector<double> y){
    return axpy(x,1,y);
}
vector<double> Utilities::axpy(vector<double> x, double a){
    vector<double> y;
    return axpy(x,a,y);
}
vector<double> Utilities::axpy(vector<double> x, double a, vector<double> y){
    vector<double> out(x);
    if (a != 1){
        for(unsigned int i = 0; i<x.size();i++){
            out[i] = a*out[i];
        }
    }
    unsigned int minsize = min(x.size(),y.size());
    for(unsigned int i = 0; i<minsize;i++){
        out[i] = out[i] + y[i];
    }
    return out;
}

vector<double> Utilities::matvec(vector< vector<double> > A, vector<double> x, bool rowFirst){
    vector<double> out;
    if(rowFirst){
        unsigned int minSize = A.size();
        for(unsigned int i = 0;i<minSize;i++){
            out.push_back(dotProd(A[i],x));
        }
        return out;
    }
    unsigned int minSize = A.back().size();
    for(unsigned int i = 0;i<minSize;i++){
        out.push_back(0);
    }
    minSize = min(A.size(),x.size());
    for(unsigned int i = 0; i<minSize;i++){
        out = axpy(A[i],x[i],out);
    }
    return out;
}

vector< vector<double> > Utilities::transpose(vector< vector<double> > A){
    vector< vector<double> > At;
    unsigned int numcols = A[0].size();
    unsigned int numrows = A.size();
    
    for(unsigned int i = 0; i < numcols; i++){
        vector<double>  temp;
        for(unsigned int j = 0; j < numrows; j++){
            temp.push_back(0);
        }
        At.push_back(temp);
    }
    
    for(unsigned int i = 0; i < numcols; i++){
        for(unsigned int j = 0; j < numrows; j++){
            At[i][j] = A[j][i];
        }
    }
    
    return At;
}

double Utilities::dotProd(vector<double> x, vector<double> y){
    double sum = 0;
    unsigned int minSize = min(x.size(),y.size());
    for(unsigned int i = 0; i<minSize;i++){
        sum += x[i]*y[i];
    }
    return sum;
}

vector<double> Utilities::backSub(vector< vector<double> > A, vector<double> b,bool rowFirst){
    if(!rowFirst){
        A = transpose(A);
    }
    
    vector<double> x(b);
    for(int i = A.size()-1; i>=0;i--){
        for(unsigned int j = i+1; j < A.size(); j++){
            x[i] = x[i] - A[i][j]*x[j];
        }
        x[i] = x[i]/A[i][i];
    }
    
    if(!rowFirst){
        A = transpose(A);
    }
    
    return x;
}

//NOTE THAT THIS RETURNS Q-transpose, NOT Q
pair<vector< vector<double> >, vector< vector<double> > > Utilities::mgs(vector< vector<double> > mat){
    unsigned int numcols = mat.back().size();
    vector< vector<double> > Q = zeros(numcols,mat.size());
    vector< vector<double> > R = zeros(numcols,numcols);
    
    return mgs(mat, numcols,R,Q);
}

pair<vector< vector<double> >, vector< vector<double> > > Utilities::mgs(vector< vector<double> > mat, unsigned int numcols, vector<vector<double> > R, vector<vector<double> > Q){
    vector< vector<double> > At = transpose(mat);
    
    for(unsigned int i = 0; i<numcols;i++){
        R[i][i] = twoNorm(At[i]);
        if (R[i][i]==0){
			Q[i] = At[i];
        }
        else{
            Q[i] = axpy(At[i],1.0/R[i][i]);
        }
        for(unsigned int j = i+1; j<numcols;j++){
            R[i][j] = dotProd(Q[i],At[j]);
            At[j] = axpy(axpy(Q[i],R[i][j]),-1,At[j]);
        }
    }
    
    //NOTE THAT THIS RETURNS Q-transpose, NOT Q
    pair<vector< vector<double> >, vector< vector<double> > > output (Q,R);
    return output;
}


//*************************Fixed array size functions****************************************
//************************************************************************************************
//************************************************************************************************

//NOTE THAT THIS RETURNS Q-transpose, NOT Q
void Utilities::mgs(double At[NUMCOLS][BLOCK_SIZE2], unsigned int numcols, double R[NUMCOLS][NUMCOLS], double Q[NUMCOLS][BLOCK_SIZE2]){
    for(unsigned int i = 0; i<numcols;i++){
        R[i][i] = twoNorm(At[i]);
        if (R[i][i]==0){
			axpy(At[i],1,Q[i]);
        }
        else{
            axpy(At[i],1.0/R[i][i],Q[i]);
        }
        for(unsigned int j = i+1; j<numcols;j++){
            R[i][j] = dotProd(Q[i],At[j]);
            axpy(Q[i],-1*R[i][j],At[j],At[j]);
        }
    }
}
/*
//NOTE THAT THIS RETURNS Q-transpose, NOT Q
void Utilities::mgs(vector<vector<double> > &At, unsigned int numcols, double Q[NUMCOLS][BLOCK_SIZE2], unsigned int ind){
    for(unsigned int i = 0; i<numcols;i++){
        for(unsigned int k = 0; k<BLOCK_SIZE2;k++){
            At[i][i] += At[i][k+ind]*At[i][k+ind];
        }
        At[i][i] = sqrt(At[i][i]);
        if (At[i][i]==0){
			for(unsigned int k = 0; k<BLOCK_SIZE2;k++){
                Q[i][k] = 0;
            }
        }
        else{
            for(unsigned int k = 0; k<BLOCK_SIZE2;k++){
                Q[i][k] = At[i][k+ind]/At[i][i];
            }
        }
        for(unsigned int j = i+1; j<numcols;j++){
            for(unsigned int k = 0;  k<BLOCK_SIZE2;k++){
                At[j][i] += Q[i][k]*At[j][k+ind];
            }
            for(unsigned int k = 0; k<BLOCK_SIZE2;k++){
                At[j][k+ind] -= Q[i][k]*At[j][i];
            }
        }
    }
}*/

double Utilities::twoNorm(double x[BLOCK_SIZE2]){
    return sqrt(dotProd(x,x));
}

double Utilities::dotProd(double x[BLOCK_SIZE2], double y[BLOCK_SIZE2]){
    double sum = 0;
    for(unsigned int i = 0; i<BLOCK_SIZE2;i++){
        sum += x[i]*y[i];
    }
    return sum;
}

void Utilities::axpy(double x[BLOCK_SIZE2], double y[BLOCK_SIZE2], double out[BLOCK_SIZE2]){
    return axpy(x,1,y);
}
void Utilities::axpy(double x[BLOCK_SIZE2], double a, double out[BLOCK_SIZE2]){
    if (a != 1){
        for(unsigned int i = 0; i<BLOCK_SIZE2;i++){
            out[i] = a*x[i];
        }
    }else{
        for(unsigned int i = 0; i<BLOCK_SIZE2;i++){
            out[i] = x[i];
        }
    }
}
void Utilities::axpy(double x[BLOCK_SIZE2], double a, double y[BLOCK_SIZE2],double out[BLOCK_SIZE2]){
    if (a != 1){
        for(unsigned int i = 0; i<BLOCK_SIZE2;i++){
            out[i] = a*x[i] + y[i];
        }
    } else{
        for(unsigned int i = 0; i<BLOCK_SIZE2;i++){
            out[i] = x[i] + y[i];
        }
    }
}

void Utilities::RAtoAi(const vector<vector<double> >  &A, vector<vector<double> >  R, double Ai[NUMCOLS][BLOCK_SIZE2], double ind1){
    for(unsigned int j = 0; j<NUMCOLS;j++){
        for(unsigned int i = 0; i< NUMCOLS;i++){
            Ai[i][j] = R[j][i];
        }
    }
    for(unsigned int j = 0; j<BLOCK_SIZE;j++){
        for(unsigned int i = 0; i< NUMCOLS;i++){
            Ai[i][j+NUMCOLS] = A[j+ind1][i];
        }
    }
}
void Utilities::RAtoAi(const vector<vector<double> >  &A, double R[NUMCOLS][NUMCOLS], double Ai[NUMCOLS][BLOCK_SIZE2], double ind1){
    for(unsigned int j = 0; j<NUMCOLS;j++){
        for(unsigned int i = 0; i< NUMCOLS;i++){
            Ai[j][i] = R[i][j];
        }
    }
    for(unsigned int j = 0; j<BLOCK_SIZE;j++){
        for(unsigned int i = 0; i< NUMCOLS;i++){
            Ai[i][j+NUMCOLS] = A[j+ind1][i];
        }
    }
}

vector < vector <double> > Utilities::tsQR_fixed(vector < vector < double> > A){
    unsigned int numrow = A.size();
    unsigned int numblk = numrow/BLOCK_SIZE;
    
    vector < vector <double> > Ai = subMatrix(A, make_pair(0,BLOCK_SIZE), make_pair(0,NUMCOLS));
    pair<vector< vector<double> >, vector< vector<double> > >  QR = mgs(Ai);
    vector < vector <double> > R = QR.second;
    
    
    // i = 2, since Q isn't the proper size yet
    Ai = subMatrix(A, make_pair(BLOCK_SIZE,2*BLOCK_SIZE), make_pair(0,NUMCOLS));
    Ai = stackMat(R,Ai);
    QR = mgs(Ai);
    R = QR.second;
    /*
    unsigned int start = clock();
    vector < vector<double> > At = transpose(A);
    cout<<"transpose took "<<clock()-start<<endl;
    
    for(unsigned int i = 0; i<NUMCOLS;i++){
        for(unsigned int j = 0; j<NUMCOLS;j++){
            At[i][j] = R[j][i];
        }
    }
    double Q_arr[NUMCOLS][BLOCK_SIZE2];
    cout<<"numblk is "<<numblk<<endl;
    for(unsigned int i = 2; i<numblk;i++){
       // start = clock();
        mgs(At,NUMCOLS,Q_arr,i*BLOCK_SIZE-NUMCOLS);
        //cout<<"iteration took "<<clock()-start<<endl;
    }
    
    
    for(unsigned int i = 0; i<NUMCOLS;i++){
        for(unsigned int j = 0; j<NUMCOLS;j++){
            R[j][i] = At[i][j];
        }
    }*/
    
    // i = 3, R isn't a matrix yet
    double Ai_arr[NUMCOLS][BLOCK_SIZE2];
    RAtoAi(A,R,Ai_arr,2*BLOCK_SIZE);
    double R_arr[NUMCOLS][NUMCOLS];
    double Q_arr[NUMCOLS][BLOCK_SIZE2];
    mgs(Ai_arr,NUMCOLS,R_arr,Q_arr);

    // Q and R are the right format i = 4 and beyond
    for(unsigned int i=3;i<numblk;i++){
        //unsigned int start = clock();
        RAtoAi(A,R_arr,Ai_arr,i*BLOCK_SIZE); // RAtoAi(A,R_arr,Ai_arr,(i-1)*numblk);
        //cout<<"RAtoAi took "<<clock()-start<<endl;
        //start = clock();
        mgs(Ai_arr,NUMCOLS,R_arr,Q_arr);
        //cout<<"mgs took "<<clock()-start<<endl;
    }
    
    //Convert R back to a vector:
    for(unsigned int i = 0; i<NUMCOLS;i++){
        for(unsigned int j = 0; j<NUMCOLS;j++){
            R[i][j] = R_arr[i][j];
        }
    }
    //Convert Q back to a vector:
    for(unsigned int i = 0; i<NUMCOLS;i++){
        for(unsigned int j = 0; j<BLOCK_SIZE2;j++){
            QR.first[i][j] = Q_arr[i][j];
        }
    }
    
    
    //If there's an odd block out:
    if (numrow % BLOCK_SIZE != 0){
        Ai = subMatrix(A, make_pair((numblk-1)*BLOCK_SIZE,numrow), make_pair(0,NUMCOLS));
        Ai = stackMat(R,Ai);
        QR = mgs(Ai,NUMCOLS,R,QR.first);
        R = QR.second;
    }
    
    
    return R;
}
//************************************************************************************************
//************************************************************************************************
//************************************************************************************************


vector< vector<double> > Utilities::subMatrix(vector< vector<double> > A, pair<unsigned int,unsigned int> ind1, pair<unsigned int,unsigned int> ind2){
    vector< vector<double> > Aout;
    
    unsigned int i = ind1.first;
    while(i<A.size() && i < ind1.second){
        unsigned int j = ind2.first;
        vector<double> temp;
        while(j < A[i].size()){
            temp.push_back(A[i][j]);
            j++;
        }
        while(j<ind2.second){
            temp.push_back(0);
            j++;
        }
        Aout.push_back(temp);
        i++;
    }
    while(i<ind1.second){
        vector<double> temp;
        for(unsigned int j = ind1.first; j<A[0].size() && j<ind2.second; j++){
            temp.push_back(0);
        }
        Aout.push_back(temp);
        i++;
    }
    return Aout;
}

vector< vector<double> > Utilities::stackMat(vector<vector<double> > A, vector<vector<double> > B, bool rowFirst){
    if(rowFirst){
        A = transpose(A);
        B = transpose(B);
    }
    vector< vector<double> > out(A);
    
    unsigned int minSize = min(out.size(),B.size());
    
    for(unsigned int i = 0; i<minSize;i++){
        out[i].insert(out[i].end(),B[i].begin(),B[i].end());
    }
    
    if(rowFirst){
        A = transpose(A);
        B = transpose(B);
        out = transpose(out);
    }
    
    return out;
}

vector<double> Utilities::leastSquaresNormal(vector< vector<double> > A, vector<double> y){
    vector<double> x;
    
    //Stuff goes here
    
    return x;
}

vector<double> Utilities::leastSquaresQR(vector< vector<double> > A, vector<double> y){
    pair<vector< vector<double> >, vector< vector<double> > >  QR = mgs(A);
    vector<double> x = matvec(QR.first,y);
    x = backSub(QR.second,x);
    return x;
}

vector < vector <double> > Utilities::tsQR(vector < vector < double> > A,unsigned int blksiz){
    unsigned int numrow = A.size();
    unsigned int numcol = A[0].size();
    unsigned int numblk = numrow/blksiz;
    
    vector < vector <double> > Ai = subMatrix(A, make_pair(0,blksiz), make_pair(0,numcol));
    pair<vector< vector<double> >, vector< vector<double> > >  QR = mgs(Ai);
    vector < vector <double> > R = QR.second;
    
    // i = 2, since Q isn't the proper size yet
    Ai = subMatrix(A, make_pair(blksiz,2*blksiz), make_pair(0,numcol));
    Ai = stackMat(R,Ai);
    QR = mgs(Ai);
    R = QR.second;
    
    // Q and R are the right size i = 3 and beyond
    for(unsigned int i=3;i<=numblk;i++){
        Ai = subMatrix(A, make_pair((i-1)*blksiz,i*blksiz), make_pair(0,numcol));
        
        Ai = stackMat(R,Ai);
        
        QR = mgs(Ai,numcol,R,QR.first);
        R = QR.second;
    }
    //If there's an odd block out:
    if (numrow % blksiz != 0){
        Ai = subMatrix(A, make_pair((numblk-1)*blksiz,numrow), make_pair(0,numcol));
        Ai = stackMat(R,Ai);
        QR = mgs(Ai,numcol,R,QR.first);
        R = QR.second;
    }
        
    return R;
}

//Classical GMRES algorithm:
struct GMRES_sol Utilities::classicalGMRES(SparseMat* A, vector<double> b, double tol, unsigned int max_it){
    return classicalGMRES(A,b,zeros(b.size()),tol,max_it);
}
struct GMRES_sol Utilities::classicalGMRES(SparseMat* A, vector<double> b, vector<double> x, double tol, unsigned int max_it){
    GMRES_sol sol;
    sol.converged = false;
    vector< vector<double> > v;
    vector<double> res;
    vector<double> y;
    
    
    vector< vector<double> > h;
    
    vector<double> r = axpy(A->smvp(x),-1,b);
    double beta = twoNorm(r);
    res.push_back(beta);
    if(beta == 0){
        sol.converged = true;
        sol.num_its = 0;
        sol.x = x;
        sol.res =res;
        return sol;
    }
    v.push_back(axpy(r,1./beta));
    
    vector<double> x_0(x);
    vector<double> e_1;
    e_1.push_back(beta);
    unsigned int j = 0;
    while(j<max_it){
        cout<<"Running iteration "<<j+1<<endl;
        v.push_back(A->smvp(v[j]));
        h = expandMat(h,j+2,j+1);
        for(unsigned int i = 0; i<=j;i++){
            h[i][j] = dotProd(v[j+1],v[i]);
            v[j+1] = axpy(v[i],-h[i][j],v[j+1]);
        }
        h[j+1][j] = twoNorm(v[j+1]);
        if(h[j+1][j] == 0){
            sol.converged = true;
            res.push_back(0);
            break;
        }
        v[j+1] = axpy(v[j+1],1.0/h[j+1][j]);
        
        e_1 = expandVec(e_1,j+2);
        y = leastSquaresQR(h,e_1);
        x = axpy(x_0,matvec(v,y,false));
        res.push_back(twoNorm(axpy(A->smvp(x),-1,b)));
        
        j++;
        if(res.back()<=tol){
            sol.converged = true;
            break;
        }
    }
    
    
    sol.num_its = j;
    sol.x = x;
    sol.res =res;
    return sol;
}


double Utilities::infNorm(vector<double> x){
    double max = abs(x[0]);
    
    for(unsigned int i=1; i<x.size(); i++){
        double val=abs(x[i]);
        if(val>max){
            max = val;
        }
    }
    return max;
}

double Utilities::twoNorm(vector<double> x){
    return sqrt(dotProd(x,x));
}