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

void Utilities::printDVector(const vector<double> &vec){
    cout<<"[ ";
    for(unsigned int i = 0; i < vec.size(); i++){
        cout<<vec[i]<<" ";
    }
    cout<<"]"<<endl;
}

void Utilities::printFullMatrix(const vector< vector<double> > &mat, bool rowFirst){
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

vector<unsigned int> Utilities::unsignedint_zeros(unsigned int rows){
    vector<unsigned int> out;
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

void Utilities::expandVec(vector<double> &vec, unsigned int rows){
    unsigned int oldNumRows = vec.size();
    if(oldNumRows == 0){
        vec = zeros(rows);
    }
    for(unsigned int i = oldNumRows; i<rows; i++){
        vec.push_back(0);
    }
}
void Utilities::expandMat(vector< vector<double> > &mat, unsigned int rows, unsigned int cols){
    unsigned int oldNumRows = mat.size();
    if(oldNumRows == 0){
        mat = zeros(rows,cols);
    }
    unsigned int oldNumCols = mat[0].size();
    if(oldNumCols == 0){
        mat = zeros(rows,cols);
    }
    for(unsigned int i = oldNumRows; i<rows; i++){
        mat.push_back(zeros(cols));
    }
    for(unsigned int i = 0; i<oldNumRows; i++){
        for(unsigned int j = oldNumCols;j<cols;j++){
            mat[i].push_back(0);
        }
    }
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

vector<double> Utilities::axpy(const vector<double> &x,const vector<double> &y){
    return axpy(x,1,y);
}
vector<double> Utilities::axpy(const vector<double> &x, double a){
    vector<double> y;
    return axpy(x,a,y);
}
vector<double> Utilities::axpy(const vector<double> &x, double a, const vector<double> &y){
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

vector<double> Utilities::matvec(const vector< vector<double> > &A, const vector<double> &x, bool rowFirst){
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

vector< vector<double> > Utilities::transpose(const vector< vector<double> > &A){
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

vector<double> Utilities::backSub(const vector< vector<double> > &A,const vector<double> &b,bool rowFirst){
    vector<double> x(b);
    for(int i = A.size()-1; i>=0;i--){
        for(unsigned int j = i+1; j < A.size(); j++){
            x[i] = x[i] - A[i][j]*x[j];
        }
        x[i] = x[i]/A[i][i];
    }
    
    return x;
}

//NOTE THAT THIS RETURNS Q-transpose, NOT Q
pair<vector< vector<double> >, vector< vector<double> > > Utilities::mgs(const vector< vector<double> > &mat){
    unsigned int numcols = mat.back().size();
    vector< vector<double> > Q = zeros(numcols,mat.size());
    vector< vector<double> > R = zeros(numcols,numcols);
    
    return mgs(mat, numcols,R,Q);
}

pair<vector< vector<double> >, vector< vector<double> > > Utilities::mgs(const vector< vector<double> > &mat, unsigned int numcols, vector<vector<double> > R, vector<vector<double> > Q){
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

void Utilities::mgs(const vector< vector<double> > &mat, vector<vector<double> > &R, vector<vector<double> > &Q){
    vector< vector<double> > At = transpose(mat);
    
    for(unsigned int i = 0; i<s;i++){
        R[i][i] = twoNorm(At[i]);
        if (R[i][i]==0){
			Q[i] = At[i];
        }
        else{
            Q[i] = axpy(At[i],1.0/R[i][i]);
        }
        for(unsigned int j = i+1; j<s;j++){
            R[i][j] = dotProd(Q[i],At[j]);
            At[j] = axpy(axpy(Q[i],R[i][j]),-1,At[j]);
        }
    }
}


pair<vector< vector<double> >, vector< vector<double> > > Utilities::mgs_col(const vector< vector<double> > &mat){
    unsigned int numcols = mat.size();
    vector< vector<double> > Q = zeros(numcols,mat[0].size());
    vector< vector<double> > R = zeros(numcols,numcols);
    
    return mgs_col(mat, numcols,R,Q);
}

pair<vector< vector<double> >, vector< vector<double> > > Utilities::mgs_col(vector< vector<double> > At, unsigned int numcols, vector<vector<double> > R, vector<vector<double> > Q){
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
            At[j] = axpy(Q[i],-1*R[i][j],At[j]);
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
void Utilities::mgs(double At[s][BLOCK_SIZE2], double R[s][s], double Q[s][BLOCK_SIZE2]){
    for(unsigned int i = 0; i<s;i++){
        R[i][i] = twoNorm(At[i]);
        if (R[i][i]==0){
			axpy(At[i],1,Q[i]);
        }
        else{
            axpy(At[i],1.0/R[i][i],Q[i]);
        }
        for(unsigned int j = i+1; j<s;j++){
            R[i][j] = dotProd(Q[i],At[j]);
            axpy(Q[i],-1*R[i][j],At[j],At[j]);
        }
    }
}


/*void Utilities::mgs(const vector<vector<double> > &mat,  double R[s][s], double Q[s][A_SIZE]){
    double At[s][A_SIZE];
    for(unsigned int j = 0; j<A_SIZE;j++){
        for(unsigned int i = 0;i<s;i++){
                At[i][j] = mat[j][i];
        }
    }
    
    for(unsigned int i = 0; i<s;i++){
        R[i][i] = twoNorm(At[i],A_SIZE);
        if (R[i][i]==0){
			axpy(At[i],1,Q[i],A_SIZE);
        }
        else{
            axpy(At[i],1.0/R[i][i],Q[i],A_SIZE);
        }
        for(unsigned int j = i+1; j<s;j++){
            R[i][j] = dotProd(Q[i],At[j],A_SIZE);
            axpy(Q[i],-1*R[i][j],At[j],At[j],A_SIZE);
        }
    }
}*/



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

double Utilities::twoNorm(double x[],unsigned int len){
    return sqrt(dotProd(x,x,len));
}

double Utilities::dotProd(double x[], double y[],unsigned int len){
    double sum = 0;
    for(unsigned int i = 0; i<len;i++){
        sum += x[i]*y[i];
    }
    return sum;
}

void Utilities::axpy(double x[], double y[], double out[],unsigned int len){
    return axpy(x,1,y,len);
}
void Utilities::axpy(double x[], double a, double out[],unsigned int len){
    if (a != 1){
        for(unsigned int i = 0; i<len;i++){
            out[i] = a*x[i];
        }
    }else{
        for(unsigned int i = 0; i<len;i++){
            out[i] = x[i];
        }
    }
}
void Utilities::axpy(double x[], double a, double y[],double out[],unsigned int len){
    if (a != 1){
        for(unsigned int i = 0; i<len;i++){
            out[i] = a*x[i] + y[i];
        }
    } else{
        for(unsigned int i = 0; i<len;i++){
            out[i] = x[i] + y[i];
        }
    }
}

void Utilities::RAtoAi(const vector<vector<double> >  &A, vector<vector<double> >  R, double Ai[s][BLOCK_SIZE2], double ind1){
    for(unsigned int j = 0; j<s;j++){
        for(unsigned int i = 0; i< s;i++){
            Ai[i][j] = R[j][i];
        }
    }
    for(unsigned int j = 0; j<BLOCK_SIZE;j++){
        for(unsigned int i = 0; i< s;i++){
            Ai[i][j+s] = A[j+ind1][i];
        }
    }
}
void Utilities::RAtoAi(const vector<vector<double> >  &A, double R[s][s], double Ai[s][BLOCK_SIZE2], double ind1){
    for(unsigned int j = 0; j<s;j++){
        for(unsigned int i = 0; i< s;i++){
            Ai[j][i] = R[i][j];
        }
    }
    for(unsigned int j = 0; j<BLOCK_SIZE;j++){
        for(unsigned int i = 0; i< s;i++){
            Ai[i][j+s] = A[j+ind1][i];
        }
    }
}

vector < vector <double> > Utilities::tsQR_fixed(const vector < vector < double> > &A){
    unsigned int numrow = A.size();
    unsigned int numblk = numrow/BLOCK_SIZE;
    
    vector < vector <double> > Ai = subMatrix(A, make_pair(0,BLOCK_SIZE), make_pair(0,s));
    pair<vector< vector<double> >, vector< vector<double> > >  QR = mgs(Ai);
    vector < vector <double> > R = QR.second;
    
    if(numblk > 1){
        // i = 2, since Q isn't the proper size yet
        Ai = subMatrix(A, make_pair(BLOCK_SIZE,2*BLOCK_SIZE), make_pair(0,s));
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
        double Ai_arr[s][BLOCK_SIZE2];
        RAtoAi(A,R,Ai_arr,BLOCK_SIZE2);
        double R_arr[s][s];
        double Q_arr[s][BLOCK_SIZE2];
        mgs(Ai_arr,R_arr,Q_arr);

        // Q and R are the right format i = 4 and beyond
        for(unsigned int i=3;i<numblk;i++){
            //unsigned int start = clock();
            RAtoAi(A,R_arr,Ai_arr,i*BLOCK_SIZE); // RAtoAi(A,R_arr,Ai_arr,(i-1)*numblk);
            //cout<<"RAtoAi took "<<clock()-start<<endl;
            //start = clock();
            mgs(Ai_arr,R_arr,Q_arr);
            //cout<<"mgs took "<<clock()-start<<endl;
        }
        
        //Convert R back to a vector:
        for(unsigned int i = 0; i<s;i++){
            for(unsigned int j = 0; j<s;j++){
                R[i][j] = R_arr[i][j];
            }
        }
        //Convert Q back to a vector:
        for(unsigned int i = 0; i<s;i++){
            for(unsigned int j = 0; j<BLOCK_SIZE2;j++){
                QR.first[i][j] = Q_arr[i][j];
            }
        }
    }
    
    
    //If there's an odd block out:
    if (numrow % BLOCK_SIZE != 0){
        Ai = subMatrix(A, make_pair((numblk-1)*BLOCK_SIZE,numrow), make_pair(0,s));
        Ai = stackMat(R,Ai);
        QR = mgs(Ai,s,R,QR.first);
        R = QR.second;
    }
    
    
    return R;
}

//*****************************************************************************
//***************COLUMN BASED STUFF BELOW**************************************
//*****************************************************************************

void Utilities::RAtoAi_col(const vector<vector<double> >  &A, vector<vector<double> >  R, double Ai[s][BLOCK_SIZE2], double ind1){
    for(unsigned int j = 0; j<s;j++){
        for(unsigned int i = 0; i< s;i++){
            Ai[i][j] = R[j][i];
        }
    }
    for(unsigned int i = 0; i< s;i++){
        for(unsigned int j = 0; j<BLOCK_SIZE;j++){
            Ai[i][j+s] = A[i][j+ind1];
        }
    }
}
void Utilities::RAtoAi_col(const vector<vector<double> >  &A, double R[s][s], double Ai[s][BLOCK_SIZE2], double ind1){
    for(unsigned int j = 0; j<s;j++){
        for(unsigned int i = 0; i< s;i++){
            Ai[i][j] = R[j][i];
        }
    }
    for(unsigned int i = 0; i< s;i++){
        for(unsigned int j = 0; j<BLOCK_SIZE;j++){
            Ai[i][j+s] = A[i][j+ind1];
        }
    }
}

//************************************************************************************************
//************************************************************************************************
//************************************************************************************************


vector< vector<double> > Utilities::subMatrix(const vector< vector<double> > &A, pair<unsigned int,unsigned int> ind1, pair<unsigned int,unsigned int> ind2){
    vector< vector<double> > Aout;
    
    unsigned int i = ind1.first;
    while(i<A.size() && i < ind1.second){
        unsigned int j = ind2.first;
        vector<double> temp;
        while(j < A[i].size() && j<ind2.second){
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

vector<double> Utilities::leastSquaresQR(const vector< vector<double> > &A,const  vector<double> &y){
    pair<vector< vector<double> >, vector< vector<double> > >  QR = mgs(A);
    vector<double> x = matvec(QR.first,y);
    x = backSub(QR.second,x);
    return x;
}

vector < vector <double> > Utilities::tsQR(const vector < vector < double> > &A,unsigned int blksiz){
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

double Utilities::infNorm(const vector<double> &x){
    double max = abs(x[0]);
    
    for(unsigned int i=1; i<x.size(); i++){
        double val=abs(x[i]);
        if(val>max){
            max = val;
        }
    }
    return max;
}

double Utilities::twoNorm(const vector<double> &x){
    return sqrt(dotProd(x,x));
}


/*****************************FINAL PUSH***************************************/
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

/************ MATRIX-MATRIX MULTIPLY****************/
void Utilities::matmat(const vector<vector<double> > &A, const vector<vector<double> > &B, vector<vector<double> > &AB){
    unsigned int indrowAB[2] = {0,AB.size()};
    unsigned int indcolAB[2] = {0,AB[0].size()};
    unsigned int shift[2] = {0,0};
    matmat(A,B,AB,shift,shift,indrowAB,indcolAB,B.size());
    
}

void Utilities::matmat(const vector<vector<double> > &A, const vector<vector<double> > &B, vector<vector<double> > &AB, unsigned int shiftA[2], unsigned int shiftB[2], unsigned int indrowAB[2],unsigned int indcolAB[2],unsigned int m_max){
    
    for(unsigned int j = indcolAB[0]; j<indcolAB[1];j++){
        unsigned int B_col = j-indcolAB[0]+shiftB[1];
        for(unsigned int i = indrowAB[0]; i<indrowAB[1];i++){
            unsigned int A_row = i-indrowAB[0]+shiftA[0];
            
            AB[j][i] = 0;
            for(unsigned m = 0; m<m_max;m++){
                AB[j][i] += A[m+shiftA[1]][A_row]*B[B_col][m+shiftB[0]];
            }
            
        }
    }
    
}

void Utilities::matmat(const vector<vector<double> > &A, double B[s][BLOCK_SIZE2], vector<vector<double> > &AB, unsigned int shiftA[2], unsigned int shiftB[2], unsigned int indrowAB[2],unsigned int indcolAB[2],unsigned int m_max){
    
    for(unsigned int j = indcolAB[0]; j<indcolAB[1];j++){
        unsigned int B_col = j-indcolAB[0]+shiftB[1];
        for(unsigned int i = indrowAB[0]; i<indrowAB[1];i++){
            unsigned int A_row = i-indrowAB[0]+shiftA[0];
            
            AB[j][i] = 0;
            for(unsigned m = 0; m<m_max;m++){
                AB[j][i] += A[m+shiftA[1]][A_row]*B[B_col][m+shiftB[0]];
            }
            
        }
    }
    
}

vector < vector <double> > Utilities::tsQR_col(const vector < vector < double> > &A, vector< vector<double> > &Q){
    unsigned int numrow = A[0].size();
    unsigned int numblk = numrow/BLOCK_SIZE;
    
    vector < vector <double> > Ai = subMatrix(A, make_pair(0,s), make_pair(0,BLOCK_SIZE));
    pair<vector< vector<double> >, vector< vector<double> > >  QR = mgs_col(Ai);
    vector < vector <double> > R = QR.second;
    vector < vector<double> > Qtemp(Q);
    for(unsigned int j = 0;j<s;j++){
        for(unsigned int i = 0; i<BLOCK_SIZE;i++){
            Q[j][i] = QR.first[j][i];
        }
    }
    
    unsigned int indrowAB[2] = {0,BLOCK_SIZE};
    unsigned int indcolAB[2] = {0,s};
    unsigned int shiftA[2] = {0,0};
    unsigned int shiftB[2] = {0,0};
    
    if(numblk > 1){
        // i = 2, since Q isn't the proper size yet
        Ai = subMatrix(A, make_pair(0,s), make_pair(BLOCK_SIZE,2*BLOCK_SIZE));
        Ai = stackMat(transpose(R),Ai,false);
        QR = mgs_col(Ai);
        R = QR.second;
        matmat(Q,QR.first,Qtemp,shiftA,shiftB,indrowAB,indcolAB,s);
        for(unsigned int j = 0;j<s;j++){
            for(unsigned int i = 0; i <indrowAB[1];i++){
                Q[j][i] = Qtemp[j][i];
            }
            for(unsigned int i = s; i<BLOCK_SIZE2;i++){
                Q[j][i+BLOCK_SIZE-s] = QR.first[j][i];
            }
        }
    }
    if(numblk > 2){
        // i = 3, R isn't a matrix yet
        double Ai_arr[s][BLOCK_SIZE2];
        RAtoAi_col(A,R,Ai_arr,BLOCK_SIZE2);
        double R_arr[s][s];
        double Q_arr[s][BLOCK_SIZE2];
        mgs(Ai_arr,R_arr,Q_arr);
        
        //BLOCK_SIZE2 = s+BLOCK_SIZE, static constant
        indrowAB[1] = 2*BLOCK_SIZE;
        matmat(Q,Q_arr,Qtemp,shiftA,shiftB,indrowAB,indcolAB,s);
        for(unsigned int j = 0;j<s;j++){
            for(unsigned int i = 0; i <indrowAB[1];i++){
                Q[j][i] = Qtemp[j][i];
            }
            for(unsigned int i = s; i<BLOCK_SIZE2;i++){
                Q[j][i+indrowAB[1]-s] = Q_arr[j][i];
            }
        }
        
        // Q and R are the right format i = 4 and beyond
        for(unsigned int i=3;i<numblk;i++){
            RAtoAi_col(A,R_arr,Ai_arr,i*BLOCK_SIZE);
            mgs(Ai_arr,R_arr,Q_arr);
            
            indrowAB[1] = i*BLOCK_SIZE;
            matmat(Q,Q_arr,Qtemp,shiftA,shiftB,indrowAB,indcolAB,s);
            for(unsigned int j = 0;j<s;j++){
                for(unsigned int i = 0; i <indrowAB[1];i++){
                    Q[j][i] = Qtemp[j][i];
                }
                for(unsigned int i = s; i<BLOCK_SIZE2;i++){
                    Q[j][i+indrowAB[1]-s] = Q_arr[j][i];
                }
            }
        }
        
        //Convert R back to a vector:
        for(unsigned int i = 0; i<s;i++){
            for(unsigned int j = 0; j<s;j++){
                R[i][j] = R_arr[i][j];
            }
        }
        //Convert Q back to a vector:
        for(unsigned int i = 0; i<s;i++){
            for(unsigned int j = 0; j<BLOCK_SIZE2;j++){
                QR.first[i][j] = Q_arr[i][j];
            }
        }
    }
    
    
    //If there's an odd block out:
    if (numrow % BLOCK_SIZE != 0){
        Ai = subMatrix(A, make_pair(0,s), make_pair(numblk*BLOCK_SIZE,numrow));
        Ai = stackMat(transpose(R),Ai,false);
        QR = mgs_col(Ai,s,R,QR.first);
        R = QR.second;
        
        
        indrowAB[1] = (numrow/BLOCK_SIZE)*BLOCK_SIZE;
        matmat(Q,QR.first,Qtemp,shiftA,shiftB,indrowAB,indcolAB,s);
        printFullMatrix(Qtemp);
        for(unsigned int j = 0;j<s;j++){
            for(unsigned int i = 0; i <indrowAB[1];i++){
                Q[j][i] = Qtemp[j][i];
            }
            for(unsigned int i = s; i<(numrow%BLOCK_SIZE)+s;i++){
                Q[j][i+indrowAB[1]-s] = QR.first[j][i];
            }
        }
    }
    
    return R;
}


/*****************************FINAL PUSH***************************************/
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
