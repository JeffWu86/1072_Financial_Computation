#include<iostream>
#include<stdio.h>
#include<iomanip>
#include<cmath>
#include<fstream>
#include<ctime>
#include<random>
#include<algorithm>
using namespace std;
vector<vector<double> > inverse(vector<vector<double> > A){
    vector<vector<double> > invA=A;
    for(int i=0;i<A.size();i++){
        for(int j=0;j<A[0].size();j++){
            if(i==j) invA[i][j]=1;
            else invA[i][j]=0;
        }
    }

    for(int i=0;i<A.size();i++){
        for(int j=0;j<=i;j++){
            if(i==j && A[i][j]!=1){
                double num=A[i][j];
                for(int k=0;k<A[0].size();k++){
                    A[i][k]/=num;
                    invA[i][k]/=num;
                }
            }
            else if(i!=j && A[i][j]!=0){
                double num=A[i][j];
                for(int k=0;k<A[0].size();k++){
                    A[i][k]=A[i][k]-num*A[j][k];
                    invA[i][k]=invA[i][k]-num*invA[j][k];
                }
            }
        }
    }
    for(int i=0;i<A.size()-1;i++){
        for(int j=i+1;j<A.size();j++){
            if(A[i][j]!=0){
                double num=A[i][j];
                for(int k=0;k<A[0].size();k++){
                    A[i][k]=A[i][k]-num*A[j][k];
                    invA[i][k]=invA[i][k]-num*invA[j][k];
                }
            }
        }
    }
for(int i=0;i<invA.size();i++){
        for(int j=0;j<invA[0].size();j++){
            cout<<setw(12)<<A[i][j];
        }cout<<endl;
    }
    cout<<endl;
        
    for(int i=0;i<invA.size();i++){
        for(int j=0;j<invA[0].size();j++){
            cout<<setw(12)<<invA[i][j];
        }cout<<endl;
    }cout<<endl;


    return invA;
}
int main(int argc, char* argv[])
{
	
    vector<vector<double> > A;
    vector<double> B;
    B.push_back(1.01);
    B.push_back(0.07);
    B.push_back(0.08);
    B.push_back(0.09);
    B.push_back(0.010);
    A.push_back(B);
    B[0]=0.02;
    B[1]=1.04;
    B[2]=0.06;
    B[3]=0.07;
    B[4]=0.12;
    A.push_back(B);
    B[0]=0.05;
    B[1]=0.07;
    B[2]=0.99;
    B[3]=0.008;
    B[4]=0.02;
    A.push_back(B);
    B[0]=0.08;
    B[1]=0.011;
    B[2]=0.001;
    B[3]=0.996;
    B[4]=0.0007;
    A.push_back(B);
    B[0]=0.001;
    B[1]=0.04;
    B[2]=0.01;
    B[3]=0.07;
    B[4]=0.989;
    A.push_back(B);
    vector<vector<double> > invA=inverse(A);
    
    for(int i=0;i<invA.size();i++){
        for(int j=0;j<invA[0].size();j++){
            cout<<setw(12)<<A[i][j];
        }cout<<endl;
    }
        cout<<endl;
    for(int i=0;i<invA.size();i++){
        for(int j=0;j<invA[0].size();j++){
            cout<<setw(12)<<invA[i][j];
        }cout<<endl;
    }
}

