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

    return invA;
}
vector<vector<double> > cholesky(vector<vector<double> > &covmat){
    vector<vector<double> > A=covmat;
    for(int i=0;i<covmat.size();i++){
        for(int j=0;j<covmat[0].size();j++){
            A[i][j]=0;
        }
    }
    int n=covmat.size();
    //Step 1
    A[0][0]=sqrt(covmat[0][0]);
    for(int j=1;j<n;j++)
        A[0][j]=covmat[0][j]/A[0][0];
    //Step 2 and 3
    for(int i=1;i<n-1;i++){
        double temp1=0;
        for(int k=0;k<i;k++){
            temp1=temp1+A[k][i]*A[k][i];
        }
        A[i][i]=sqrt(covmat[i][i]-temp1);
        for(int j=i+1;j<n;j++){
            temp1=0;
            for(int k=0;k<i;k++){
                temp1=temp1+A[k][i]*A[k][j];
            }
            A[i][j]=(covmat[i][j]-temp1)/A[i][i];
        }

    }
    //Step 4
    double temp1=0;
    for(int k=0;k<n-1;k++) 
        temp1=temp1+A[k][n-1]*A[k][n-1];

    A[n-1][n-1]=sqrt(covmat[n-1][n-1]-temp1);

    return A;
}
double rainbow(int nrolls, int n, double K, double r, double T, vector<double> &S, vector<double> &q,\
    vector<double> &sigma, vector<vector<double> > A, int cases){
    std::random_device rd;
    std::default_random_engine generator(rd());
    std::normal_distribution<double> distribution(0.0,1.0);
    vector<vector<double> > zmat;
    vector<vector<double> > hstar=A;
    vector<vector<double> > Astar=A;
    vector<vector<double> > invAstar=A;
    vector<vector<double> > temp3=A;

    //Cases 1 for basic requirememnt, 2 for bonus 1
    if(cases==1){
        for(int i=0;i<nrolls;i++){
            vector<double> ztemp;
            for(int j=0;j<n;j++){
                ztemp.push_back(distribution(generator));
            }
            zmat.push_back(ztemp);
        }
    }
    else if(cases==2 || cases==3){
        for(int i=0;i<nrolls/2;i++){
            vector<double> ztemp;
            for(int j=0;j<n;j++){
                ztemp.push_back(distribution(generator));
            }
            zmat.push_back(ztemp);
        }
        for(int i=0;i<nrolls/2;i++){
            vector<double> ztemp;
            for(int j=0;j<n;j++){
                ztemp.push_back(-zmat[i][j]);
            }
            zmat.push_back(ztemp);
        }
        
        for(int j=0;j<n;j++){
            double meanz=0;
            for(int i=0;i<zmat.size();i++){ 
                meanz+=zmat[i][j];
            }
            meanz=meanz/zmat.size();
            double standardDeviation=0;
            for(int i=0;i<zmat.size();i++){
                standardDeviation+=(zmat[i][j]-meanz)*(zmat[i][j]-meanz);
            }
            standardDeviation=sqrt(standardDeviation/zmat.size());

            if(cases==2){
                for(int i=0;i<zmat.size();i++){
                    zmat[i][j]=zmat[i][j]/standardDeviation;
                }
            }
        
            hstar[j][j]=standardDeviation*standardDeviation;
            
        }
    }
    if(cases==3){
        for(int i=0;i<n-1;i++){
            for(int j=i+1;j<n;j++){
                double sum=0;
                for(int k=0;k<nrolls;k++){
                    sum+=zmat[k][i]*zmat[k][j];
                }
                sum/=nrolls;
                hstar[i][j]=sum;
                hstar[j][i]=sum;
            }
        }
        Astar=cholesky(hstar);
        // cout<<"hstar\n";
        // for(int i=0;i<n;i++){
        //     for(int j=0;j<n;j++){
        //         cout<<setw(12)<<hstar[i][j];
        //     }
        //     cout<<endl;
        // }
        // cout<<"Astar\n";
        // for(int i=0;i<n;i++){
        //     for(int j=0;j<n;j++){
        //         cout<<setw(12)<<Astar[i][j];
        //     }
        //     cout<<endl;
        // }
        // cout<<"Inverse Astar\n";
        invAstar=inverse(Astar);
        // for(int i=0;i<n;i++){
        //     for(int j=0;j<n;j++){
        //         cout<<setw(12)<<invAstar[i][j];
        //     }
        //     cout<<endl;
        // }
        for(int i=0;i<Astar.size();i++){
            for(int j=0;j<Astar[0].size();j++){
                double kk=0;
                for(int k=0;k<Astar.size();k++){
                    kk=kk+invAstar[i][k]*A[k][j];
                }
                temp3[i][j]=kk;
            }
        }
        // cout<<"A :\n";
        // for(int i=0;i<n;i++){
        //     for(int j=0;j<n;j++){
        //         cout<<setw(12)<<A[i][j];
        //     }
        //     cout<<endl;
        // }
        // for(int i=0;i<A.size();i++)
        //     for(int j=0;j<A[0].size();j++)
        //         A[i][j]=temp3[i][j];

        // cout<<"A become:\n";
        // for(int i=0;i<n;i++){
        //     for(int j=0;j<n;j++){
        //         cout<<setw(12)<<A[i][j];
        //     }
        //     cout<<endl;
        // }cout<<endl;
        // vector<vector<double> > rmat;
        // for(int i=0;i<nrolls;i++){
        //     vector<double> rtemp;
        //     for(int j=0;j<n;j++){
        //         double temp2=0;
        //         for(int k=0;k<n;k++){
        //             temp2+=zmat[i][k]*invAstar[k][j];
        //         }
        //         rtemp.push_back(temp2);
        //     }
        //     rmat.push_back(rtemp);
        // }
        // for(int i=0;i<n-1;i++){
        //     for(int j=i+1;j<n;j++){
        //         double sum=0;
        //         for(int k=0;k<nrolls;k++){
        //             sum+=rmat[k][i]*rmat[k][j];
        //         }
        //         sum/=nrolls;
        //         hstar[i][j]=sum;
        //         hstar[j][i]=sum;
        //     }
        // }
        // cout<<"hstar\n";
        // for(int i=0;i<n;i++){
        //     for(int j=0;j<n;j++){
        //         cout<<setw(12)<<hstar[i][j];
        //     }
        //     cout<<endl;
        // }
    }

    
    vector<vector<double> > rmat;
    for(int i=0;i<nrolls;i++){
        vector<double> rtemp;
        for(int j=0;j<n;j++){
            double temp2=0;
            for(int k=0;k<n;k++){
                if(cases!=3)
                    temp2+=zmat[i][k]*A[k][j];
                else
                    temp2+=zmat[i][k]*temp3[k][j];
            }
            rtemp.push_back(temp2);
        }
        rmat.push_back(rtemp);
    }
    
    for(int i=0;i<n;i++){
        double mu=log(S[i])+(r-q[i]-sigma[i]*sigma[i]/2)*T;
        for(int j=0;j<nrolls;j++){
            rmat[j][i]+=mu;
            rmat[j][i]=exp(rmat[j][i]);
        }
    }
    
    vector<double> payoff;
    for(int i=0;i<nrolls;i++){
        double maxi=-10000;
        for(int j=0;j<n;j++){
            if(rmat[i][j]>maxi) maxi=rmat[i][j];
        }
        if((maxi-K)>0) 
            payoff.push_back((maxi-K)*exp(-r*T));
        else payoff.push_back(0);
    }

    double total=0;
    for(int i=0;i<nrolls;i++){
        total+=payoff[i];
    }

    return total/nrolls;

}
int main(int argc, char* argv[])
{
	ifstream input1;
	input1.open(argv[1]);
	if(!input1){ 
        return 1; 
    }

    double K;
    double r;
    double T;
    int nrolls;         //number of simulations
    int repeat_times;   //number of repetitions
    int n;              //number of binomial tree
    vector<double> S;
    vector<double> q;
    vector<double> sigma;
    vector<vector<double> > rho;

    input1>>K>>r>>T>>nrolls>>repeat_times>>n;
 
    cout<<"K:"<<setw(12)<<K<<'\n';
    cout<<"r:"<<setw(12)<<r<<'\n';
    cout<<"T:"<<setw(12)<<T<<'\n';
    cout<<"# of sim:"<<setw(5)<<nrolls<<'\n';
    cout<<"# of rep:"<<setw(5)<<repeat_times<<'\n';
    cout<<"n:"<<setw(12)<<n<<'\n';
    
    S.reserve(n);
    for(int i=0;i<n;i++){
        input1>>S[i];
        cout<<"S["<<i+1<<"]: "<<setw(8)<<S[i]<<'\n';
    }
    q.reserve(n);
    for(int i=0;i<n;i++){
        input1>>q[i];
        cout<<"q["<<i+1<<"]: "<<setw(8)<<q[i]<<'\n';
    }
    sigma.reserve(n);
    for(int i=0;i<n;i++){
        input1>>sigma[i];
        cout<<"\u03C3["<<i+1<<"]: "<<setw(8)<<sigma[i]<<'\n';
    }

    //Initaite 2d array for rho
    vector<double> temp(n,0.0);
    for(int i=0;i<n;i++) rho.push_back(temp);

    for(int i=0;i<n-1;i++)
        for(int j=i+1;j<n;j++){
            input1>>rho[i][j];
            rho[j][i]=rho[i][j];
            cout<<"\u03C1"<<i+1<<j+1<<":"<<setw(10)<<rho[i][j]<<'\n';
        }
    vector<vector<double> > covmat=rho;
    for(int i=0;i<n;i++)
    for(int j=0;j<n;j++){
        if(i==j) covmat[i][j]=sigma[i]*sigma[i]*T;
        else covmat[i][j]=covmat[i][j]*sigma[i]*sigma[j]*T;
    }

    cout<<"\nVariance and covariance matrix:\n";
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            cout<<setw(7)<<covmat[i][j];
        }
        cout<<'\n';
    }

    //Generate A matrix
    vector<vector<double> > A;
    for(int i=0;i<n;i++) A.push_back(temp);
    
    A=cholesky(covmat);


    cout<<"\nA matrix:\n";
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            cout<<setw(12)<<A[i][j];
        }
        cout<<'\n';
    }

    vector<double> data,data1,data2;
    double price,total=0;
    for(int k=0;k<repeat_times;k++){
        price=rainbow(nrolls,n,K,r,T,S,q,sigma,A,1);
        data.push_back(price);
        total+=price;
    }
    double meanprice=total/repeat_times;

    double standardDeviation;
    for(int i=0;i<repeat_times;i++)
        standardDeviation+=(data[i]-meanprice)*(data[i]-meanprice);
    standardDeviation=sqrt(standardDeviation/repeat_times);

    cout<<"Basic requirement: \n";
    cout<<"Standard Deviation: "<<standardDeviation<<'\n';
    cout<<"Under "<<repeat_times<<" times of testing, 95\% of confidence interval for the call option value are: "<<'\n';
    cout<<"Lower bounds: "<<meanprice-2*standardDeviation<<'\n';
    cout<<"Mean:         "<<meanprice<<'\n';
    cout<<"Upper bounds: "<<meanprice+2*standardDeviation<<'\n';


    total=0;
    standardDeviation=0;
    for(int k=0;k<repeat_times;k++){
        price=rainbow(nrolls,n,K,r,T,S,q,sigma,A,2);
        data1.push_back(price);
        total+=price;
    }
    meanprice=total/repeat_times;

    for(int i=0;i<repeat_times;i++)
        standardDeviation+=(data1[i]-meanprice)*(data1[i]-meanprice);
    standardDeviation=sqrt(standardDeviation/repeat_times);

    cout<<"\nBonus1 : \n";
    cout<<"Standard Deviation: "<<standardDeviation<<'\n';
    cout<<"Under "<<repeat_times<<" times of testing, 95\% of confidence interval for the call option value are: "<<'\n';
    cout<<"Lower bounds: "<<meanprice-2*standardDeviation<<'\n';
    cout<<"Mean:         "<<meanprice<<'\n';
    cout<<"Upper bounds: "<<meanprice+2*standardDeviation<<'\n';

    total=0;
    standardDeviation=0;
    for(int k=0;k<repeat_times;k++){
        price=rainbow(nrolls,n,K,r,T,S,q,sigma,A,3);
        data2.push_back(price);
        total+=price;
    }
    meanprice=total/repeat_times;

    for(int i=0;i<repeat_times;i++)
        standardDeviation+=(data2[i]-meanprice)*(data2[i]-meanprice);
    standardDeviation=sqrt(standardDeviation/repeat_times);

    cout<<"\nBonus2 : \n";
    cout<<"Standard Deviation: "<<standardDeviation<<'\n';
    cout<<"Under "<<repeat_times<<" times of testing, 95\% of confidence interval for the call option value are: "<<'\n';
    cout<<"Lower bounds: "<<meanprice-2*standardDeviation<<'\n';
    cout<<"Mean:         "<<meanprice<<'\n';
    cout<<"Upper bounds: "<<meanprice+2*standardDeviation<<'\n';
}

