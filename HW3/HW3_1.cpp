#include<iostream>
#include<stdio.h>
#include<iomanip>
#include<cmath>
#include<fstream>
#include<ctime>
#include<random>
#include<algorithm>
using namespace std;
double rainbow(int nrolls, int n, double K, double r, double T, vector<double> &S, vector<double> &q,\
    vector<double> &sigma, vector<vector<double> > &A, int cases){
    std::random_device rd;
    std::default_random_engine generator(rd());
    std::normal_distribution<double> distribution(0.0,1.0);
    vector<vector<double> > zmat;
    if(cases==1){
        for(int i=0;i<nrolls;i++){
            vector<double> ztemp;
            for(int j=0;j<n;j++){
                ztemp.push_back(distribution(generator));
            }
            zmat.push_back(ztemp);
        }
    }
    else if(cases==2){
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
            for(int i=0;i<zmat.size();i++){
                zmat[i][j]=zmat[i][j]/standardDeviation;
            }
        }
    }
    
    vector<vector<double> > rmat;
    for(int i=0;i<nrolls;i++){
        vector<double> rtemp;
        for(int j=0;j<n;j++){
            double temp2=0;
            for(int k=0;k<n;k++){
                temp2+=zmat[i][k]*A[k][j];
            }
            rtemp.push_back(temp2);
        }
        rmat.push_back(rtemp);
    }
    vector<double> mean;
    for(int i=0;i<n;i++){
        double total=0;
        double mu=log(S[i])+(r-q[i]-sigma[i]*sigma[i]/2)*T;
        for(int j=0;j<nrolls;j++){
            rmat[j][i]+=mu;
            rmat[j][i]=exp(rmat[j][i]);
            total+=rmat[j][i];
        }
        mean.push_back(total/nrolls);
    }
    double maxi=0;
    for(int i=0;i<mean.size();i++){
        //cout<<mean[i]<<' '<<endl;
        if(mean[i]>maxi) maxi=mean[i];
    }
    if((maxi-K)>0) return (maxi-K)*exp(-r*T);
    else return 0;

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

    vector<vector<double> > A;
    for(int i=0;i<n;i++) A.push_back(temp);
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
    for(int k=0;k<n-1;k++) temp1=temp1+A[k][n-1]*A[k][n-1];
    A[n-1][n-1]=sqrt(covmat[n-1][n-1]-temp1);

    cout<<"\nA matrix:\n";
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            cout<<setw(12)<<A[i][j];
        }
        cout<<'\n';
    }

    vector<double> data,data1;
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


}

