#include<iostream>
#include<stdio.h>
#include<iomanip>
#include<cmath>
#include<fstream>
#include<ctime>
#include<random>
#include<algorithm>
#include<vector>
using namespace std;
double normalCDF(double x) // Phi(-∞, x) aka N(x)
{
    return erfc(-x/sqrt(2))/2;
}
double binom(int n, int k)
{
    double value = 1;

    for (int i = 0; i < k; i++)
        value = (value * (n - i)) / (i + 1);
    return value;
}
double binomial(int n, int j, double p)
{
    double value = 1;
    int pcount=n-j;
    int qcount=j;

    for (int i = 0; i < j; i++){
        value = (value * (n - i)) / (i + 1);
        while(value>1 && pcount>0){
            value*=p;
            pcount=pcount-1;
        }
        while(value>1 && qcount>0){
            value=value*(1-p);
            qcount=qcount-1;
        }
    }
    while(pcount>0){
            value*=p;
            pcount=pcount-1;
        }
    while(qcount>0){
        value=value*(1-p);
        qcount=qcount-1;
    }
    return value;
}
int main(int argc, char* argv[])
{
	ifstream input1;
	input1.open(argv[1]);
	if(!input1){ 
        return 1; 
    }

    double S;
    double K;
    double r;
    double q;
    double sigma;
    double T;
    int nrolls;     //number of experiments
    int repeat_times;
    int n;          //number of binomial tree

    input1>>S>>K>>r>>q>>sigma>>T>>nrolls>>repeat_times>>n;

    cout<<"S:"<<setw(12)<<S<<'\n';
    cout<<"K:"<<setw(12)<<K<<'\n';
    cout<<"r:"<<setw(12)<<r<<'\n';
    cout<<"q:"<<setw(12)<<q<<'\n';
    cout<<"\u03C3:"<<setw(12)<<sigma<<'\n';
    cout<<"T:"<<setw(12)<<T<<'\n';
    cout<<"# of sim:"<<setw(5)<<nrolls<<'\n';
    cout<<"# of rep:"<<setw(5)<<repeat_times<<'\n';
    cout<<"n:"<<setw(12)<<n<<'\n';
    double u=exp(sigma*sqrt(T/n));
    double d=1/u;
    double p=(exp(r*T/n-q*T/n)-d)/(u-d);

    



    //Bonus 2
    cout<<"\nBonus 2:\n";
    double cprice_euro_call_two=0,cprice_euro_put_two=0;;
    for(int j=0;j<=n;j++){
        //cprice_euro_call_two+=binom(n,j)*pow(p,n-j)*pow((1-p),j)*max(S*pow(u,n-j)*pow(d,j)-K,0.0);
        //cprice_euro_put_two+=binom(n,j)*pow(p,n-j)*pow((1-p),j)*max(K-S*pow(u,n-j)*pow(d,j),0.0);
        cprice_euro_call_two+=binomial(n,j,p)*max(S*pow(u,n-j)*pow(d,j)-K,0.0);
        cprice_euro_put_two+=binomial(n,j,p)*max(K-S*pow(u,n-j)*pow(d,j),0.0);
    }
    cout<<"Euro Call: "<<cprice_euro_call_two*exp(-r*T)<<'\n';
    cout<<"Euro Put : "<<cprice_euro_put_two*exp(-r*T)<<'\n';


}

