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
const double PI =3.141592653589793238463;
void readdata(double &S, double &K, double &r, double &q, double &T, double &market_price, int &n, double &converge, ifstream &input1);
double Bisection(double S, double K, double r, double q, double T, double market_price, double converge, int choice, int n=-1);
double Newton(double S, double K, double r, double q, double T, double market_price, double converge, int choice, int n=-1);
double bino(double S, double K, double r, double q, double sigma, double T, int n, int choice);
double normalCDF(double x);
double normalCDFdiff(double x);
//double bino_euro(double S, double K, double r, double q, double sigma, double T, int n, int choice);
//double binomial(int n, int j, double p);
int main(int argc, char* argv[])
{
	double S;
    double K;
    double r;
    double q;
    double T;
    double market_price;
    int n;          //number of binomial tree
    double converge;
    ifstream input1;
    input1.open(argv[1]);
    if(!input1){ 
        return 1; 
    }
    readdata(S,K,r,q,T,market_price,n,converge,input1);
    cout<<fixed<<setprecision(6);

    //Black-Scholes formulas
    cout<<"\nBlack-Scholes formulas:\n";
    cout<<"Bisection Method\n";
    cout<<"Euro calls \u03C3: "<<Bisection(S,K,r,q,T,market_price,converge,1)<<'\n';
    cout<<"Euro puts  \u03C3: "<<Bisection(S,K,r,q,T,market_price,converge,2)<<'\n';

    cout<<"Newton's Method\n";
    cout<<"Euro calls \u03C3: "<<Newton(S,K,r,q,T,market_price,converge,1)<<'\n';
    cout<<"Euro puts  \u03C3: "<<Newton(S,K,r,q,T,market_price,converge,2)<<'\n';
  
    //CRR binomial tree model
    cout<<"\nCRR binomial tree model:\n";
    cout<<"Bisection Method\n";
    cout<<"Euro calls \u03C3: "<<Bisection(S,K,r,q,T,market_price,converge,3,n)<<'\n';
    cout<<"Euro puts  \u03C3: "<<Bisection(S,K,r,q,T,market_price,converge,4,n)<<'\n';
    cout<<"Amer calls \u03C3: "<<Bisection(S,K,r,q,T,market_price,converge,5,n)<<'\n';
    cout<<"Amer puts  \u03C3: "<<Bisection(S,K,r,q,T,market_price,converge,6,n)<<'\n';


    cout<<"Newton's Method\n";
    cout<<"Euro calls \u03C3: "<<Newton(S,K,r,q,T,market_price,converge,3,n)<<'\n';
    cout<<"Euro puts  \u03C3: "<<Newton(S,K,r,q,T,market_price,converge,4,n)<<'\n';
    cout<<"Amer calls \u03C3: "<<Newton(S,K,r,q,T,market_price,converge,5,n)<<'\n';
    cout<<"Amer puts  \u03C3: "<<Newton(S,K,r,q,T,market_price,converge,6,n)<<'\n';
}


void readdata(double &S, double &K, double &r, double &q, double &T, double &market_price, int &n, double &converge,ifstream &input1){
    
    input1>>S>>K>>r>>q>>T>>market_price>>n>>converge;

    cout<<"S:"<<setw(13)<<S<<'\n';
    cout<<"K:"<<setw(13)<<K<<'\n';
    cout<<"r:"<<setw(13)<<r<<'\n';
    cout<<"q:"<<setw(13)<<q<<'\n';
    cout<<"T:"<<setw(13)<<T<<'\n';
    cout<<"Price:  "<<market_price<<'\n';
    cout<<"n:"<<setw(13)<<n<<'\n';
    cout<<"Converge:"<<setw(5)<<converge<<'\n';
}
double Bisection(double S, double K, double r, double q, double T, \
                    double market_price, double converge, int choice,int n){
    //Choice 1 for Euro calls, 2 for Euro puts
    double d1a,d2a;
    double d1c,d2c;
    double sigma_a=0.0001,sigma_b=2,sigma_c;
    double price_a,price_c;
    while(1){
        sigma_c=(sigma_a+sigma_b)/2;
        if(n==-1){
            d1a=(log(S/K)+(r-q+sigma_a*sigma_a/2)*T)/sigma_a/sqrt(T);
            d2a=(log(S/K)+(r-q-sigma_a*sigma_a/2)*T)/sigma_a/sqrt(T);
            d1c=(log(S/K)+(r-q+sigma_c*sigma_c/2)*T)/sigma_c/sqrt(T);
            d2c=(log(S/K)+(r-q-sigma_c*sigma_c/2)*T)/sigma_c/sqrt(T);
        }
        if(choice==1){
            price_a=S*exp(-q*T)*normalCDF(d1a)-K*exp(-r*T)*normalCDF(d2a)-market_price;
            price_c=S*exp(-q*T)*normalCDF(d1c)-K*exp(-r*T)*normalCDF(d2c)-market_price;
        }
        else if(choice==2){
            price_a=K*exp(-r*T)*normalCDF(-d2a)-S*exp(-q*T)*normalCDF(-d1a)-market_price;
            price_c=K*exp(-r*T)*normalCDF(-d2c)-S*exp(-q*T)*normalCDF(-d1c)-market_price;
        }
        else if(choice==3){
            price_a=bino(S,K,r,q,sigma_a,T,n,1)-market_price;
            price_c=bino(S,K,r,q,sigma_c,T,n,1)-market_price;
        }
        else if(choice==4){
            price_a=bino(S,K,r,q,sigma_a,T,n,2)-market_price;
            price_c=bino(S,K,r,q,sigma_c,T,n,2)-market_price;
        }
        else if(choice==5){
            price_a=bino(S,K,r,q,sigma_a,T,n,3)-market_price;
            price_c=bino(S,K,r,q,sigma_c,T,n,3)-market_price;
        }
        else if(choice==6){
            price_a=bino(S,K,r,q,sigma_a,T,n,4)-market_price;
            price_c=bino(S,K,r,q,sigma_c,T,n,4)-market_price;
        }

        if(price_a*price_c<0) sigma_b=sigma_c;
        else sigma_a=sigma_c;
        if(abs(price_c-price_a)<converge)
            break;
    }
    return sigma_c;

}
double Newton(double S, double K, double r, double q, double T, \
                 double market_price, double converge, int choice, int n){
    //Choice 1 for Euro calls, 2 for Euro puts
    double d1a,d2a;
    double sigma_n1,sigma_n0=0.5;
    while(1){
        if(n==-1){
            d1a=(log(S/K)+(r-q+sigma_n0*sigma_n0/2)*T)/sigma_n0/sqrt(T);
            d2a=(log(S/K)+(r-q-sigma_n0*sigma_n0/2)*T)/sigma_n0/sqrt(T);
        }
        if(choice==1)
            sigma_n1=sigma_n0-(S*exp(-q*T)*normalCDF(d1a)-K*exp(-r*T)*normalCDF(d2a)-market_price)/ \
            (S*sqrt(T)*normalCDFdiff(d1a)*exp(-q*T));
        else if(choice==2)
            sigma_n1=sigma_n0-(K*exp(-r*T)*normalCDF(-d2a)-S*exp(-q*T)*normalCDF(-d1a)-market_price)/ \
            (S*sqrt(T)*normalCDFdiff(d1a)*exp(-q*T));
        else if(choice==3){
            double f2=bino(S,K,r,q,sigma_n0+0.00000001,T,n,1);
            double f1=bino(S,K,r,q,sigma_n0,T,n,1);
            sigma_n1=sigma_n0-(f1-market_price)/((f2-f1)/0.00000001);
        }
        else if(choice==4){
            double f2=bino(S,K,r,q,sigma_n0+0.00000001,T,n,2);
            double f1=bino(S,K,r,q,sigma_n0,T,n,2);
            sigma_n1=sigma_n0-(f1-market_price)/((f2-f1)/0.00000001);
        }
        else if(choice==5){
            double f2=bino(S,K,r,q,sigma_n0+0.00000001,T,n,3);
            double f1=bino(S,K,r,q,sigma_n0,T,n,3);
            sigma_n1=sigma_n0-(f1-market_price)/((f2-f1)/0.00000001);
        }
        else if(choice==6){
            double f2=bino(S,K,r,q,sigma_n0+0.00000001,T,n,4);
            double f1=bino(S,K,r,q,sigma_n0,T,n,4);
            sigma_n1=sigma_n0-(f1-market_price)/((f2-f1)/0.00000001);
        }
        if(abs(sigma_n1-sigma_n0)<converge)
            break;
        else sigma_n0=sigma_n1;
    }
    return sigma_n1;
}
double bino(double S, double K, double r, double q, double sigma, double T, int n, int choice){
    double u=exp(sigma*sqrt(T/n));
    double d=1/u;
    double p=(exp(r*T/n-q*T/n)-d)/(u-d);
    double price=0;
    vector<double> prices(n+1,0);
    //Choice 1 for euro call, 2 for euro put, 3 for amer call, 4 for amer put
    for(int i=0;i<=n;i++){
        if(choice==1 || choice==3)
            prices[i]=max(S*pow(u,n-i)*pow(d,i)-K,0.0);
        else if(choice==2 || choice==4)
            prices[i]=max(K-S*pow(u,n-i)*pow(d,i),0.0);
    }
    for(int i=n-1;i>=0;i--){
        for(int j=0;j<=i;j++){
            if(choice==1)
                prices[j]=exp(-r*T/n)*(p*prices[j]+(1-p)*prices[j+1]);
            else if(choice==2)
                prices[j]=exp(-r*T/n)*(p*prices[j]+(1-p)*prices[j+1]);
            else if(choice==3)
                prices[j]=max(exp(-r*T/n)*(p*prices[j]+(1-p)*prices[j+1]),S*pow(u,i-j)*pow(d,j)-K);
            else if(choice==4)
                prices[j]=max(exp(-r*T/n)*(p*prices[j]+(1-p)*prices[j+1]),K-S*pow(u,i-j)*pow(d,j));
        }
    }
    return prices[0];
}
double normalCDF(double x) // Phi(-∞, x) aka N(x)
{
    return erfc(-x/sqrt(2))/2;
}
double normalCDFdiff(double x){
    return exp(-(x*x/2))/sqrt(2*PI);
}
/*
//For bigger case, not necessary
double bino_euro(double S, double K, double r, double q, double sigma, double T, int n, int choice){
    double u=exp(sigma*sqrt(T/n));
    double d=1/u;
    double p=(exp(r*T/n-q*T/n)-d)/(u-d);
    double price=0;
    for(int j=0;j<=n;j++){
        if(choice==1)
            price+=binomial(n,j,p)*max(S*pow(u,n-j)*pow(d,j)-K,0.0);
        else
            price+=binomial(n,j,p)*max(K-S*pow(u,n-j)*pow(d,j),0.0);
    }
    return price*exp(-r*T);
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
}*/
