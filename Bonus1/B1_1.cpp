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
double normalCDF(double x) // Phi(-∞, x) aka N(x)
{
    return erfc(-x/sqrt(2))/2;
}
double normalCDFdiff(double x){
    return exp(-(x*x/2))/sqrt(2*PI);
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
    //double sigma;
    double T;
    double market_price;
    int n;          //number of binomial tree
    double converge;
    

    input1>>S>>K>>r>>q>>T>>market_price>>n>>converge;

    cout<<"S:"<<setw(13)<<S<<'\n';
    cout<<"K:"<<setw(13)<<K<<'\n';
    cout<<"r:"<<setw(13)<<r<<'\n';
    cout<<"q:"<<setw(13)<<q<<'\n';
    //cout<<"\u03C3:"<<setw(12)<<sigma<<'\n';
    cout<<"T:"<<setw(13)<<T<<'\n';
    cout<<"Price:  "<<market_price<<'\n';
    cout<<"n:"<<setw(13)<<n<<'\n';
    cout<<"Converge:"<<setw(5)<<converge<<'\n';

    //Black-Scholes formulas
    cout<<"\nBlack-Scholes formulas:\n";
    double d1a,d1c;//=(log(S/K)+(r-q+sigma*sigma/2)*T)/sigma/sqrt(T);
    double d2a,d2c;//=(log(S/K)+(r-q-sigma*sigma/2)*T)/sigma/sqrt(T);

    cout<<"Bisection Method\n";
    double sigma_a=0.0001,sigma_b=2,sigma_c;
    double price_a,price_c;
    while(1){
        sigma_c=(sigma_a+sigma_b)/2;
        d1a=(log(S/K)+(r-q+sigma_a*sigma_a/2)*T)/sigma_a/sqrt(T);
        d2a=(log(S/K)+(r-q-sigma_a*sigma_a/2)*T)/sigma_a/sqrt(T);
        d1c=(log(S/K)+(r-q+sigma_c*sigma_c/2)*T)/sigma_c/sqrt(T);
        d2c=(log(S/K)+(r-q-sigma_c*sigma_c/2)*T)/sigma_c/sqrt(T);
        price_a=S*exp(-q*T)*normalCDF(d1a)-K*exp(-r*T)*normalCDF(d2a)-market_price;
        price_c=S*exp(-q*T)*normalCDF(d1c)-K*exp(-r*T)*normalCDF(d2c)-market_price;
        if(price_a*price_c<0) sigma_b=sigma_c;
        else sigma_a=sigma_c;
        if(abs(price_c-price_a)<converge)
            break;
    }
    cout<<"Euro calls \u03C3: "<<fixed<<setprecision(6)<<sigma_c<<'\n';
    sigma_a=0.0001,sigma_b=2;
    while(1){
        sigma_c=(sigma_a+sigma_b)/2;
        d1a=(log(S/K)+(r-q+sigma_a*sigma_a/2)*T)/sigma_a/sqrt(T);
        d2a=(log(S/K)+(r-q-sigma_a*sigma_a/2)*T)/sigma_a/sqrt(T);
        d1c=(log(S/K)+(r-q+sigma_c*sigma_c/2)*T)/sigma_c/sqrt(T);
        d2c=(log(S/K)+(r-q-sigma_c*sigma_c/2)*T)/sigma_c/sqrt(T);
        price_a=K*exp(-r*T)*normalCDF(-d2a)-S*exp(-q*T)*normalCDF(-d1a)-market_price;
        price_c=K*exp(-r*T)*normalCDF(-d2c)-S*exp(-q*T)*normalCDF(-d1c)-market_price;
        if(price_a*price_c<0) sigma_b=sigma_c;
        else sigma_a=sigma_c;
        if(abs(price_c-price_a)<converge)
            break;
    }
    cout<<"Euro puts  \u03C3: "<<fixed<<setprecision(6)<<sigma_c<<'\n';


    cout<<"Newton's Method\n";
    double sigma_n1,sigma_n0=0.5;
    for(int i=0;i<5;i++){
        d1a=(log(S/K)+(r-q+sigma_n0*sigma_n0/2)*T)/sigma_n0/sqrt(T);
        d2a=(log(S/K)+(r-q-sigma_n0*sigma_n0/2)*T)/sigma_n0/sqrt(T);
        sigma_n1=sigma_n0-(S*exp(-q*T)*normalCDF(d1a)-K*exp(-r*T)*normalCDF(d2a)-market_price)/ \
        (S*sqrt(T)*normalCDFdiff(d1a)*exp(-q*T));
        if(abs(sigma_n1-sigma_n0)<converge)
            break;
        else sigma_n0=sigma_n1;
    }
    cout<<"Euro calls \u03C3: "<<fixed<<setprecision(6)<<sigma_n1<<'\n';
    sigma_n0=0.5;
    for(int i=0;i<5;i++){
        d1a=(log(S/K)+(r-q+sigma_n0*sigma_n0/2)*T)/sigma_n0/sqrt(T);
        d2a=(log(S/K)+(r-q-sigma_n0*sigma_n0/2)*T)/sigma_n0/sqrt(T);
        sigma_n1=sigma_n0-(K*exp(-r*T)*normalCDF(-d2a)-S*exp(-q*T)*normalCDF(-d1a)-market_price)/ \
        (S*sqrt(T)*normalCDFdiff(d1a)*exp(-q*T));
        if(abs(sigma_n1-sigma_n0)<converge)
            break;
        else sigma_n0=sigma_n1;
    }
    cout<<"Euro puts  \u03C3: "<<fixed<<setprecision(6)<<sigma_n1<<'\n';
  /*  cout<<"Call:"<<setw(9)<<S*exp(-q*T)*normalCDF(d1)-K*exp(-r*T)*normalCDF(d2)<<'\n';
    cout<<"Put:"<<setw(10)<<K*exp(-r*T)*normalCDF(-d2)-S*exp(-q*T)*normalCDF(-d1)<<'\n';
    

    

    //CRR binomial tree model
    cout<<"\nCRR binomial tree model:\n";
    double u=exp(sigma*sqrt(T/n));
    double d=1/u;
    double p=(exp(r*T/n-q*T/n)-d)/(u-d);
    vector<double> zeroo(n+1,0);
    vector<vector<double> > sprice,cprice_amer_call,cprice_euro_call,cprice_amer_put,cprice_euro_put;
    for(int i=0;i<n+1;i++){
        sprice.push_back(zeroo);
        cprice_amer_call.push_back(zeroo);
        cprice_euro_call.push_back(zeroo);
        cprice_amer_put.push_back(zeroo);
        cprice_euro_put.push_back(zeroo);
    }
    for(int i=0;i<=n;i++){
        for(int j=0;j<=i;j++){
            sprice[i][j]=S*pow(u,i-j)*pow(d,j);
        }
    }
    for(int i=0;i<n;i++){
        cprice_euro_call[n-1][i]=exp(-r*T/n)*(p*max(sprice[n][i]-K,0.0)+(1-p)*max(sprice[n][i+1]-K,0.0));
        cprice_euro_put[n-1][i]=exp(-r*T/n)*(p*max(K-sprice[n][i],0.0)+(1-p)*max(K-sprice[n][i+1],0.0));
        cprice_amer_call[n-1][i]=exp(-r*T/n)*(p*max(sprice[n][i]-K,0.0)+(1-p)*max(sprice[n][i+1]-K,0.0));
        cprice_amer_put[n-1][i]=exp(-r*T/n)*(p*max(K-sprice[n][i],0.0)+(1-p)*max(K-sprice[n][i+1],0.0));
    }
    for(int i=n-2;i>=0;i--){
        for(int j=0;j<=i;j++){
            cprice_amer_call[i][j]=exp(-r*T/n)*(p*max(cprice_amer_call[i+1][j],sprice[i+1][j]-K)+(1-p)*max(cprice_amer_call[i+1][j+1],sprice[i+1][j+1]-K));
            cprice_amer_put[i][j]=exp(-r*T/n)*(p*max(cprice_amer_put[i+1][j],K-sprice[i+1][j])+(1-p)*max(cprice_amer_put[i+1][j+1],K-sprice[i+1][j+1]));
            cprice_euro_call[i][j]=exp(-r*T/n)*(p*cprice_euro_call[i+1][j]+(1-p)*cprice_euro_call[i+1][j+1]);
            cprice_euro_put[i][j]=exp(-r*T/n)*(p*cprice_euro_put[i+1][j]+(1-p)*cprice_euro_put[i+1][j+1]);
        }
    }
    cout<<"Euro Call: "<<cprice_euro_call[0][0]<<'\n';
    cout<<"Euro Put : "<<cprice_euro_put[0][0]<<'\n';
    cout<<"Amer Call: "<<cprice_amer_call[0][0]<<'\n';
    cout<<"Amer Put : "<<cprice_amer_put[0][0]<<'\n';

    //Bonus 1
    cout<<"\nBonus 1:\n";
    vector<double> cprice_euro_call_one(n+1,0),cprice_euro_put_one(n+1,0);
    vector<double> cprice_amer_call_one(n+1,0),cprice_amer_put_one(n+1,0);
    for(int i=0;i<=n;i++){
        cprice_euro_call_one[i]=max(S*pow(u,n-i)*pow(d,i)-K,0.0);
        cprice_euro_put_one[i]=max(K-S*pow(u,n-i)*pow(d,i),0.0);
        cprice_amer_call_one[i]=max(S*pow(u,n-i)*pow(d,i)-K,0.0);
        cprice_amer_put_one[i]=max(K-S*pow(u,n-i)*pow(d,i),0.0);
    }
    for(int i=n-1;i>=0;i--){
        for(int j=0;j<=i;j++){
            cprice_euro_call_one[j]=exp(-r*T/n)*(p*cprice_euro_call_one[j]+(1-p)*cprice_euro_call_one[j+1]);
            cprice_euro_put_one[j]=exp(-r*T/n)*(p*cprice_euro_put_one[j]+(1-p)*cprice_euro_put_one[j+1]);
            cprice_amer_call_one[j]=max(exp(-r*T/n)*(p*cprice_amer_call_one[j]+(1-p)*cprice_amer_call_one[j+1]),S*pow(u,i-j)*pow(d,j)-K);
            cprice_amer_put_one[j]=max(exp(-r*T/n)*(p*cprice_amer_put_one[j]+(1-p)*cprice_amer_put_one[j+1]),K-S*pow(u,i-j)*pow(d,j));
        }
    }
    cout<<"Euro Call: "<<cprice_euro_call_one[0]<<'\n';
    cout<<"Euro Put : "<<cprice_euro_put_one[0]<<'\n';
    cout<<"Amer Call: "<<cprice_amer_call_one[0]<<'\n';
    cout<<"Amer Put : "<<cprice_amer_put_one[0]<<'\n';

    //Bonus 2
    cout<<"\nBonus 2:\n";
    double cprice_euro_call_two=0,cprice_euro_put_two=0;;
    for(int j=0;j<=n;j++){
        cprice_euro_call_two+=binomial(n,j,p)*max(S*pow(u,n-j)*pow(d,j)-K,0.0);
        cprice_euro_put_two+=binomial(n,j,p)*max(K-S*pow(u,n-j)*pow(d,j),0.0);
    }
    cout<<"Euro Call: "<<cprice_euro_call_two*exp(-r*T)<<'\n';
    cout<<"Euro Put : "<<cprice_euro_put_two*exp(-r*T)<<'\n';

*/
}

