#include<iostream>
#include<stdio.h>
#include<iomanip>
#include<cmath>
#include<fstream>
using namespace std;
#include<ctime>
#include<random>
#include<algorithm>
double normalCDF(double x) // Phi(-∞, x) aka N(x)
{
    return erfc(-x/sqrt(2))/2;
}
double bin(int n, int m)
{
if(m == 0 || n == m)
return 1;
else 
return (bin(n-1, m) + bin(n-1, m-1));
}
double binom(int n, int k)
{
    double value = 1;
    // need to be careful here - can't just use *= due to integer division
    for (int i = 0; i < k; i++)
        value = (value * (n - i)) / (i + 1);
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

    double mean,var;
    mean=log(S)+(r-q-0.5*sigma*sigma)*T;
    var=sigma*sqrt(T);

    double test;

    std::random_device rd;
	std::default_random_engine generator(rd());
 	std::lognormal_distribution<double> distribution(mean,var);


    cout<<"S:"<<setw(12)<<S<<'\n';
    cout<<"K:"<<setw(12)<<K<<'\n';
    cout<<"r:"<<setw(12)<<r<<'\n';
    cout<<"q:"<<setw(12)<<q<<'\n';
    cout<<"\u03C3:"<<setw(12)<<sigma<<'\n';
    cout<<"T:"<<setw(12)<<T<<'\n';
    cout<<"# of sim:"<<setw(5)<<nrolls<<'\n';
    cout<<"# of rep:"<<setw(5)<<repeat_times<<'\n';
    cout<<"n:"<<setw(12)<<n<<'\n';

    //Black-Scholes formulas
    cout<<"\nBlack-Scholes formulas:\n";
    double d1=(log(S/K)+(r-q+sigma*sigma/2)*T)/sigma/sqrt(T);
    double d2=(log(S/K)+(r-q-sigma*sigma/2)*T)/sigma/sqrt(T);

    cout<<"Call:"<<setw(9)<<S*exp(-q*T)*normalCDF(d1)-K*exp(-r*T)*normalCDF(d2)<<'\n';
    cout<<"Put:"<<setw(10)<<K*exp(-r*T)*normalCDF(-d2)-S*exp(-q*T)*normalCDF(-d1)<<'\n';

    //Monte Carlo Simulation
    cout<<"\nMonte Carlo simulation:\n";
    double Call=0,Put=0;

    for(int i=0;i<nrolls;i++){
    	test=distribution(generator);
    	Call+=max(test-K,0.0);
    	Put+=max(K-test,0.0);
    }

    double avgpresentcall=Call/nrolls*exp(-r*T);
    double avgpresentput=Put/nrolls*exp(-r*T);
    cout<<setiosflags(ios::fixed);
    cout<<"Call:"<<setw(9)<<avgpresentcall<<'\n';
    cout<<"Put:"<<setw(10)<<avgpresentput<<'\n';


    vector <double> data,data1;
    double total=0;
    double avgpresentvalue;
    double value=0;

    for(int k=0;k<repeat_times;k++){
        value=0;
        for(int i=0;i<nrolls;i++){
            test=distribution(generator);
            value+=max(test-K,0.0);
        }
        avgpresentvalue=value/nrolls*exp(-r*T);
        data.push_back(avgpresentvalue);
        total+=avgpresentvalue;
    }
    double meanpresentvalue=total/repeat_times;

    double standardDeviation;
    for(int i=0;i<repeat_times;i++)
        standardDeviation+=(data[i]-meanpresentvalue)*(data[i]-meanpresentvalue);
    standardDeviation=sqrt(standardDeviation/repeat_times);

    cout<<"Under "<<repeat_times<<" times of testing, 95\% of confidence interval for the call option value are: "<<'\n';
    cout<<"Lower bounds: "<<meanpresentvalue-2*standardDeviation<<'\n';
    cout<<"Mean:         "<<meanpresentvalue<<'\n';
    cout<<"Upper bounds: "<<meanpresentvalue+2*standardDeviation<<'\n';

    total=0;
    for(int k=0;k<repeat_times;k++){
        value=0;
        for(int i=0;i<nrolls;i++){
            test=distribution(generator);
            value+=max(K-test,0.0);
        }
        avgpresentvalue=value/nrolls*exp(-r*T);
        data1.push_back(avgpresentvalue);
        total+=avgpresentvalue;
    }
    meanpresentvalue=total/repeat_times;
    standardDeviation=0;
    for(int i=0;i<repeat_times;i++)
        standardDeviation+=(data1[i]-meanpresentvalue)*(data1[i]-meanpresentvalue);
    standardDeviation=sqrt(standardDeviation/repeat_times);
    cout<<"standardDeviation"<<standardDeviation<<'\n';
    cout<<"var"<<var<<'\n';
    cout<<"Under "<<repeat_times<<" times of testing, 95\% of confidence interval for the put option value are: "<<'\n';
    cout<<"Lower bounds: "<<meanpresentvalue-2*standardDeviation<<'\n';
    cout<<"Mean:         "<<meanpresentvalue<<'\n';
    cout<<"Upper bounds: "<<meanpresentvalue+2*standardDeviation<<'\n';
    

    //CRR binomial tree model
    cout<<"\nCRR binomial tree model:\n";
    double u=exp(sigma*sqrt(T/n));
    double d=1/u;
    double p=(exp(r*T/n-q*T/n)-d)/(u-d);
    double sprice[n+1][n+1];
    for(int i=0;i<=n;i++){
        for(int j=0;j<=i;j++){
            sprice[i][j]=S*pow(u,i-j)*pow(d,j);
        }
    }
    double cprice_amer_call[n][n],cprice_euro_call[n][n];
    double cprice_amer_put[n][n],cprice_euro_put[n][n];
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
    double cprice_euro_call_one[n+1],cprice_euro_put_one[n+1];
    double cprice_amer_call_one[n+1],cprice_amer_put_one[n+1];
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
        cprice_euro_call_two+=binom(n,j)*pow(p,n-j)*pow((1-p),j)*max(S*pow(u,n-j)*pow(d,j)-K,0.0);
        cprice_euro_put_two+=binom(n,j)*pow(p,n-j)*pow((1-p),j)*max(K-S*pow(u,n-j)*pow(d,j),0.0);
    }
    cout<<"Euro Call: "<<cprice_euro_call_two*exp(-r*T)<<'\n';
    cout<<"Euro Put : "<<cprice_euro_put_two*exp(-r*T)<<'\n';

    cout<<"C(100 50):\n";
    cout<<setprecision(0)<<binom(100,50)<<'\n';


}

