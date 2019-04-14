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

int main(int argc, char* argv[])
{
	ifstream input1;
	input1.open(argv[1]);
	if(!input1){ 
        return 1; 
    }
    const int nrolls=10000;  // number of experiments

    double S;
    double r;
    double q;
    double sigma;
    double T;
    double K1;
    double K2;
    double K3;
    double K4;

    input1>>S>>r>>q>>sigma>>T>>K1>>K2>>K3>>K4;

    double mean,var;
    mean=log(S)+(r-q-0.5*sigma*sigma)*T;
    var=sigma*sqrt(T);

    double test;

    std::random_device rd;
	std::default_random_engine generator(rd());
 	std::lognormal_distribution<double> distribution(mean,var);

    cout<<"S: "<<setw(5)<<S<<endl;
    cout<<"r: "<<setw(5)<<r<<endl;
    cout<<"q: "<<setw(5)<<q<<endl;
    cout<<"\u03C3: "<<setw(5)<<sigma<<endl;
    cout<<"T: "<<setw(5)<<T<<endl;
    cout<<"K1:"<<setw(5)<<K1<<endl;
    cout<<"K2:"<<setw(5)<<K2<<endl;
    cout<<"K3:"<<setw(5)<<K3<<endl;
    cout<<"K4:"<<setw(5)<<K4<<endl;

    double value=0;

    for(int i=0;i<nrolls;i++){
    	test=distribution(generator);
        if(test>=K1 && test<=K2)
            value+=(test-K1);
        else if(test>K2 && test<K3)
            value+=(K2-K1);
        else if(test>=K3 && test<=K4)
            value+=((K2-K1)-(test-K3)*(K2-K1)/(K4-K3));
    }

    double avgpresentvalue=value/nrolls*exp(-r*T);
    cout<<"Value:"<<setw(5)<<fixed<<setprecision(3)<<avgpresentvalue<<endl;

    vector <double> data;
    double total=0;
    int repeat_times=20;

    for(int k=0;k<repeat_times;k++){
        value=0;
        for(int i=0;i<nrolls;i++){
            test=distribution(generator);
            if(test>=K1 && test<=K2)
                value+=(test-K1);
            else if(test>K2 && test<K3)
                value+=(K2-K1);
            else if(test>=K3 && test<=K4)
                value+=((K2-K1)-(test-K3)*(K2-K1)/(K4-K3));
        }
        avgpresentvalue=value/nrolls*exp(-r*T);
        data.push_back(avgpresentvalue);
        total+=avgpresentvalue;
        //cout<<"No. "<<k<<'\t'<<"Value:"<<setw(5)<<avgpresentvalue<<endl;  
    }
    const double meanpresentvalue=total/repeat_times;

    double standardDeviation;
    for(int i=0;i<repeat_times;i++)
        standardDeviation += pow(data[i] - meanpresentvalue, 2);
    standardDeviation=sqrt(standardDeviation/repeat_times);

    cout<<"Under "<<repeat_times<<" times of testing, 95\% of confidence interval for the option value are: "<<endl;
    cout<<"Lower bounds: "<<meanpresentvalue-2*standardDeviation<<endl;
    cout<<"Mean:         "<<meanpresentvalue<<endl;
    cout<<"Upper bounds: "<<meanpresentvalue+2*standardDeviation<<endl;

}

