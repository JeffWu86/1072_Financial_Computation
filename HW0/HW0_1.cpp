#include <iostream>
#include <stdio.h>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <ctime>
#include <random>
#include <vector>
#include <algorithm>
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
    double K;
    double r;
    double sigma;
    double T;

    input1>>S>>K>>r>>sigma>>T;

    double mean,var;
    mean=log(S)+(r-0.5*sigma*sigma)*T;
    var=sigma*sqrt(T);

    double test;

    std::random_device rd;
	std::default_random_engine generator(rd());
 	std::lognormal_distribution<double> distribution(mean,var);


    cout<<"S:"<<setw(5)<<S<<endl;
    cout<<"K:"<<setw(5)<<K<<endl;
    cout<<"t:"<<setw(5)<<r<<endl;
    cout<<"\u03C3:"<<setw(5)<<sigma<<endl;
    cout<<"T:"<<setw(5)<<T<<endl;

    double Call=0,Put=0;

    for(int i=0;i<nrolls;i++){
    	test=distribution(generator);
    	Call+=max(test-K,0.0);
    	Put+=max(K-test,0.0);
    }

    double avgpresentcall=Call/nrolls*exp(-r*T);
    double avgpresentput=Put/nrolls*exp(-r*T);
    cout<<"Call:"<<setw(5)<<avgpresentcall<<endl;
    cout<<"Put :"<<setw(5)<<avgpresentput<<endl;

}

