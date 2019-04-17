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
double normalCDF(double x) // Phi(-∞, x) aka N(x)
{
    return erfc(-x/sqrt(2))/2;
}


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

    cout<<"S: "<<setw(5)<<S<<endl;
    cout<<"r: "<<setw(5)<<r<<endl;
    cout<<"q: "<<setw(5)<<q<<endl;
    cout<<"\u03C3: "<<setw(5)<<sigma<<endl;
    cout<<"T: "<<setw(5)<<T<<endl;
    cout<<"K1:"<<setw(5)<<K1<<endl;
    cout<<"K2:"<<setw(5)<<K2<<endl;
    cout<<"K3:"<<setw(5)<<K3<<endl;
    cout<<"K4:"<<setw(5)<<K4<<endl;

    double var1=(log(S/K1)+(r-q-sigma*sigma/2)*T)/sigma/sqrt(T);
    double var2=(log(S/K2)+(r-q-sigma*sigma/2)*T)/sigma/sqrt(T);
    double var3=(log(S/K3)+(r-q-sigma*sigma/2)*T)/sigma/sqrt(T);
    double var4=(log(S/K4)+(r-q-sigma*sigma/2)*T)/sigma/sqrt(T);
    double var5=(log(S/K1)+(r-q+sigma*sigma/2)*T)/sigma/sqrt(T);
    double var6=(log(S/K2)+(r-q+sigma*sigma/2)*T)/sigma/sqrt(T);
    double var7=(log(S/K3)+(r-q+sigma*sigma/2)*T)/sigma/sqrt(T);
    double var8=(log(S/K4)+(r-q+sigma*sigma/2)*T)/sigma/sqrt(T);



    double avgpresentvalue=exp(-r*T)*(-K1*normalCDF(var1)+K2*normalCDF(var2)+K3*(K2-K1)/(K4-K3)*normalCDF(var3)\
           -((K2-K1)+K3*(K2-K1)/(K4-K3))*normalCDF(var4)+S*exp(r*T-q*T)*(normalCDF(var5)-normalCDF(var6)\
           -(K2-K1)/(K4-K3)*normalCDF(var7)+(K2-K1)/(K4-K3)*normalCDF(var8)));


    cout<<"Value:"<<setw(5)<<fixed<<setprecision(5)<<avgpresentvalue<<endl;


}

