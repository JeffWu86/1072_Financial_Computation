function [monte_price] = MonteCarlo(St,K,r,q,sigma,T_t,n,Save_t,passing_time,nrolls)
%MONTECARLO Summary of this function goes here
%   Detailed explanation goes here
    dT=(T_t)/n;
    npass=n*passing_time/T_t;
    
    sim(1:nrolls,1:n+1)=nan;
    sim(1:nrolls,1)=St;

    for i=1:n
        for j=1:nrolls
            mu=log(sim(j,i))+(r-q-sigma*sigma/2)*dT;
            sigma1=sigma*sqrt(dT);
            sim(j,i+1)=exp(sigma1*randn()+mu);
        end
    end
    
    price(1:nrolls,1)=nan;
    for i=1:nrolls
        price(i)= max( 0, (sum(sim(i,:))+(npass+1)*Save_t-St)/(npass+n+1) -K );
    end
    
    monte_price=mean(price)*exp(-r*T_t);
end

