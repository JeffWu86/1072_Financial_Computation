function [price] = MonteCarlo(St,r,q,sigma,t,T,Smax_t,n,nrolls)
%MONTECARLO Summary of this function goes here
%   Detailed explanation goes here
    price(1:10000,1:n+1)=nan;
    price(1:10000,1)=St;

    dT=(T-t)/n;
    
    for i=1:n
        for j=1:nrolls
            mu=log(price(j,i))+(r-q-sigma*sigma/2)*dT;
            sigma1=sigma*sqrt(dT);
            price(j,i+1)=exp(sigma1*randn()+mu);
        end
    end
    
    out(1:10000)=nan;
    for i=1:nrolls
        out(i)=max( max(max(price(i,:)),Smax_t) - price(i,n+1), 0);
    end
%     This is for call with K=Smax_t (Check)    
%     for i=1:nrolls
%         out(i)=max(price(i,n+1)-Smax_t,0);
%     end
%     mean(zmat1(:,n+1))
    price=mean(out)*exp(-r*T);
end

