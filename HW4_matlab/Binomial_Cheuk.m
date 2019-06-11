function [price] = Binomial_Cheuck(St,r,q,sigma,t,T,n,optType)
%BINOMIAL_CHEUK Summary of this function goes here
%   Detailed explanation goes here
    dT=(T-t)/n;
    u=exp(sigma*sqrt(dT));
    d=1/u;
    mu=exp(r*dT-q*dT);
    pu=(mu*u-1)/(mu*(u-d));
    pd=1-pu;
    tree(1:n+1,1:n+1)=nan;
    for i=1:n+1
        tree(i,n+1)=max(power(u,n-i+1)-1,0);
    end

    for i=n:-1:1
        for j=n+2-i:n
            tree(j,i)=(pd*tree(j-1,i+1)+pu*tree(j+1,i+1))*exp(-r*dT)*mu;
            if(optType=='AMER')
                tree(j,i)=max(tree(j,i),tree(j,n+1));
            end
        end
        tree(n+1,i)=(pd*tree(n,i+1)+pu*tree(n+1,i+1))*exp(-r*dT)*mu;
    end
    price=tree(n+1,1)*St;
end


