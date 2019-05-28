function [Explicit_Price] = Explicit(S0,K,r,q,sigma,T,Smin,Smax,m,n,optType,AMER)
    %EXLICIT Summary of this function goes here
    %   Detailed explanation goes here
    if(nargin==11)
        AMER=false;
    end
    %tic;
    
    dS=(Smax-Smin)/m;
    dT=T/n;
    %Boundary condition for price
    price(1:m+1,1:n+1)=nan;
    switch optType
        case 'CALL'
            price(1,:)=m*dS-K;
            price(end,:)=max(Smin-K,0);
            %Stock price in time T
            for i=1:m+1
                price(i,n+1)=max((m-i+1)*dS-K,0);
            end
        case 'PUT'
            price(1,:)=max(K-m*dS,0);
            price(end,:)=K-Smin;
            %Stock price in time T
            for i=1:m+1
                price(i,n+1)=max(K-(m-i+1)*dS,0);
            end
    end

    % coef of aj, bj, cj
    j=1:m-1;
    aj=(-(r-q)/2*j*dT+0.5*power(sigma*j,2)*dT)/(1+r*dT);
    bj=(1-power(sigma*j,2)*dT)/(1+r*dT);
    cj=((r-q)/2*j*dT+0.5*power(sigma*j,2)*dT)/(1+r*dT);

    % A matrix
    expmat=zeros(m-1,m+1);
    for i=1:m-1
        expmat(i,i)=cj(1,m-i);
        expmat(i,i+1)=bj(1,m-i);
        expmat(i,i+2)=aj(1,m-i);
    end

    for idx=n:-1:1
        price(2:m,idx)=expmat*price(:,idx+1);
        for idy=2:m
            if(price(idy,idx)<0)
                price(idy,idx)=0;
            end
            if(AMER)
                price(idy,idx)=max(price(idy,idx),price(idy,n+1));
            end
        end
    end
    %pre_t=toc;
    %fprintf('Compute %f (s)\n',pre_t);

    Explicit_Price=interp1(Smax:-dS:Smin,price(:,1),S0);
end