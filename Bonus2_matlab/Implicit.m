function [Implicit_Price] = Implicit(S0,K,r,q,sigma,T,Smin,Smax,m,n,optType,AMER)
    %IMPLICIT Summary of this function goes here
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
    aj=(r-q)/2*j*dT-0.5*power(sigma*j,2)*dT;
    bj=1+power(sigma*j,2)*dT+r*dT;
    cj=-(r-q)/2*j*dT-0.5*power(sigma*j,2)*dT;
    
    % A matrix
    impmat=zeros(m-1,m-1);
    for i=1:m-1
        impmat(i,i)=bj(1,m-i);
        if(i~=1)
            impmat(i,i-1)=cj(1,m-i);
        end
        if(i~=m-1)
            impmat(i,i+1)=aj(1,m-i);
        end
    end

    [L,U] = lu(impmat);
    for idx=n:-1:1
        b=price(2:m,idx+1);
        b(end)=b(end)-aj(1)*price(m+1,idx); %Alwayw zero
        b(1)=b(1)-cj(end)*price(1,idx);
        %x=impmat\b;
        x=U\(L\b);
        for j=1:m-1
            if(x(j,1)<0) 
                x(j,1)=0;
            end
        end
        price(2:m,idx)=x;
        for idy=2:m
            if(AMER)
                price(idy,idx)=max(price(idy,idx),price(idy,n+1));
            end
        end
    end
    %pre_t=toc;
    %fprintf('Compute %f (s)\n',pre_t);

    Implicit_Price=interp1(Smax:-dS:Smin,price(:,1),S0);
end

