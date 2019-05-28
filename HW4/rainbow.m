function [price1,price2,price3] = rainbow(nrolls,n,K,r,T,S,q,sigma,A,cases)
    %RAINBOW Summary of this function goes here
    %   Detailed explanation goes here
    zmat1 = randn(nrolls,n);
    zmat2 = zmat1;

    %Bonus 1
    for i=1:nrolls/2
        for j=1:n
            zmat2(i+nrolls/2,j)=-zmat2(i,j);
        end
    end
    %Bonus 2
    zmat3=zmat2;
    %Bonus 1 
    for j=1:n
        std2=std(zmat2(:,j));
        for i=1:nrolls
            zmat2(i,j)=zmat2(i,j)/std2;
        end
    end
    %Bonus 2
    hmat3=zeros(n,n);
    for i=1:n
        for j=1:n
            s1=cov(zmat3(:,i),zmat3(:,j));
            hmat3(i,j)=s1(1,2);
        end
    end
    [Astar]=choles(hmat3,n);
    A3=Astar\A;


    rmat1 = zmat1*A;
    rmat2 = zmat2*A;
    rmat3 = zmat3*A3;

    for i=1:n
        mu=log(S(i))+(r-q(i)-sigma(i)*sigma(i)/2)*T;
        for j=1:nrolls
            rmat1(j,i)=rmat1(j,i)+mu;
            rmat1(j,i)=exp(rmat1(j,i));
            rmat2(j,i)=rmat2(j,i)+mu;
            rmat2(j,i)=exp(rmat2(j,i));
            rmat3(j,i)=rmat3(j,i)+mu;
            rmat3(j,i)=exp(rmat3(j,i));
        end
    end
    
    price1=0;
    price2=0;
    price3=0;
    for i=1:nrolls
        maxi1=-10000;
        maxi2=-10000;
        maxi3=-10000;
        for j=1:n
            if(rmat1(i,j)>maxi1) 
                maxi1=rmat1(i,j);
            end
            if(rmat2(i,j)>maxi2) 
                maxi2=rmat2(i,j);
            end
            if(rmat3(i,j)>maxi3) 
                maxi3=rmat3(i,j);
            end
        end
        if((maxi1-K)>0) 
            price1=price1+(maxi1-K);
        end
        if((maxi2-K)>0) 
            price2=price2+(maxi2-K);
        end
        if((maxi3-K)>0) 
            price3=price3+(maxi3-K);
        end
    end
    price1=price1/nrolls*exp(-r*T);
    price2=price2/nrolls*exp(-r*T);
    price3=price3/nrolls*exp(-r*T);
    
end

