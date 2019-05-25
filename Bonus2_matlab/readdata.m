function [K,r,T,nrolls,repeat_times,n,S,q,sigma,rho,covmat]=readdata(path)
    %READDATA Summary of this function goes here
    %   Detailed explanation goes here
    fileID = fopen(path,'r');
    data = fscanf(fileID,'%f');
    K=data(1);
    r=data(2);
    T=data(3);
    nrolls=data(4);
    repeat_times=data(5);
    n=data(6);
    count=6;
    S=zeros(n);
    q=zeros(n);
    sigma=zeros(n);
    for i = 1:n
        S(i)=data(count+i);
        q(i)=data(count+i+n);
        sigma(i)=data(count+i+2*n);
    end
    count=count+3*n+1;
    rho=ones(n,n);
    for i = 1:n-1
        for j = i+1:n
            rho(i,j)=data(count);
            rho(j,i)=data(count);
            count=count+1;
        end
    end
    covmat=zeros(n,n);
    for i = 1:n
        for j = 1:n
            covmat(i,j)=rho(i,j)*sigma(i)*sigma(j)*T;
        end
    end
    fclose(fileID);
end

