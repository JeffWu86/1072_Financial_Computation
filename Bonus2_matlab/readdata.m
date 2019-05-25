function [S0,K,r,q,sigma,T,Smin,Smax,m,n]=readdata(path)
    %READDATA Summary of this function goes here
    %   Detailed explanation goes here
    fileID = fopen(path,'r');
    data = fscanf(fileID,'%f');
    S0=data(1);
    K=data(2);
    r=data(3);
    q=data(4);
    sigma=data(5);
    T=data(6);
    Smin=data(7);
    Smax=data(8);
    m=data(9);
    n=data(10);
    
    fclose(fileID);
end

