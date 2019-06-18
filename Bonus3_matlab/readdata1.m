function [S0,K,r,q,sigma,T,n,nrolls,num_of_rep]=readdata1(path)
    %READDATA1 Summary of this function goes here
    %   Detailed explanation goes here
    fileID = fopen(path,'r');
    data = fscanf(fileID,'%f');
    S0=data(1);
    K=data(2);
    r=data(3);
    q=data(4);
    sigma=data(5);
    T=data(6);
    n=data(7);
    nrolls=data(8);
    num_of_rep=data(9);
    fclose(fileID);
end

