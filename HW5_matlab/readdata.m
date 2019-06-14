function [St,K,r,q,sigma,T_t,M,n,Save_t,passing_time,nrolls,num_of_rep]=readdata(path)
    %READDATA Summary of this function goes here
    %   Detailed explanation goes here
    fileID = fopen(path,'r');
    data = fscanf(fileID,'%f');
    St=data(1);
    K=data(2);
    r=data(3);
    q=data(4);
    sigma=data(5);
    T_t=data(6);
    M=data(7);
    n=data(8);
    Save_t=data(9);
    passing_time=data(10);
    nrolls=data(11);
    num_of_rep=data(12);
    
    fclose(fileID);
end

