function [St,r,q,sigma,t,T,Smax_t,n,nrolls,num_of_rep]=readdata(path)
    %READDATA Summary of this function goes here
    %   Detailed explanation goes here
    fileID = fopen(path,'r');
    data = fscanf(fileID,'%f');
    St=data(1);
    r=data(2);
    q=data(3);
    sigma=data(4);
    t=data(5);
    T=data(6);
    Smax_t=data(7);
    % n number of tree
    n=data(8);
    nrolls=data(9);
    num_of_rep=data(10);
    
    fclose(fileID);
end

