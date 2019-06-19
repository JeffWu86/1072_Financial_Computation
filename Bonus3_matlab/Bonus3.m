clear
% Plain Vanilla
[S0,K,r,q,sigma,T,n,nrolls,num_of_rep]=readdata1('input_1_6.txt');
PlainVanilla=zeros(num_of_rep,1);
for i=1:num_of_rep
    [price]=Plain_Vanilla(S0,K,r,q,sigma,T,n,nrolls);
    PlainVanilla(i)=price;
end

fprintf('Plain Vanilla: \nStandard Deviation: %f\n',std(PlainVanilla));
fprintf('Under %d times of testing, 95%% of confidence interval for the call option value are: \n',num_of_rep);
fprintf('Lower bounds: %f\n',mean(PlainVanilla)-2*std(PlainVanilla));
fprintf('Mean          %f\n',mean(PlainVanilla));
fprintf('Upper bounds: %f\n',mean(PlainVanilla)+2*std(PlainVanilla));
fprintf('\n');
clear;

% Lookback Put
[St,r,q,sigma,t,T,n,Smax_t,nrolls,num_of_rep]=readdata2('input_2_7.txt');
LookbackPuts=zeros(num_of_rep,1);
for i=1:num_of_rep
    [price,pay]=LookbackPut(St,r,q,sigma,t,T,n,Smax_t,nrolls);
    LookbackPuts(i)=price;
end
fprintf('Lookback Put: \nStandard Deviation: %f\n',std(LookbackPuts));
fprintf('Under %d times of testing, 95%% of confidence interval for the call option value are: \n',num_of_rep);
fprintf('Lower bounds: %f\n',mean(LookbackPuts)-2*std(LookbackPuts));
fprintf('Mean          %f\n',mean(LookbackPuts));
fprintf('Upper bounds: %f\n',mean(LookbackPuts)+2*std(LookbackPuts));
fprintf('\n');