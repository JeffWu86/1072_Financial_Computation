clear
% Plain Vanilla
[S0,K,r,q,sigma,T,n,nrolls,repeat_times]=readdata1('input_1_4.txt');
PlainVanilla=zeros(repeat_times,1);
for i=1:repeat_times
    [price]=Plain_Vanilla(S0,K,r,q,sigma,T,n,nrolls);
    PlainVanilla(i)=price;
end

fprintf('Plain Vanilla: \nStandard Deviation: %f\n',std(PlainVanilla));
fprintf('Under %d times of testing, 95%% of confidence interval for the call option value are: \n',repeat_times);
fprintf('Lower bounds: %f\n',mean(PlainVanilla)-2*std(PlainVanilla));
fprintf('Mean          %f\n',mean(PlainVanilla));
fprintf('Upper bounds: %f\n',mean(PlainVanilla)+2*std(PlainVanilla));
fprintf('\n');