
[K,r,T,nrolls,repeat_times,n,S,q,sigma,rho,covmat]=readdata('input_test3.txt');
%chol(covmat)
[A]=choles(covmat,n);

basic1=zeros(repeat_times,1);
bonus1=zeros(repeat_times,1);
bonus2=zeros(repeat_times,1);
for i=1:repeat_times
    [price1,price2,price3]=rainbow(nrolls,n,K,r,T,S,q,sigma,A);
    basic1(i)=price1;
    bonus1(i)=price2;
    bonus2(i)=price3;
end
fprintf('Basic requirement: \nStandard Deviation: %f\n',std(basic1));
fprintf('Under %d times of testing, 95%% of confidence interval for the call option value are: \n',repeat_times);
fprintf('Lower bounds: %f\n',mean(basic1)-2*std(basic1));
fprintf('Mean          %f\n',mean(basic1));
fprintf('Upper bounds: %f\n',mean(basic1)+2*std(basic1));
fprintf('\n');
fprintf('Bonus1: \nStandard Deviation: %f\n',std(bonus1));
fprintf('Under %d times of testing, 95%% of confidence interval for the call option value are: \n',repeat_times);
fprintf('Lower bounds: %f\n',mean(bonus1)-2*std(bonus1));
fprintf('Mean          %f\n',mean(bonus1));
fprintf('Upper bounds: %f\n',mean(bonus1)+2*std(bonus1));
fprintf('Bonus2: \nStandard Deviation: %f\n',std(bonus2));
fprintf('\n');
fprintf('Under %d times of testing, 95%% of confidence interval for the call option value are: \n',repeat_times);
fprintf('Lower bounds: %f\n',mean(bonus2)-2*std(bonus2));
fprintf('Mean          %f\n',mean(bonus2));
fprintf('Upper bounds: %f\n',mean(bonus2)+2*std(bonus2));