clear

[St,K,r,q,sigma,T_t,M,n,Save_t,passing_time,nrolls,num_of_rep]=readdata('input_1_1.txt');

% Binomial
tic;
[bino_euro,bino_amer]=Binomial(St,K,r,q,sigma,T_t,M,n,passing_time);
tt=toc;
fprintf('Basic requirement(i):\n');
fprintf('Binomial euro arithmetic average call : %f\n',bino_euro);
fprintf('Binomial euro arithmetic average call : %f\n',bino_amer);
fprintf('Time %f(s)\n\n',tt);

% Monte Carlo
monte_euro(1:num_of_rep)=nan;
for i=1:num_of_rep
    monte_euro(i)=MonteCarlo(St,K,r,q,sigma,T_t,n,Save_t,passing_time,nrolls);
end

fprintf('Monte Carlo: \nStandard Deviation: %f\n',std(monte_euro));
fprintf('Under %d times of testing, 95%% of confidence interval for the call option value are: \n',num_of_rep);
fprintf('Lower bounds: %f\n',mean(monte_euro)-2*std(monte_euro));
fprintf('Mean          %f\n',mean(monte_euro));
fprintf('Upper bounds: %f\n',mean(monte_euro)+2*std(monte_euro));