clear
[St,r,q,sigma,t,T,Smax_t,n,nrolls,num_of_rep]=readdata('input_2.txt');

% Basic 
% Binomial
[bino_euro,bino_amer]=Binomial(St,r,q,sigma,t,T,n);
fprintf('Basic requirement(i):\n');
fprintf('Binomial euro lookback options : %f\n',bino_euro);
fprintf('Binomial amer lookback options : %f\n',bino_amer);

% Monte Carlo
monte_euro(1:num_of_rep)=nan;
for i=1:num_of_rep
    monte_euro(i)=MonteCarlo(St,r,q,sigma,t,T,Smax_t,n,nrolls);
end

fprintf('\nBasic requirement(ii): \nStandard Deviation: %f\n',std(monte_euro));
fprintf('Under %d times of testing, 95%% of confidence interval for the call option value are: \n',num_of_rep);
fprintf('Lower bounds: %f\n',mean(monte_euro)-2*std(monte_euro));
fprintf('Mean          %f\n',mean(monte_euro));
fprintf('Upper bounds: %f\n',mean(monte_euro)+2*std(monte_euro));

% Bonus 2
bonus2_euro=Binomial_Cheuk(St,r,q,sigma,t,T,n,'EURO');
bonus2_amer=Binomial_Cheuk(St,r,q,sigma,t,T,n,'AMER');
fprintf('\nBonus 2:\n');
fprintf('European lookback puts %f dollars.\n',bonus2_euro);
fprintf('American lookback puts %f dollars.\n',bonus2_amer);

