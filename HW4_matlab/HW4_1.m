
[St,r,q,sigma,t,T,Smax_t,n,nrolls,num_of_rep]=readdata('input_2.txt');

% Binomial
dT=(T-t)/n;
u=exp(sigma*sqrt(T/n));
d=1/u;
p=(exp(r*T/n-q*T/n)-d)/(u-d);
pricelist(1:2*n+1)=St;
for i=1:n
   pricelist(i)=St*power(d,n-i+1);
   pricelist(2*n-i+2)=St*power(u,n-i+1);
end

s=struct('name',[],'dollar',10,'i',[],'j',[]);
s1(1:10,1:10)=struct(s);


% Monte Carlo
price(1:num_of_rep)=nan;
for i=1:num_of_rep
    price(i)=MonteCarlo(St,r,q,sigma,t,T,Smax_t,n,nrolls);
end

fprintf('Basic requirement(ii): \nStandard Deviation: %f\n',std(price));
fprintf('Under %d times of testing, 95%% of confidence interval for the call option value are: \n',num_of_rep);
fprintf('Lower bounds: %f\n',mean(price)-2*std(price));
fprintf('Mean          %f\n',mean(price));
fprintf('Upper bounds: %f\n',mean(price)+2*std(price));

% Bonus 2
bonus2_euro=Binomial_Cheuk(St,r,q,sigma,t,T,n,'EURO');
bonus2_amer=Binomial_Cheuk(St,r,q,sigma,t,T,n,'AMER');
fprintf('\nBonus 2:\n');
fprintf('European lookback puts %f dollars.\n',bonus2_euro);
fprintf('American lookback puts %f dollars.\n',bonus2_amer);

