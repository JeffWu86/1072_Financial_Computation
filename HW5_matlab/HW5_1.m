clear

[St,K,r,q,sigma,T_t,M,n,Save_t,passing_time,nrolls,num_of_rep]=readdata('input_1_4.txt');

% Binomial
tic;
[bino_euro,bino_amer]=Binomial(St,K,r,q,sigma,T_t,M,n,Save_t,passing_time);  
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

fprintf('Basic requirement(ii):\n');
fprintf('Monte Carlo: \nStandard Deviation: %f\n',std(monte_euro));
fprintf('Under %d times of testing, 95%% of confidence interval for the call option value are: \n',num_of_rep);
fprintf('Lower bounds: %f\n',mean(monte_euro)-2*std(monte_euro));
fprintf('Mean          %f\n',mean(monte_euro));
fprintf('Upper bounds: %f\n\n',mean(monte_euro)+2*std(monte_euro));


% Bonus 1 log
tic;
[binolog_euro,binolog_amer]=Binomial_log(St,K,r,q,sigma,T_t,M,n,Save_t,passing_time);
tt=toc;
fprintf('Bonus 1:\n');
fprintf('Binomial euro arithmetic average call : %f\n',binolog_euro);
fprintf('Binomial euro arithmetic average call : %f\n',binolog_amer);
fprintf('Time %f(s)\n\n',tt);

start=5;
jump=5;
nn=10;
eurodata=nan(3,nn);
amerdata=nan(3,nn);
for M1=start:jump:(start+jump*(nn-1))
    tic;
    [bino_euro,bino_amer]=Binomial(St,K,r,q,sigma,T_t,M1,n,Save_t,passing_time);
    [binolog_euro,binolog_amer]=Binomial_log(St,K,r,q,sigma,T_t,M1,n,Save_t,passing_time);
    tt=toc;
    eurodata(1,M1/jump)=M1;
    eurodata(2,M1/jump)=bino_euro;
    eurodata(3,M1/jump)=binolog_euro;
    amerdata(1,M1/jump)=M1;
    amerdata(2,M1/jump)=bino_amer;
    amerdata(3,M1/jump)=binolog_amer;
    fprintf('M1 is now %d , compute time: %f(s)\n',M1,tt);
end
figure(1);
subplot(2,1,1);
plot(eurodata(1,:),eurodata(2,:),'magenta',eurodata(1,:),eurodata(3,:),'blue');
title('Euro price');
xlabel('M');
ylabel('price');
legend('linear','log');  

subplot(2,1,2);
plot(amerdata(1,:),amerdata(2,:),'magenta',amerdata(1,:),amerdata(3,:),'blue');
title('Amer price');
xlabel('M');
ylabel('price');
legend('linear','log');  
