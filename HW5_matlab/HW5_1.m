clear

[St,K,r,q,sigma,T_t,M,n,Save_t,passing_time,nrolls,num_of_rep]=readdata('input_1_2.txt');

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

fprintf('Basic requirement(ii):\n');
fprintf('Monte Carlo: \nStandard Deviation: %f\n',std(monte_euro));
fprintf('Under %d times of testing, 95%% of confidence interval for the call option value are: \n',num_of_rep);
fprintf('Lower bounds: %f\n',mean(monte_euro)-2*std(monte_euro));
fprintf('Mean          %f\n',mean(monte_euro));
fprintf('Upper bounds: %f\n\n',mean(monte_euro)+2*std(monte_euro));


% Bonus 1 log
tic;
[binolog_euro,binolog_amer]=Binomial_log(St,K,r,q,sigma,T_t,M,n,passing_time);
tt=toc;
fprintf('Bonus 1:\n');
fprintf('Binomial euro arithmetic average call : %f\n',binolog_euro);
fprintf('Binomial euro arithmetic average call : %f\n',binolog_amer);
fprintf('Time %f(s)\n\n',tt);


eurodata=nan(3,10);
amerdata=nan(3,10);
for M1=50:50:500
    [bino_euro,bino_amer]=Binomial(St,K,r,q,sigma,T_t,M1,n,passing_time);
    [binolog_euro,binolog_amer]=Binomial_log(St,K,r,q,sigma,T_t,M1,n,passing_time);
    eurodata(1,M1/50)=M1;
    eurodata(2,M1/50)=bino_euro;
    eurodata(3,M1/50)=binolog_euro;
    amerdata(1,M1/50)=M1;
    amerdata(2,M1/50)=bino_amer;
    amerdata(3,M1/50)=binolog_amer;
    M1
end
figure(1);
subplot(2,1,1);
plot(eurodata(1,:),eurodata(2,:),'magenta',eurodata(1,:),eurodata(3,:),'blue');
title('Euro price');
xlabel('M');
ylabel('price');

subplot(2,1,2);
plot(amerdata(1,:),amerdata(2,:),'magenta',amerdata(1,:),amerdata(3,:),'blue');
title('Amer price');
xlabel('M');
ylabel('price');
