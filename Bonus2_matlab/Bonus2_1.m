
[S0,K,r,q,sigma,T,Smin,Smax,m,n]=readdata('input_1.txt');
tic;
%Implicit methond
[Im_AC, Im_AP, Im_EC, Im_EP] = Implicit(S0,K,r,q,sigma,T,Smin,Smax,m,n);

deltaS=(Smax-Smin)/m;
deltaT=T/n;
price=zeros(m+1,n+1);
for i=1:n+1
    price(1,i)=m*deltaS-K;
    price(n+1,i)=max(Smin-K,0);
end
for i=1:m+1
    price(i,n+1)=max((m-i+1)*deltaS-K,0);
end
avec=zeros(1,m-1);
bvec=zeros(1,m-1);
cvec=zeros(1,m-1);
for i=1:m-1
    avec(1,i)=(r-q)/2*i*deltaT-0.5*power(sigma*i,2)*deltaT;
    bvec(1,i)=1+power(sigma*i,2)*deltaT+r*deltaT;
    cvec(1,i)=-(r-q)/2*i*deltaT-0.5*power(sigma*i,2)*deltaT;
end
impmat=zeros(m-1,m-1);
for i=1:m-1
    impmat(i,i)=bvec(1,i);
    if(i~=1)
        impmat(i,i-1)=cvec(1,i);
    end
    if(i~=m-1)
        impmat(i,i+1)=avec(1,i);
    end
end
invimpmat=inv(impmat);
for i=n:-1:1
    b=price(2:m,i+1);
    b(m-1,1)=b(m-1,1)-avec(1,1)*price(m+1,i);
    b(1,1)=b(1,1)-cvec(1,m-1)*price(1,i);
    %x=impmat\b;
    x=invimpmat*b;
    for j=1:m-1
        if(b(j,1)<0) 
            b(j,1)=0;
        end
    end
    price(2:m,i)=x;
end
pre_t=toc;
fprintf('Compute %f (s)\n',pre_t);
euro_call=price(m+1-S0/deltaS,1);
fprintf('Euro call is %f dollars\n',euro_call);
% fprintf('Basic requirement: \nStandard Deviation: %f\n',std(basic1));
% fprintf('Under %d times of testing, 95%% of confidence interval for the call option value are: \n',repeat_times);
% fprintf('Lower bounds: %f\n',mean(basic1)-2*std(basic1));
% fprintf('Mean          %f\n',mean(basic1));
% fprintf('Upper bounds: %f\n',mean(basic1)+2*std(basic1));
% fprintf('\n');
% fprintf('Bonus1: \nStandard Deviation: %f\n',std(bonus1));
% fprintf('Under %d times of testing, 95%% of confidence interval for the call option value are: \n',repeat_times);
% fprintf('Lower bounds: %f\n',mean(bonus1)-2*std(bonus1));
% fprintf('Mean          %f\n',mean(bonus1));
% fprintf('Upper bounds: %f\n',mean(bonus1)+2*std(bonus1));
% fprintf('Bonus2: \nStandard Deviation: %f\n',std(bonus2));
% fprintf('\n');
% fprintf('Under %d times of testing, 95%% of confidence interval for the call option value are: \n',repeat_times);
% fprintf('Lower bounds: %f\n',mean(bonus2)-2*std(bonus2));
% fprintf('Mean          %f\n',mean(bonus2));
% fprintf('Upper bounds: %f\n',mean(bonus2)+2*std(bonus2));