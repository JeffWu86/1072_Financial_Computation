
[S0,K,r,q,sigma,T,Smin,Smax,m,n]=readdata('input_3.txt');
tic;
%Implicit methond
[Im_AC, Im_AP, Im_EC, Im_EP] = Implicit(S0,K,r,q,sigma,T,Smin,Smax,m,n);
optType='EURO_CALL';
dS=(Smax-Smin)/m;
dT=T/n;
%Boundary condition for price
price(1:m+1,1:n+1)=nan;
switch optType
    case 'EURO_CALL'
        price(1,:)=m*dS-K;
        price(end,:)=max(Smin-K,0);
        %Stock price in time T
        for i=1:m+1
            price(i,n+1)=max((m-i+1)*dS-K,0);
        end
    case 'EURO_PUT'
        price(1,:)=max(K-m*dS,0);
        price(end,:)=K-Smin;
        %Stock price in time T
        for i=1:m+1
            price(i,n+1)=max(K-(m-i+1)*dS,0);
        end
end

% coef of aj, bj, cj
i=1:m-1;
aj=(r-q)/2*i*dT-0.5*power(sigma*i,2)*dT;
bj=1+power(sigma*i,2)*dT+r*dT;
cj=-(r-q)/2*i*dT-0.5*power(sigma*i,2)*dT;

% A matrix
impmat=zeros(m-1,m-1);
for i=1:m-1
    impmat(i,i)=bj(1,m-i);
    if(i~=1)
        impmat(i,i-1)=cj(1,m-i);
    end
    if(i~=m-1)
        impmat(i,i+1)=aj(1,m-i);
    end
end

[L,U] = lu(impmat);
for idx=n:-1:1
    b=price(2:m,idx+1);
    b(end)=b(end)-aj(1)*price(m+1,idx); %Alwayw zero
    b(1)=b(1)-cj(end)*price(1,idx);
    %x=impmat\b;
    x=U\(L\b);
    for j=1:m-1
        if(x(j,1)<0) 
            x(j,1)=0;
        end
    end
    price(2:m,idx)=x;
end
pre_t=toc;
fprintf('Compute %f (s)\n',pre_t);

euro_call=interp1(Smax:-dS:Smin,price(:,1),S0);
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