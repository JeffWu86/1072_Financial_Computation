
[S0,K,r,q,sigma,T,Smin,Smax,m,n]=readdata('input_6.txt');
%Implicit methond
Im_EC= Implicit(S0,K,r,q,sigma,T,Smin,Smax,m,n,'CALL');
Im_EP= Implicit(S0,K,r,q,sigma,T,Smin,Smax,m,n,'PUT');
Im_AC= Implicit(S0,K,r,q,sigma,T,Smin,Smax,m,n,'CALL',true);
Im_AP= Implicit(S0,K,r,q,sigma,T,Smin,Smax,m,n,'PUT',true);

fprintf('Implicit method\n');
fprintf('Euro call is %f dollars\n',Im_EC);
fprintf('Euro put  is %f dollars\n',Im_EP);
fprintf('Amer call is %f dollars\n',Im_AC);
fprintf('Amer put  is %f dollars\n',Im_AP);
fprintf('\n');

%Explicit method
Ex_EC= Explicit(S0,K,r,q,sigma,T,Smin,Smax,m,n,'CALL');
Ex_EP= Explicit(S0,K,r,q,sigma,T,Smin,Smax,m,n,'PUT');
Ex_AC= Explicit(S0,K,r,q,sigma,T,Smin,Smax,m,n,'CALL',true);
Ex_AP= Explicit(S0,K,r,q,sigma,T,Smin,Smax,m,n,'PUT',true);

fprintf('Explicit method\n');
fprintf('Euro call is %f dollars\n',Ex_EC);
fprintf('Euro put  is %f dollars\n',Ex_EP);
fprintf('Amer call is %f dollars\n',Ex_AC);
fprintf('Amer put  is %f dollars\n',Ex_AP);

% a=zeros(1,4);
% a(1)=finDiffImplicit(60,50,0.05,0.2,0:1:100,0:0.001:1,'CALL');
% a(2)=finDiffImplicit(60,50,0.05,0.2,0:1:100,0:0.001:1,'PUT');
% a(3)=finDiffExplicit(60,50,0.05,0.2,0:1:100,0:0.001:1,'CALL');
% a(4)=finDiffExplicit(60,50,0.05,0.2,0:1:100,0:0.001:1,'PUT');
% fprintf('%f %f %f %f\n',a);


