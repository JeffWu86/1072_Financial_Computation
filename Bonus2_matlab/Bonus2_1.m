
[S0,K,r,q,sigma,T,Smin,Smax,m,n]=readdata('input_3.txt');
%Implicit methond
Im_EC= Implicit(S0,K,r,q,sigma,T,Smin,Smax,m,n,'EURO_CALL');
Im_EP= Implicit(S0,K,r,q,sigma,T,Smin,Smax,m,n,'EURO_PUT');


fprintf('Euro call is %f dollars\n',Im_EC);
fprintf('Euro put  is %f dollars\n',Im_EP);

%Explicit method
Ex_EC= Explicit(S0,K,r,q,sigma,T,Smin,Smax,m,n,'EURO_CALL');
Ex_EP= Explicit(S0,K,r,q,sigma,T,Smin,Smax,m,n,'EURO_PUT');
fprintf('Euro call is %f dollars\n',Ex_EC);
fprintf('Euro put  is %f dollars\n',Ex_EP);

a=zeros(1,4);
a(1)=finDiffImplicit(60,50,0.05,0.2,0:1:100,0:0.001:1,'CALL');
a(2)=finDiffImplicit(60,50,0.05,0.2,0:1:100,0:0.001:1,'PUT');
a(3)=finDiffExplicit(60,50,0.05,0.2,0:1:100,0:0.001:1,'CALL');
a(4)=finDiffExplicit(60,50,0.05,0.2,0:1:100,0:0.001:1,'PUT');
fprintf('%f %f %f %f\n',a);


