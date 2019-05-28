
[St,r,q,sigma,t,T,Smax_t,n,nrolls,num_of_rep]=readdata('input_1.txt');

zmat1=randn(nrolls,n);
zmat2(1:nrolls,1:n)=nan;

dT=(T-t)/n;
for i=1:n
    mu=log(St)+(r-q-sigma*sigma/2)*dT*i;
    sigma1=sigma*dT*i;
    for j=i:nrolls
        zmat1(j,i)=lognrnd(mu,sigma1);
        zmat2(j,i)=exp(zmat1(j,i));
    end
end