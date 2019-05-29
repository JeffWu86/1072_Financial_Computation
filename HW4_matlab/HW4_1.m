
[St,r,q,sigma,t,T,Smax_t,n,nrolls,num_of_rep]=readdata('input_2.txt');

zmat1=randn(nrolls,n);
zmat2(1:nrolls,1:n)=nan;

u=exp(sigma*sqrt(T/n));
d=1/u;
p=(exp(r*T/n-q*T/n)-d)/(u-d);

dT=(T-t)/n;
for i=1:n
    mu=log(St)+(r-q-sigma*sigma/2)*dT*i;
    sigma1=sigma*dT*i;
    for j=i:nrolls
        zmat1(j,i)=lognrnd(mu,sigma1);
        zmat2(j,i)=exp(zmat1(j,i));
    end
end

% Binomail
pricelist(1:2*n+1)=St;
for i=1:n
   pricelist(i)=St*power(d,n-i+1);
   pricelist(2*n-i+2)=St*power(u,n-i+1);
end

s=struct('name',[],'dollar',10,'i',[],'j',[]);
s1(1:10,1:10)=struct(s);