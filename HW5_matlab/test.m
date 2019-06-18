S = 50; E = 50; sigma = 0.4; r = 0.1; T = 0.5;
Dt = 0.05; N = T/Dt; M = 1e3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%% Geom Asian exact mean  %%%%%%%%%%%%
sigsqT= sigma^2*T*(N+1)*(2*N+1)/(6*N*N);
muT = 0.5*sigsqT + (r - 0.5*sigma^2)*T*(N+1)/(2*N);

d1 = (log(S/E) + (muT + 0.5*sigsqT))/(sqrt(sigsqT));
d2 = d1 - sqrt(sigsqT);

N1 = 0.5*(1+erf(d1/sqrt(2)));
N2 = 0.5*(1+erf(d2/sqrt(2)));

geo =  exp(-r*T)*( S*exp(muT)*N1 - E*N2 );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Spath = S*cumprod(exp((r-0.5*sigma^2)*Dt+sigma*sqrt(Dt)*randn(M,N)),2);

% Standard Monte Carlo
arithave = mean(Spath,2);
Parith = exp(-r*T)*max(arithave-E,0);  % payoffs 
Pmean = mean(Parith)
Pstd = std(Parith);
confmc = [Pmean-1.96*Pstd/sqrt(M), Pmean+1.96*Pstd/sqrt(M)]
