function oPrice = finDiffExplicit(X,S0,r,sig,Svec,tvec,oType)
% Function to calculate the price of a vanilla European
% Put or Call option using the explicit finite difference method
%
% oPrice = finDiffExplicit(X,r,sig,Svec,tvec,oType)
%
% Inputs: X - strike
%       : S0 - stock price
%       : r - risk free interest rate
%       : sig - volatility
%       : Svec - Vector of stock prices (i.e. grid points)
%       : tvec - Vector of times (i.e. grid points)
%       : oType - must be 'PUT' or 'CALL'.
%
% Output: oPrice - the option price
%
% Notes: This code focuses on details of the implementation of the
%        explicit finite difference scheme.
%        It does not contain any programatic essentials such as error
%        checking.
%        It does not allow for optional/default input arguments.
%        It is not optimized for memory efficiency, speed or
%        use of sparse matrces.

% Author: Phil Goddard (phil@goddardconsulting.ca)
% Date  : Q4, 2007

% Get the number of grid points
M = length(Svec)-1;
N = length(tvec)-1;
% Get the grid sizes (assuming equi-spaced points)
dt = tvec(2)-tvec(1);

% Calculate the coefficients
% To do this we need a vector of j points
j = 1:M-1;
sig2 = sig*sig;
j2 = j.*j;
aj = 0.5*dt*(sig2*j2-r*j);
bj = 1-dt*(sig2*j2+r);
cj = 0.5*dt*(sig2*j2+r*j);

% Pre-allocate the output
price(1:M+1,1:N+1) = nan;

% Specify the boundary conditions
switch oType
    case 'CALL'
        % Specify the expiry time boundary condition
        price(:,end) = max(Svec-X,0);
        % Put in the minimum and maximum price boundary conditions
        % assuming that the largest value in the Svec is
        % chosen so that the following is true for all time
        price(1,:) = 0;
        price(end,:) = (Svec(end)-X)*exp(-r*tvec(end:-1:1));
    case 'PUT'
        % Specify the expiry time boundary condition
        price(:,end) = max(X-Svec,0);
        % Put in the minimum and maximum price boundary conditions
        % assuming that the largest value in the Svec is
        % chosen so that the following is true for all time
        price(1,:) = (X-Svec(end))*exp(-r*tvec(end:-1:1));
        price(end,:) = 0;
end

% Form the tridiagonal matrix
A = diag(bj);  % Diagonal terms
A(2:M:end) = aj(2:end); % terms below the diagonal
A(M:M:end) = cj(1:end-1); % terms above the diagonal

% Calculate the price at all interior nodes
offsetConstants = [aj(1); cj(end)];
for i = N:-1:1
    price(2:end-1,i) = A*price(2:end-1,i+1);
    % Offset the first and last terms
    price([2 end-1],i) = price([2 end-1],i) + ...
        offsetConstants.*price([1 end],i+1);
end

% Calculate the option price
oPrice = interp1(Svec,price(:,1),S0);