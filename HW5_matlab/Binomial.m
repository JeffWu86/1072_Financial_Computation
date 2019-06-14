function [bino_euro,bino_amer] = Binomial(St,K,r,q,sigma,T_t,M,n,passing_time);
%BINOMIAL Summary of this function goes here
%   Detailed explanation goes here
    dT=(T_t)/n;
    npass=n*passing_time/T_t;
    u=exp(sigma*sqrt(dT));
    d=1/u;
    p=(exp(r*dT-q*dT)-d)/(u-d);

    % Initiate Smax table
    s=struct('Amax',[],'Amin',[],'price',[]);
    tree(1:n+1,1:n+1)=struct(s);

    for i=1:n+1
        for j=1:i
            tree(j,i).Amax=St*( npass + (1-power(u,i-j+1))/(1-u) + ...
                           power(u,i-j)*d*(1-power(d,j-1))/(1-d) ) / (i+npass);
            tree(j,i).Amin=St*( npass + (1-power(d,j))/(1-d) + ...
                           power(d,j-1)*u*(1-power(u,i-j))/(1-u) ) / (i+npass);
        end
    end

    
    for j=1:n+1
        tree(j,n+1).price=nan(M+1,2);
        for k=1:M+1
            tree(j,n+1).price(k,1) = ( (M-k+1)*tree(j,n+1).Amax+(k-1)*tree(j,n+1).Amin ) / M;
            tree(j,n+1).price(k,2)=max( tree(j,n+1).price(k,1) - K , 0 );
        end
    end
    
    for i=n:-1:1
        for j=1:i
            tree(j,i).price=nan(M+1,2);
            for k=1:M+1
                if(tree(j,i).Amax~=tree(j,i).Amin)
                    tree(j,i).price(k,1) = ( (M-k+1)*tree(j,i).Amax+(k-1)*tree(j,i).Amin ) / M;
                else
                    tree(j,i).price(k,1) = tree(j,i).Amax;
                end
                
                Au = ((i+npass)*tree(j,i).price(k,1)+St*power(u,i-j+1)*power(d,j-1)) / (i+1+npass);
                if( abs(Au-tree(j,i+1).price(1,1))<1e-5 )
                    Cu=tree(j,i+1).price(1,2);
                elseif( abs(Au-tree(j,i+1).price(M+1,1))<1e-5 )
                    Cu=tree(j,i+1).price(M+1,2);
                else
                    Cu=interp1(tree(j,i+1).price(:,1), tree(j,i+1).price(:,2), Au, 'linear');
                end
                
                Ad = ((i+npass)*tree(j,i).price(k,1)+St*power(u,i-j)*power(d,j)) / (i+1+npass);
                if( abs(Ad-tree(j+1,i+1).price(1,1))<1e-5 )
                    Cd=tree(j+1,i+1).price(1,2);
                elseif( abs(Au-tree(j+1,i+1).price(M+1,1))<1e-5 )
                    Cd=tree(j+1,i+1).price(M+1,2);
                else
                    Cd=interp1(tree(j+1,i+1).price(:,1), tree(j+1,i+1).price(:,2), Ad, 'linear');
                end
                tree(j,i).price(k,2) = ( p*Cu + (1-p)*Cd ) * exp(-r*dT);
            end
        end
    end
    bino_euro=tree(1,1).price(1,2);
    bino_amer=0;
end
