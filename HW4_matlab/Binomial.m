function [euro_price] = Binomial(St,r,q,sigma,t,T,n,optType)
    %BINOMIAL Summary of this function goes here
    %   Detailed explanation goes here
% Binomial
dT=(T-t)/n;
u=exp(sigma*sqrt(T/n));
d=1/u;
p=(exp(r*T/n-q*T/n)-d)/(u-d);

%pricelist is Smax struct
pricelist(1:2*n+1)=struct('price',[],'bool',false,'value',0);
pricelist(n+1).price=St;
for i=1:n
   pricelist(i).price=St*power(u,n-i+1);
   pricelist(2*n-i+2).price=St*power(d,n-i+1);
end

% Initiate Smax table
s=struct('Price',[],'Smax',[]);
tree(1:2*n+1,1:n+1)=struct(s);
tree(n+1,1).Price=St;
tree(n+1,1).Smax=pricelist;
tree(n+1,1).Smax(n+1).bool=true;
for i=2:n+1
    for j=n+2-i:2:n+i
        tree(j,i).Price=pricelist(j).price;
        tree(j,i).Smax=pricelist;
        % Not top of tree, inherit from up
        if(j~=n+2-i)
            for k=1:2*n+1
                if(tree(j-1,i-1).Smax(k).bool==true)
                    if(tree(j-1,i-1).Smax(k).price >= tree(j,i).Price )
                        tree(j,i).Smax(k).bool=true;
                    else
                        tree(j,i).Smax(j).bool=true;
                    end
                end
            end
        end
        % Not bottom of tree, inherit from down
        if(j~=n+i)
            for k=1:2*n+1
                if(tree(j+1,i-1).Smax(k).bool==true)
                    if(tree(j+1,i-1).Smax(k).price >= tree(j,i).Price )
                        tree(j,i).Smax(k).bool=true;
                    else
                        tree(j,i).Smax(j).bool=true;
                    end
                end
            end
        end
    end
end

% Last price
for j=1:2:2*n+1
    for k=1:2*n+1
        if(tree(j,n+1).Smax(k).bool==true)
            tree(j,n+1).Smax(k).value=max(tree(j,n+1).Smax(k).price-tree(j,n+1).Price,0.0);
        end
    end
end

%Backward induction
for i=n:-1:1
    for j=n+2-i:2:n+i
        for k=1:2*n+1
            if(tree(j,i).Smax(k).bool==true)
                if(tree(j,i).Smax(k).price>tree(j,i).Price)
                    tree(j,i).Smax(k).value=(p*tree(j-1,i+1).Smax(k).value+...
                        (1-p)*tree(j+1,i+1).Smax(k).value)*exp(-r*dT);
                else
                    tree(j,i).Smax(k).value=(p*tree(j-1,i+1).Smax(k-1).value+...
                        (1-p)*tree(j+1,i+1).Smax(k).value)*exp(-r*dT);
                end
                if(optType=='AMER')
                    tree(j,i).Smax(k).value=max(tree(j,i).Smax(k).value,tree(j,i).Smax(k).price-tree(j,i).Price);
                end
            end
        end
    end
end

euro_price=tree(n+1,1).Smax(n+1).value;
end

