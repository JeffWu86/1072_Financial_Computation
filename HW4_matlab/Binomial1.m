function [euro_price,amer_price,c,d] = Binomial1(St,r,q,sigma,t,T,Smax_t,n)
    %BINOMIAL1 Summary of this function goes here
    %   Detailed explanation goes here

    % Binomial
    dT=(T-t)/n;
    u=exp(sigma*sqrt(T/n));
    d=1/u;
    p=(exp(r*T/n-q*T/n)-d)/(u-d);

    %Smaxlist is Smax struct
    Smaxlist(1:n+1)=struct('price',[],'bool',false,'eurovalue',0,'amervalue',0);
    
    level=0;
    for i=1:n
        if(St*power(u,n-i+1)>=Smax_t)
            Smaxlist(i).price=St*power(u,n-i+1);
        else
            Smaxlist(i).price=Smax_t;
            level=i;
            break;
        end
    end
    if(level==0)
        Smaxlist(n+1).price=max(Smax_t,St);
    end
    
    pricelist(1:2*n+1)=nan;
    pricelist(n+1)=St;
    for i=1:n
       pricelist(i)=St*power(u,n-i+1);
       pricelist(2*n-i+2)=St*power(d,n-i+1);
    end
    % Initiate Smax table
    s=struct('Price',[],'Smax',[]);
    tree(1:2*n+1,1:n+1)=struct(s);
    tree(n+1,1).Price=St;
    tree(n+1,1).Smax=Smaxlist;
    if(level==0)
        tree(n+1,1).Smax(n+1).bool=true;
    else
        tree(n+1,1).Smax(level).bool=true;
    end
    
    for i=2:n+1
        for j=n+2-i:2:n+i
            tree(j,i).Price=pricelist(j);
            tree(j,i).Smax=Smaxlist;
            % Not top of tree, inherit from up
            if(j~=n+2-i)
                for k=1:n+1
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
                for k=1:n+1
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
        for k=1:n+1
            if(tree(j,n+1).Smax(k).bool==true)
                tree(j,n+1).Smax(k).eurovalue=max(tree(j,n+1).Smax(k).price-tree(j,n+1).Price,0.0);
                tree(j,n+1).Smax(k).amervalue=tree(j,n+1).Smax(k).eurovalue;
            end
        end
    end

    %Backward induction
    for i=n:-1:1
        for j=n+2-i:2:n+i
            for k=1:n+1
                if(tree(j,i).Smax(k).bool==true)
                    if(tree(j,i).Smax(k).price>=tree(j,i).Price*u)
                        tree(j,i).Smax(k).eurovalue=(p*tree(j-1,i+1).Smax(k).eurovalue+...
                            (1-p)*tree(j+1,i+1).Smax(k).eurovalue)*exp(-r*dT);
                        tree(j,i).Smax(k).amervalue=(p*tree(j-1,i+1).Smax(k).amervalue+...
                            (1-p)*tree(j+1,i+1).Smax(k).amervalue)*exp(-r*dT);
                    else
                        tree(j,i).Smax(k).eurovalue=(p*tree(j-1,i+1).Smax(k-1).eurovalue+...
                            (1-p)*tree(j+1,i+1).Smax(k).eurovalue)*exp(-r*dT);
                        tree(j,i).Smax(k).amervalue=(p*tree(j-1,i+1).Smax(k-1).amervalue+...
                            (1-p)*tree(j+1,i+1).Smax(k).amervalue)*exp(-r*dT);
                    end
                    tree(j,i).Smax(k).amervalue=max(tree(j,i).Smax(k).amervalue,...
                        tree(j,i).Smax(k).price-tree(j,i).Price);
                end
            end
        end
    end

    if(level==0)
        euro_price=tree(n+1,1).Smax(n+1).eurovalue;
        amer_price=tree(n+1,1).Smax(n+1).amervalue;
    else
        euro_price=tree(n+1,1).Smax(level).eurovalue;
%         euro_price=0;
%         amer_price=0;
        amer_price=tree(n+1,1).Smax(level).amervalue;
    end
    
    c=tree;
    d=Smaxlist;
end

