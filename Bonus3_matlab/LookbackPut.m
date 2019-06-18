function [output,pay] = LookbackPut(St,r,q,sigma,t,T,n,Smax_t,nrolls)
%LOOKBACKPUT Summary of this function goes here
%   Detailed explanation goes here
    price(1:nrolls,1:n+1)=nan;
    price(1:nrolls,1)=St;

    dT=(T-t)/n;
    for i=1:n
        for j=1:nrolls
            mu=log(price(j,i))+(r-q-sigma*sigma/2)*dT;
            sigma1=sigma*sqrt(dT);
            price(j,i+1)=exp(sigma1*randn()+mu);
        end
    end
    
    payoff=nan(nrolls,n);
    for j=1:nrolls
        payoff(j,n)=max(max(max(price(j,:)),Smax_t) - price(j,n+1), 0);
    end
    for i=n:-1:2
        exervalue=nan(nrolls,1);
        for j=1:nrolls
            exervalue(j,1)=max(max(max(price(j,1:i)),Smax_t) - price(j,i), 0);
        end

        % Compute how many EV not equal to zero
        count=0;
        for j=1:nrolls
            if(exervalue(j)~=0)
                count=count+1;
            end
        end
        % For those EV not zero, computer expected holding value
        if(count~=0)
            holdingvalue=nan(count,3);
            counter=1;
            for j=1:nrolls
                if(exervalue(j)~=0)
                    holdingvalue(counter,1)=payoff(j,i)*exp(-r*dT);
                    holdingvalue(counter,2)=price(j,i);
                    holdingvalue(counter,3)=power(price(j,i),2);
                    counter=counter+1;
                end
            end
            b=regress(holdingvalue(:,1),[ones(size(holdingvalue(:,1))),holdingvalue(:,2),holdingvalue(:,3)]);
            Expectholding=b(1)*1+b(2)*holdingvalue(:,2)+b(3)*holdingvalue(:,3);
        end

        % Renew payoff
        counter=1;
        for j=1:nrolls
            if(exervalue(j)==0)
                payoff(j,i-1)=payoff(j,i)*exp(-r*dT);
            elseif (Expectholding(counter)<exervalue(j))
                payoff(j,i-1)=exervalue(j);
                payoff(j,i:end)=0;
                counter=counter+1;
            else
                payoff(j,i-1)=holdingvalue(counter,1);
                counter=counter+1;
            end
        end
    end
    pay=payoff;
    output=max(max(Smax_t-St,sum(payoff(:,1))*exp(-r*dT)/nrolls),0);     
end

