function [A] = choles(covmat,n)
    %CHOLES Summary of this function goes here
    %   Detailed explanation goes here
    
    %Step 1
    A=zeros(n,n);
    A(1,1)=sqrt(covmat(1,1));
    
    for j =2:n
        A(1,j)=covmat(1,j)/A(1,1);
    end
    %Step 2 and 3
    for i=2:n-1
        temp1=0;
        for k=1:i-1
            temp1=temp1+A(k,i)*A(k,i);
        end
        A(i,i)=sqrt(covmat(i,i)-temp1);
        for j=i+1:n
            temp1=0;
            for k=1:i-1
                temp1=temp1+A(k,i)*A(k,j);
            end
            if(covmat(i,j)~=temp1)
                A(i,j)=(covmat(i,j)-temp1)/A(i,i);
            else
                A(i,j)=0;
            end
        end
    end
    %Step 4
    temp1=0;
    for k=1:n
        temp1=temp1+A(k,n)*A(k,n);
    end
    A(n,n)=sqrt(covmat(n,n)-temp1);
end

