function [ y ] = Tosz( x,D )
%×Óº¯Êý Tosz
 y=zeros(D,1);
for i=1:D
    if x(i)~=0
        xj=log(abs(x(i)));
    else 
        xj=0;
    end
    if x(i)<0
        signx=-1;
    elseif x(i)==0
        signx=0;
    else
        signx=1;
    end
    if x(i)>0
        c1=10;
    else
        c1=5.5;
    end
    if x(i)>0
        c2=7.9;
    else
        c2=3.1;
    end
    y(i)=signx*exp(xj+0.049*(sin(c1*xj)+sin(c2*xj)));
end



end

