function [ y ] = Tasy( x,beita,D )
%�Ӻ��� Tasy
y=zeros(D,1);
for i=1:D
    if x(i)>0
        y(i)=x(i)^(1+beita*(i-1)*sqrt(x(i))/(D-1));
    end
end


end

