function [ y ] = Dia( arerfa,D )
%�Ӻ��� �ԽǾ���
y=zeros(D,D);
for i=1:D
    y(i,i)=arerfa^((i-1)/(2*(D-1)));
end


end

