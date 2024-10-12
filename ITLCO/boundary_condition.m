function [output] = boundary_condition(X,Lb,Ub,dim)
for i = 1:dim
    if X(i) > Ub(i)
        X(i) = Lb(i) + rand.*(Ub(i)-Lb(i));
    end
    if X(i) < Lb(i)
        X(i) = Lb(i) + rand*(Ub(i)-Lb(i));
    end
end
output = X;
end