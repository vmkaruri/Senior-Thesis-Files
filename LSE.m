function [a,b] = LSE(x,y)
    N=size(x,2);
    a=[N*x*y'-sum(x)*sum(y)]/[N*sum(x.^2)-sum(x)^2];
    b=[sum(x.^2)*sum(y)-x*y'*sum(x)]/[N*sum(x.^2)-sum(x)^2];
end