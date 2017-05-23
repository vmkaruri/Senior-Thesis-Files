function [CD] = correlation_dimension(x,varargin)
    if nargin == 1
        StepR = 0.5;
    else
        StepR = varargin{1};
    end
    N = length(x);
    a = zeros(N-1,N);
    for i = 1:N
        j = 1:N;
        j(i) = [];
        a(:,i) = abs(x(i)-x(j));
    end
    m = log10(min(min(a)));
    if abs(log10(min(min(a)))) == inf
        m = -3;
    end
    R = 10.^(m:StepR:log10(max(max(a)))-1);
    
    C = zeros(length(R));
    for k = 1:length(R)
        C(k) = log10(sum(sum(heaviside(R(k)-a)))/(N*(N-1)));
    end
    IX = find(C == -inf);
    C(IX) = [];
    R(IX) = [];
    [CD,b] = LSE(log10(R),C);
    if CD > 1
        CD = 1;
    end
end