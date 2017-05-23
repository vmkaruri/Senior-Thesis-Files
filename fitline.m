function [p,perr] = fitline(x,y)

%  FITLINE(X,Y)  Fitting straight line trough datapoints specified by
%     vectors X and Y.
%     Simplified and slightly changed version of POLYFIT
%     Returns parameters of the fit and their estimated errors
%     instead of Cholesky factor matrix

%  Kirill Pankratov,  Jan. 5, 1994.

x = x(:); y = y(:);  % Make input column vectors

% Construct Vandermonde matrix.
n = 1;
V(:,n+1) = ones(length(x),1);
for j = n:-1:1
    V(:,j) = x.*V(:,j+1);
end

% Solve least squares problem.
[Q,R] = qr(V,0);
p = R\(Q'*y);   % Same as p = V\y;
r = y - V*p;
p = p';         % Polynomial coefficients are row vectors by convention.

df = length(y) - (n+1);       % Degrees of freedom
perr = norm(r)*diag(inv(R))'; % Errors of the fit
