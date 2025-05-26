function [out] = GCmy(x, y, p)
%GCmy: compute the Granger causality from x to y with order p VAR.
% x: dx*T time series with T time points.
% y: dy*T time series with T time points.
% p: the order of VAR model to estimate GC.
%Use the multivariate regression.

if (nargin<3 || p<1), p = 1; end
[dy, T] = size(y); dx = size(x, 1);
Y = y(:, (p+1):end)'; %(T-p)*dy matrix

X0 = zeros(T-p, dy*p); % dy*p columns
for i = 1:p
    X0(:, ((i-1)*dy+1):(i*dy)) = y(:, (p+1-i):(end-i))';
end
X1 = zeros(T-p, dx*p); % dx*p columns
for i = 1:p
    X1(:, ((i-1)*dx+1):(i*dx)) = x(:, (p+1-i):(end-i))';
end

%H1 model
[~, ~, e1] = mvregress([X0, X1], Y);
%H0 model
[~, ~, e0] = mvregress(X0, Y);
%GC index
out = -log(mean(var(e1, 0, 1))/mean(var(e0, 0, 1)));
end

