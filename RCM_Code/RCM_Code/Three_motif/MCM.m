function out = MCM(x, y, p, Thei)
%ECmy: compute the embedding causality from x to y with order p.
% x: dx*T time series with T time points.
% y: dy*T time series with T time points.
% p: the order of model. The embedding dimension is dy*(p+1) for Y and dx*p for X^*, which should be >2*innerdimension.
% Thei: half the length of Theiler correction window. [-Thei, Thei] around point i. Thei>=p.
%EC uses Y to predict X^*, *means t-1, t-2, ... t-p
% Covergent cross mapping (CCM) uses rho = |corr(\hat{X^*}, X^*)| as the index.

if (nargin<3 || p<2), p = 2; end
if (nargin<4 || Thei<p), Thei = p; end

[dy, T] = size(y); dx = size(x, 1);

X = zeros(T-p, dx*(p+1)); 
X(:, 1:dx) = x(:, (p+1):end)';
for i = 1:p
    X(:, (i*dx+1):((i+1)*dx)) = x(:, (p+1-i):(end-i))'; %(T-p)*dx matrix
end

Y = zeros(T-p, dy*(p+1)); %dy*(p+1) columns
Y(:, 1:dy) = y(:, (p+1):end)';
for i = 1:p
    Y(:, (i*dy+1):((i+1)*dy)) =  y(:, (p+1-i):(end-i))'; %(T-p)*dy matrix
end

% kNN in Y space. Use dy*(p+1)+1 points in dy*(p+1) space.
hatX = zeros(size(X));
for i = 1:(T-p)
    % Theiler correction. Points around i are excluded.
    ids = max(1, i-Thei); idf = min(i+Thei, T-p); idx = ones(T-p, 1); idx(ids:idf) = 0; idx = logical(idx);
    tempY = Y(idx, :); tempX = X(idx, :);
    [idy, D] = knnsearch(tempY, Y(i, :), 'K', dy*(p+1)+1);
    U = exp(-D./D(1)); U = U./sum(U); %calculate the weight of dy*(p+1)+1 nearest neighbors.
    hatX(i, :) = sum(tempX(idy, :).*U', 1);
end

%MCM index
MCM = zeros(1, size(X, 2));
for i = 1:size(X, 2)
    MCM(i) = abs(corr(X(:, i), hatX(:, i))); %correlation of prediction hatX and real X, for each component.
end
out = mean(MCM);
end
