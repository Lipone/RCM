function out = EE(x, y, p, k, Thei)
%EE: compute the embedding entropy causality from x to y with order p.
% x: dx*T time series with T time points.
% y: dy*T time series with T time points.
% p: the order of model to estimate causality.
% k: the k-th nearest number to use in calculating entropy, at least 2.
% Thei: half the length of Theiler correction window. [-Thei, Thei] around point i. Thei>=p.
%Use the kNN method to estimate the MI(X, YpNN), where X = x_*, YpNN = (dy*(p+1)+1) NN around [y_t, y_*], *means t-1, t-2, ... t-p.
%kNN formular MI(X, YpNN) = psi(k) - <psi(nX+1) + psi(nYpNN+1) - psi(nN+1)>.

if (nargin<3 || p<1), p = 1; end
if (nargin<4 || k<2), k = 2; end
if (nargin<5 || Thei<p), Thei = p; end

[dy, T] = size(y); dx = size(x, 1);

X = zeros(T-p, dx*p); % dx*p columns
for i = 1:p
    X(:, ((i-1)*dx+1):(i*dx)) = x(:, (p+1-i):(end-i))';
end

Y = zeros(T-p, dy*(p+1)); %dy*(p+1) columns
Y(:, 1:dy) = y(:, (p+1):end)';
for i = 1:p
    Y(:, (i*dy+1):((i+1)*dy)) =  y(:, (p+1-i):(end-i))'; %dy*p columns
end

% the dy*(p+1)+1 NN of Y_t on Manifold M_Y.
YpNN = zeros(T-p, size(Y, 2)*(dy*(p+1)+1));
for i = 1:(T-p)
    % Theiler correction. Points around i are excluded.
    ids = max(1, i-Thei); idf = min(i+Thei, T-p); idx = ones(T-p, 1); idx(ids:idf) = 0; idx = logical(idx);
    %dy*(p+1)+1 NN in Y_t space
    tempY1 = Y(idx, :);
    [pnnidx, ~] = knnsearch(tempY1, Y(i, :), 'K', dy*(p+1)+1, 'Distance', 'euclidean'); %dy*(p+1)+1 th NN
    tempY2 = tempY1(pnnidx, :)';
    YpNN(i, :) = tempY2(:)';
end

%Calculate mutual information of prediction
% kNN in (X, YpNN) space
nX = zeros(T-p, 1); nYpNN = zeros(T-p, 1); nN = zeros(T-p, 1);
for i = 1:(T-p)
     %Theiler correction. Points around i are excluded.
    ids = max(1, i-Thei); idf = min(i+Thei, T-p); idx = ones(T-p, 1); idx(ids:idf) = 0; idx = logical(idx);
    %kNN in (X, YpNN) space
    [~, D] = knnsearch([X(idx, :), YpNN(idx, :)], [X(i, :), YpNN(i, :)], 'K', k, 'Distance', 'chebychev'); %kth NN
    halfepsilonXYkNN = D(k); %the distance of kNN = epsilon/2
    nX(i) = sum(pdist2(X(idx, :), X(i, :), 'chebychev')<halfepsilonXYkNN);
    nYpNN(i) = sum(pdist2(YpNN(idx, :), YpNN(i, :), 'chebychev')<halfepsilonXYkNN);
    nN(i) = sum(idx);
end

%tauDC index
idN = (nX~=0)&(nYpNN~=0);
out = abs(psi(k) - mean(psi(nX(idN)+1)) - mean(psi(nYpNN(idN)+1)) + mean(psi(nN(idN)+1)));
end
