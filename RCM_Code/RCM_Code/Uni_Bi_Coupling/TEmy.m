function out = TEmy(x, y, p, k, Thei)
%TEmy: compute the Transfer Entropy causality from x to y with order p.
% x: dx*T time series with T time points.
% y: dy*T time series with T time points.
% p: the order of model to estimate TE.
% k: the k-th nearest number to use, at least 2.
% Thei: half the length of Theiler correction window. [-Thei, Thei] around point i. Thei>=p.
%Use the kNN method to estimate the CMI(X, Y|Z), where X = x_*, Y = y_t, Z =
%y_*, *means t-1, t-2, ... t-p.
%kNN formular CMI(X, Y|Z) = psi(k) + <psi(nZ+1) - psi(nYZ+1) - psi(nXZ+1)>.

if (nargin<3 || p<1), p = 1; end
if (nargin<4 || k<2), k = 2; end
if (nargin<5 || Thei<p), Thei = p; end

[dy, T] = size(y); dx = size(x, 1);
Y = y(:, (p+1):end)'; %(T-p)*dy matrix

X = zeros(T-p, dx*p); % dx*p columns
for i = 1:p
    X(:, ((i-1)*dx+1):(i*dx)) = x(:, (p+1-i):(end-i))';
end
Z = zeros(T-p, dy*p); % dy*p columns
for i = 1:p
    Z(:, ((i-1)*dy+1):(i*dy)) = y(:, (p+1-i):(end-i))';
end

nZ = zeros(T-p, 1); nYZ = zeros(T-p, 1); nXZ = zeros(T-p, 1);
for i = 1:(T-p)
     %Theiler correction. Points around i are excluded.
    ids = max(1, i-Thei); idf = min(i+Thei, T-p); idx = ones(T-p, 1); idx(ids:idf) = 0; idx = logical(idx);
    %kNN in (X, Y, Z) space
    [~, D] = knnsearch([X(idx, :), Y(idx, :), Z(idx, :)], [X(i, :), Y(i, :), Z(i, :)], 'K', k, 'Distance', 'chebychev'); %kth NN
    halfepsilonXYZ = D(k); %the distance of kNN = epsilon/2
    nZ(i) = sum(pdist2(Z(idx, :), Z(i, :), 'chebychev')<halfepsilonXYZ);
    nYZ(i) = sum(pdist2([Y(idx, :), Z(idx, :)], [Y(i, :), Z(i, :)], 'chebychev')<halfepsilonXYZ);
    nXZ(i) = sum(pdist2([X(idx, :), Z(idx, :)], [X(i, :), Z(i, :)], 'chebychev')<halfepsilonXYZ);
end
%TE index
idN = (nZ~=0)&(nYZ~=0)&(nXZ~=0);
out = abs(psi(k) + mean(psi(nZ(idN)+1)) - mean(psi(nXZ(idN)+1)) - mean(psi(nYZ(idN)+1)));
end
