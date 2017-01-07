function d = distmat(x,y)

% distmat - compute pairwise L^2 distance matrix
%
%   d = distmat(x,y);
%
%   d_ij = |x(:,i)-y(:,j)|
%
%   Copyright (c) 2015 Gabriel Peyre

if nargin<2
    y = x;
end

N = size(x,2);
a = sum(x.^2);
b = sum(y.^2);
d = repmat(b,N,1) + repmat(a',1,N) - 2*x'*y;
d = sqrt(max(d,0));
        
end