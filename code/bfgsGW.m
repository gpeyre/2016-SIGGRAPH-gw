function [gamma,optimTimes,optimVals] = bfgsGW(D0,D,alpha)
alpha
n = size(D0,1);

assert(n == size(D,1)); % only going to make this work for simplest case
assert(norm(D-D','fro')<1e-8);
assert(norm(D0-D0','fro')<1e-8);

idx = reshape((1:(n*n))',n,n);

% Slow way to construct constraint matrices, not counting toward timing in the paper
constraintMtx = sparse(0,n*n);
constraintRhs = [];
c = 1;

for i=1:n
    constraintMatrix(c,idx(i,:)) = 1; % row sum
    constraintRhs(c,1) = n;
    c = c + 1;
    
    constraintMatrix(c,idx(:,i)) = 1; % column sum
    constraintRhs(c,1) = n;
    c = c + 1;
end

constraintMatrix = sparse(constraintMatrix);

f = @(x) gwObjective(x,D0,D,n,alpha);
x0 = ones(n*n,1);

% algorithm = interior-point, active-set, sqp (latter two are not large scale)
options = optimoptions('fmincon',...
                       'GradObj','on',...
                       'display','none',...
                       'algorithm','interior-point',...
                       'hessian','lbfgs',...
                       'derivativecheck','off');

global optimTimes optimVals

optimVals = [];
optimTimes = []; % in seconds
f(x0); % starts timer...

warning('off','MATLAB:nearlySingularMatrix');
result = fmincon(f, x0,...
                 [], [],... % inequality constraints
                 constraintMatrix,constraintRhs,... % equality constraints
                 zeros(n*n,1),[],... % bounds
                 [], ... nonlinear constraints
                 options);
warning('on','MATLAB:nearlySingularMatrix');
             
gamma = reshape(result,n,n);
optimTimes = optimTimes - optimTimes(1); % shift to start at t=0

function [f,g] = gwObjective(gamma,D0,D,n,alpha)

gamma = reshape(gamma,n,n);

o = ones(n,1);
prod = D0*gamma*D-.5*(D0.^2)*gamma*o*o'-.5*o*o'*gamma*(D.^2);
logGamma = log(gamma);

f = -dot(gamma(:),prod(:))/n^4 + (alpha/(n*n))*dot(gamma(:),logGamma(:));

global optimTimes optimVals % every time we compute f, keep track
optimTimes = [optimTimes; cputime];
optimVals = [optimVals; f];

if nargout > 1 % gradient required
    grad = -2*prod/n^4 + (alpha/(n*n))*(logGamma+1);
    g = grad(:);
end
