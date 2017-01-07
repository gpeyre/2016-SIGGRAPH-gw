function [gamma,optimTimes,optimVals] = minimalGW(D0,D,alpha)
% Minimalistic implementation of GW optimization for use in timing tests
% (only fair, since BFGS implementation works on simple problems that don't
% need some of the extra matrices).

n0 = size(D0,1); % Size of metric matrices
n = size(D,1);

assert(n0==n);
assert(norm(D-D','fro')<1e-8);
assert(norm(D0-D0','fro')<1e-8);

gamma = ones(n,n);
sinkhornTol = 1e-6;
maxSinkhornIter = 500; % Max inner loop

optimTimes = cputime;
optimVals = gwObjective(gamma,D0,D,n,alpha);

for iteration=1:100 % Gromov-Wasserstein iteration
    target = D0*gamma*D/(n^2); % Will go in exponent
    K = exp(target/alpha); % Evaluate f_alpha 
    
    if iteration == 1 % Warm start in future iterations
        v = ones(n0,1);
        w = ones(n,1);
    end
    
    Kw = K*(w/n);
    for i=1:maxSinkhornIter % check this...
        v = 1./Kw;
        Kv = K'*(v/n);
        w = 1./Kv;
        Kw = K*(w/n);

        % Stop when margins are close to 1
        if norm(v.*Kw-1,1)<n*sinkhornTol && norm(w.*Kv-1,1)<n*sinkhornTol
            break;
        end
    end
    gamma = diag(v)*K*diag(w);
    
    optimTimes = [optimTimes;cputime];
    optimVals = [optimVals;gwObjective(gamma,D0,D,n,alpha)];
end

optimTimes = optimTimes - optimTimes(1); % start at t=0

function f = gwObjective(gamma,D0,D,n,alpha)

o = ones(n,1);
prod = D0*gamma*D-.5*(D0.^2)*gamma*o*o'-.5*o*o'*gamma*(D.^2);
logGamma = log(gamma);
f = -dot(gamma(:),prod(:))/n^4 + (alpha/(n*n))*dot(gamma(:),logGamma(:));