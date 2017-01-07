function gamma = weightedGW(D0,W0,D,W,alpha,maxIter)
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
eta = .025;

for iteration=1:maxIter % Gromov-Wasserstein iteration
%     target = D0*gamma*D/(n^2); % Will go in exponent

    t1 = -.5*((D0.^2).*W0)*gamma*W*0;
    t2 = (D0.*W0)*gamma*(D.*W);
    t3 = -.5*W0*gamma*((D.^2).*W)*0;
    target = (t1+t2+t3) / n^2;
    
    K = exp(target/alpha); % Evaluate f_alpha 
    K = (K.^eta).*(gamma.^(1-eta));
    
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
    
    imagesc(gamma); axis equal; axis off;
    title(sprintf('Iteration %d...',iteration));
    colorbar;
    drawnow;
end
