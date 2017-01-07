function gamma = regularizedPartialGW(D0,D,g,gradG,options)
% Assuming uniform area weights for simplicity --- assumes the target is
% partial, so gamma*1 = 1.

n0 = size(D0,1); % Size of metric matrices
n = size(D,1);

regularizer = getoptions(options, 'regularizer', .005); % alpha
maxIter = getoptions(options, 'maxIter', 40); % Maximum number of GW iterations
GWTol = getoptions(options,'GWTol',1e-7); % Tolerance for Gromov-Wasserstein iteration
innerTol = getoptions(options,'innerTol',1e-4); % Tolerance for Gromov-Wasserstein iteration
maxInnerIter = getoptions(options,'maxInnerIter',500); % Max inner loop
normalizeDistances = getoptions(options,'normalizeDistances',0); % Divide distances by max
gamma = getoptions(options,'initialGuess',ones(n0,n)); % Initial guess of matching
mu = getoptions(options,'mu',100);
fixedSource = getoptions(options,'fixedSource',[]);
fixedTarget = getoptions(options,'fixedTarget',[]);

gamma(fixedSource,:) = fixedTarget; % Input data

% To give units to the regularizer
if normalizeDistances
    D = D / max(D(:));
    D0 = D0 / max(D0(:));
end

o0 = ones(n0,1); % Vectors of ones -- useful later
o = ones(n,1);

for iteration=1:maxIter % Gromov-Wasserstein iteration
    fprintf('Iteration %d...\n',iteration);
    
    oldGamma = gamma; % For convergence criterion
    
    % From expanding the square from the optimal transportation objective
    target = (.5*o0*o0'*gamma*(D.^2) + .5*(D0.^2)*gamma*o*o' - D0*gamma*D)/ (n0*n);
    
    % Kernel we'll need to project onto singly-stochastic matrices (with regularization)
    K = exp(-target/regularizer);
    
    gamma = K; % Reasonable starting point for proximal gradient
    gamma = n*bsxfun(@rdivide,gamma,gamma*o); % Normalize the starting point
    
    gamma(fixedSource,:) = fixedTarget; % Build in fixed matches
        
    initialG = g(gamma'*o0/n0); % Track regularizer g(x) to make sure it decreases
    fprintf('\tg before = %g\n',initialG);
    
    expK = K.^(1/(1+mu)); % Precomputation needed in each iteration
    for innerIteration=1:maxInnerIter
        curG = g(gamma'*o0/n0);

        if curG > initialG % just in case the step is too large -- conservative criterion
            mu = mu*2;
            fprintf('\tChanged step to %g.\n',mu);
        end
        
        marginal = gamma'*o0/n0; % The free marginal
        oldInnerGamma = gamma; % Last iteration (for convergence criterion)

        % Convenience substeps for proximal gradient formula
        P = gradG(marginal)'/n0;
        expP = exp(-n0*n*P/(1+mu));
        prod = expK.*(gamma.^(mu/(1+mu)));
        target = bsxfun(@times,expP,prod);

        % Project onto the fixed marginal constraint
        gamma = n*bsxfun(@rdivide,target,target*o);
        gamma(fixedSource,:) = fixedTarget; % Keep building in fixed matches

        % Convergence criteria
        innerChange = norm(gamma-oldInnerGamma,'fro') / norm(gamma,'fro');
        if innerIteration > 1 && innerChange < innerTol
            fprintf('\tProximal gradient took %d iterations.\n',innerIteration);
            break;
        end
        
        if mod(innerIteration,10) == 0 % Display matrix/marginals for entertainment
            subplot(1,2,2); plot(marginal); title('Marginal');
            axis([1 n 0 max(marginal)]);
            drawnow;
        end
    end
    fprintf('\tg after = %g\n',g(gamma'*o0/n0));
    
    subplot(1,2,1);imagesc(gamma); title('Matching'); drawnow;
    
    % Check for convergence
    if norm(gamma-oldGamma,'fro')/norm(gamma,'fro') < GWTol
        break
    end
end

function v = getoptions(options, name, v)

if isfield(options, name)
    v = eval(['options.' name ';']);
else
    fprintf('Setting option %s to ',name);
    disp(v(1));
    fprintf('\b');
end 