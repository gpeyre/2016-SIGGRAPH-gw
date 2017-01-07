function [gamma,objectives] = gromovWassersteinDistance(D0,D,options)

n0 = size(D0,1); % Size of metric matrices
n = size(D,1);

mu0 = getoptions(options, 'mu0', ones(n0,1)/n0); % Measure on the source domain
mu = getoptions(options, 'mu', ones(n,1)/n); % Measure on the target domain

regularizer = getoptions(options, 'regularizer', .005); % alpha
initialRegularizer = getoptions(options, 'initialRegularizer', regularizer); % for multiscale
regularizerChangeRate = getoptions(options, 'regularizerChangeRate', 1); % how fast to approach alpha

maxIter = getoptions(options, 'maxIter', 40); % Maximum number of GW iterations
display = getoptions(options, 'display', 0); % Show map matrix in each iteration
relativeRegularizer = getoptions(options, 'relativeRegularizer', 0); % Alpha relative to matrix
GWTol = getoptions(options,'GWTol',1e-7); % Tolerance for Gromov-Wasserstein iteration
sinkhornTol = getoptions(options,'sinkhornTol',1e-6); % Tolerance for Sinkhorn iteration
maxSinkhornIter = getoptions(options,'maxSinkhornIter',500); % Max inner loop
toSubtract = getoptions(options,'toSubtract',sparse(n0,n)); % Experimental for symmetry detection
normalizeDistances = getoptions(options,'normalizeDistances',1); % Divide distances by max
plotObjective = getoptions(options,'plotObjective',0); % Plot objective in each iteration
redrawFrequency = getoptions(options,'redrawFrequency',1);
gamma = getoptions(options,'initialGuess',ones(n0,n)); % Initial guess of matching
eta = getoptions(options,'eta',1);

stencil = getoptions(options,'stencil',ones(n0,n));

% Partial matching
partialSource = getoptions(options,'partialSource',0);
partialTarget = getoptions(options,'partialTarget',0);

% Rescale initial guess
if ~partialSource && ~partialTarget
    gamma = bsxfun(@rdivide,gamma,gamma*mu);
end

% Some mesh code gives sparse area weights, so here I'll make them dense.
% Also, scale to integrate to 1.
mu0 = full(mu0);  mu0 = mu0 / sum(mu0);
mu = full(mu);    mu = mu / sum(mu);

% Diagional matrix of area weights for convenience
dmu0 = spdiags(mu0,0,n0,n0);
dmu = spdiags(mu,0,n,n);

% To give units to the regularizer
if normalizeDistances
    scale = max(max(D(:)),max(D0(:)));
    D = D / scale;
    D0 = D0 / scale;
end

% Array of objective values, useful for illustrating convergence
if nargout > 1 || plotObjective
    objectives = [];
end

% Create figure windows if desired
if display && plotObjective
    displayFigure = @() subplot(1,2,1);
    objectiveFigure = @() subplot(1,2,2);
elseif display, df = figure; displayFigure = @() figure(df); 
elseif plotObjective, of = figure; objectiveFigure = @() figure(of); objectives = [];
end

% Precompute some useful matrices
leftProd = D0*dmu0;
rightProd = dmu*D;

if partialSource
    psLeft = .5*ones(n0,1)*mu0';
    psRight = dmu*D.^2;
end

if partialTarget
    ptLeft = .5*D0.^2*dmu0;
    ptRight = mu*ones(1,n);
end

alpha = initialRegularizer;
for iteration=1:maxIter % Gromov-Wasserstein iteration
    target = leftProd*gamma*rightProd; % Will go in exponent
    oldGamma = gamma; % For convegence criterion
    
    % For partial matching, the kernel matrix is different
    if partialSource, target=target-psLeft*gamma*psRight; end    
    if partialTarget, target=target-ptLeft*gamma*ptRight; end
    
    % Scale regularizer by mean of target so that it's relative
    if relativeRegularizer, reg = mean(target(:))*alpha;
    else reg = alpha; end
    
    alpha = regularizerChangeRate*regularizer + (1-regularizerChangeRate)*alpha;
    
    % Subtracting min(target) doesn't affect Sinkhorn but improves numerics
    K = exp((target - toSubtract - max(target(:)))/reg); % Evaluate f_alpha 
    
    K = K.*stencil; % to zero out bad matches
    
    if eta ~= 1
        K = (K.^eta).*(gamma.^(1-eta));
    end
    
    % Add in ground truth
%    K(groundTruthRowIndices,:) = groundTruthRows;
%    K(:,groundTruthColIndices) = groundTruthCols;
    
    if ~partialSource && ~partialTarget % Sinkhorn projection
        if iteration == 1 % Warm start in future iterations
            v = ones(n0,1);
            w = ones(n,1);
        end
        Kw = K*(mu.*w);
        for i=1:maxSinkhornIter % check this...
            v = 1./Kw;
            Kv = K'*(mu0.*v);
            w = 1./Kv;
            Kw = K*(mu.*w);

            % Stop when margins are close to 1
            if norm(mu0.*(v.*Kw-1),1)<sinkhornTol && norm(mu.*(w.*Kv-1),1)<sinkhornTol
                break;
            end
        end
        gamma = diag(v)*K*diag(w);
    elseif partialSource && ~partialTarget % Partial-to-full matching
        gamma = bsxfun(@rdivide,K,K*mu);
    elseif ~partialSource && partialTarget
        gamma = bsxfun(@rdivide,K,mu0'*K);
    else gamma = K / (mu0'*K*mu);
    end
    
    % Add in ground truth
%    gamma(groundTruthRowIndices,:) = groundTruthRows;
%    gamma(:,groundTruthColIndices) = groundTruthCols;
    
    % Show mapping matrix for entertainment
    if display 
        displayFigure();
        imagesc(gamma); colorbar;
        title(sprintf('Iteration %d',iteration));
        if mod(iteration,redrawFrequency) == 0, drawnow; end 
    end
    
    % Evaluate objective function
    if plotObjective || nargout > 1
        % Do the full objective even though some of it is constant
        target = leftProd*gamma*rightProd - .5*ones(n0,1)*mu0'*gamma*dmu*D.^2 ...
            - .5*D0.^2*dmu0*gamma*mu*ones(1,n);
        gammaBar = exp(target/reg);
        
        % Generalized KL divergence
        objective = sum(sum(dmu0*(gamma.*log(gamma./gammaBar) - gamma + gammaBar)*dmu));
        objectives(end+1) = reg*objective; % Technically it's alpha * KL(gamma|F(gamma))
    end
    
    % Generate plot of the objective function
    if plotObjective
        objectiveFigure();
        plot(objectives); title('GW Convergence');
        xlabel('Iteration'); ylabel('Objective');
        if mod(iteration,redrawFrequency) == 0, drawnow; end
    end
     
    % Stopping conditions
    change = norm(oldGamma-gamma,'fro')/norm(gamma,'fro');
    if iteration > 2 && change < GWTol && abs(alpha-regularizer) < 1e-5
        fprintf('\tGW took %g iterations.\n', iteration);
        break;
    end
end

function v = getoptions(options, name, v)

if isfield(options, name)
    v = eval(['options.' name ';']);
else
    fprintf('Setting option %s to ',name);
    if ~isempty(v)
        disp(v(1));
        fprintf('\b');
    else
        fprintf('[]\n');
    end
end 