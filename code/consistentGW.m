function [G,A,otherdata] = consistentGW(d, options)
% Consistent Gromov-Wasserstein distance
% Assumes d{} is a set of distance matrices, all of the same size.

if nargin < 2
    options = [];
end

p = length(d); % number of domains
m = size(d{1},1); % number of samples per domain

targetRank = getoptions(options,'targetRank',m);
A = ones(m*p,targetRank);
G = ones(m*p,m*p);

gwRegularizer = getoptions(options,'gwRegularizer',.005);
eta = getoptions(options,'eta',1);

beta = getoptions(options,'lowRankRegularizer',1);
delta = 1/(1+beta);

maxIter = getoptions(options,'maxIter',100);
maxGWIter = getoptions(options,'maxGWIter',500);
maxSinkhornIter = getoptions(options,'maxSinkhornIter',500);
sinkhornTol = getoptions(options,'sinkhornTol',1e-7);
GWTol = getoptions(options,'GWTol',1e-7);
display = getoptions(options,'display',1);

if display
    matrixFig = figure;
end

for iteration=1:maxIter
    for i=1:p % update blocks of G using Gromov-Wasserstein
        for j=1:p
            idxI = (i-1)*m+1;
            idxJ = (j-1)*m+1;
            
            % Current block from the nonnegative factored maps
            pairwiseI = A(idxI:(idxI+m-1),:);
            pairwiseJ = A(idxJ:(idxJ+m-1),:);
            
            V = pairwiseI*pairwiseJ'; % takes distances on J, outputs distances on I
            gamma = G(idxI:(idxI+m-1),idxJ:(idxJ+m-1)); % Output of GW iteration
            
            for subIteration=1:maxGWIter % Block projection
                oldGamma = gamma;
                
                % TODO:  FIGURE OUT RIGHT DATA
                if iteration > 1
                    K = (V.^(1-delta)).*exp(delta*d{i}*gamma*d{j}/(gwRegularizer*m*m));
                else
                    K = exp(d{i}*gamma*d{j}/(gwRegularizer*m*m));
                end
                K = (K.^eta).*(gamma.^(1-eta));
                
                % gamma = d_v * K * d_w
                v = ones(m,1);
                w = ones(m,1);
                Kw = K*w;
                for subsubIteration=1:maxSinkhornIter % Sinkhorn iteration
                    v = 1./Kw;
                    Kv = K'*v;
                    w = 1./Kv;
                    Kw = K*w;
                    
                    % Stop when margins are close to 1
                    if norm(v.*Kw-1,1)/m<sinkhornTol && norm(w.*Kv-1,1)/m<sinkhornTol
                        break;
                    end
                end
                gamma = diag(v)*K*diag(w)*m; % update Gamma 
                
                change = norm(oldGamma-gamma,'fro')/norm(gamma,'fro');
                if subIteration > 2 && change < GWTol
                    fprintf('\tGW took %g iterations.\n',subIteration);
                    break;
                end
            end
            
            % Store corresponding block of inconsistent map matrix
            G(idxI:(idxI+m-1),idxJ:(idxJ+m-1)) = gamma;
        end
    end
    
    G(isnan(G)) = 1e-8;
    G(G<1e-8) = 1e-8;
    
    if iteration == 1 % no consistency, might as well save for comparison
        fprintf('Saving maps.\n');
        otherdata.noConsistencyMaps = G;
    end
    
    % Update Gamma using pLSA-style nonnegative matrix factorization
    factorizationOptions = options;
    if iteration > 1 % Uniform matrix guess isn't good!
        factorizationOptions.initialGuess = A;
    end
    
    if display
        figure(matrixFig);
        
        if nargout == 3
            subplot(1,3,1);
            imagesc(otherdata.noConsistencyMaps); axis equal; axis off;
            title('Original maps');
        end
        
        subplot(1,3,2);
        imagesc(G); axis equal; axis off;
        title(sprintf('Inconsistent maps %d',iteration));
        drawnow;
    end
    
    A = symmetricKLFactorization((G+G')/2,targetRank,factorizationOptions);
    
    A(isnan(A)) = 1e-8;
    A(A<1e-8) = 1e-8;
    
    % Display some entertaining data
    if display
        subplot(1,3,3);
        imagesc(A*A'); axis equal; axis off;
        title(sprintf('Consistent maps %d',iteration));
        drawnow;
    end
end

function v = getoptions(options, name, v) % from Gabriel's code

if isfield(options, name)
    v = eval(['options.' name ';']);
else
    fprintf('Setting option %s to ',name);
    disp(v);
    fprintf('\b');
end 
