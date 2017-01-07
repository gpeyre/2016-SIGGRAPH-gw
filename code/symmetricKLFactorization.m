function A = symmetricKLFactorization(B,rank,options)
% Finds A so that B \approx A A^T

if nargin < 3
    options = [];
end

n = size(B,1);

A = getoptions(options,'initialGuess',rand(n,rank)+.5);
display = getoptions(options,'display',1);
maxFactorizationIter = getoptions(options,'maxFactorizationIter',1000);
factorizationTol = getoptions(options,'factorizationTol',1e-5);

assert(norm(B-B','fro')/norm(B,'fro') < 1e-9); % Only works when B is symmetric!

oldKL = inf; % to keep track of convergence

for iteration = 1:maxFactorizationIter
    P = A*A';
    
    if mod(iteration,50) == 0
        KL = sum(sum(B.*log(B./P)-B+P));
        change = norm(oldKL-KL,'fro')/norm(KL,'fro');
        oldKL = KL; fprintf('\tKL %d of %d = %g\n',iteration,maxFactorizationIter,KL);
        if change < factorizationTol && iteration > 5 && display
            fprintf('\tNNMF took %d iterations.\n',iteration);
            break;
        end
    end

    U = A .* ((B./P)*A);
    A = bsxfun(@rdivide,U,sqrt(sum(U)));
end

function v = getoptions(options, name, v) % from Gabriel's code

if isfield(options, name)
    v = eval(['options.' name ';']);
end 