clear;
n = 30;
k = 5;
A = rand(n,n); A = A+A';
close all;

% for i=1:100 % Sinkhorn -- algorithm works with and without this
%     A = bsxfun(@rdivide,A,sum(A,1));
%     A = bsxfun(@rdivide,A,sum(A,2));
% end

U = rand(n,k);
objs = [];
const = [];

o = ones(n,1);

for i=1:3000
%     U = U .* ((A./(U*U'))*U) ./ (o*o'*U);
%     U = sqrt(sum(sum(A)))* U / norm(U'*o);

    W = U .* ((A./(U*U'))*U);
    U = W * diag(1 ./ sqrt(W'*o));

    subplot(1,3,1); imagesc(U); axis equal; axis off;
    
    title('Nonnegative factor');
    
    objs(end+1) = sum(sum(A.*log(A./(U*U')))) + norm(U'*o)^2;
    subplot(1,3,2); plot(objs);
    title('Objective value');
    
    subplot(1,3,3);
    const(end+1) = norm(sum(U*U',2)-sum(A,2))/norm(sum(A,2));
    semilogy(const);
    title('Constraint violation');
    
    if mod(i,20)==0, drawnow; end
end