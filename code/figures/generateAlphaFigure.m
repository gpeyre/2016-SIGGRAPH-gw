%% Load mesh, put into convenient data structure

[X,T] = readOff('../data/meshes/octopus1.off');
n0 = size(X,1);
M0 = getMeshData(X,T);

[X,T] = readOff('../data/meshes/octopus2.off');
n = size(X,1);
M = getMeshData(X,T);

%% Compute full distance matrix

fprintf('Computing pairwise geodesic distance...\n');

D0 = zeros(n0,n0);
for i=1:n0
    D0(:,i) = perform_fast_marching_mesh(M0.vertices,double(M0.triangles),i);
end
D0 = D0 + D0'; % Symmetrize -- fast marching may not be symmetric (I think)

D = zeros(n,n);
for i=1:n
    D(:,i) = perform_fast_marching_mesh(M.vertices,double(M.triangles),i);
end
D = D + D';

%% Compute Gromov-Wasserstein

fprintf('Optimizing regularized Gromov-Wasserstein...\n');

options = [];

options.mu0 = M0.areaWeights;
options.mu = M.areaWeights;
options.display = 0;
options.plotObjective = 0;
options.maxIter = 200;
options.GWTol = 0; % run all iterations

options.eta = 1;

alphas = logspace(log(.0007)/log(10),log(1)/log(10),100);

gammas = [];
for i=1:length(alphas)
    fprintf('alpha %d of %d...\n',i,length(alphas));
    options.regularizer = alphas(i);
    gammas{i} = gromovWassersteinDistance(D0,D,options);
end

save alphatest.mat

%% Plot convergence with different etas

x = alphas(1:length(gammas));
y = zeros(1,length(gammas));
entropies = zeros(1,length(gammas));

o = ones(n,1);
o0 = ones(n0,1);
dmu0 = spdiags(options.mu0,0,n0,n0);
dmu = spdiags(options.mu,0,n,n);

for p=1:length(gammas)
    fprintf('p=%d of %d...\n',p,length(gammas));
    L1 = .5*(D0.^2)*dmu0*gammas{p}*options.mu*o';
    L2 = -D0*dmu0*gammas{p}*dmu*D;
    L3 = .5*o0*options.mu0'*gammas{p}*dmu*(D.^2);
    Lambda = L1+L2+L3;
    y(p) = sum(sum(bsxfun(@times,bsxfun(@times,Lambda.*gammas{p},options.mu0),options.mu')));
    entropies(p) = sum(sum(bsxfun(@times,bsxfun(@times,log(gammas{p}).*gammas{p},options.mu0),options.mu')));
end

%%

loglog(alphas,y);
hold on
loglog(alphas,entropies)

legend('GW','Negative entropy');
xlabel('alpha');
ylabel('objective');

%% Write out in a TikZ-friendly way

fid = fopen('alphatest.txt','w');
for i=1:length(alphas)
    fprintf(fid,'%g,%g,%g\n',alphas(i),full(y(i)),full(entropies(i)));
end
fclose(fid);
