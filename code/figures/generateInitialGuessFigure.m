%% Load mesh, put into convenient data structure

[X,T] = readOff('../data/meshes/bull.off');
n0 = size(X,1);
M0 = getMeshData(X,T);

[X,T] = readOff('../data/meshes/dog.off');
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

% super fuzzy map at .1
% sharp map at .0005

options.mu0 = M0.areaWeights;
options.mu = M.areaWeights;
options.display = 1;
options.plotObjective = 1;
options.eta = 1;
options.maxIter = 100;
options.GWTol = 0; % force running all the iterations

o = ones(n,1);
o0 = ones(n0,1);
dmu0 = spdiags(options.mu0,0,n0,n0);
dmu = spdiags(options.mu,0,n,n);

nTests = 500;

alphas = zeros(nTests,1);
values = zeros(nTests,1);
minAlpha = .0005;
maxAlpha = .1;
for i=1:nTests
    fprintf('Test %d of %d...\n',i,nTests);
    t = rand();
    alphas(i) = minAlpha + exp((1-t)*log(minAlpha) + t*log(maxAlpha));
    options.regularizer = alphas(i);
    
    options.initialGuess = rand(n0,n)+1e-5;
    [gamma,objectives] = gromovWassersteinDistance(D0,D,options);
    
    L1 = .5*(D0.^2)*dmu0*gamma*options.mu*o';
    L2 = -D0*dmu0*gamma*dmu*D;
    L3 = .5*o0*options.mu0'*gamma*dmu*(D.^2);
    Lambda = L1+L2+L3;
    values(i) = sum(sum(bsxfun(@times,bsxfun(@times,Lambda.*gamma,options.mu0),options.mu')));
end

save initialGuessTest.mat

%% Plot experiments

semilogx(alphas,values,'.');
title('Experiments');
xlabel('alpha');
ylabel('objective');

%% Write out in a TikZ-friendly way

fid = fopen('initialGuess.txt','w');
for i=1:nTests
    fprintf(fid,'%g,%g\n',alphas(i),values(i));
end
fclose(fid);

%% Illustrate map in a corny way -- just select out some points randomly

nPoints = 3; %419 = back foot, 349 = midsection, 4 = horn
points = [419,349,4];%randi(M0.numVertices,[nPoints 1]);

for i=1:nPoints
    p = points(i);
    figure;
    
    fig = subplot(1,2,1);
    f = zeros(n0,1); f(p) = 1;
    showDescriptor(M0,f,[],[],[],fig); colorbar off; hold on;
    plot3(M0.vertices(p,1),M0.vertices(p,2),M0.vertices(p,3),'.','markersize',50,'markeredgecolor',[1 0 0]);
    title('Source');
    
    fig = subplot(1,2,2);
    showDescriptor(M,gamma(p,:)',[],[],[],fig); colorbar off;
    title('Target');
    colorbar;
end

%% Write out example meshes

options = [];
options.mu0 = M0.areaWeights;
options.mu = M.areaWeights;
options.display = 1;
options.plotObjective = 1;
options.eta = 1;
options.maxIter = 100;
options.regularizer = minAlpha;

[gamma,~] = gromovWassersteinDistance(D0,D,options);

generateSoftMapMeshes(M0,M,'minAlphaSource.obj','minAlphaTarget.obj',...
    [],[],points,gamma(points,:)',.1);

options.regularizer = .01;
[gamma,~] = gromovWassersteinDistance(D0,D,options);

generateSoftMapMeshes(M0,M,'maxAlphaSource.obj','maxAlphaTarget.obj',...
    [],[],points,gamma(points,:)',.1);