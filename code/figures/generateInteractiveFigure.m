%% Load mesh, put into convenient data structure

clear
close all

[X,T] = readOff('../data/meshes/star_subdivided.off');
X = [X(:,1) X(:,3) X(:,2)]; % easier camera angle...
% X(:,3) = X(:,3)*-1; 
X(:,2) = X(:,2)*-1;
n0 = size(X,1);
M0 = getMeshData(X,T);

[X,T] = readOff('../data/meshes/trim-star-simplified.off');
n = size(X,1);
M = getMeshData(X,T);

%% Compute full distance matrix

fprintf('Computing pairwise geodesic distance...\n');

D0 = distmat(M0.vertices');
D = distmat(M.vertices');

%% Compute Gromov-Wasserstein

fprintf('Optimizing regularized Gromov-Wasserstein...\n');

options = [];

options.mu0 = M0.areaWeights;
options.mu = M.areaWeights;
options.display = 1;
options.regularizer = .008;
options.plotObjective = 1;
options.maxIter = 50;

gamma = gromovWassersteinDistance(D0,D,options);
close all;

%% Show the user the highest entropy point, ask for a target

entropies = -gamma.*log(gamma);
entropyIntegrals = entropies*options.mu;
[~,sourcePoint] = max(entropyIntegrals);

f = zeros(n0,1); f(sourcePoint) = 1;
showDescriptor(M0,f); colorbar off; hold on;
plot3(M0.vertices(sourcePoint,1),M0.vertices(sourcePoint,2),M0.vertices(sourcePoint,3),'.','markersize',50,'markeredgecolor',[1 0 0]);
title('Source of largest entropy');
    
fig = showDescriptor(M,gamma(sourcePoint,:)'); colorbar off;
title('Target of largest entropy');

targetPoint = selectPoint(fig,M);
close all;

%% Write meshes

gg = gamma(sourcePoint,:);
generateSoftMapMeshes(M0,M,'initial_map_source.obj','initial_map_target.obj',...
    [],[],sourcePoint,gg',1);

%% Now, re-map with a stencil

fprintf('Optimizing regularized Gromov-Wasserstein...\n');

options = [];

options.mu0 = M0.areaWeights;
options.mu = M.areaWeights;
options.display = 1;
options.regularizer = .0005;
options.plotObjective = 1;
options.maxIter = 50;

stencil = ones(n0,n);
stencil(sourcePoint,:) = 1e-5;
stencil(:,targetPoint) = 1e-5;
options.stencil = stencil;

gamma2 = gromovWassersteinDistance(D0,D,options);

%% Write meshes

gg = gamma2(points,:);
generateSoftMapMeshes(M0,M,'point_map_source.obj','point_map_target.obj',...
    sourcePoint,targetPoint,points,sqrt(gg)',1);

%% Illustrate map in a corny way -- just select out some points randomly

close all

% 5 = left tip
% 179 = mid left top
points = [5 179 348 352 51];%randperm(M0.numVertices,nPoints);
nPoints = length(points);

for i=1:nPoints
    p = points(i);
    figure;
    
    fig = subplot(1,2,1);
    f = zeros(n0,1); f(p) = 1;
    showDescriptor(M0,f,[],[],[],fig); colorbar off; hold on;
    plot3(M0.vertices(p,1),M0.vertices(p,2),M0.vertices(p,3),'.','markersize',50,'markeredgecolor',[1 0 0]);
    title('Source');
    
    fig = subplot(1,2,2);
    showDescriptor(M,gamma2(p,:)',[],[],[],fig); colorbar off;
    title('Target');
end

fprintf('Press any key to continue...\n');
pause
close all

%% Now get rid of symmetry plane

front0 = M0.vertices(:,3) > mean(M0.vertices(:,3));
front1 = M.vertices(:,3) > mean(M.vertices(:,3));
stencil2 = stencil;
stencil2(front0,~front1) = 1e-8;

%% Now, re-map with a stencil

fprintf('Optimizing regularized Gromov-Wasserstein...\n');

options = [];

options.mu0 = M0.areaWeights;
options.mu = M.areaWeights;
options.display = 1;
options.regularizer = .00025;
options.plotObjective = 1;
options.maxIter = 20;
options.stencil = stencil2;

gamma3 = gromovWassersteinDistance(D0,D,options);

%% Illustrate map in a corny way -- just select out some points randomly

close all

for i=1:nPoints
    p = points(i);
    figure;
    
    fig = subplot(1,2,1);
    f = zeros(n0,1); f(p) = 1;
    showDescriptor(M0,f,[],[],[],fig); colorbar off; hold on;
    plot3(M0.vertices(p,1),M0.vertices(p,2),M0.vertices(p,3),'.','markersize',50,'markeredgecolor',[1 0 0]);
    title('Source');
    
    fig = subplot(1,2,2);
    showDescriptor(M,gamma3(p,:)',[],[],[],fig); colorbar off;
    title('Target');
end

gg = gamma3(points,:);
generateSoftMapMeshes(M0,M,'reflection_map_source.obj','reflection_map_target.obj',...
    sourcePoint,targetPoint,points,sqrt(gg)',1);