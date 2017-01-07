clear

shrecdir = 'C:\Users\Justin\Documents\Stanford\research\brain_geometry\data\CorrsBenchmark\Data\watertight_shrec07\Meshes\';

% mesh1 = 381; % cow --> camel
% mesh2 = 393;

mesh1 = 386; % bull --> giraffe
mesh2 = 390;

% mesh1 = 393; % camel --> pig
% mesh2 = 383;

% mesh1 = 396;
% mesh2 = 390;

[X,T] = readOff(sprintf('%s\\%d.off',shrecdir,mesh1));
X = X/max(abs(X(:)));
Mfull{1} = getMeshData(X,T);

[X,T] = readOff(sprintf('%s\\%d.off',shrecdir,mesh2));
X = [X(:,3) X(:,2) X(:,1)];
X = X/max(abs(X(:)));
Mfull{2} = getMeshData(X,T);

map = 1+dlmread(sprintf('../data/bim/%d_%d.map',mesh1,mesh2));

showDescriptor(Mfull{2},Mfull{2}.vertices(:,3));
f = Mfull{2}.vertices(map,3);
showDescriptor(Mfull{1},f);

%% Simplified

[X,T] = readOff('../data/meshes/bull.off');
n0 = size(X,1);
M0 = getMeshData(X,T);

[X,T] = readOff('../data/meshes/giraffe.off');
n = size(X,1);
M = getMeshData(X,T);

%% Compute full distance matrix

fprintf('Computing pairwise geodesic distance...\n');

D0 = zeros(n0,n0);
for i=1:n0
    if mod(i,100) == 0
        fprintf('Vertex %d of %d...\n',i,n0);
    end
    D0(:,i) = perform_fast_marching_mesh(M0.vertices,double(M0.triangles),i);
end
D0 = D0 + D0'; % Symmetrize -- fast marching may not be symmetric (I think)

D = zeros(n,n);
for i=1:n
    if mod(i,100) == 0
        fprintf('Vertex %d of %d...\n',i,n);
    end
    D(:,i) = perform_fast_marching_mesh(M.vertices,double(M.triangles),i);
end
D = D + D';

%% Compute Gromov-Wasserstein

fprintf('Optimizing regularized Gromov-Wasserstein...\n');

options = [];

options.mu0 = M0.areaWeights;
options.mu = M.areaWeights;
options.display = 1;
options.regularizer = .0005; % for all examples
options.plotObjective = 1;
options.maxIter = 50;

gamma = gromovWassersteinDistance(D0,D,options);

%% Illustrate map in a corny way -- just select out some points randomly

nPoints = 6;
points = randperm(M0.numVertices,nPoints);

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
end

%% Generate figure

showDescriptor(Mfull{2},Mfull{2}.vertices(:,3));
f = Mfull{2}.vertices(map,3);
showDescriptor(Mfull{1},f);

%% Choose source points

sparseSource = [327 168 411 490 212]; % bottom of leg to the top
denseSource = zeros(size(sparseSource));
sparseTarget = [];
denseTarget = [];

for i=1:length(sparseSource)
    p = M0.vertices(sparseSource(i),:);
    
    d = bsxfun(@minus,Mfull{1}.vertices,p);
    d = sum(d.^2,2);
    [~,denseSource(i)] = min(d);
    
    [~,sparseTarget(i)] = max(gamma(sparseSource(i),:));
    denseTarget(i) = map(denseSource(i));
    
    showDescriptor(M,gamma(sparseSource(i),:)');
end

%% Show meshes

close all

[Xout,Tout,index] = addSpheres(M0.vertices,M0.triangles,sparseSource,.1*ones(size(sparseSource)));
MM.vertices = Xout;
MM.triangles = Tout;
MM.numVertices = size(Xout,1);
showDescriptor(MM,index);
writeTexturedObj('bimSource.obj', MM, index);

[Xout,Tout,index] = addSpheres(M.vertices,M.triangles,sparseTarget,.1*ones(size(sparseSource)));
MM.vertices = Xout;
MM.triangles = Tout;
MM.numVertices = size(Xout,1);
showDescriptor(MM,index);
writeTexturedObj('gwTarget.obj', MM, index);

[Xout,Tout,index] = addSpheres(Mfull{1}.vertices,Mfull{1}.triangles,denseSource,.1*ones(size(denseSource)));
MM.vertices = Xout;
MM.triangles = Tout;
MM.numVertices = size(Xout,1);
showDescriptor(MM,index);

[Xout,Tout,index] = addSpheres(Mfull{2}.vertices,Mfull{2}.triangles,denseTarget,.1*ones(size(denseSource)));
MM.vertices = Xout;
MM.triangles = Tout;
MM.numVertices = size(Xout,1);
showDescriptor(MM,index);
writeTexturedObj('bimTarget.obj', MM, index);