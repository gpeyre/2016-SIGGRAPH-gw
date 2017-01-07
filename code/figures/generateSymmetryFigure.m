%% Load mesh, put into convenient data structure

[X,T] = readOff('../data/meshes/scape2_remeshed.off');
n0 = size(X,1);
M0 = getMeshData(X,T);

[X,T] = readOff('../data/meshes/3_1024.off');
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
options.regularizer = .0005; % for humans
options.plotObjective = 1;
options.maxIter = 50;

flipGamma = gromovWassersteinDistance(D0,D,options);

%% Compute symmetric Gromov-Wasserstein

options2 = options;
options2.toSubtract = flipGamma*.000001;
gamma = gromovWassersteinDistance(D0,D,options2);

%% Compute bad Gromov-Wasserstein

% options3 = options;
% options3.toSubtract = flipGamma;
% gammaBad = gromovWassersteinDistance(D0,D,options3);

%% Illustrate map in a corny way -- just select out some points randomly

nPoints = 6;

% For humans
% 1479 = forehead
% 599 = right hand
% 1042 = left elbow
% 947 = side
% 391 = right thigh
% 167 = left ankle
% 1185 = chest
points = [1479 599 1042 947 391 167];% 1185];

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

%% Write meshes

gg = gamma(points,:);
generateSoftMapMeshes(M0,M,'source_symmetric_humans.obj','target_symmetric_humans.obj',...
    [],[],points,gg',.1);

gg = flipGamma(points,:);
gg = bsxfun(@rdivide,gg,max(gg')');
generateSoftMapMeshes(M0,M,'source_symmetric_humans_orig.obj','target_symmetric_humans_orig.obj',...
    [],[],points,gg',.1);
