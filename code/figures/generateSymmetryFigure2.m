%% Load mesh, put into convenient data structure

[X,T] = readOff('../data/meshes/mouse2.off');
n0 = size(X,1);
M0 = getMeshData(X,T);

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

%% Compute Gromov-Wasserstein

fprintf('Optimizing regularized Gromov-Wasserstein...\n');

options = [];

options.mu0 = M0.areaWeights;
options.mu = M0.areaWeights;
options.display = 1;
options.regularizer = .001; % for humans
options.plotObjective = 1;
options.maxIter = 50;
options.toSubtract = eye(n0)*.001;

gamma = gromovWassersteinDistance(D0,D0,options);

%% Illustrate map in a corny way -- just select out some points randomly

nPoints = 6;

% 243 = ear
% 1208 = palm
% 1505 = elbow
% 1621 = belly
% 1434 = knee
% 1264 = foot (flipper?)
points = [243 1208 1505 1621 1434 1264];

for i=1:nPoints
    p = points(i);
    figure;
    
    fig = subplot(1,2,1);
    f = zeros(n0,1); f(p) = 1;
    showDescriptor(M0,f,[],[],[],fig); colorbar off; hold on;
    plot3(M0.vertices(p,1),M0.vertices(p,2),M0.vertices(p,3),'.','markersize',50,'markeredgecolor',[1 0 0]);
    title('Source');
    
    fig = subplot(1,2,2);
    showDescriptor(M0,gamma(p,:)',[],[],[],fig); colorbar off;
    title('Target');
end

%% Write meshes

gg = gamma(points,:);
generateSoftMapMeshes(M0,M0,'source_symmetric_mouse.obj','target_symmetric_mouse.obj',...
    [],[],points,gg',.1);
