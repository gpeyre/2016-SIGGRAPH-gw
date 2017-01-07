%% Load mesh, put into convenient data structure

clear
[X,T] = readOff('../data/meshes/alligator.off');
n = size(X,1);
M = getMeshData(X,T);

%% Compute full distance matrix

fprintf('Computing pairwise geodesic distance...\n');

D = zeros(n,n);
for i=1:n
    if mod(i,100) == 0
        fprintf('Vertex %d of %d...\n',i,n);
    end
    D(:,i) = perform_fast_marching_mesh(M.vertices,double(M.triangles),i);
end
D = D + D'; % Symmetrize -- fast marching may not be symmetric (I think)

%% Design a graph

edges = [1 3;2 3; 3 4; 4 5; 5 8; 6 7; 7 8; 8 9; 9 10; 8 11; 11 12; 12 13; ...
         13 14; 14 15; 15 18; 16 17; 17 18; 18 19; 19 20; 18 21; 21 22; 22 23; ...
         23 24; 24 25; 25 26; 26 27];
     
n0 = max(edges(:));

adj = zeros(n0, n0);
for i=1:size(edges,1)
    adj(edges(i,1),edges(i,2)) = 1;
    adj(edges(i,2),edges(i,1)) = 1;
end
D0 = graphallshortestpaths(sparse(adj));

%% Compute Gromov-Wasserstein

fprintf('Optimizing regularized Gromov-Wasserstein...\n');

options = [];

options.mu = M.areaWeights;
options.display = 1;
options.regularizer = .00018;
options.plotObjective = 1;
options.maxIter = 50;

gamma = gromovWassersteinDistance(D0,D,options);

%% Illustrate map in a corny way -- just select out some points randomly

close all;
points = [1 6 13 20 23 27];
nPoints = length(points);

for i=1:nPoints
    p = points(i);

    showDescriptor(M,gamma(p,:)); colorbar off;
    title(sprintf('Target %d', i));
end

%% Write meshes

gg = gamma(points,:);
generateSoftMapMeshes(M,M,'alligator_source.obj','alligator_target.obj',...
    [],[],points,gg',.1);
