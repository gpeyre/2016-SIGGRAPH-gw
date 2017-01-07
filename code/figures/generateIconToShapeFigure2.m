im = double(imread('../data/icon2.png')==0);

% Number off the points
n = 0;
r = size(im,1);
c = size(im,2);
for i=1:r
    for j=1:c
        if im(i,j)
            n = n+1;
            im(i,j) = n;
        end
    end
end

% As a hack for distances, I'll just triangulate this thing
X = [];
T = [];
for i=1:(r-1)
    for j=1:(c-1)
        if im(i,j)
            X(im(i,j),:) = [i j];
            if im(i+1,j) && im(i,j+1)
                T = [T ; im(i,j) im(i,j+1) im(i+1,j)];
            end
        end
        
        if im(i,j+1) && im(i+1,j) && im(i+1,j+1)
            T = [T; im(i,j+1) im(i+1,j) im(i+1,j+1)];
        end
    end
end
M0 = getMeshData([X zeros(size(X,1),1)],T);
n0 = size(X,1);

%% Compute distances on source

D0 = zeros(n0,n0);
for i=1:n0
    if mod(i,100) == 0
        fprintf('Vertex %d of %d...\n',i,n0);
    end
    D0(:,i) = perform_fast_marching_mesh(M0.vertices,double(M0.triangles),i);
end
D0 = D0 + D0'; % Symmetrize -- fast marching may not be symmetric (I think)

%% Load mesh, put into convenient data structure

[X,T] = readOff('../data/meshes/skirt.off');
n = size(X,1);
M = getMeshData(X,T);

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

options.mu = M.areaWeights;
options.display = 1;
options.regularizer = .00002;
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

%% Write out maps --- in reverse direction

% 934 = hand
% 1827 = forehead
% 538 = skirt corner
% 214 = leg
% 33 = foot
% 1524 = shoulder
points = [934 1827 538 214 33 1524];

gg = gamma(:,points)';
generateSoftMapMeshes(M,M0,'source_icon.obj','target_icon2.obj',...
    [],[],points,gg',.1);