%% Load mesh, put into convenient data structure

% [X,T] = readOff('../data/meshes/chair1_2k.off');
%[X,T] = readOff('../data/meshes/scape2_remeshed.off');
%[X,T] = readOff('../data/meshes/18_1024.off');
%[X,T] = readOff('../data/meshes/mug1.off');
% [X,T] = readOff('../data/meshes/hand1.off');
% [X,T] = readOff('../data/meshes/skirt.off');
[X,T] = readOff('../data/meshes/skirt.off');
n0 = size(X,1);
M0 = getMeshData(X,T);

% [X,T] = readOff('../data/meshes/chair2_2k.off');
%[X,T] = readOff('../data/meshes/3_1024.off');
%[X,T] = readOff('../data/meshes/stick_human2.off');
%[X,T] = readOff('../data/meshes/mug3.off');
% [X,T] = readOff('../data/meshes/hand2.off');
% [X,T] = readOff('../data/meshes/robot.off');
[X,T] = readOff('../data/meshes/18_1024.off');
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

% For chair
% 194 = front of seat
% 64 = mid front leg
% 1647 = back seat
% 878 = top left
% 1417 = mid seat
% 1723 = mid right back
% points = [194 64 1647 878 1417 1723];

% For humans
% 1479 = forehead
% 599 = right hand
% 1042 = left elbow
% 947 = side
% 391 = right thigh
% 167 = left ankle
% 1185 = chest
% points = [1479 599 1042 947 391 167];% 1185];

% Human --> wood thing
% 563 = hand
% 187 = shoulder
% 393 = side
% 847 = leg
% 122 = neck
% 907 = ankle
% points = [563 187 393 847 122 907];

% Mugs
% 795 = upper handle
% 877 = lower handle
% 1784 = inside
% 743 = lip
% 300 = lower front
% 495 = closer to handle
% points = [795 877 1784 743 300 495];

% Hand
% 122 = mid thumb
% 415 = palm
% 520 = index finger tip
% 697 = middle finger base
% 334 = mid 4th finger
% 746 = space between fingers
% points = [122 415 520 697 334 746];

% Skirt
% 934 = hand
% 1827 = forehead
% 538 = skirt corner
% 214 = leg
% 33 = foot
% 1524 = shoulder
points = [934 1827 538 214 33 1524];

% points = randperm(M0.numVertices,nPoints);

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
generateSoftMapMeshes(M0,M,'source_skirt2.obj','target_skirt2.obj',...
    [],[],points,gg',.1);
