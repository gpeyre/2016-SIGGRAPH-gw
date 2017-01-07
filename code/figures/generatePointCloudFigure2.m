%% Load mesh, put into convenient data structure

[X,T] = readOff('../data/meshes/skirt.off');
n = size(X,1);
M = getMeshData(X,T);
D = distmat(X');

%% Read mesh

X0 = X(randperm(size(X,1),1000),:);
X0 = X0+randn(size(X0))*.02;
n0 = size(X0,1);
D0 = distmat(X0');

%% Compute Gromov-Wasserstein

fprintf('Optimizing regularized Gromov-Wasserstein...\n');

options = [];

options.mu = M.areaWeights;
options.display = 1;
options.regularizer = .0005;
options.plotObjective = 1;
options.maxIter = 50;

gamma = gromovWassersteinDistance(D0,D,options);

%% Illustrate map in a corny way -- just select out some points randomly

nPoints = 6;
xyz = X0;
points = randperm(size(xyz,1),nPoints);

for i=1:nPoints
    p = points(i);
    figure;
    
    fig = subplot(1,2,1);
    plot3(xyz(:,1),xyz(:,2),xyz(:,3),'k.','markersize',10); hold on;
    plot3(xyz(p,1),xyz(p,2),xyz(p,3),'r.','markersize',20);
    axis equal; axis off;
    title('Source');
    
    fig = subplot(1,2,2);
    showDescriptor(M,gamma(p,:)',[],[],[],fig); colorbar off;
    title('Target');
end

%% Try to generate maps /into/ the point cloud since they'll be easier to see

% 934 = hand
% 1827 = forehead
% 538 = skirt corner
% 214 = leg
% 33 = foot
% 1524 = shoulder
points = [934 1827 538 214 33 1524];
nPoints = length(points);

dists = gamma(:,points);
dists = bsxfun(@rdivide,dists,max(dists));

colors = [192 0 0; 237 125 49; 255 192 0;68 114 196;255 0 255;112 48 160];
colors = colors/255;

vtxColors = zeros(n0,3);
for i=1:n0
    [t,idx] = max(dists(i,:));
    vtxColors(i,:) = colors(idx,:)*t + [242 242 242]*(1-t)/255;
end

figure('color',[0 0 0]);
scatter3(xyz(:,1),xyz(:,2),xyz(:,3),40,vtxColors,'filled');
axis equal;
axis off;
cameratoolbar;
