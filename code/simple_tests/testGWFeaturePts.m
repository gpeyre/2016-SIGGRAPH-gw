dir = '/Users/vokim/Research/DBShapes3D/DBShapesCode/bin/AnalyzedDBs/GWTestDir/';
mesh = 'mesh000.off';
mesh0 = 'mesh001.off';
type = 'EvenFast_0.1';
pts= [dir 'PerModel/' mesh '/' type '.pts'];
pts0 = [dir 'PerModel/' mesh0 '/' type '.pts'];
dists= [dir 'PerModel/' mesh '/' type '.dists'];
dists0 = [dir 'PerModel/' mesh0 '/' type '.dists'];
gwReg = 0.005;
gwIter = 200;
mapFile = [dir 'Pairwise/' mesh '_' mesh0 '/' type '.gw'];
[gamma, D, D0] = runGW(dists, dists0, mapFile, gwReg, gwIter, 1);



if false
%% Load mesh, put into convenient data structure

% [X,T] = readOff('../data/meshes/moomoo_s0.off');
%[X,T] = readOff('../data/meshes/octopus1.off');
%pointsType = 'feature';
%pointsType = 'even';
%pointsType = 'even5';
pointsType = 'even16';
dir = '/Users/vokim/Research/DBShapes3D/DBShapesCode/bin/AnalyzedDBs/GWTestDir/';
datadir = '/Users/vokim/Research/DBShapes3D/Data/CorrsBenchmark_BIM/Data/SCAPE/Meshes/';
mname0 = 'mesh000.off';
mname = 'mesh001.off';
[X,T] = readOff([datadir, mname0]);
[pts0, vtx0] = readPts([dir 'PerModel/' mname0 '/' pointsType '.pts'], T);
n0 = size(X,1);
M0 = getMeshData(X,T);

% [X,T] = readOff('../data/meshes/moomoo_s0_d.off');
%[X,T] = readOff('../data/meshes/octopus2.off');
[X,T] = readOff([datadir mname]);
[pts, vtx] = readPts([dir 'PerModel/' mname '/' pointsType '.pts'], T);

n = size(X,1);
M = getMeshData(X,T);

fprintf('Computing pairwise geodesic distance...\n');
%% Compute full distance matrix
%D0 = zeros(n0,n0);
%for i=1:n0
%    D0(:,i) = perform_fast_marching_mesh(M0.vertices,double(M0.triangles),i);
%end
%D0 = D0 + D0'; % Symmetrize -- fast marching may not be symmetric (I think)

%D = zeros(n,n);
%for i=1:n
%    D(:,i) = perform_fast_marching_mesh(M.vertices,double(M.triangles),i);
%end
%D = D + D';
%pointsDists = false;

%% Compute distance matrix for feature points
nvtx0 = size(vtx0, 2);
nvtx = size(vtx, 2);
pointsDists = true;
recomputeDists = true;

if recomputeDists
    D0 = zeros(nvtx0, nvtx0);
    D = zeros(nvtx, nvtx);
    for i=1:nvtx0
        row = perform_fast_marching_mesh(M0.vertices,double(M0.triangles), vtx0(i));
        D0(i,:) = row(vtx0);
    end
    for i=1:nvtx
        row = perform_fast_marching_mesh(M.vertices,double(M.triangles), vtx(i));
        D(i,:) = row(vtx);
    end

    % Write or distances
    writeMtx(D, [dir 'PerModel/' mname '/geodesic_' pointsType '.dists']);
    writeMtx(D0, [dir 'PerModel/' mname0 '/geodesic_' pointsType '.dists']);
else % read distances
    D = readMtx([dir 'PerModel/' mname '/geodesic_' pointsType '.dists']);
    D0 = readMtx([dir 'PerModel/' mname0 '/geodesic_' pointsType '.dists']);
end
%size(D)
%size(D0)


%% Compute Gromov-Wasserstein

fprintf('Optimizing regularized Gromov-Wasserstein...\n');

options = [];
if pointsDists
    options.mu0 = ones(nvtx0,1) / nvtx0;
    options.mu = ones(nvtx,1) / nvtx;
else
    options.mu0 = M0.areaWeights;
    options.mu = M.areaWeights;
end
options.display = 1;
options.regularizer = .005; % .0007 matches single tentacle, .005 matches some, .008 middle
options.plotObjective = 1;
options.eta = 1;
options.maxIter = 50;

% Partial-to-full test -- broken!
% options.regularizer = .1;
% options.partialSource = 1;

[gamma,objectives] = gromovWassersteinDistance(D0,D,options);

%% Illustrate map in a corny way -- just select out some points randomly

%nPoints = 2;
%points = [277 97];

%for i=1:nPoints
%    p = points(i);
%    figure;
%    
%    fig = subplot(1,2,1);
%    f = zeros(n0,1); f(p) = 1;
%    showDescriptor(M0,f,[],[],[],fig); colorbar off; hold on;
%    plot3(M0.vertices(p,1),M0.vertices(p,2),M0.vertices(p,3),'.','markersize',50,'markeredgecolor',[1 0 0]);
%    title('Source');
%    
%    fig = subplot(1,2,2);
%    showDescriptor(M,gamma(p,:)',[],[],[],fig); colorbar off;
%    title('Target');
%end

%% Write GW distances
writeMtx(gamma, [dir 'Pairwise/' mname0 '_' mname '/' pointsType '.gw']);


%% Write meshes

%gg = gamma(points,:);
%gg(1,:) = 0;
%generateSoftMapMeshes(M0,M,'source_octopus.obj',sprintf('target_octopus1_%g.obj',options.regularizer),...
%    [],[],points,gg',.1);

%gg = gamma(points,:);
%gg(2,:) = 0;
%generateSoftMapMeshes(M0,M,'source_octopus.obj',sprintf('target_octopus2_%g.obj',options.regularizer),...
%    [],[],points,gg',.1);
end

close all