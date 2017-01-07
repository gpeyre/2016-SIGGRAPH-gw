%% Load mesh, put into convenient data structure

clear

p = 20;

%dir = 'C:\\Users\\Justin\\Documents\\Stanford\\research\\brain_geometry\\data\\CorrsBenchmark\\Data\\watertight_shrec07\\Meshes\\';
dir = '/var/scratch/justinms/MeshsegBenchmark-1.0/data/off/';

for i=0:19
    [X,T] = readOff(sprintf('%s%d.off',dir,120+i)); % 120 = octopus, 1 = human, 20 = cup
    M{i+1} = getMeshData(X,T);
end

%% Compute distances

k = 50;

samples = cell(p,1);
D = cell(p,1);

for i=1:p
    samples{i} = perform_farthest_point_sampling_mesh(M{i}.vertices, ...
        double(M{i}.triangles), [], k);
    
    D{i} = zeros(k);
    for j=1:k
        d = perform_fast_marching_mesh(M{i}.vertices,double(M{i}.triangles),samples{i}(j));
        D{i}(:,j) = d(samples{i});
    end
    D{i} = (D{i} + D{i}')/2;
end

%% Compute Gromov-Wasserstein

fprintf('Optimizing consistent Gromov-Wasserstein...\n');

options = [];
options.gwRegularizer = .03; % .0075 before
options.maxIter = 30;
options.eta = 1;
options.lowRankRegularizer = .5;
options.maxFactorizationIter = 1000;
options.sinkhornTol = 1e-5;
options.GWTol = 1e-5;

[G,A,otherdata] = consistentGW(D, options);
save /var/scratch/justinms/octopi.mat

%% Illustrate map in a corny way -- just select out some points randomly

nPoints = 10;
close all;
m = size(A,2);

for i=1:nPoints
    figure;
    
    fig = subplot(3,p+1,1);
    f = zeros(M{1}.numVertices,1);
    showDescriptor(M{1},f,[],[],[],fig); colorbar off; hold on;
    
    pt = M{1}.vertices(samples{1}(i),:);
    scatter3(pt(1),pt(2),pt(3),50,1,'filled'); title('Source');
    
    src = A(1:m,:);
    
    for j=1:p
        f = zeros(M{j}.numVertices,1);
        pt = M{j}.vertices(samples{j},:);
        
        idx = 1+(j-1)*m;
        tgt = A(idx:(idx+m-1),:);
        
        map = src*tgt';
        map2 = otherdata.noConsistencyMaps(1:m,idx:(idx+m-1));
        map3 = G(1:m,idx:(idx+m-1));
        
        fig = subplot(3,p+1,j+1);
        showDescriptor(M{j},f,[],[],[],fig); colorbar off; hold on;
        scatter3(pt(:,1),pt(:,2),pt(:,3),50,map(i,:),'filled');
        
        fig = subplot(3,p+1,j+1+(p+1));
        showDescriptor(M{j},f,[],[],[],fig); colorbar off; hold on;
        scatter3(pt(:,1),pt(:,2),pt(:,3),50,map3(i,:),'filled');
        
        fig = subplot(3,p+1,j+1+2*(p+1));
        showDescriptor(M{j},f,[],[],[],fig); colorbar off; hold on;
        scatter3(pt(:,1),pt(:,2),pt(:,3),50,map2(i,:),'filled');
    end
end