clear

path = 'C:\\Users\\Justin\\Documents\\Stanford\\research\\brain_geometry\\data\\CorrsBenchmark\\Data\\watertight_shrec07\\Meshes';

for i=1:20
    [X,T] = readOff(sprintf('%s\\%d.off',path,i+40));
    M{i} = getMeshData(X,T);
end

%% Sample points

k = 100;

D = cell(length(M),1);
samples = cell(length(M),1);

for i=1:length(M)
    fprintf('Processing shape %d...\n',i);
     pts = randi(M{i}.numVertices);
     
     D{i} = inf*ones(M{i}.numVertices,k);
     D{i}(:,1) = perform_fast_marching_mesh(M{i}.vertices,double(M{i}.triangles),pts);
     
     for j=2:k
         minDists = min(D{i},[],2);
         [~,pts(j)] = max(minDists);
         D{i}(:,j) = perform_fast_marching_mesh(M{i}.vertices,double(M{i}.triangles),pts(j));
     end
     D{i} = D{i}(pts,:);
     D{i} = .5*(D{i}+D{i}');
     samples{i} = pts;
end

%% Normalize

for i=1:length(M)
    D{i} = D{i} / max(D{i}(:));
end

%% Consistent maps

options = [];
options.gwRegularizer = .001;%was .002
options.maxIter = 20;
options.eta = 1;
options.lowRankRegularizer = 1;
options.maxFactorizationIter = 100;
options.sinkhornTol = 1e-5;
options.GWTol = 1e-5;
options.targetRank = 15; % left-right symmetry? (15 worked well in one iter)
options.maxGWIter = getoptions(options,'maxSinkhornIter',50);

[G,A,otherdata] = consistentGW(D, options);

%% Illustrate map in a corny way -- just select out some points randomly

nPoints = 1;
close all;

source = 1;
points =  2;%randi(k,nPoints,1);%floor(points/k+1);

for p=1:nPoints
    figure;
    i = points(p);
    
    fig = subplot(3,length(M)+1,1);
    f = zeros(M{source}.numVertices,1);
    showDescriptor(M{source},f,[],[],[],fig); colorbar off; hold on;
    
    pt = M{source}.vertices(samples{source}(i),:);
    scatter3(pt(1),pt(2),pt(3),50,1,'filled'); title('Source');
    
    sourceRange = ((source-1)*k+1):source*k;
    src = A(sourceRange,:);
    
    for j=1:length(M)
        f = zeros(M{j}.numVertices,1);
        pt = M{j}.vertices(samples{j},:);
        
        idx = 1+(j-1)*k;
        tgt = A(idx:(idx+k-1),:);
        
        map = src*tgt';
        %map = map';
        map2 = otherdata.noConsistencyMaps(sourceRange,idx:(idx+k-1));
        map3 = G(sourceRange,idx:(idx+k-1));
        
        fig = subplot(3,length(M)+1,1+j);
        showDescriptor(M{j},f,[],[],[],fig); colorbar off; hold on;
        scatter3(pt(:,1),pt(:,2),pt(:,3),50,map(i,:),'filled');
        
        fig = subplot(3,length(M)+1,2+j+length(M));
        showDescriptor(M{j},f,[],[],[],fig); colorbar off; hold on;
        scatter3(pt(:,1),pt(:,2),pt(:,3),50,map3(i,:),'filled');
        
        fig = subplot(3,length(M)+1,3+j+2*length(M));
        showDescriptor(M{j},f,[],[],[],fig); colorbar off; hold on;
        scatter3(pt(:,1),pt(:,2),pt(:,3),50,map2(i,:),'filled');
    end
end

%% Merge the meshes

XX = [];
TT = [];
kk = 0;
shift = [1.7 0 0];
shift2 = [0 1.25 0];
fn = [];
fn2 = [];
pts = [];
whichVerts = [];

for i=1:length(M)
    curXX = M{i}.vertices;
    curXX = bsxfun(@minus,curXX,mean(curXX));
    axes = pca(curXX);
    if det(axes) < 0
        axes = axes*-1;
    end
    curXX = curXX*axes';
    
    if i==16
        curXX = [curXX(:,2) curXX(:,1) curXX(:,3)];
        curXX = curXX*-1;
    end
    
    curXX = bsxfun(@plus,curXX,shift*mod(i-1,5) + shift2*floor((i-1)/5));
    
    XX = [XX; curXX];
    TT = [TT; M{i}.triangles+kk];
    
    idx = 1+(i-1)*k;
    tgt = A(idx:(idx+k-1),:);
    map = src*tgt';
    
    map2 = otherdata.noConsistencyMaps(sourceRange,idx:(idx+k-1));
    
    fn2 = [fn2;map2(2,:)'/max(map2(2,:))];
    fn = [fn;map(2,:)'/max(map(2,:))];
    pts = [pts; curXX(samples{i},:)];
    
    whichVerts = [whichVerts; kk+samples{i}'];
    
    kk = kk+M{i}.numVertices;
end

% showDescriptor(MM,MM.vertices(:,1)*0);
% hold on
% scatter3(pts(:,1),pts(:,2),pts(:,3),10,fn,'filled');

[Xout,Tout,index] = addSpheres(XX,TT,whichVerts,ones(length(whichVerts),1)*.075);

MM.vertices = Xout;
MM.triangles = Tout;
MM.numVertices = size(Xout,1);
MM.numTriangles = size(Tout,1);

ff = zeros(MM.numVertices,1);
ff2 = zeros(MM.numVertices,1);
for i=1:length(fn)
    ff(index==i) = fn(i);
    ff2(index==i) = fn2(i);
end
showDescriptor(MM,ff);
showDescriptor(MM,ff2);

%%

% writeTexturedObj('consistentTarget.obj', MM, ff);
% writeTexturedObj('inconsistentTarget.obj', MM, ff2);

generateSoftMapMeshes(M{1},MM,'source_inconsistent.obj','target_inconsistent.obj',...
    [],[],samples{1}(2),ff2,.1);

generateSoftMapMeshes(M{1},MM,'source_inconsistent.obj','target_consistent.obj',...
    [],[],samples{1}(2),ff,.1);