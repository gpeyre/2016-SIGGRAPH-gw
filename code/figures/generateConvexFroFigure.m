%% Load data

close all
class = [2 5];
which = [1 6];
dir = '../data/2d_shapes/plane_data';

n = 75;

minSize = inf;
for i=1:2
    load(sprintf('%s/Class%d_Sample%d.mat',dir,class(i),which(i)));
    
    div = ceil(size(x,1)/n);
    x = x(1:div:end,:);
    
    minSize = min(minSize,size(x,1));
    X{i} = x;
end

for i=1:2
    if size(X{i},1) > minSize
        k = size(X{i},1)-minSize;
        toRemove = randperm(size(X{i},1),k);
        X{i}(toRemove,:) = [];
    end
    
    D{i} = distmat(X{i}');
    figure;
    imagesc(D{i});
    axis equal; axis off;
end

X{1} = X{1}*-1;

for i=1:2
    figure; plot(X{i}(:,1),X{i}(:,2),'.'); axis equal; axis off;
end

%% Kimmel solution

cvx_begin
    cvx_solver sdpt3
    cvx_precision high
    variable map(size(X{1},1),size(X{2},1))
    
    minimize norm(D{1}*map - map*D{2},'fro')
    subject to 
        map >= 0
        sum(map) == 1
        sum(map') == 1
cvx_end

%% Compute Gromov-Wasserstein

fprintf('Optimizing regularized Gromov-Wasserstein...\n');

options = [];

options.display = 1;
options.regularizer = .0004;
options.plotObjective = 1;
options.maxIter = 50;

gamma = gromovWassersteinDistance(D{1},D{2},options);

%% Render map

points = 1:13:size(D{1},1);
nPoints = length(points);

for i=1:nPoints
    p = points(i);
    figure;
    
    figure('color',[1 1 1]);%fig = subplot(1,3,1);
    plot(X{1}(:,1),X{1}(:,2),'k.','markersize',30); hold on;
    plot(X{1}(p,1),X{1}(p,2),'r.','markersize',40);
    axis equal; axis off;
    saveas(gcf,sprintf('source%d.pdf',i));
    close all
    
    figure('color',[1 1 1]);%fig = subplot(1,3,2);
    scatter(X{2}(:,1),X{2}(:,2),60,map(p,:)','filled');
    axis equal; axis off;
    saveas(gcf,sprintf('fro%d.pdf',i));
    close all
    
    figure('color',[1 1 1]);%fig = subplot(1,3,3);
    scatter(X{2}(:,1),X{2}(:,2),60,gamma(p,:)','filled');
    axis equal; axis off;
    saveas(gcf,sprintf('gw%d.pdf',i));
    close all
end
