%% Load data

clear

close all
class = 3; % 2 for hearts
dir = '../data/2d_shapes/bicego_data';
load(sprintf('%s/Class%d_Sample%d.mat',dir,class,1));

X0 = x;
n0 = size(X0,1);
D0 = distmat(X0');

k = 19;
for i=1:k
    load(sprintf('%s/Class%d_Sample%d.mat',dir,class,i+1));
    X{i} = x;
    n{i} = size(x,1);
    D{i} = distmat(x');
    figure; plot(x(:,1),x(:,2),'.'); axis equal; axis off;
end

%% Compute Gromov-Wasserstein

close all
fprintf('Optimizing regularized Gromov-Wasserstein...\n');

options = [];

options.display = 1;
options.regularizer = .0006;%00075
options.plotObjective = 1;
options.maxIter = 50;

gamma = cell(k,1);
for i=1:k
    gamma{i} = gromovWassersteinDistance(D0,D{i},options);
end

%% Round to a point-to-point map

gammaRounded = cell(k,1);
for i=1:k
    cvx_begin
        cvx_solver mosek
        variable g(size(gamma{i}))
        maximize sum(g(:).*gamma{i}(:))
        subject to
            sum(g) == sum(gamma{i})
            sum(g') == sum(gamma{i}')
            g >= 0
    cvx_end
    gammaRounded{i} = g;
end

%% Assign colors to source

colors = zeros(n0,3);

XX = bsxfun(@minus,X0,mean(X0));
theta = atan2(XX(:,1),XX(:,2));
c = (theta+pi)/(2*pi);

for i=1:n0
    colors(i,:) = hsv2rgb([c(i) 1 1]);
end

scatter(X0(:,1),X0(:,2),[],colors,'filled');

%% Colored results

figure('color',[1 1 1]);

% fig = subplot(1,length(gamma)+1,1);
scatter(X0(:,1),X0(:,2),20,colors,'filled');
axis equal; axis off;

saveas(gcf,'source_glass.png');
close

for j=2:(length(gamma)+1)
    curColors = zeros(n{j-1},3);
    for k=1:n{j-1}
        [~,p] = max(gammaRounded{j-1}(:,k));
        curColors(k,:) = colors(p,:);
    end
    
%     fig = subplot(1,length(gamma)+1,j);
    figure('color',[1 1 1])
    scatter(X{j-1}(:,1),X{j-1}(:,2),20,curColors,'filled');
    axis equal; axis off;
    
    saveas(gcf,sprintf('target_glass%d.png',j-1));
    close
end