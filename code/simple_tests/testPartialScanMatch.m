%% Load data

clear;close all;

X0 = dlmread('../data/simulated_scans/mesh000.off/mesh000.xyz');
X = dlmread('../data/simulated_scans/mesh000.off/View2.xyz');

X0 = X0(randi(size(X0,1),[500 1]),:);
[~,order] = sort(X0(:,1));
X0 = X0(order,:);

% X = X0(1:500,:);
X = X(randi(size(X,1),[250 1]),:);
[~,order] = sort(X(:,1));
X = X(order,:);

n0 = size(X0,1);
n = size(X,1);

plot3(X0(:,1),X0(:,2),X0(:,3),'rx');
hold on;
plot3(X(:,1),X(:,2),X(:,3),'b.');
axis equal; axis off; cameratoolbar;

%% Ground truth match

[~,fixedSource] = min(X(:,1));
fixedTarget = knnsearch(X0,X(fixedSource,:));

p = [X0(fixedTarget,:) ; X(fixedSource,:)];
plot3(p(:,1),p(:,2),p(:,3),'g-','linewidth',3);

fixedTargetDist = zeros(n0,1)+1e-5;
fixedTargetDist(fixedTarget) = n0;

%% Compute distances -- for now just extrinsic

D0 = distmat(X0',X0');
D = distmat(X',X');

%% Compute Gromov-Wasserstein

fprintf('Optimizing regularized Gromov-Wasserstein...\n');

options = [];
options.regularizer = .0003;
options.maxIter = 1000;
% options.mu = 1000;
options.eta = .1;
options.maxInnerIter = 500;
options.GWTol = 1e-5;

options.partialTarget = 1;

options.display = 1;
options.plotObjective = 1;

% options.initialGuess = zeros(n0,n);
% options.initialGuess(1:n,1:n) = eye(n)*n0;
% options.initialGuess = options.initialGuess + 1e-8;
% options.initialGuess = options.initialGuess';

options.normalizeDistances = 1; % Important!!

options.groundTruthColIndices = fixedSource;
options.groundTruthCols = fixedTargetDist;

gamma = gromovWassersteinDistance(D0,D,options)';

% uniformMarginal = 1;
% coveringFrac = 1/8; 
% regularizerStrength = 10;
% 
% a = uniformMarginal/coveringFrac;
% g = @(x) sum((x.^2).*(x-a).^2)/((.5*a^4)*n) * regularizerStrength;
% gradG = @(x) 2*x.*(a^2-3*a*x+2*x.^2)/((.5*a^4)*n) * regularizerStrength;
% 
% options.fixedSource = fixedSource;
% options.fixedTarget = fixedTargetDist;
% 
% close all;
% gamma = regularizedPartialGW(D,D0,g,gradG,options);

%% Visualize output

nPoints = 10;
close all

for i=1:nPoints
    figure;
    
    p = randi([1 n]);
    
    subplot(1,3,1);
    scatter3(X0(:,1),X0(:,2),X0(:,3),[],sum(gamma,1),'.'); axis equal; axis off;
    title('Marginals');
    
    subplot(1,3,2);
    hold on;
    scatter3(X(:,1),X(:,2),X(:,3),'.'); axis equal; axis off;
    plot3(X(p,1),X(p,2),X(p,3),'r.','markersize',50); axis equal; axis off;
    cameratoolbar('SetCoordSys','x');
    title('Source');
    
    subplot(1,3,3);
    scatter3(X0(:,1),X0(:,2),X0(:,3),[],gamma(p,:)','.'); axis equal; axis off;
    axis equal;
    title('Target');
    axis off;
    
    cameratoolbar
end