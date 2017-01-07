%% Load mesh, put into convenient data structure

addpath('../');

% Arbitrary parameterized shape without any symmetries
f0 = @(x) [ cos(x)+.2*sin(5*x).^5 ; sin(x)+.2*cos(4*x) ];

% Sample the shape
n0 = 400;
t0 = linspace(0,2*pi,n0+1);
t0(end) = [];
X0 = f0(t0)';

n=200; % Take a segment for partial matching
X = X0(1:n,:);
X(:,2) = X(:,2) + 1;

%% Compute full distance matrix

D0 = zeros(n0,n0);
for i=1:n0
    for j=1:n0
        D0(i,j) = norm(X0(i,:)-X0(j,:));
    end
end

D = zeros(n,n);
for i=1:n
    for j=1:n
        D(i,j) = norm(X(i,:)-X(j,:));
    end
end

%% Compute Gromov-Wasserstein

fprintf('Optimizing regularized Gromov-Wasserstein...\n');

options = [];

options.display = 1;
options.regularizer = .001; % 0003, 005 are bad, 0005 is good
options.plotObjective = 1;
options.maxIter = 500;
options.GWTol = 1e-9;

% options.initialGuess = eye(n0); % Uncomment me and it'll reach the right optimum!
% options.initialGuess(:,(n+1):end) = [];
options.initialRegularizer = .1;
options.regularizerChangeRate = .05;

% Partial-to-full test
options.partialSource = 0;
options.partialTarget = 1;

gamma = gromovWassersteinDistance(D0,D,options);
close all;

%% Visualize 2d map by connecting each point to its most likely match

nPlots = 6;
step = ceil(n0/nPlots);
pos = 1;
for i=1:step:n0
    subplot(2,nPlots,pos);
    scatter(X0(:,1),X0(:,2),[],sum(gamma,2),'filled');
    axis equal; axis off;  title('Row sums');
    if i+step > n0
        colorbar; caxis([0 max(sum(gamma,2))]); 
    end
    hold on;
    plot(X0(i,1),X0(i,2),'r.','markersize',50);
    
    subplot(2,nPlots,pos+nPlots);
    scatter(X(:,1),X(:,2),[],gamma(i,:),'filled');
    axis equal; axis off; title(sprintf('Row %d',i));
    
    pos = pos+1;
end
