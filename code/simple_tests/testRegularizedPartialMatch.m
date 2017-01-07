%% Arbitrary parameterized shape without any symmetries
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
options.regularizer = .0005;
options.maxIter = 100;
options.mu = 100;
options.maxInnerIter = 1000;

% g = @(x) .5*sum(x.^2)/n*100; % L2
% gradG = @(x) x/n*100;

% g = @(x) 0%sum(.5*(x-1).^2.*double(x>1))*.01; % Hinge
% gradG = @(x) ones(size(x));%double(x>1).*(x-1)*.01;

uniformMarginal = 1;
coveringFrac = .5; % expect this much out of the marginal
regularizerStrength = 1;

a = uniformMarginal/coveringFrac;
g = @(x) sum((x.^2).*(x-a).^2)/((.5*a^4)*n) * regularizerStrength;
gradG = @(x) 2*x.*(a^2-3*a*x+2*x.^2)/((.5*a^4)*n) * regularizerStrength;

close all;
gamma = regularizedPartialGW(D,D0,g,gradG,options);

%% Visualize 2d map by connecting each point to its most likely match

nPlots = 7;
step = ceil(n0/nPlots);
pos = 1;
for i=1:step:n0
    subplot(2,nPlots,pos);
    scatter(X0(:,1),X0(:,2),[],sum(gamma,1),'filled');
    axis equal; axis off;  title('Row sums');
    if i+step > n0
        colorbar; caxis([0 max(sum(gamma,2))]); 
    end
    hold on;
    plot(X0(i,1),X0(i,2),'r.','markersize',50);
    
    subplot(2,nPlots,pos+nPlots);
    scatter(X(:,1),X(:,2),[],gamma(:,i),'filled');
    axis equal; axis off; title(sprintf('Row %d',i));
    
    pos = pos+1;
end
