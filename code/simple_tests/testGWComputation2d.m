%% Load mesh, put into convenient data structure

addpath('../');

% Arbitrary parameterized shape without any symmetries
f0 = @(x) [ cos(x)+.1*sin(3*x) ; .5*sin(2*x).^3+.1*cos(5*x) ];

% Sample the shape
n0 = 300;
t0 = linspace(0,2*pi,n0+1);
t0(end) = [];
X0 = f0(t0)';

% Sample a rotated version of the same shape
n = 400;
%R = rotx(10);
%R = R(2:3,2:3);
theta = 10;
R = [cos(theta), sin(theta); -sin(theta), cos(theta)];
f = @(x) R*f0(x);
t = linspace(0,2*pi,n+1);
t(end) = [];
X = f(t)';

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
options.regularizer = .001;
options.plotObjective = 1;
options.maxIter = 500;


% Partial-to-full test
options.partialSource = 0;
options.partialTarget = 0;

gamma = gromovWassersteinDistance(D0,D,options);

%% Visualize 2d map by connecting each point to its most likely match

plot(X0(:,1),X0(:,2),'b'); axis equal; hold on;
plot(X(:,1),X(:,2),'r'); axis equal;
axis off;

for i=1:n0
    [~,p] = max(gamma(i,:));
    plot([X0(i,1) X(p,1)],[X0(i,2) X(p,2)],'k');
end

figure;
plot(1:n0,sum(gamma,2),1:n,sum(gamma,1));
legend('Row sums','Column sums');
