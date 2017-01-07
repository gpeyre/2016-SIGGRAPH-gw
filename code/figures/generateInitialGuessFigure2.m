%% Generate points in the unit square

X0 = poissonSamples([100 100],5,20)/100;%rand(n,2);
X0((1-X0(:,2))<X0(:,1),:) = []; % sample from a triangle
n0 = size(X0,1);

X = poissonSamples([100 100],5,20)/100;
X((1-X(:,2))<X(:,1),:) = []; % sample from a triangle
n = size(X,1);

close all
plot(X0(:,1),X0(:,2),'r.'); axis equal; axis off; hold on;
plot(X(:,1),X(:,2),'b.');

%% Compute full distance matrix

fprintf('Computing pairwise distance...\n');

D0 = distmat(X0');
D = distmat(X');

%% Compute Gromov-Wasserstein

fprintf('Optimizing regularized Gromov-Wasserstein...\n');

options = [];

options.mu0 = ones(n0,1)/n0;
options.mu = ones(n,1)/n;
options.display = 0;
options.plotObjective = 0;
options.eta = 1;
options.maxIter = 100;
% options.GWTol = 1e-8; % force running all the iterations

o = ones(n,1);
o0 = ones(n0,1);
dmu0 = spdiags(options.mu0,0,n0,n0);
dmu = spdiags(options.mu,0,n,n);

nTests = 500;

alphas = zeros(nTests,1);
values = zeros(nTests,1);
minAlpha = .0003;%.0003;
maxAlpha = .0015;%.0006;%.1;

for i=1:nTests
    fprintf('Test %d of %d...\n',i,nTests);
    t = rand();
    alphas(i) = minAlpha + exp((1-t)*log(minAlpha) + t*log(maxAlpha));
    options.regularizer = alphas(i);
    
    options.initialGuess = rand(n0,n)+1e-5;
    [gamma,objectives] = gromovWassersteinDistance(D0,D,options);
    
    L1 = .5*(D0.^2)*dmu0*gamma*options.mu*o';
    L2 = -D0*dmu0*gamma*dmu*D;
    L3 = .5*o0*options.mu0'*gamma*dmu*(D.^2);
    Lambda = L1+L2+L3;
    values(i) = sum(sum(bsxfun(@times,bsxfun(@times,Lambda.*gamma,options.mu0),options.mu')));
    
    clf;
    semilogx(alphas,values,'.');
    title('Experiments');
    xlabel('alpha');
    ylabel('objective');
    xlim([min(alphas) max(alphas)]);
    ylim([min(values) max(values)]);
    drawnow;
end

save initialGuessTest2.mat

%% Write out in a TikZ-friendly way

fid = fopen('initialGuess2.txt','w');
for i=1:nTests
    fprintf(fid,'%g,%g\n',alphas(i),values(i));
end
fclose(fid);

%% Generate two maps with same alpha, different initial guesses

p = 13;%randi(n0,[1 1]);

options = [];
options.mu0 = ones(n0,1)/n0;
options.mu = ones(n,1)/n;
options.display = 0;
options.plotObjective = 0;
options.eta = 1;
options.maxIter = 100;
options.regularizer = .001;

options.initialGuess = rand(n0,n)+1e-5;
gamma1 = gromovWassersteinDistance(D0,D,options);
options.initialGuess = rand(n0,n)+1e-5;
gamma2 = gromovWassersteinDistance(D0,D,options);

%%

figure('color',[1 1 1]);
plot(X0(:,1),X0(:,2),'k.','markersize',30); hold on;
plot(X0(p,1),X0(p,2),'r.','markersize',50);
axis equal; axis off;
    
figure('color',[1 1 1]);
scatter(X(:,1),X(:,2),60,gamma1(p,:),'filled');
axis equal; axis off;

figure('color',[1 1 1]);
scatter(X(:,1),X(:,2),60,gamma2(p,:),'filled');
axis equal; axis off;

%%

options.initialGuess = rand(n0,n)+1e-5;
options.regularizer = .02;
gamma3 = gromovWassersteinDistance(D0,D,options);

figure('color',[1 1 1]);
scatter(X(:,1),X(:,2),60,gamma3(p,:),'filled');
axis equal; axis off;