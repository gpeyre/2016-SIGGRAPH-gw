%% Generate points in the unit square

clear
density = 5;

X0 = poissonSamples([100 100],density,20)/100;%rand(n,2);
X0((1-X0(:,2))<X0(:,1),:) = []; % sample from a triangle

X = poissonSamples([100 100],density,20)/100;
X((1-X(:,2))<X(:,1),:) = []; % sample from a triangle

n = min(size(X,1),size(X0,1));
n0 = n;

X = X(1:n,:);
X0 = X0(1:n,:);

close all
% plot(X0(:,1),X0(:,2),'r.'); axis equal; axis off; hold on;
% plot(X(:,1),X(:,2),'b.');

%% Compute full distance matrix

fprintf('Computing pairwise distance...\n');

D0 = distmat(X0');
D = distmat(X');

%% Compute Gromov-Wasserstein

alphas = [.01 .001];

gwGamma = cell(length(alphas),1);
gwOptimTimes = cell(length(alphas),1);
gwOptimVals = cell(length(alphas),1);

bfgsGamma = cell(length(alphas),1);
bfgsOptimTimes = cell(length(alphas),1);
bfgsOptimVals = cell(length(alphas),1);

for i=1:length(alphas)
    options.regularizer = alphas(i);
    [gwGamma{i},gwOptimTimes{i},gwOptimVals{i}] = minimalGW(D0,D,options.regularizer);
    [bfgsGamma{i},bfgsOptimTimes{i},bfgsOptimVals{i}] = bfgsGW(D0,D,options.regularizer);
end

%% Plot results

for i=1:length(alphas)
    figure
    plot(gwOptimTimes{i},gwOptimVals{i},'-',bfgsOptimTimes{i},bfgsOptimVals{i},'-');
    legend('GW','BFGS');
    title(sprintf('alpha=%g',alphas(i)));
end

%% Write out files for TikZ plots

for i=1:length(alphas)
    fid = fopen(sprintf('gw%dalpha%g.txt',n,alphas(i)),'w');
    for j=1:length(gwOptimTimes{i})
        fprintf(fid,'%g,%g\n',gwOptimTimes{i}(j),gwOptimVals{i}(j));
    end
    fclose(fid);
    
    fid = fopen(sprintf('bfgs%dalpha%g.txt',n,alphas(i)),'w');
    for j=1:length(bfgsOptimTimes{i})
        fprintf(fid,'%g,%g\n',bfgsOptimTimes{i}(j),bfgsOptimVals{i}(j));
    end
    fclose(fid);
end