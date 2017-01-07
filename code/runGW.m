%% Load mesh, put into convenient data structure
function [gamma, D, D0] = runGW(dists, dists0, outfile, gwReg, gwIter, vis)
setupPaths
D = readMtx(dists);
D0 = readMtx(dists0);
nvtx0 = size(D0, 1);
nvtx = size(D, 1);
if ~exist('vis','var') || isempty(vis)
  vis=0;
end

%% Compute Gromov-Wasserstein
fprintf('Optimizing regularized Gromov-Wasserstein...\n');

options = [];
options.mu0 = ones(nvtx0,1) / nvtx0;
options.mu = ones(nvtx,1) / nvtx;
options.display = vis;
options.regularizer = gwReg; % .0007 matches single tentacle, .005 matches some, .008 middle
options.plotObjective = vis;
options.eta = 1;
options.maxIter = gwIter;

[gamma,objectives] = gromovWassersteinDistance(D0,D,options);
% size(gamma)

%% Write GW distances
Lambda = 0.5*D0.^2 * diag(options.mu0)*gamma*options.mu*ones(1,nvtx);
Lambda = Lambda - D0 * diag(options.mu0)*gamma*diag(options.mu)*D;
Lambda = Lambda + 0.5*ones(nvtx0,1)*options.mu0'*gamma*diag(options.mu)*D.^2;
gwDist = sum(sum(diag(options.mu0)*Lambda.*gamma*diag(options.mu)));
writeMtx(gamma', outfile);

fid = fopen([outfile '.txt'],'w');
fwrite(fid, ['dist=' num2str(gwDist)]);
fclose(fid);

close all