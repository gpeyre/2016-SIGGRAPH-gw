pfx = '/Volumes/VovaAdobe5T/DATA/2016/GWTestDir/';
roiFile=  [pfx '/PerModel/mesh004.off/highlight_bent_knee.ids'];
gammaFile= [pfx '/Pairwise/mesh004.off_mesh062.off/EvenFast_0.1_manifold.0.005.gw'];
distsFile = [pfx '/PerModel/mesh004.off/EvenFast_0.1_manifold.dists'];
dists0File = [pfx '/PerModel/mesh062.off/EvenFast_0.1_manifold.dists'];

D = readMtx(dists);
D0 = readMtx(dists0);
gamma = readMtx(gammaFile);

mu0 = ones(length(D0),1) / length(D0);
mu = ones(length(D),1) / length(D);

Lambda = 0.5*D0.^2 * diag(mu0)*gamma*mu*ones(1,nvtx);
Lambda = Lambda - D0 * diag(mu0)*gamma*diag(mu)*D;
Lambda = Lambda + 0.5*ones(nvtx0,1)*mu0'*gamma*diag(mu)*D.^2;
gwDist = sum(sum(diag(mu0)*Lambda.*gamma*diag(mu)));

gwDist
