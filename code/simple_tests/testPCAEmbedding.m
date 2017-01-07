models = {};
labels = {};
% 48 -- horse gallop
% 64 -- sphere
% 71 -- scape
% SHREC
%for i=1:64
%for i=1:71
%    if i==51
%        continue
%    end
%     if i<10
         %models{length(models)+1} = ['horse-gallop-0' num2str(i) '.obj'];
%         models{length(models)+1} = ['mesh00' num2str(i) '.off'];
%     else
         %models{length(models)+1} = ['horse-gallop-' num2str(i) '.obj'];
%         models{length(models)+1} = ['mesh0' num2str(i) '.off'];
%     end
%    models{length(models)+1} = ['deform_sph_' num2str(i) '.off'];
%    labels{length(models)} = num2str(i);
%end

for i=1:20
    models{i} = [num2str(i) '.off'];
    models{20*i} = [num2str(160+i) '.off'];
    models{40*i} = [num2str(280+i) '.off'];
    models{60*i} = [num2str(380+i) '.off'];
end

%reference = 'horse-gallop-reference.obj';
%reference = 'sphereuniflarge.off';
reference = 'mesh000.off';

dir = '/Users/vokim/Research/DBShapes3D/DBShapesCode/bin/AnalyzedDBs/';
dir = [dir 'GWTestDir/Pairwise/'];
mapFile = 'EvenFast_0.1_manifold.0.005.gw';

Dall = [];
for i=1:length(models)
    mapfile = [dir models{i} '_' reference '/' mapFile];
    gamma = readMtx(mapfile);
    D = gamma'*gamma;
    ds = size(D);
    Dvec = reshape(D, [1,ds(1)*ds(2)]);
    if i==1
        Dall = Dvec;
    else
        Dall = [Dall;Dvec];
    end 
end
[coeff,score,latent,tsquared,explained,mu] = pca(Dall');
x = coeff(:, 1);
y = coeff(:, 2);
x = (x-min(x))/(max(x)-min(x));
y = (y-min(y))/(max(y)-min(y));
labelpoints(x, y, labels, 'adjust_axes', 1);