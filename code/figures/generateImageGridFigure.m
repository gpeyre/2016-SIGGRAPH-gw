%% Load images

clear
load ../data/images.mat

w = 14;
n = w*w;
images = images(1:n);

sz = size(images{1},1)*size(images{1},2);
avLab = zeros(length(images),3);
for i=1:length(images)
    avLab(i,:) = squeeze(sum(sum(rgb2lab(images{i}/255))))'/sz;
end
% avLab = avLab(:,2:3); % get rid of L

%% Generate grid

r = repmat((1:w)',1,w);
c = r';

gridCoord = [r(:) c(:)];

%% Compute distances

D0 = distmat(avLab');
D = distmat(gridCoord');

%% Try to align the two

%% Compute Gromov-Wasserstein

fprintf('Optimizing regularized Gromov-Wasserstein...\n');

options = [];

options.display = 1;
options.regularizer = .00002;
options.plotObjective = 1;
options.maxIter = 300;

gamma = gromovWassersteinDistance(D0,D,options);

%% Round using a linear assignment problem

cvx_begin
    variable gammaRounded(n,n)
    
    cvx_solver mosek
    cvx_precision best
    
    maximize sum(sum(gamma.*gammaRounded))
    subject to
        sum(gamma) == sum(gammaRounded)
        sum(gamma') == sum(gammaRounded')
        gammaRounded >= 0
cvx_end

%% Reconstruct image

pos = 1;

im = [];
for col = 1:w
    curCol = [];
    for row = 1:w
        [~,match] = max(gammaRounded(:,pos)); % might need to flip
        curCol = [curCol ; images{match}];
        pos = pos+1;
    end
    im = [im curCol];
end

figure;imshow(im);

%% Save image

imwrite(im,'imageLayout.png');