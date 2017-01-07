%% Load pattern and construct a graph

pattern = rgb2gray(imread('../data/yo.png'));

pattern = imresize(pattern,[14,21]);
pattern = double(pattern>0);

% Pad with zeros to simplify code a bit
pattern = padarray(pattern,1);
pattern = padarray(pattern',1)';

% Label each nonzero pixel
pos = 1;
for r=1:size(pattern,1)
    for c=1:size(pattern,2)
        if pattern(r,c) ~= 0
            pattern(r,c) = pos;
            pos = pos+1;
        end
    end
end

n = pos;

adj = zeros(n,n);
for r=1:size(pattern,1)
    for c=1:size(pattern,2)
        if pattern(r,c) ~= 0
            i = pattern(r,c);
            if pattern(r+1,c) ~= 0, adj(i,pattern(r+1,c)) = 1; end
            if pattern(r-1,c) ~= 0, adj(i,pattern(r-1,c)) = 1; end
            if pattern(r,c+1) ~= 0, adj(i,pattern(r,c+1)) = 1; end
            if pattern(r,c-1) ~= 0, adj(i,pattern(r,c-1)) = 1; end
        end
    end
end

D = graphallshortestpaths(sparse(adj));

W = double(~isinf(D));
D(isinf(D)) = 20; % just to avoid unhappy math

%% Load images

load ../data/images.mat

images = images(1:n);

sz = size(images{1},1)*size(images{1},2);
avLab = zeros(length(images),3);
for i=1:length(images)
    avLab(i,:) = squeeze(sum(sum(rgb2lab(images{i}/255))))'/sz;
end

D0 = distmat(avLab');
W0 = ones(size(D0));

%% Compute Gromov-Wasserstein

fprintf('Optimizing regularized Gromov-Wasserstein...\n');

gamma = weightedGW(D0,W0,D,W,.003,300);%.00002

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

for c=1:size(pattern,2)
    curCol = [];
    for r=1:size(pattern,1)
        if pattern(r,c) ~= 0
            pos = pattern(r,c);
            [~,match] = max(gammaRounded(:,pos)); 
            curCol = [curCol ; images{match}];
        else
            curCol = [curCol ; 255*ones(size(images{1}))];
        end
    end
    im = [im curCol];
end

figure;imshow(im);

%% Save image

imwrite(im,'weightedImageLayout.png');