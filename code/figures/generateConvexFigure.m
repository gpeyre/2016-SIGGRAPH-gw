clear

n0 = 7;
n = n0;

theta0 = linspace(0,2*pi,1+n0)';
theta0 = theta0(1:n);
X0{1} = rand(n0,2);
X0{1} = X0{1} / max(X0{1}(:));

X0{2} = [sin(theta0) 2*cos(theta0)];
X0{2}(:,2) = X0{2}(:,2)*-1;
X0{2} = X0{2} / max(X0{2}(:));

X = X0;

for i=1:2
    X{i}(:,1) = X{i}(:,1) + 1.2; % shift for visualization
    D0{i} = distmat(X0{i}');
    D{i} = D0{i}; % same point set
end

% plot(X0(:,1),X0(:,2),'bx');
% hold on; axis equal;
% plot(X(:,1),X(:,2),'rx');

idx = reshape(1:(n0*n),n0,n);

%% CVX optimization

maps = [];
outerprods = [];
for iter = 1:length(X0)
    cvx_begin
        variable map(n0,n)
        variable outerprod(n0*n,n0*n)
        cvx_solver mosek

        objective = [];
        for i=1:n0
            i
        for j=1:n0
        for k=1:n0
        for l=1:n0
            curTerm = (D0{iter}(i,j)-D{iter}(k,l))^2 * outerprod(idx(i,k),idx(j,l));
            if isempty(objective)
                objective = curTerm;
            else
                objective = objective+curTerm;
            end
        end
        end
        end
        end

        minimize objective
        subject to
            sum(map,1) == 1
            sum(map,2) == 1
            map >= 0
            [1 map(:)' ; map(:) outerprod] == semidefinite(n0*n+1)
            trace(outerprod) == n0
            sum(outerprod(:)) == n0*n0
            outerprod >= 0

            for q=1:n0
                q
            for r=1:n0
            for s=1:n0
            for t=1:n0
                idx1 = idx(q,r);
                idx2 = idx(s,t);
                outerprod(idx1,idx2) <= map(q,r)
                outerprod(idx1,idx2) <= map(s,t)
                if q == s && r ~= t || r == t && q ~= s
                    outerprod(idx1,idx2) == 0
                end
            end
            end
            end
            end
    cvx_end
    
    maps{iter} = map;
    outerprods{iter} = outerprod;
end

%% Compute Gromov-Wasserstein

close all
fprintf('Optimizing regularized Gromov-Wasserstein...\n');

options = [];

options.display = 1;
options.regularizer = .001;
options.plotObjective = 1;
options.maxIter = 20;

for i=1:length(X0)
    maps{end+1} = gromovWassersteinDistance(D0{i},D{i},options);
    return
end

%% Show matches

for i=1:length(maps)
    figure('color',[1 1 1]); hold on;
    
    curX0 = X0{mod(i-1,length(X0))+1};
    curX = X{mod(i-1,length(X))+1}+.1;
    
    mtx = maps{i};
    mtx = bsxfun(@rdivide,mtx,sum(mtx,2));
    mtx = 1-mtx;
    
    [~,order] = sort(mtx(:));
    order = order(end:-1:1);
    for l=1:(n0*n0)
        j = mod(order(l)-1,n0)+1;
        k = 1+floor((order(l)-1)/n0);
        plot([curX0(j,1) curX(k,1)],[curX0(j,2) curX(k,2)],'color',[1 1 1]*mtx(j,k),'linewidth',2);
    end
    
    plot(curX0(:,1),curX0(:,2),'k.','markersize',20);
    plot(curX(:,1),curX(:,2),'k.','markersize',20);
    
    minY = min(curX0(:,2));
    maxY = max(curX(:,2));
    
    l = max(curX0(:,1));
    r = min(curX(:,1));
    c = (l+r)/2;
    plot([c c],[minY-.1 maxY+.1],'k:');
    
    axis equal;
    axis off;
end