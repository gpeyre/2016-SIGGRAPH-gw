D = dlmread('shrec-pairwise.dat');
D = .5*(D+D');

figure
p = mdscale(D,2);
plot(p(1:20,1),p(1:20,2),'r.',...
     p(21:40,1),p(21:40,2),'g.',...
     p(41:60,1),p(41:60,2),'b.',...
     p(61:80,1),p(61:80,2),'k.','markersize',25);
 
legend('Humans','Teddies','Armadillo','Four-legged');
title('MDS embedding 2D');
axis equal;

figure
p = mdscale(D,3);
plot3(p(1:20,1),p(1:20,2),p(1:20,3),'r.',...
     p(21:40,1),p(21:40,2),p(21:40,3),'g.',...
     p(41:60,1),p(41:60,2),p(41:60,3),'b.',...
     p(61:80,1),p(61:80,2),p(61:80,3),'k.','markersize',25);
 
legend('Humans','Teddies','Armadillo','Four-legged');
title('MDS embedding 3D');
axis equal; cameratoolbar;

dlmwrite('humans.dat',p(1:20,:));
dlmwrite('teddies.dat',p(21:40,:));
dlmwrite('armadillo.dat',p(41:60,:));
dlmwrite('fourlegged.dat',p(61:80,:));