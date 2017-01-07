function [pts, vtx] = readPts(filename, T)

fid = fopen(filename,'r');
if( fid==-1 )
    error('Cannot open the file.');
    return;
end

%str = fgets(fid);   % -1 if eof
[pts, cnt] = fscanf(fid,'%lf %lf %lf %lf %lf %lf %lf\n', [7, inf]);
fclose(fid);

cnt = cnt/7;

pts = pts';
vtx = zeros(1, cnt);
for i=1:cnt
    if pts(i,2)>=pts(i,3) && pts(i,2)>=pts(i,4)
        vtx(i) = T(pts(i,1),1);
    elseif pts(i,3)>=pts(i,4)
        vtx(i) = T(pts(i,1),2);
    else
        vtx(i) = T(pts(i,1),3);
    end
end
