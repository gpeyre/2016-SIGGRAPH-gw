function D = readMtx(filename)
fid = fopen(filename,'r');
if( fid==-1 )
    error('Cannot open the file.');
    return;
end
sizeD = fread(fid, [1, 2], 'int');
D = fread(fid, sizeD, 'double');
fclose(fid);

