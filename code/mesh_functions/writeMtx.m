function  writMtx(D, filename)
fid = fopen(filename,'w');
if( fid==-1 )
    error('Cannot open the file.');
    return;
end
fwrite(fid, size(D), 'int');
fwrite(fid, D, 'double');
fclose(fid);