function nBytes = fileSize(filename)
%FILESIZE returns the size of the file at the specified path
%   
% N = FILESIZE(F) returns the file size in bytes of the file specified
% by filename F. F may also be a file ID.

if ischar(filename) || isa(filename, 'string')
    fid = fopen(filename, 'r');
    isFile = false;
else
    fid = filename;
    isFile = true;
    pos = ftell(fid);
end

fseek(fid, 0, 'eof');
nBytes = ftell(fid);

if isFile
    fseek(fid, pos, 'bof');
else
    fclose(fid);
end

end

