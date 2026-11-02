function v = fgetln(fin,n)

% function v = fgetln(fin,n)
% reads at most n lines from a file
%
% input parameters
%    fin        : file name or file identifier
%    n          : number of lines    (optional; default: Inf)
%
% output parameters
%    v          : cell array of strings, one per line, not including the newline character

%% A. Parameter check
if ~exist('n','var'), n = Inf; end
if ischar(fin), fid = fopen(fin,'r'); else fid = fin; end

%% B. Determine the number of rows to be read
i = 0;
while ~feof(fid) & (i<n)
  i = i+1;
  v{i} = fgetl(fid);
end

%% C. Close file
if ischar(fin), fclose(fid); end
