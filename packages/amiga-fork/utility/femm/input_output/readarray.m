function v = readarray(fin,numrow,numcol,format)

% function v = readarray(fin,numrow,numcol,format)
% reads an array from file
%
% input parameters
%    fin        : file name or file identifier
%    numrow     : number of rows    (optional; default: [])
%               :   numrow==[]  reads the number of rows from the first line
%               :   numrow==Inf reads until the end of the file
%    numcol     : number of columns (optional; default: [])
%               :   numcol==[]  reads all columns
%    format     : number format     (optional; default: '%f')
%
% output parameters
%    v          : numrow-by-numcol array
%
% Example 1 with file format 1
%    4
%    1 2 3 4
%    2 3 4 5
%    3 4 5 6
%    4 5 6 7
%    5 6 7 8
% v = readarray(fin,[],3)
% will result in
% v = [
%    1 2 3
%    2 3 4
%    3 4 5
%    4 5 6
% ];
%
% Example 2 with file format 2
%    1 2 3 4
%    2 3 4 5
%    3 4 5 6
%    4 5 6 7
%    5 6 7 8
% v = readarray(fin,3,[])
% will result in
% v = [
%    1 2 3 4
%    2 3 4 5
%    3 4 5 6
% ];

if 0
  % test 1
  fid = fopen('test1.txt','wt');
  fprintf(fid,'4\n    1 2 3 4\n    2 3 4 5\n    3 4 5 6\n    4 5 6 7\n    5 6 7 8\n');
  fclose(fid);
  v = readarray('test1.txt',[],5)  % also test with 5 instead of 3
  % test 2
  fid = fopen('test2.txt','wt');
  fprintf(fid,'    1 2 3 4\n    2 3 4 5\n    3 4 5 6\n    4 5 6 7\n    5 6 7 8\n');
  fclose(fid);
  v = readarray('test2.txt',3,[])
end

%% A. Parameter check
if ~exist('numrow','var'), numrow = [];   end
if ~exist('numcol','var'), numcol = [];   end
if ~exist('format','var'), format = '%f'; end
if ischar(fin), fid = fopen(fin,'r'); else fid = fin; end

%% B. Determine the number of rows to be read
if isempty(numrow)
  numrow = fscanf(fid,'%d',1);          % read the first number
  ln = fgetl(fid);                      % discard the line termination
end

%% C. Determine the number of columns to be read, scan all lines and combine into the right format
if numrow~=0
  ln = fgetl(fid);                        % read the first line to detect the number of columns
  [v1,numcol1] = sscanf(ln,format);
  v = [ v1'; fscanf(fid,format,[numcol1 numrow-1])' ];
  ln = fgetl(fid);                      % discard the last line termination
  if ~isempty(numcol)
    v = [ v(:,1:min(numcol,numcol1)) zeros(size(v,1),max(0,numcol-numcol1)) ];
  end
else
  if ~isempty(numcol) & ~isinf(numcol)
    v = zeros(0,numcol);
  else
    v = zeros(0,10);   % 10 is set arbitrarily
  end
end

%% D. Close file
if ischar(fin), fclose(fid); end
