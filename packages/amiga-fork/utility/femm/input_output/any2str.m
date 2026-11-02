function str = any2str(a)

% function str = any2str(a)
% puts a value of any type into a string
%
% input parameters
%    a    : value (of any type)
%
% output parameters
%    str  : string

if isnumeric(a)
  if prod(size(a))>1,  str = mat2str(a);
  elseif isinteger(a), str = int2str(a);
  elseif isfloat(a),   str = num2str(a);
  else,                str = evalc('disp a');
  end
elseif ischar(a),      str = a;
elseif iscell(a),
  allstr = cellfun(@any2str,a,'UniformOutput',false);
  switch length(allstr)
    case 0,            str = '{}';
    case 1,            str = sprintf('{%s}',allstr{1});
    otherwise,         str = [ '{' allstr{1} sprintf(',%s',allstr{2:length(allstr)}) '}' ];
  end
else,                  str = evalc('disp a');
end
