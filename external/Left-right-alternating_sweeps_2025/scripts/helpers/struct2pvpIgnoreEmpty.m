function c = struct2pvpIgnoreEmpty(s, fieldList)
%Adapted from struct2arglist
%
%Ignores empty-valued parameters

if nargin < 2, fieldList = fieldnames(s); end

valid = arrayfun(@(field) ~isempty(s.(field)), fieldList);
f = fieldList(valid);

c = cell(1,2*length(f));
for i = 1:length(f)
  c{2*i-1} = f{i};
  c{2*i} = s.(f{i});
end
