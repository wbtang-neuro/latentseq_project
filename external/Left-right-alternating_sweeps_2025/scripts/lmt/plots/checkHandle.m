function [tf, ax] = checkHandle(hstruct, name, ind)
if isempty(hstruct), hstruct = struct(); end
hasfield = isfield(hstruct, name);
if ~hasfield
    tf = false;
else
    hvalid = isgraphics(hstruct.(name));
    tf = numel(hvalid) >= ind && hvalid(ind);
end

if tf
    ax = hstruct.(name)(ind);
else
    ax = [];
end
end