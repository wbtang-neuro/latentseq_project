function [pSel, inds] = extractParams(p, fieldIndsStruct)
fields = fieldnames(fieldIndsStruct);
inds = [];
for f = 1:numel(fields)
    fd = fields{f};
    inds = [inds; fieldIndsStruct.(fd)(:)];
end
inds = sort(inds);
pSel = p(inds);
end