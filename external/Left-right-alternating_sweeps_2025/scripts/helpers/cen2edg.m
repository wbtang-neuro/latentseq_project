function edg = cen2edg(cen, bw)
%CEN2EDG convert bin edges
if nargin < 2 || isempty(bw)
    bw = cen(2)-cen(1);
end
dim = find(size(cen)>1);
edg = cat(dim, cen-bw/2, cen(end)+bw/2);
end