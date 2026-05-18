function cen = edg2cen(edg)
% d = edg(2)-edg(1);
d = diff(edg);
cen = edg(1:end-1)+d/2;
end