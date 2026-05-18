function inds = gridnnFast(xg, x, mode)
% Find indices of nearest neighbours in a uniform grid

ngrid = numel(xg);
ngridI = uint32(ngrid);
m = (xg(end)-xg(1)) / (ngrid-1);

switch lower (mode)
    case 'round'
        a = 0;
    case 'ceil'
        a = 0.5;
    case 'floor'
        a = -0.5;
end

inds = uint32(a + (x - xg(1))./m) + 1;
% v0 = inds > 0;
% vG = inds <= ngridI;
% valid = inds > 0 & inds <= ngridI;
% valid = v0 & vG;
% inds(~v0) = 1;
inds(inds > ngridI) = ngridI;
% inds(inds < 1) = 1;
% inds(inds > ngridI) = ngridI;

end