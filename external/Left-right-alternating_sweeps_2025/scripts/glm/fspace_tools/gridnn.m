function [inds, valid] = gridnn(xg, x, xgStep)
% Find indices of nearest neighbours in a uniform grid

ngrid = numel(xg);
ngridI = uint32(ngrid);
% m = (xg(end)-xg(1)) / (ngrid-1);

% switch lower (mode)
%     case 'round'
%         fcn = @round;
%     case 'ceil'
%         fcn = @ceil;
%     case 'floor'
%         fcn = @floor;
% end

% inds = uint32(1 + fcn((x - xg(1))/m));
% if nargout==2
%     valid = inds > 0 & inds <= ngridI;
% end
% inds(inds < 1) = 1;
% inds(inds > ngridI) = ngridI;

inds = uint32(1 + floor((x - xg(1))/xgStep));
if nargout==2
    valid = inds > 0 & inds <= ngridI;
end
inds = max(inds, 1);
inds = min(inds, ngridI);
% inds(inds < 1) = 1;
% inds(inds > ngridI) = ngridI;

end