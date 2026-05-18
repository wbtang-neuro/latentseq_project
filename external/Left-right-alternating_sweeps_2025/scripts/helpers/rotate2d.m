function varargout = rotate2d(varargin)

if nargin == 2
    % Two inputs: prot = rotate2d(p, alpha)
    p = varargin{1};
    x = p(:, 1);
    y = p(:, 2);
    alpha = varargin{2};
    [xr, yr] = rot2dlocal(x, y, alpha);
elseif nargin == 3
    % Three inputs: [xrot, yrot] = rotate2d(x, y, alpha)
    x = varargin{1};
    y = varargin{2};
    alpha = varargin{3};
    [xr, yr] = rot2dlocal(x, y, alpha);
else
    error("Invalid number of arguments");
end

if nargout == 1
    varargout{1} = [xr(:), yr(:)];
elseif nargout == 2
    varargout{1} = xr;
    varargout{2} = yr;
else
    error("Invalid number of outputs")
end

end


function [xr, yr] = rot2dlocal(x, y, alpha)

[theta, rho] = cart2pol(x, y);
theta = theta + alpha;
[xr, yr] = pol2cart(theta, rho);

end