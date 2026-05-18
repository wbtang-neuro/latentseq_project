function ysm = gsmoothcircdomain(y, sigma)
% linear interpolation for a function whose domain is circular-valued

y = y(:);
yrep = repmat(y, 1, 3);
ysm = gsmooth(yrep(:), sigma);

% Reshape and extract the middle tiling
ysm = reshape(ysm, numel(y), 3);
ysm = ysm(:, 2);

end