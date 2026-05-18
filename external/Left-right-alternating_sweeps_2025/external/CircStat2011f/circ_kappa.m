function kappa = circ_kappa(alpha,w)
%
% kappa = circ_kappa(alpha,[w])
%   Computes an approximation to the ML estimate of the concentration 
%   parameter kappa of the von Mises distribution.
%
%   Input:
%     alpha   angles in radians OR alpha is length resultant
%     [w      number of incidences in case of binned angle data]
%
%   Output:
%     kappa   estimated value of kappa
%
%   References:
%     Statistical analysis of circular data, Fisher, equation p. 88
%
% Circular Statistics Toolbox for Matlab

% By Philipp Berens, 2009
% berens@tuebingen.mpg.de - www.kyb.mpg.de/~berens/circStat.html


alpha = alpha(:);

if nargin<2
  % if no specific weighting has been specified
  % assume no binning has taken place
	w = ones(size(alpha));
end

N = length(alpha);

if N>1
  R = circ_r(alpha,w);
else
  R = alpha;
end

v1 = R < 0.53;
v2 = R >= 0.53 & R < 0.85;
v3 = R >= 0.85;

kappa = nan(1, size(w, 2));
kappa(v1) = 2.*R(v1) + R(v1).^3 + 5.*R(v1).^5/6;
kappa(v2) = -.4 + 1.39.*R(v2) + 0.43./(1-R(v2));
kappa(v3) = 1./(R(v3).^3 - 4.*R(v3).^2 + 3.*R(v3));

% if R < 0.53
%   kappa = 2*R + R^3 + 5*R^5/6;
% elseif R>=0.53 && R<0.85
%   kappa = -.4 + 1.39*R + 0.43/(1-R);
% else
%   kappa = 1/(R^3 - 4*R^2 + 3*R);
% end

if N<15 && N>1
  if kappa < 2
    kappa = max(kappa-2*(N*kappa)^-1,0);    
  else
    kappa = (N-1)^3*kappa/(N^3+N);
  end
end
