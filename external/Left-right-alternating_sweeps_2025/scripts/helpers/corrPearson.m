function r = corrPearson(x, y)
n = size(x, 1);
x = x - sum(x,1)/n;  % Remove mean
y = y - sum(y,1)/n;  % Remove mean
r = x' * y; % 1/(n-1) doesn't matter, renormalizing anyway
dx = vecnorm(x,2,1);
dy = vecnorm(y,2,1);
r = r./dx'; r = r./dy; % coef = coef ./ dx'*dy;
end