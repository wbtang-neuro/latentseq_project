function [BBwfun, BBwTfun] = prior_kernel_sp(rhoxx,lenxx,tgrid)
% generate functions for transforming between X and U space

% rhoxx affects only the scale of ddB
dt = abs(diff(tgrid));
a = exp(-dt/lenxx); % for bidiagonal matrix Az
vvec = [rhoxx; rhoxx-rhoxx*a.^2]; % for diagonal variance for each
ddB = 1./sqrt(vvec);
ssB = -[a;0].*circshift(ddB,-1);

BBwfun = @(xx,invflag) BBwfun_AR(xx, ddB, ssB, invflag);
BBwTfun = @(xx,invflag) BBwTfun_AR(xx, ddB, ssB, invflag);

end

function BBwxx = BBwfun_AR(xx,ddB,ssB,invflag)

if invflag
    BBwxx = bidiagonal_low_multiply_matrix(ddB,ssB,xx);
else
    BBwxx = tryMLDivide(xx, ssB, ddB, 0);
end

end

function BBwTxx = BBwTfun_AR(xx,ddB,ssB,invflag)

if invflag
    % (not used)
    BBwTxx = bidiagonal_up_multiply_matrix(ddB,ssB,xx);
else
    BBwTxx = tryMLDivide(xx, ssB, ddB, 1);
end

end

function C = tryMLDivide(X, d2, d1, tflag)
% C = D'\X

% N.B. creating a new bidiagonal sparse matrix with each call is slower
% than keeping it in memory and reassigning the new values.

persistent D inds1 inds2 lastd1 lastd2 cache i j

% cast up to double do allow sparse matrix-vector division
n = size(X, 1);
X0 = X;
X = double(gather(X));

if isempty(cache)
    cache = SimpleCache();
end

if isempty(inds1) || ~isequal(numel(inds1), n-1)
    inds1 = (n+2):(n+1):n^2;
    inds2 = 2:(n+1):n^2;
    [i, j] = ind2sub([n, n], [1, inds1, inds2]');
    lastd1 = 0;
    lastd2 = 0;
end

newd1 = ~isequal(d1, lastd1);
newd2 = ~isequal(d2, lastd2);

if newd1 || newd2
    % rather than reassign elements of D, look for cached verion
    key = [n; d1; d2];
    if cache.iskey(key)
        D = cache.get(key);
    else
        vals = [d1; d2(1:n-1)];
        tic()
        D = sparse(i, j, vals, n, n);
        cache.put(key, D);
    end
    lastd1 = d1;
    lastd2 = d2;
end

if tflag
    C = D'\X;
else
    C = D\X;
end

C = cast(C, 'like', X0);

end
