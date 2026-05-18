function g = gsmooth(y, sig)
% gsmooth - smooth vector or matrix by filtering with a Gaussian 
%
% g = gsmooth(x, sig);
%
% Inputs:
%     y [MxN] - matrix or vector (if matrix, operates along columns only)  
%   sig [1x1] - stdev of smoothing Gaussian 
%
% Output:
%     g [MxN] - smoothed signal
%
% Note: normalizes filter to have unit-norm 
%
% See also: bcsmooth

if (sig <= 0)  % Return original if no smoothing width 
    g = y;
else
    [len,wid] = size(y);
    if (len == 1)  % Flip to column vector
        y = y';
        flipped=1;
        len = wid;
    else
        flipped=0;
    end
    [~,nx] = size(y);
    
%     nflt = max(min(len, sig*5),3); % Does not ensure odd-sized kernel

    nflt = max(min(len, ceil(sig*5)),3);
    %make sure nflt is odd
    x = (-nflt:nflt)';
    x = cast(x, "like", y);
    gfilt = normpdf(x, 0, sig);
    gfilt = gfilt./norm(gfilt);
    
    % quantify edge effects by convolving vector ones with kernel
    o = conv(ones(size(y,1),1, "like", y), gfilt, 'same');
    g = convn(y, gfilt, "same") ./ o;
    
    if flipped
        g = g';
    end
end
