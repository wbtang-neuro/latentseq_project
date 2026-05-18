function y = gsmoothcirc(y, sigma, method)
if nargin < 3 || isempty(method), method = "decompose"; end
if strcmpi(method, "decompose")
    % Use a sine/cosine decomposition, which gives us two 1D variables wich
    % we can smooth separately, before converting back into angles. This
    % method will perform best when the input signal is noisy. A caveat is
    % that the output may sometimes in the opposite direction to the 
    % input data, if the dynamics of the signal are faster than the
    % smoothing window.
    [u, v] = pol2cart(y, 1);
    u = gsmooth(u, sigma);
    v = gsmooth(v, sigma);
    y = atan2(v, u);
elseif strcmpi(method, "unwrap")
    % Unwrap the angles to give us a 1D variable which we can smooth as a
    % non-angular 1D variable. This method is only suitable for low-noise
    % signals which can be unwrapped reliably, but for such signals it will
    % follow the dynamics more reliably than the "decompose" method.
    yu = unwrap(y);
    y = gsmooth(yu, sigma);
    y = wrapToPi(y);
else
    error("'%s' is not a valid method");
end

end

% function g = gsmooth(y, sig)
% % gsmooth - smooth vector or matrix by filtering with a Gaussian 
% %
% % g = gsmooth(x, sig);
% %
% % Inputs:
% %     y [MxN] - matrix or vector (if matrix, operates along columns only)  
% %   sig [1x1] - stdev of smoothing Gaussian 
% %
% % Output:
% %     g [MxN] - smoothed signal
% %
% % Note: normalizes filter to have unit-norm 
% %
% % See also: bcsmooth
% 
% if (sig <= 0)  % Return original if no smoothing width 
%     g = y;
% else
%     [len,wid] = size(y);
%     if (len == 1)  % Flip to column vector
%         y = y';
%         flipped=1;
%         len = wid;
%     else
%         flipped=0;
%     end
%     [~,nx] = size(y);
%     
% %     nflt = max(min(len, sig*5),3); % Does not ensure odd-sized kernel
% 
%     nflt = max(min(len, ceil(sig*5)),3);
%     %make sure nflt is odd
%     x = (-nflt:nflt)';
%     x = cast(x, "like", y);
%     gfilt = normpdf(x, 0, sig);
%     gfilt(nflt+1)=0;
%     gfilt = gfilt./norm(gfilt);
%     
%     % quantify edge effects by convolving vector ones with kernel
%     o = conv(ones(size(y,1),1, "like", y), gfilt, 'same');
%     g = convn(y, gfilt, "same") ./ o;
%     
%     if flipped
%         g = g';
%     end
% end
% end
