function [C,numberOfOverlapPixels] = normxcorr2_general_stack_helper(varargin)
% Adapted to handle stacks of 2D images

[T, A, requiredNumberOfOverlapPixels] = ParseInputs(varargin{:});

sizeA = size(A);
sizeT = size(T);

% Find the number of pixels used for the calculation as the two images are
% correlated.  The size of this image will be the same as the correlation
% image.  
numberOfOverlapPixels = local_sum(ones(sizeA(1:2)),sizeT(1),sizeT(2));

local_sum_A = local_sum(A,sizeT(1),sizeT(2));
local_sum_A2 = local_sum(A.*A,sizeT(1),sizeT(2));

% Note: diff_local_sums should be nonnegative, but it may have negative
% values due to round off errors. Below, we use max to ensure the radicand
% is nonnegative.
diff_local_sums_A = ( local_sum_A2 - (local_sum_A.^2)./ numberOfOverlapPixels );
clear local_sum_A2;
denom_A = max(diff_local_sums_A,0); 
clear diff_local_sums_A;

% Flip T in both dimensions so that its correlation can be more easily
% handled.
rotatedT = rot90(T,2);
local_sum_T = local_sum(rotatedT,sizeA(1),sizeA(2));
local_sum_T2 = local_sum(rotatedT.*rotatedT,sizeA(1),sizeA(2));

diff_local_sums_T = ( local_sum_T2 - (local_sum_T.^2)./ numberOfOverlapPixels );
denom_T = max(diff_local_sums_T,0); 
denom = sqrt(denom_T .* denom_A);


T_size = size(T);
A_size = size(A);
outsize = A_size(1:2) + T_size - 1;
xcorr_TA = freqxcorr(T, A, outsize); % this is much faster for 3D stacks
% xcorr_TA = convn(rot90(T,2),A);

numerator = xcorr_TA - local_sum_A .* local_sum_T ./ numberOfOverlapPixels;
clear xcorr_TA local_sum_A;

% denom is the sqrt of the product of positive numbers so it must be
% positive or zero.  Therefore, the only danger in dividing the numerator
% by the denominator is when dividing by zero. We know denom_T~=0 from
% input parsing; so denom is only zero where denom_A is zero, and in these
% locations, C is also zero.
% C = zeros(size(numerator), 'like', numerator);
tol = 1000*eps( max(abs(denom(:))) );
% i_nonzero = denom > tol;
% C(i_nonzero) = numerator(i_nonzero) ./ denom(i_nonzero);
C = numerator ./ denom;
C(denom<tol) = nan;
clear numerator denom;

% Remove the border values since they result from calculations using very
% few pixels and are thus statistically unstable.
% By default, requiredNumberOfOverlapPixels = 0, so C is not modified.
if( requiredNumberOfOverlapPixels > max(numberOfOverlapPixels(:)) )
    error(['ERROR: requiredNumberOfOverlapPixels ' num2str(requiredNumberOfOverlapPixels) ...
    ' must not be greater than the maximum number of overlap pixels ' ...
    num2str(max(numberOfOverlapPixels(:))) '.']);
end
C(numberOfOverlapPixels < requiredNumberOfOverlapPixels, :) = 0;

end


%-------------------------------
% Function  local_sum
%
function local_sum_A = local_sum(A,m,n)

% This algorithm depends on precomputing running sums.

% If m,n are equal to the size of A, a faster method can be used for
% calculating the local sum.  Otherwise, the slower but more general method
% can be used.  The faster method is more than twice as fast and is also
% less memory intensive. 
if( m == size(A,1) && n == size(A,2) )
    s = cumsum(A,1);
    c = [s; s(end,:,:) - s(1:end-1,:,:)];
    s = cumsum(c,2);
    local_sum_A = [s, s(:,end, :) - s(:,1:end-1, :)];
else
    % Break the padding into parts to save on memory. Clear variables as
    % soon as they're no longer needed.
    B = zeros(size(A,1)+2*m,size(A,2),size(A,3),'like',A);
    B(m+1:m+size(A,1),:,:) = A;
    s = cumsum(B,1);
    clear B
    c = s(1+m:end-1,:,:)-s(1:end-m-1,:,:);
    d = zeros(size(c,1),size(c,2)+2*n,size(A,3),'like',A);
    d(:,n+1:n+size(c,2), :) = c;
    clear c
    s = cumsum(d,2);
    clear d
    local_sum_A = s(:,1+n:end-1,:)-s(:,1:end-n-1,:);
end

% local_sum_A = gather(local_sum_A);

end

%-----------------------------------------------------------------------------
function [T, A, requiredNumberOfOverlapPixels] = ParseInputs(varargin)

if( nargin < 2 || nargin > 3 )
    error('ERROR: The number of arguments must be either 2 or 3.  Please see the documentation for details.');
end

T = varargin{1};
A = varargin{2};

if( nargin == 3 )
    requiredNumberOfOverlapPixels =  varargin{3};
else
    requiredNumberOfOverlapPixels = 0;
end


checkSizesTandA(T,A)

% See geck 342320. If either A or T has a minimum value which is negative, we
% need to shift the array so all values are positive to ensure numerically
% robust results for the normalized cross-correlation.
A = shiftData(A);
T = shiftData(T);

checkIfFlat(T);

end

%-----------------------------------------------------------------------------
function B = shiftData(A)

B = A;
% B = double(A);

is_unsigned = isa(A,'uint8') || isa(A,'uint16') || isa(A,'uint32');
if ~is_unsigned
    
    min_B = min(B(:)); 
    
    if min_B < 0
        B = B - min_B;
    end
    
end

end

%-----------------------------------------------------------------------------
function checkSizesTandA(T,A)

if numel(T) < 2
    eid = sprintf('Images:%s:invalidTemplate',mfilename);
    msg = 'TEMPLATE must contain at least 2 elements.';
    error(eid,'%s',msg);
end

end

%-----------------------------------------------------------------------------
function checkIfFlat(T)

if std(T(:)) == 0
    eid = sprintf('Images:%s:sameElementsInTemplate',mfilename);
    msg = 'The values of TEMPLATE cannot all be the same.';
    error(eid,'%s',msg);
end

end

function xcorr_ab = freqxcorr(a,b,outsize)
  
% Find the next largest size that is a multiple of a combination of 2, 3,
% and/or 5.  This makes the FFT calculation much faster.
optimalSize(1) = FindClosestValidDimension(outsize(1));
optimalSize(2) = FindClosestValidDimension(outsize(2));

% Calculate correlation in frequency domain
Fa = fft2(rot90(a,2),optimalSize(1),optimalSize(2));
Fb = fft2(b,optimalSize(1),optimalSize(2));
xcorr_ab = real(ifft2(Fa .* Fb));
% xcorr_ab = gather(xcorr_ab);
xcorr_ab = xcorr_ab(1:outsize(1),1:outsize(2),:);
end


%-----------------------------------------------------------------------------
function [newNumber] = FindClosestValidDimension(n)

% Find the closest valid dimension above the desired dimension.  This
% will be a combination of 2s, 3s, and 5s.

% Incrementally add 1 to the size until
% we reach a size that can be properly factored.
newNumber = n;
result = 0;
newNumber = newNumber - 1;
while( result ~= 1 )
    newNumber = newNumber + 1;
    result = FactorizeNumber(newNumber);
end

end

%-----------------------------------------------------------------------------
function [n] = FactorizeNumber(n)

for ifac = [2 3 5]
    while( rem(n,ifac) == 0 )
        n = n/ifac;
    end
end

end