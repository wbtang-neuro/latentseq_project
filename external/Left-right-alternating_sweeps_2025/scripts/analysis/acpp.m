function varargout = acpp(varargin)
%ACPP autocorrelation histogram for point-process data.
% Usage: acpp(eventTimes, binSize, nBins, normalize)
[varargout{1:nargout}] = aCorrPointProcess_safe(varargin{:});
% [varargout{1:nargout}] = aCorrPointProcess(varargin{:}); % DON'T USE THIS VERSION ON MAC
end