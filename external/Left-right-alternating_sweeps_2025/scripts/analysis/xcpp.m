function varargout = xcpp(varargin)
[varargout{1:nargout}] = xCorrPointProcess_safe(varargin{:});
% [varargout{1:nargout}] = xCorrPointProcess(varargin{:}); % DON'T USE THIS VERSION ON MAC
end