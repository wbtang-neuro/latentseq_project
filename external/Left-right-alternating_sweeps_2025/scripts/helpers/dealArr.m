function varargout = dealArr(A)
% Like deal(), but works elementwise on array input
A = num2cell(A);
[varargout{1:nargout}] = deal(A{:});
end