function A = reshapeSquare(A)
% reshape array A into a square matrix
n = sqrt(numel(A));
A = reshape(A, n, n);
end