function [jumps,maxJump] = compute_sequenceJumps(X)

% jumps = diag(squareform(pdist(X(~isnan(X(:,1)),:))),1);

if size(X,2)==2
    jumps = sqrt(diff(X(:,1)).^2 + diff(X(:,2)).^2);
else
    jumps = abs(diff(X(:,1)));
end

maxJump = max(jumps);
if isempty(maxJump)==1
    maxJump = NaN;
end

