function straightness = movingStraightness(pos,window)
%MOVINGSTRAIGHTNESS moving-window quantifiaction of a 2D trajectory's
%directness.

% Get distances
dpos = diff(pos);
dists = hypot(dpos(:,1), dpos(:,2));
travelDist = movsum(dists, window);
dshortest = circshift(pos,-window(1))-circshift(pos, window(2));
shortestDist = hypot(dshortest(:,1), dshortest(:,2));
straightness = travelDist./shortestDist(1:end-1);

straightness(1:window(1)) = nan;       % assuming window(1)
straightness(end-window(1):end) = nan;
straightness(end+1)=nan;
end

