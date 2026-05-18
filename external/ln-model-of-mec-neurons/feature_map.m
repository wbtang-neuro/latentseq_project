function [feature_grid,featureVec] = feature_map(feature,nbins)

max_feature = max(feature); 
min_feature = min(feature);
range_feature = max_feature-min_feature;
% featureVec = min_feature+range_feature/nbins/2:range_feature/nbins:max_feature-range_feature/nbins/2;
featureVec = linspace(min_feature + range_feature/nbins/2, ...
                      max_feature - range_feature/nbins/2, ...
                      nbins);
feature_grid = zeros(numel(feature),numel(featureVec));

for i = 1:numel(feature)
    % figure out the feature index
    [~, idx] = min(abs(feature(i)-featureVec));
    feature_grid(i,idx) = 1;
end
