function [posgrid, bins_x,bins_y] = feauture2D_map(pos, nbins)

% take the histogram
maxposx = max(pos(:,1));
minposx = min(pos(:,1));
range_posx = maxposx-minposx;
bins_x = minposx + range_posx/nbins/2:range_posx/nbins:maxposx-range_posx/nbins/2;


% take the histogram
maxposy = max(pos(:,2));
minposy = min(pos(:,2));
range_posy = maxposy-minposy;
bins_y = minposy + range_posy/nbins/2:range_posy/nbins:maxposy-range_posy/nbins/2;

% store grid
posgrid = zeros(length(pos), nbins^2);

% loop over positions
for idx = 1:length(pos)
    
    % figure out the position index
    [~, xcoor] = min(abs(pos(idx,1)-bins_x));
    [~, ycoor] = min(abs(pos(idx,2)-bins_y));
    
    bin_idx = sub2ind([nbins, nbins], nbins - ycoor + 1, xcoor);
    
    posgrid(idx, bin_idx) = 1;
    
end

end