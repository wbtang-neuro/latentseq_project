function [outputArg1,outputArg2] = plotShank2(p)
%PLOTSHANK Summary of this function goes here
arguments
    p.width = [12, 0];
    p.lims = []
    p.x = 0;
end
width = p.width;
lims = p.lims;
bottom = [-width(1)+p.x,0]; top = [-width(1)+p.x,.5e4];
top = [-width(1)+p.x,1e4];
tip = [0,100];
nshanks = 1;
% 
% Plot shank

vertices = [top-width; top+width; bottom + width; bottom - tip; bottom-width];
% shank = patch(vertices(:, 1), vertices(:, 2), [.8, .8, .8], 'EdgeColor', 'none');
shank = patch(vertices(:, 1), vertices(:, 2), [.8, .8, .8]+.1, 'EdgeColor', [.3,.3,.3], 'LineWidt', .1);

if~isempty(lims)
    % Highlight sampled portion of shank
    bottom = [-width(1)+p.x,lims(1)-20]; top = [-width(1)+p.x,lims(2)+20];
    vertices = [top-width; top+width; bottom + width; bottom; bottom-width];
    shank_recorded = patch(vertices(:, 1), vertices(:, 2), [.3, .3, .3], 'EdgeColor', 'none');
end
end

