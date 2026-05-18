function hpatch = plotShankOverview(units)
%PLOTSHANKOVERVIEW Summary of this function goes here
%   Detailed explanation goes here
%%
shanks = unique([units.shank]);
edges = 0:50:5000;
centers = edg2cen(edges);
S = SweepsSettings;
% clf
nshank = 4;
if numel(shanks)==1
    nshank = 1;
end
xpos = [0:250:750];
for s = 1:nshank
    
    us = units([units.shank]==s);
    
    lims = minmax([us.shankPos]);
  
    idx = [us.isId];
    hpatch(1) = patchCellHistogram(units = us, idx = idx & ~[us.isGrid], edges=edges, color=S.col_id, xpos=xpos(s));
    hpatch(2) = patchCellHistogram(units = us, idx = idx&[us.isGrid], edges=edges, color=S.col_conj, xpos=xpos(s));
    idx = [us.isGrid]&~idx;
    hpatch(3) = patchCellHistogram(units = us, idx = idx, edges=edges, color=[0,0,1], xpos=xpos(s));
    plotShank2(lims = lims, x=xpos(s));
end
legend(hpatch, ["ID", "Conj", "Pure grid"])
end

