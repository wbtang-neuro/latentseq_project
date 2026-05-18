function [outputArg1,outputArg2] = plotOnePair(units, xcorr)
%PLOTONEPAIR Summary of this function goes here
%   Detailed explanation goes here
S = SweepsSettings;


% Plot xcorr
tlags = linspace(-50,50, 101)'+.2; % for plotting xcorr spike
bar(nexttile, tlags, xcorr, 1, 'FaceColor', [1,1,1]*.4); axis square; xline(0, 'r'), xticks([]); yticks([0, max(xcorr)]);
xlim([-20,20]); ylim([0, max(xcorr)]);
title(units(1).id+" to "+units(2).id, 'FontWeight','normal');

% Plot dir tuning curves
rms = arrayfun(@(u) u.rmf.id, units);
cols = [S.("col_"+units(1).cellType);S.("col_"+units(2).cellType)];
nexttile; plotRateMap(rms, "colors", cols);
set(gca, "ThetaTickLabels",[])

meandir = @(units, vname) arrayfun(@(u) mod(circ_mean(S.gv.angular, u.rmf.(vname).z), 2*pi), units);
if ismember(units(1).cellType, ["id", "conjunctive"]) && ismember(units(2).cellType, ["id", "conjunctive"])
%     str = sprintf("Pref dir (%s): %.1f%s", units(1).cellType,rad2deg(meandir(units(1), "id")), char(176));
%     str = str+sprintf("\nPref dir (%s): %.1f%s", units(2).cellType,rad2deg(meandir(units(2), "id")), char(176));
    str = sprintf("Pref dir: %.1f%s",rad2deg(meandir(units(1), "id")), char(176));
    str = str+sprintf("\nPref dir (): %.1f%s", rad2deg(meandir(units(2), "id")), char(176));
elseif ismember(units(1).cellType, ["id", "conjunctive"])
%     str = sprintf("Pref dir (%s): %.1f%s", units(1).cellType,rad2deg(meandir(units(1), "id")), char(176));
    str = sprintf("Pref dir (): %.1f%s",rad2deg(meandir(units(1), "id")), char(176));
else
    str = "";
end
title(str, 'FontWeight','normal');


% Plot overlayed spatial tuning
vbins = restrictq(S.gv.pos_of_fine_extended, [min(S.gv.pos_of_fine)-.01, max(S.gv.pos_of_fine)+.01]);
gv = S.gv.pos_of_fine_extended(vbins);
plotcontour.id = false;plotcontour.conjunctive = true; plotcontour.grid = true;

rms = arrayfun(@(u) u.rmf.pos_id_shift, units);
rms = arrayfun(@(u) u.rmf.pos, units);
rms1 = rms(:, 1);
rms2 = rms(:, 2);

z1 = rms1.z(vbins, vbins); z1 = regionfill(z1, isnan(z1));%z1(isnan(z1))=0;
z2 = rms2.z(vbins, vbins); z2 = regionfill(z2, isnan(z2));%z2(isnan(z2))=0;

ax = nexttile; 
%         plotRatemapColor(z1', [1,0,0]); plotRatemapColor(z2', [.3,0.3,0.3]);
if units(1).cellType==units(2).cellType
    plotRatemapColor(z1', S.("col_"+units(1).cellType), plotcontour.(units(1).cellType)); 
    plotRatemapColor(z2', [.5,.5,.5], plotcontour.(units(2).cellType));
else

    plotRatemapColor(z1', S.("col_"+units(1).cellType), plotcontour.(units(1).cellType)); 
    plotRatemapColor(z2', S.("col_"+units(2).cellType), plotcontour.(units(2).cellType));
end
ax.YDir = "normal";% title('Superimosed rate maps'); ax.FontSize = 12;

if units(2).cellType == "grid"
    [phaseoffset, shiftdir, shiftdist] = getPhaseOffset(units(1), units(2));
    shiftdir = wrapToPi(shiftdir-pi);
    shiftdir = mod(shiftdir, 2*pi);
    title(sprintf("Grid phase offset: %.1f", rad2deg(shiftdir)), 'FontWeight','normal')
end
%( Plot spatial xcorr)
% if units(1).cellType=="conjunctive" && units(2).celltype=="grid";
   
% end
end

function  [phaseoffset, shiftdir, shiftdist] = getPhaseOffset(unit1, unit2)
% Todo: Can interpolate xcorr for better resolution
S = SweepsSettings;
% Should be fine to include a bit more here. Maybe the whole thing? No,
% some junk along the edges.
vbins = restrictq(S.gv.pos_of_fine_extended, [-.85, .85]);

vbins = restrictq(S.gv.pos_of_fine_extended, [min(S.gv.pos_of_fine)-.01, max(S.gv.pos_of_fine)+.01]);
% vbins(:)=true;
rms1 = unit1.rmf.pos_id_shift;
rms2 = unit2.rmf.pos_id_shift;
      
z1 = rms1.z(vbins, vbins); z1 = regionfill(z1, isnan(z1));%z1(isnan(z1))=0;
z2 = rms2.z(vbins, vbins); z2 = regionfill(z2, isnan(z2));%z2(isnan(z2))=0;

% z1 = imgaussfilt(z1, 2);
% z2 = imgaussfilt(z2, 2);
        
% xc = xcorr2(z2', z1');
xc = normxcorr2_general(z2', z1');
% %     peaks = fastPeakFind(xc(30:90, 30:90));
bw = imregionalmax(xc);
s = regionprops(bw);
peaks = cat(1, s.Centroid);
origo = size(xc)./2;
pkrel = peaks-origo;
if ~isempty(peaks)
    [~,minIdx] = min(hypot(pkrel(:, 1), pkrel(:, 2)));
end
nearestPk = peaks(minIdx, :)-origo;


[shiftdir, shiftdist] = cart2pol(nearestPk(1), nearestPk(2));
phaseoffset = nearestPk;

% imagesc(nexttile, xc(50:150, 50:150))
% imagesc(nexttile, xc); axis square
% xline(origo, 'w'); yline(origo, 'w')
% scatter(peaks(minIdx, 1), peaks(minIdx, 2), 'r')
end
