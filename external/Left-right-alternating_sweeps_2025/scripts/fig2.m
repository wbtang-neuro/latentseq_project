%% 2a Plot example ID cells (hd tuning, pos tuning and acorr)
sweepsSetup
units = getExampleIdCells(process=0, load=1);
plotExampleIdCells(units, plotPos=true);
%% 2b Plot shank overview with grid. and theta dir cells from example rat
fname = fullfile(S.codeRoot, "results", "unit_data_25953.mat");
tmp = load(fname);
units = tmp.units;
clf
plotShankOverview(units);
xlim([-50,900])
ylim([-100,5000])
ax = gca;
ax.XAxis.Visible="off";
ax.YTick = [0, 5000];
text([0:250:750]-75, repmat(5150, 1, 4), ["1", "2", "3", "4"], "FontSize",8);
pbaspect([1,2,1])
plot([0,0,200,200], [25,0,0,25]-70, 'k')
text(-5,-130, "0"); text(150,-130, "100")
text(120, -250, "Fraction of total (%)", 'HorizontalAlignment','center')
ylabel("Distance from shank tip (microns)")

%% 2c Plot example ID raster
sweepsSetup
fname = fullfile(S.dataRoot, "navigation", "of", "25843_1.mat");
tmp = load(fname);
D = tmp.Dsession;
res_id = runPvIdDecoding(recs="25843_1");
res_pos = runPvPosDecoding(recs="25843_1");
%%
D.thetaChunks = res_id;
D.sweeps = res_pos.sweeps;
fontSz = 10;
trngPad = [-10, 20];
trngShort = [9910, 9914];
[vtShort, t] = restrictq(D.t, trngShort);
hdMean = circ_mean(D.hd(vtShort));
trng = trngShort + trngPad;
[vt, t] = restrictq(D.t, trng);
clear axs
clf
ax = subplot(3, 4, [1:4]);
axs(1) = ax;
units = D.units.mec;
mvl = @(units, vname) arrayfun(@(u) u.rmf.(vname).mvl, units);
units = units(mvl(units, "hd")>.3&mvl(units, "theta")>.3);
meandir = @(units, vname) arrayfun(@(u) circ_mean(D.gv.hd, u.rmf.(vname).z), units);
prefdir = meandir(units, "hd");
prefdir = mod(prefdir-hdMean+0*pi, 2*pi);
y = D.hd(vt);
y = mod(y-hdMean+pi*0, 2*pi);
ydiff = circ_dist(y, circshift(y,1));
y(ydiff>1)=nan;

raster2(units, timeRange=trng, yvar=prefdir, yHeight=.05, color=[0,0,0], tickWidth=.5)
plot(ax, t, y, 'blue')
formatCircAxes(ax, "y", "Direction");
margin = 45;
degRange = [-margin, 360+margin]+180;
for i = 1:2 %remove samples outside range
    h = ax.Children(i);
    vt = restrictq(h.YData, degRange);
    h.YData(~vt)=nan;
end
ax.XLim = ax.XLim;
ax.YLim = ax.YLim;
ax.XTick = [];
ax.FontSize = fontSz;

xline(trngShort, 'r')
xline(5)
ax.YLim = degRange;
ax.YTick = [0,360]+180; 

ax.YTickLabel = ["0", "360"];
ax.XLim = [trng];
ax.Clipping = "off";
scale_bar = plot(ax, trng(1)+[0, 5], ax.YLim(1) - 20 + [0, 0], "k", "lineWidth", 1);
text(trng(1)+2.5, ax.YLim(1) - 45, "5s", 'HorizontalAlignment','center', 'FontSize',fontSz)
plot([trngShort(1),trng(1), nan, trngShort(2),trng(2)], [0, -150, nan, 0, -150]+ax.YLim(1), 'r')
%%%%%%%%%%% EXAMPLE ALTERNATING PERIODS %%%%%%%%%%%%%%
trng = trngShort;
[vt, t] = restrictq(D.t, trng);
t0 = t(1);

ax = subplot(3, 4, [5:8]);
axs(2) = ax;
y = D.hd(vt);
y = mod(y-hdMean+pi, 2*pi);
s = plot(ax, t, y, "blue");
raster2(units, timeRange=trng, yvar=prefdir+pi, yHeight=.13, color=[0,0,0], tickWidth=.5)
% Add decoded ID for the second half
[vcyc, cyct] = restrictq(D.thetaChunks.tCen, trng);
x = cyct; y = (D.thetaChunks.id(vcyc));
idx = cyct>(cyct(1)+.5*diff(minmax(cyct')));
col =[23, 163, 60]./255;
y = mod(y-hdMean+pi, 2*pi);
hid = plot(ax, x(idx), y(idx),  '-o',"color",col, "lineWidth", 1); % hollow
formatCircAxes(ax, "y", "Direction");
margin = 45;
degRange = [-margin, 360+margin]+180*0;
for i = 1:3 %remove samples outside range
    h = ax.Children(i);
    vt = restrictq(h.YData, degRange);
    h.YData(~vt)=nan;
end
x0 = ax.XLim(1);
y0 = ax.YLim(1);

ax.XLim = ax.XLim;
ax.YLim = ax.YLim;
ax.XTick = [];
ax.FontSize = fontSz;

xline(5)

ax.YLim = degRange;
ax.YTick = [0,360]; 

ax.YTickLabel = ["0", "360"];
ax.XLim = [trng];
plotShadedThetaCycles(gca, D.thetaChunks.tStart, trng);
ax.Clipping = "off";
ax.XAxis.Visible = "on";
scale_bar = plot(ax, trng(1)+[0, .5], ax.YLim(1) - 20 + [0, 0], "k", "lineWidth", 1);
text(trng(1)+.25, ax.YLim(1) - 45, "0.5s", 'HorizontalAlignment','center', 'FontSize',fontSz)

ax = subplot(3, 4, [9:11]);
ax.FontSize = fontSz;
ratSpeed = .8;
stepSize = ratSpeed*.01;
x = 0:stepSize:round(numel(D.t)*stepSize);
chk = D.thetaChunks;
chk.id=circ_dist(chk.id, chk.hd);
sweeps = D.sweeps;
chk.sweepdir = circ_dist([sweeps.hpfPosDirection]', chk.hd);
trng = [1.1498, 1.1502]*1e4;
trng=trng+[2.255,-.45];
h = sweepIdTrajectoryPlotShaded(ax, D.t, x, x*0, trng, sweeps=chk, id=chk, col = [.6,.6,.8]); % This func can probably be used to plot
axis on
ax.XAxis.Visible = "off";
ax.YLim = [-.1,.1]; ax.YTick = [-.1,.1]; ax.YTickLabel = ["Right", "Left"];
ax.Clipping = "off";
egoego_heatmap(x=chk.id, y=chk.sweepdir, ax=subplot(3,4,12), lims=70, nsmooth=1)

%% Run/load PV directional decoding
sweepsSetup
recs = runPvIdDecoding(process=1, save=0);
res_pos = runPvPosDecoding(process=1, save=0);
%%
clf
clear nopposite egohists
egohists = struct("all", [],"right", [],"left", []);
nrecs = numel(recs);
for r = 1:nrecs
    dec = recs(r);
    id = dec.id;
    hd = dec.hd;
    egoid = circ_dist(id, hd);
    egoid(dec.speed<S.minSpeed) = nan;
    
    [isRight,isLeft, prevRight, prevLeft] = egoRightLeft(egoid);
    
    egohists.right(:,r) = histcounts(egoid(prevLeft), linspace(-pi, pi, 101));
    egohists.left(:,r) = histcounts(egoid(prevRight), linspace(-pi, pi, 101));

    [~, acorrs(r, :)] = circAlternationAcorrAdjacent(circ_diff(egoid), 7);
    [recs(r).prcAltern, recs(r).pAltern] = computeAlternationPercent(egoid);
    [recs(r).modes, lefthist(r, :), righthist(r, :)] = computeModesKsd(egoid);
    [recs(r).hdoffset, recs(r).hdcorr, recs(r).hdpval, recs(r).absHdoffset] = computeDirAlignment(id, hd, 1);
    recs(r).ndir = dec.ndir;
    recs(r).sweeps = res_pos(r).sweeps;
    if r == 5
        egodir = egoid;
    end
end
%% Print stats
modes = rad2deg(circ_mean(cat(1, recs.modes)));
fprintf("Modes: %.2f, %.2f\n", modes(1), modes(2))
modes = (cat(1, recs.modes));
modes(:, 1) = modes(:, 1)*-1;
fprintf("Modes: %.2f, %.2f\n", rad2deg(circ_mean(modes(:))), rad2deg(circ_std(modes(:)))./sqrt(numel(recs)))
fprintf("Alternation: %.2f, %.2f, max p=%.3f\n", 100*mean([recs.prcAltern]), 100*std([recs.prcAltern])./sqrt(numel(recs)),max([recs.pAltern]) );
fprintf("ID vs HD: mean=%.1f, %.2f, rho=%.2f, %.3f, max p=%.4f  \n", rad2deg(mean(abs([recs.hdoffset]'))), rad2deg(std(abs([recs.hdoffset]')))./sqrt(numel(recs)), mean([recs.hdcorr]'), std([recs.hdcorr]')./sqrt(numel(recs)), max([recs.hdpval]'))

%% d-e plots
figure("WindowStyle","normal")
tl = tiledlayout(2,2, 'TileIndexing','rowmajor');
plotEgoHistPolar(tl = tl, egodir=egodir,alpha=1)
pax = polaraxes(tl);
pax.ThetaZeroLocation = "top";
pax.Layout.Tile = 2;
for r = 1:nrecs
    h = egohists.left(:, r);
    h = gsmooth(h, .6);
    polarplot(pax, linspace(-pi,pi, 100), h./sum(h), 'Color',[.2,.2,1, .5]);
    h = egohists.right(:, r);
    h = gsmooth(h, .6);
    polarplot(pax, linspace(-pi,pi, 100), h./sum(h), 'Color',[1,.2,.2, .5]);
end
pax.ThetaZeroLocation = "top"; 
pax.ThetaTick = [0, 90, 180, 270]; 
pax.ThetaTickLabel = ["Front", "Left", "Back", "Right"]; 
pax.RTick = []; pax.FontSize = 10;
scatterEgodirAcorrs(acorrs = acorrs, plotHalf=0, ax=nexttile([1,2]));ylim([-1,1])
set(gca, "FontSize", 10)

%% Plot ID vs sweeps, panel g
mean_center=0;
clf
dirType = "hpfPosDirection";
clear muoffset absoffsets
for r = 1:nrecs
    D = recs(r);
    dec = recs(r);
    id = dec.id;
    hd = dec.hd;
    egoid = circ_dist(id, hd);
    egoid(dec.speed<S.minSpeed) = nan;
    
    % Get sweep dir (ego allo)
    sweeps = D.sweeps;
    vswp = D.vswp;
    sd = [sweeps.(dirType)]';
    sd(~vswp)=nan;
    egosd = circ_dist(sd, hd);
    
    % compute coalternation
    [coflickering(r), pCoflicker(r)] = computeCoAlternation([egosd, egoid], 0);
    coflickering(r) = coflickering(r)*100;
    
    if r==1 
        dirs = [egosd, egoid]; % Keep for example plot
    end
    
end

%%
clf
tiledlayout(1,4);
egoego_heatmap("ax",nexttile([1,3]), "computeStats",true, "x",dirs(:, 1), "y",dirs(:, 2), "gve",linspace(-pi, pi, 90), 'lims', 60, "nsmooth",1);
swarmchart(nexttile, coflickering*0, coflickering, 6,[1,1,1]*.7)
ylim([0,100])
yline(50)
xlim([-1,1])
errorbar(0, mean(coflickering), std(coflickering), 'k')
scatter(0, mean(coflickering), 15, 'red', 'MarkerEdgeColor', 'k')
xticks([])
yticks([0,50,100])
ax = gca;
ax.FontSize = 12;


%% Functions
function [res] = getExampleIdCells(p)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
arguments
    p.datafld = [];
    p.process = 1;
    p.load = 1;
    p.save = 0;
end
S = SweepsSettings;
if isempty(p.datafld)
    p.datafld = fullfile(S.dataRoot, "navigation", "of");
end

if p.process
    i = 1; res = struct;
    recs = S.example_idcells_main;
    recNames = char(fields(recs)');
    
    % acorr settings
    binSize = 5e-3;
    nBins = 200;
    norm = 0;
    for r = 1:size(recNames,1)
        recName = string(recNames(r, 4:end));
        disp(recName+"...")
        fname = fullfile(p.datafld, recName+".mat");
        tmp = load(fname);
        D = tmp.Dsession;
        units = D.units.mec;
        uids = recs.(recNames(r, :));
        for u = 1:numel(uids)
            unit = units([units.id]==uids(u));
            spkspeed = D.speed(unit.spikeInds);
            spkt = unit.spikeTimes(spkspeed>.2);
            [rate, lags] = acpp(spkt, binSize, nBins);
            res(i).recName = recName;
            res(i).id = unit.id;
            res(i).hdTuning = unit.rmf.hd;
            res(i).posTuning = unit.rmf.pos;
            res(i).tempAcorr.rate = rate;
            res(i).tempAcorr.lags = lags;
            i = i+1;
        end
        disp("Done")
    end
    if p.save
        fld = "results";
        if~isempty(dir("*code"))
            fld = fullfile("code", "results");
        end
        save(fullfile(fld, "exampleIdCells.mat"), "res", "-v7", "-nocompression");
    end
end

if p.load
    fld = "results";
        if~isempty(dir("*code"))
            fld = fullfile("code", "results");
        end
   tmp = load(fullfile(fld, "exampleIdCells.mat"));
   res = tmp.res
end
end

