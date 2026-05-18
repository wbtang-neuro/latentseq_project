%% Run/load PV pos decoding
sweepsSetup
S.dataRoot = '/Users/wt248/Downloads/Left-right-alternating_sweeps_2025/sample_data/';
% S.dataRoot = S.dataRoot_ % Set data directory
decs = runPvPosDecoding(process=1, save=1);

%% Append useful info
recnames = S.recs_of_mec;
datafld = fullfile(S.dataRoot_, "navigation", "of");
for r = 1:numel(recnames)
    rec = recnames(r);
    fname = fullfile(datafld, rec+".mat");
    tmp = load(fname);
    D = tmp.Dsession;
    D.nu = numel(D.units.mec);
    D.ngrid = sum([D.units.mec.isGrid]);
    D = rmfield(D, "units");
    D.thetaChunks = decs(r).chk;
    D.dec = decs(r);
    recs(r) = D;
    disp(rec)
end
clear decs
%% Compute stats
clf
clear nopposite egohists
egohists = struct("all", [],"right", [],"left", []);
nrecs = numel(recs);
for r = 1:nrecs
    dec = recs(r).dec;
    chk = dec.chk;
    sweeps = dec.sweeps;
    sd = [sweeps.hpfPosDirection]'; 
    hd = chk.hd;
    egosd = circ_dist(sd, hd);
    vswp = chk.speed>S.minSpeed & [sweeps.nvalid]'>3 & [sweeps.straight]'>.5;
    egosd(~vswp) = nan;
    
    [isRight,isLeft, prevRight, prevLeft] = egoRightLeft(egosd);
    
    egohists.right(:,r) = histcounts(egosd(prevLeft), linspace(-pi, pi, 101));
    egohists.left(:,r) = histcounts(egosd(prevRight), linspace(-pi, pi, 101));

    [~, acorrs(r, :)] = circAlternationAcorrAdjacent(circ_diff(egosd), 7);
    [recs(r).prcAltern, recs(r).pAltern] = computeAlternationPercent(egosd);
    [recs(r).modes, lefthist(r, :), righthist(r, :)] = computeModesKsd(egosd);
    [recs(r).hdoffset, recs(r).hdcorr, recs(r).hdpval, recs(r).absHdoffset] = computeDirAlignment(sd, hd, 1);
    recs(r).vswp = vswp;
    recs(r).prcswp = 100*sum(vswp)./sum(chk.speed>S.minSpeed);
    recs(r).mulen = mean([dec.sweeps(vswp).length]);
    if r == 1
        egodir = egosd;
    end
end

%% Print stats
fprintf("Sweep length: %.2f, %.2f\n", mean([recs.mulen]*100), std([recs.mulen]*100)./sqrt(numel(recs)))
fprintf("Alternation: %.2f, %.2f\n", 100*mean([recs.prcAltern]), 100*std([recs.prcAltern])./sqrt(numel(recs)));
modes = rad2deg(circ_mean(cat(1, recs.modes)));
fprintf("Modes: %.2f, %.2f\n", modes(1), modes(2))

%% c plot
figure("WindowStyle","normal")
tl = tiledlayout(1,2, 'TileIndexing','rowmajor');
plotEgoHistPolar(tl = tl, egodir=egodir,alpha=1);
scatterEgodirAcorrs(acorrs = acorrs, plotHalf=0, ax=nexttile());ylim([-1,1])
set(gca, "FontSize", 10)

%% Plot averaged sweeps for all animals
clf
clear cols;
fds = ["right", "left"];
cols.right = [1,0,0, .4];
cols.left = [0,0,1, .4];
clf

ax = nexttile; 
xline(0); yline(0);
xlim([-.2,.2]);
ylim([-.1,.3]);
axis square off
plot([.1,.1], [0,.1]+.05, 'k','LineWidth',2)
text(.11, 0.1, "10 cm")

title(ax, "Averaged sweeps")
nrecs = numel(recs);
for r = 1:nrecs
    D = recs(r);
    dec = D.dec;
    [tmp] = plotAvgSweeps(D,dec);
    for fd = fds   
        plot(ax, tmp.(fd)(:, 1), tmp.(fd)(:, 2), '-',Color=cols.(fd))
        scatter(ax, tmp.(fd)(end, 1), tmp.(fd)(end, 2), 20,cols.(fd)(1:3), 'MarkerFaceAlpha',.7)
    end
end

%% Plot sweep dir modes for top 10 animals with highest fractions of swps
[sorted, isort] = sort([recs.prcswp], 'descend');
tl = tiledlayout(1,10);
tbl = table();
for i = 1:10
    r = isort(i);
    modes = recs(r).modes;
    gve = linspace(-pi,pi, 101);
    pax = polaraxes(tl);
    pax.Layout.Tile = i;
    [pax]=plotEgoHistPolar(ax=pax, egohist=egohists.left(:, r), color = [1,0,0], gve=gve);
    [pax]=plotEgoHistPolar(ax=pax, egohist=egohists.right(:, r), color = [0,0,1], gve=gve);
    pax.ThetaTick=[];
    polarplot(pax, [0, 0], [-.2,1], 'k')
    polarplot(pax, [-pi/2, pi/2], [.6,.6], 'k')
    polarplot(pax, [modes(1),0,0,modes(2)], [1,0,0,1], 'k')
    modes = rad2deg(modes);
    title(sprintf("%s, %.1f%s", S.recs_of_mec(r), recs(r).prcswp, "%"))
end

%% Plot example trajectory with sweeps
D = recs(1);
trng = [9910.4, 9913.6];
vt = restrictq(D.t, trng);
inds = find(vt);
pos = [D.x, D.y];
possm = D.dec.possm;
sweeps = D.dec.sweeps;
vswp = restrictq([sweeps.tStart], trng+[0, 0]);
sweeps = sweeps(vswp);

clf;
plot([0,1.5,1.5,0,0]-.75, [0,0,1.5,1.5,0]-.75, 'k')
ptracking=plot(pos(vt, 1), pos(vt, 2), 'Color',[1,.5,.5]*.3, 'LineWidth',.1);
plpass=plot(possm(vt, 1), possm(vt, 2), 'Color',[.5,1,.5]*.3, 'LineWidth',.1, 'LineStyle','--');
decpos = gsmooth(D.dec.decpos, 0.8);
for s = 1:(numel(sweeps)-1)
    swpinds = sweeps(s).iStart:(sweeps(s).iStop);
    swp = decpos(swpinds, :);
    if rem(s, 2)==0
        col = S.col_cyc_odd;
    else
        col = S.col_cyc_even;
    end
    hswp(s) = plot(swp(:, 1), swp(:, 2), 'Color', col);
    plot(swp(end, 1), swp(end, 2), '.', 'markerSize', 5, 'Color', col)
end
hstart = scatter(pos(inds(1), 1), pos(inds(1), 2), 20, 'r');
xlim([-1,1]); ylim([-1,1])
axis square off
ax = gca;
ax.YDir = "reverse";
legend([hswp(1:2), ptracking, plpass, hstart], ["sweep odd", "sweep even", "tracked pos", "low-pass dec", "start pos"], "Location","eastoutside")

%% Sweep expression correlates with cell yield and spiking activity
clf
tl = tiledlayout(1, 3, 'TileSpacing','compact', 'Padding','compact');
minspd = S.minSpeed;
fds = ["nu", "ngrid"];
for fd = fds
    nexttile
    scatter([recs.(fd)], [recs.prcswp], 10, 'k')
    r = refline; r.Color = 'r';
    [rho, pval] = corr([recs.(fd)]', [recs.prcswp]');
    ax = gca;
    xmx = ax.XLim(2);
    xmn = ax.XLim(1);
    ylim([0, 100])
    yticks([0:25:100])
    axis square;
    ax.FontSize=12;
    title(fd)
    xlabel(fd)
    ylabel("% Sweeps")
    text(xmn, 90, sprintf(" r=%.2f\n p=%.3f", rho, pval), 'HorizontalAlignment','left')
end

cols = gray(26).*1.2;
cols(cols>1)=1;
nexttile
for r = 1:nrecs
    sweeps = recs(r).dec.sweeps;
    nspk = [sweeps.nSpikesTotal]';
    vswp = cat(1, recs(r).vswp);
    speed = interp1(recs(r).t, recs(r).speed, recs(r).thetaChunks.tStart, 'linear');
    potswp = (speed>minspd);
    edgs = [0,prctile(nspk, 10:10:99)];
    nspkinds = discretize(nspk, edgs);
    bins = movmean(edgs, [0, 1]);
    bins = bins(1:end-1);
    clear swprate
    for i = 1:numel(bins)
        idx = nspkinds==i;
        swprate(i) = 100.*sum(vswp(idx))./sum(potswp(idx));
    end
    [rhos(r), pvals(r)] = corr(bins(:), swprate(:));
    plot(bins, swprate, '.-', 'Color', cols(r, :))

end
fprintf("\nMean corr: %.2f, %.5f\n", mean(rhos), std(rhos)./sqrt(nrecs));
xlabel("N spikes in theta cycle");
ylabel("% sweeps");
ax = gca; ax.FontSize = 12;
axis square
yticks(0:25:100)
xticks(0:100:600)

%% -----------  Single module decoding --------------------
sweepsSetup
% S.dataRoot = S.dataRoot_ % Set data directory
% decs = runPvPosDecodingModules(process=1, save=1, recs=S.recs_of_mec(:));

res = runPvPosDecodingModules(process=0, save=0, load=1, recs=S.recs_of_mec(:));

%% Plot ego-ego plots for between modules 1 vs 2 and 2 vs 3
r = 1;
dec = res(r);
mods = res(r).mods(1:3);
chk = res(r).chk;
% S.figure;
tl = tiledlayout(1,2,'TileSpacing','compact');
dirName = "hpfPosDirection";
for ipair = 1:2
    clear sd;
    for s = 1:2
        mod = mods(s+ipair-1);
        sweeps = mod.sweeps;
        vswp = [sweeps.nvalid]'>3 & [sweeps.straight]'>.5 & chk.speed>.15;
        tmp=[sweeps.(dirName)]';
        tmp(~vswp, 1)=nan;
        hd = chk.hd;      
        egosd = circ_dist(tmp, hd);
        sd(:, s) = egosd;
    end
    egoego_heatmap("ax",nexttile, "computeStats",true, "x",sd(:, 1), "y", sd(:, 2), "lims",70, "nsmooth",.7);
    title(sprintf("Module %u vs Module %u", ipair, ipair+1))
end
%% Plot histogram of sweep lengths for example rat
tl = tiledlayout(2,1,'TileSpacing','compact');
nexttile
r = 1;
dec = res(r);
recName = res(r).recName;
chk = res(r).chk;
cols = ["r", "g", "b", "k"];
i = 1;
for m = 1:3
        mod = dec.mods(m);
        mod.dec = dec;
        mod.chk = chk;

        sweeps = mod.sweeps;
        vswp = [sweeps.nvalid]'>3 & [sweeps.straight]'>.5 & chk.speed>.15;
        sweeplengths(i) = mean([sweeps(vswp).hpfLength], 'omitnan');
        histogram([sweeps(vswp).length].*mod.gridSpacingMeters, linspace(0,1,100), 'FaceColor',cols(m), 'FaceAlpha',.5)
        i = i+1;
end
ax = gca;
ax.FontSize = 12;
ylim([0, 2000])
yticks([0, 2000]);
xlim([0, .8])
xticks(0:.2:.8)
xlabel("Sweep length (m)")
ylabel("N sweeps")
%% sweep length vs spacing 
nrecs = numel(res);
i = 1;
clear sweeplengths modspacings cols nu recidx
for r = 1:nrecs
    dec = res(r);
    recName = res(r).recName;
    chk = res(r).chk;
    nmods = numel(dec.mods);
    
    for m = 1:nmods
        mod = dec.mods(m);
        mod.dec = dec;
        mod.chk = chk;

        sweeps = mod.sweeps;
        vswp = [sweeps.nvalid]'>3 & [sweeps.straight]'>.5 & chk.speed>.15;
        sweeplengths(i) = median([sweeps(vswp).length], 'omitnan').*mod.gridSpacingMeters;
        maxlen(i) = prctile([sweeps(vswp).length],99);
        modspacings(i) = mod.gridSpacingMeters;
        nu(i) =mod.nUnits;
        nOverHalf(i) = sum([sweeps(vswp).length]>.5*modspacings(i))./sum(vswp);
        
        recinds(i) = r;
        i = i+1;
    end
    
end
% Plot length figure

idx = nu>40;
x = double(modspacings(idx));
y = double(sweeplengths(idx));
[rho,pval] = corrcoef(x,y);
dlm = fitlm(x,y,'Intercept',false);
nexttile;
cols = distinguishable_colors(nrecs);
a = dlm.Coefficients.Estimate;
plot([0, 2], a(1)*[0,2], 'k')
scatter(modspacings(idx), sweeplengths(idx), 10, cols(recinds(idx), :))
%
xlabel("Module spacing (m)")
ylabel("Sweep length (m)")
ax = gca; ax.FontSize = 12;

xlim(minmax(modspacings(idx)))
xlim([0,max(modspacings(idx))])
ylim([0,.4])

%% Print length stats
normlen = sweeplengths(idx)./modspacings(idx);
mulen = 100*mean(normlen);
semlen = 100*std(normlen)./sqrt(sum(idx));
fprintf("\nSweep lengths: %.1f + %.1f%s of mod spacing\n", mulen, semlen, "%")

muoverhalf = 100*mean(nOverHalf(idx));
semoverhalf = 100*std(nOverHalf(idx))./sqrt(sum(idx));
fprintf("N over half: %.1f + %.1f\n", muoverhalf, semoverhalf)

%% Sweep modes and alternation
clf
i = 1;
clear modes flickering
for r = 1:nrecs
    dec = res(r);
    recName = res(r).recName;
    chk = res(r).chk;
    nmods = numel(dec.mods);

    for m = 1:nmods
        mod = dec.mods(m);
        if mod.nUnits<40;continue;end        
            mod.dec = dec;
            mod.chk = chk;
            sweeps = mod.sweeps;
    
            vswp = [sweeps.nvalid]'>3 & [sweeps.straight]'>.5 & chk.speed>.15;
            sd = [sweeps.hpfPosDirection];
            sd(~vswp)=nan;
            hd = chk.hd;
            egosd = circ_dist(sd(:), hd(:));
            sd = sd(:);
    
            flickering(i) = computeAlternationPercent(egosd);     
            [modes(i, :),lefthist(i, :), righthist(i, :)] = computeModesKsd(egosd);
            
            %shuffled alternation
            egosdnonan = egosd(~isnan(egosd));
            n = numel(egosdnonan);
            prcAltern = zeros(n, 1);
            for s = 1:1000
                prcAltern(s) = computeAlternationPercent(egosdnonan(randperm(n)));
            end
            nhigher(i) = sum(prcAltern>flickering(i));
            i = i+1;
            disp(i)
    end
end
% %%
%%
muflicker = mean(flickering*100);
semflicker = std(flickering*100)./sqrt(numel(flickering));
fprintf("\nLeft-right-alternation: %.1f, %.1f\n", muflicker, semflicker)
fprintf("Greater than 99th percentile of shuffled values in %.1u of %.1u modules\n", sum(nhigher<10), numel(nhigher))
mumodes = rad2deg(circ_mean(modes));
fprintf("Modes: %.2f, %.2f\n", mumodes(1), mumodes(2))

%% Compute alignment
i = 1;
clf
dirName = "hpfPosDirection";
pb = ProgressBar;
clear pvals rhos coflickering pCoflicker muoffset
for r = 1:nrecs
    dec = res(r);
    recName = res(r).recName;
    chk = res(r).chk;
    vmods = find([dec.mods.nUnits]>40);

    if numel(vmods)>1
    combs = nchoosek(vmods, 2);
    for c = 1:size(combs, 1)
        clear sd sdallo
        for s = 1:2
            mod = dec.mods(combs(c,s));
            sweeps = mod.sweeps;
            vswp = [sweeps.nvalid]'>3 & [sweeps.straight]'>.5 & chk.speed>.15;
            
            tmp=[sweeps.(dirName)]';
            tmp(~vswp, 1)=nan;
            hd = chk.hd;
    
            egosd = circ_dist(tmp, hd);
            sd(:, s) = egosd;
            sdallo(:, s) = tmp;
            
        end
        [coflickering(i, 1), pCoflicker(i, 1)] = computeCoAlternation(sd);
        invalid = isnan(sum(sd, 2));
        [rhos(i, 1), pvals(i, 1)] = circ_corrcc_no_rotation(sdallo(~invalid, 1), sdallo(~invalid, 2));
        recNames(i) = recName;
         i = i+1;
    end
    end
    % iterate through pairs of mods
    pb.update(r/nrecs)
end
%% print
idx = true(size(coflickering));
fprintf("\nCoflickering: %.2f, %.3f, max p=%.4f\n", 100*mean(coflickering(idx)), 100*std(coflickering(idx))./sqrt(sum(idx)), max(pCoflicker))
fprintf("Correlation: %.2f, %.3f, max p=%.4f\n", mean(rhos(idx)), std(rhos(idx))./sqrt(sum(idx)), max(pvals))

%%
coflickering(isnan(coflickering))=[];
swarmchart(nexttile, coflickering*0, coflickering*100, 6,[1,1,1]*.7)
ylim([0,100])
yline(50)
xlim([-1,1])
errorbar(0, mean(coflickering)*100, std(coflickering)*100, 'k')
scatter(0, mean(coflickering)*100, 15, 'red', 'MarkerEdgeColor', 'k')
xticks([])
yticks([0,50,100])
ylabel("% aligned sweeps")

%% Plot averaged sweeps for all mods
clear cols;
fds = ["right", "left"];
cols.right = [1,0,0, .4];
cols.left = [0,0,1, .4];
clf

ax = nexttile; 
xline(0); yline(0);
xlim([-.2,.2]);
ylim([-.1,.3]);
axis square off
plot([.1,.1]+.1, [0,.25]+.02, 'k','LineWidth',2)
text(.11, 0.1, "1/4 spacing")

title(ax, "Sweeps fom grid cells")
%
% tl = tiledlayout('flow')

nrecs = numel(res);
for r = 1:nrecs
    D = res(r);
    dec = D;
    recName = res(r).recName;
    chk = res(r).chk;
    chk.egoid = circ_dist(chk.id, chk.hd);
    D.thetaChunks = chk;
    nmods = numel(dec.mods);
    for m = 1:nmods
        mod = dec.mods(m);
        if mod.nUnits<40;continue;end
        [tmp] = plotAvgSweepsModules(D,mod);   
        
        for fd = fds
      
            % plot(ax, tmp.(fd)(:, 1), tmp.(fd)(:, 2), Color=cols.(fd))
            % scatter(ax, tmp.(fd)(end, 1), tmp.(fd)(end, 2), 20,cols.(fd)(1:3), 'MarkerFaceAlpha',.7)
            col = mapcolors(mod.gridSpacingMeters, [.5,1.4], "turbo");
            plot(ax, tmp.(fd)(:, 1), tmp.(fd)(:, 2), Color=col)
            scatter(ax, tmp.(fd)(end, 1), tmp.(fd)(end, 2), 20,col, 'MarkerFaceAlpha',.7)
        
        end
    end
end
%% ------------------ Place cell decoding ----------------------------
sweepsSetup
% S.dataRoot = S.dataRoot_ % Set data directory

%% Append useful info
recnames = S.recs_of_mec_hc(:);
res_hc = runPvPosDecoding(process=1, save=1, recs=recnames, brainregion="hc");
res_mec = runPvPosDecoding(process=1, save=1, recs=recnames, brainregion="mec");
datafld = fullfile(S.dataRoot_, "navigation", "of");
%%
for r = 1:numel(recnames)
    rec = recnames(r);
    fname = fullfile(datafld, rec+".mat");
    tmp = load(fname);
    D = tmp.Dsession;
    D = rmfield(D, "units");
    D.thetaChunks = res_mec(r).chk;
    D.dec.hc = res_hc(r);
    D.dec.mec = res_mec(r);
    D.sweeps.hc = res_hc(r).sweeps;
    D.sweeps.mec = res_mec(r).sweeps;
    D.recName = recnames(r);
    recs(r) = D;
    disp(rec)
end
clear decs
%% Example sweeps from MEC and HC
D = recs([recs.recName]=="28258_4");
brs = ["mec", "hc"];
trng=[0,5]+113;

[vt, t] = restrictq(D.t, trng);
clf;
inds = find(vt);
pos = [D.x, D.y];
vspeed = D.speed>.15;

clear col swps
col.hc = [107,205,127]./255;
col.mec = [44,127,184]./255;

dec = D.dec.hc;
sweeps = D.sweeps.hc;
vswp = restrictq([sweeps.tStart], trng+[0, 0]);
allsweeps = sweeps;
sweeps = sweeps(vswp);
swps.mec = D.sweeps.mec(vswp);
swps.hc = D.sweeps.hc(vswp);

tl = tiledlayout(1,5, 'TileSpacing','none');
pos = rotate2d(pos, -pi/2);
for s = 7:11
    nexttile
    plot(([0,1,1,0,0]-.5), ([0,0,1,1,0]-.5), 'k')
    posinds = sweeps(7).iStart+(0:300);
    
    plot(pos(posinds, 1), pos(posinds, 2), 'Color',[1,1,1]*0.3, 'LineWidth',.1)
    ratPatch("color",[.7,.3,.3], "edgeColor",'k', 'orientation',D.hd(sweeps(s).iStart)-pi/2,'position',pos(sweeps(s).iStart, :), 'sizeMeters',.08);% 

    for br = brs
        dec = D.dec.(br);
        decpos = rotate2d(gsmooth(dec.decpos, 0.5), -pi/2);
        swpinds = swps.(br)(s).iStart:(swps.(br)(s).iStop);
        swp = decpos(swpinds, :);
        plot(swp(:, 1), swp(:, 2), 'Color', col.(br))
        plot(swp(end, 1), swp(end, 2), '.', 'markerSize', 8, 'Color', col.(br))
    end
    xlim([-1,1]); ylim([-1,1])
    xlim([-1,1]*.5); ylim([-1,1]*.5)
    axis square off
    ax = gca;
    ax.YDir = "reverse";
end
%% Check co-alternation
brs = ["mec", "hc"];
clf
refDir = "hd";
dirType = "hpfPosDirection";
clear egosd coflickering nvalid muoffset rhos pvals muabsoffset
tiledlayout(1,2);
for r = 1:numel(recs)
    D = recs(r);
    chk = D.thetaChunks;
    speed = chk.speed;
    for br = brs
        sweeps = D.sweeps.(br);
        vswp = speed>.15 & [sweeps.nvalid]'>3 & [sweeps.straight]'>.5;
        recs(r).dec.(br).vswp = vswp;
        sd = [sweeps.(dirType)]';
        sd(~vswp)=nan;
        hd = chk.hd;

        tmp = circ_dist(sd, hd);
        egosd.(br) = circ_dist(sd, hd);
    end
    invalid = (isnan(egosd.mec)|isnan(egosd.hc));
    nvalid(r) =sum(~invalid);
    [coflickering(r), pCoflicker(r)] = computeCoAlternation([egosd.mec, egosd.hc]);
    coflickering(r) = coflickering(r)*100;
    [rhos(r), pvals(r)] = circ_corrcc(egosd.mec(~invalid), egosd.hc(~invalid));
    sdoffset = circ_dist(egosd.hc(~invalid), egosd.mec(~invalid));
    muoffset(r, 1) = circ_mean(sdoffset);
    if r == 4
        egoego_heatmap("ax",nexttile, "computeStats",true, "x",egosd.mec, "y",egosd.hc, "gve",linspace(-pi, pi, 70), 'lims', 90, "nsmooth",1);
        xlabel("MEC sweep dir"); ylabel("HC sweep dir");
        title("Head-centered sweep directions (MEC vs HC)")
    end

end
idx = nvalid>0;
fprintf("\n\nCoflickering: %.2f, %.3f, max p=%.4f", mean(coflickering(idx)), std(coflickering(idx))./sqrt(sum(idx)), max(pCoflicker))
fprintf("\nCorrelation: %.2f, %.3f", mean(rhos(idx)), std(rhos(idx))./sqrt(sum(idx)))
fprintf("\nAlignment: %.2f, %.3f\n", rad2deg(mean(abs(muoffset(idx)))), rad2deg(circ_std(muoffset(idx))./sqrt(sum(idx))))

%
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
ylabel("% aligned sweeps")
pbaspect([1,3,1])

%% Check temporal delay
clf
brs = ["mec", "hc"];
cols = pink(16);
clear pklag decpos

for i = 1:numel(recs)
    D = recs(i);
    for br = brs
        poshpf = D.dec.(br).decpos-D.dec.(br).possm;
        posr = rotate2d(poshpf, -D.hd);
        decpos.(br) = posr(:, 1);
    end
    idx = ~isnan(decpos.mec) & ~isnan(decpos.hc) & D.speed>.15;
    [r, lags] = xcorr(decpos.hc(idx), decpos.mec(idx), 8, "coef");
    lags = lags.*10;
    lagsfine = linspace(lags(1), lags(end), 100);
    r = interp1(lags, r, lagsfine, 'cubic');
    r = r./max(r);
    plot(lagsfine, r, "Color", cols(i, :))
    [~,imx] = max(r);
    pklag(i) = lagsfine(imx);
end
xline(0)
ax = gca;
xlim([-120,120]*.5)
ylim([-.5,1.25])
xticks([-120,0,120]*.5)
yticks([-.5, 1])
mulag = mean(pklag);
stdlag = std(pklag);
errorbar(mulag,1.2,stdlag,'horizontal', 'Color', 'k')
scatter(mulag, 1.2, 25, 'r', 'MarkerEdgeColor','k')
ax.Clipping="off";
xlabel('Time lag (ms)')
ylabel("Correlation (norm.)")
ax.FontSize = 12;
%% Theta triggered sweeps all animals for ED (KEEP)
clear col
col.hc = [107,205,127]./255;
col.mec = [44,127,184]./255;
clf
tl = tiledlayout(2, 1, 'TileSpacing','compact');

for r = 1:numel(recs)
%     nexttile
    D = recs(r);
    D.chk = D.thetaChunks;
    [swpt, swpx] = triggedSweepDelay2(D,"mec");
    plot(nexttile(1), swpt, swpx, 'Color', col.mec)
    scatter(nexttile(1), swpt(end), swpx(end), 10, col.mec)
    [swpt, swpx] = triggedSweepDelay2(D,"hc");
    plot(nexttile(2), swpt, swpx, 'Color', col.hc)
    scatter(nexttile(2), swpt(end), swpx(end), 10, col.hc)
end
for a = 1:2
    ax = nexttile(a);
    xline(0, 'r')
    ax.XLim = [-20, 125];
    ax.YLim = [-.1, .125];
    ax.YLim = [-.05, .125];
    ax.FontSize = 12;
    ax.FontName = "myriad pro";
    xlabel(ax, "Time from theta trough (ms)")
    ylabel(ax, "Decoding offset (m)")
    yticks(ax, [-.1, 0, .1])
    yline(0, 'Color', [1,1,1]*.5)
end
title(nexttile(1), "MEC-parasubiculum")
title(nexttile(2), "Hippocampus")
%% Theta triggered sweeps for example animal (Keep)
clf
clear col
col.hc = [107,205,127]./255;
col.mec = [44,127,184]./255;
clf
nrecs = numel(recs);
nrecs = 1;
tl = tiledlayout(1,nrecs, 'TileSpacing','compact');

for r = 1:nrecs
%     nexttile
    D = recs(r);
    D.chk = D.thetaChunks;
    [swpt, swpx, iends] = triggedSweepDelay2(D,"mec", 1);
    hplt(1)=plot(nexttile(r), swpt, swpx, 'Color', col.mec);
    scatter(nexttile(r), swpt(iends), swpx(iends), 10, col.mec)
    [swpt, swpx, iends] = triggedSweepDelay2(D,"hc", 1);
    hplt(2)=plot(nexttile(r), swpt, swpx, 'Color', col.hc);
    scatter(nexttile(r), swpt(iends), swpx(iends), 10, col.hc)
end
for a = 1:tilenum(gca)
    ax = nexttile(a);
    xline(0, 'r')
    ax.XLim = [-250, 250];
    ax.YLim = [-.1, .125];
    ax.YLim = [-.05, .125];
    ax.FontSize = 12;
    ax.FontName = "myriad pro";
    
    yticks(ax, [-.1, 0, .1])
    yline(0, 'Color', [1,1,1]*.5)
end
xlabel(nexttile(1), "Time from theta trough (ms)")
ylabel(nexttile(1), "Decoding offset (m)")
legend(hplt, ["MEC-parasubiculum", "Hippocampus"])

