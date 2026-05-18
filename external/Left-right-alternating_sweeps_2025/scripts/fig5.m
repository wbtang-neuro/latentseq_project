%{
NB: Currently not compatible with format of released dataset and may have
code dependencies outside this repository.
%}

%% Run sleep decoding
process = 0;
if process
    run_grid_decoding3
end

%% Process data
processSleepDec

%% Load data
sweepsSetup
recs = findSweepsRecs("mec", "sleep_best");
for r = 1:height(recs)
    recName = getRecName(recs(r, :));
    disp(recName)
    tmp = S.load(S.filepath("sleepData", "withsweeps2", recName+".mat"));
    res(r) = tmp;
    disp("done")
end
recs = res;
%% Sleep stats
res = struct();
clear res epochs epocs
epochs.rem = nan;
epochs.sws = nan;
for r = 1:numel(recs)
    recName = recs(r).recName;
    tmp = recs(r).sleepData;
    % Epoch durations
    epochs.rem = [epochs.rem;diff(tmp.times.rem, 1,2)];
    epochs.sws = [epochs.sws;diff(tmp.times.sws, 1,2)];
    res(r) = tmp;
end
totaltime = [res.totalTime];
sleep = table([recs.recName]', [totaltime.rem]'./60, [totaltime.sws]'./60, ([totaltime.rem]'+[totaltime.sws]')./60,'VariableNames',["Rec", "REM", "SWS", "Total"]);
fprintf("\nTotal sleep: %.2f, %.3f\n", mean(sleep.Total), std(sleep.Total))
fprintf("REM sleep: %.2f, %.2f%s\n", 100*mean(sleep.REM./sleep.Total), 100*std(sleep.REM./sleep.Total), "%")
fprintf("SWS sleep: %.2f, %.2f%s\n", 100*mean(sleep.SWS./sleep.Total), 100*std(sleep.SWS./sleep.Total), "%")
fprintf("REM epoch length: %.0f %.0f sec\n",mean(epochs.rem, 'omitnan'), std(epochs.rem, 'omitnan'))
fprintf("SWS epoch length: %.0f %.0f sec\n",mean(epochs.sws, 'omitnan'), std(epochs.sws, 'omitnan'))

%% Show that activity is rhythmic (for ED 11)
figure
tl = tiledlayout(1,3, "TileSpacing","compact", "TileIndexing","rowmajor");
fd = "id_lmt";
states = ["run", "rem", "sws"];
cols = pink(16);
for r = 1:numel(recs)
    res = recs(r);      
    for s = 1:numel(states)
        state = states(s);
        sc = gsmooth(res.idspk, 2);
        [pks, iloc] = findpeaks(sc);
        trigIdx = intersect(iloc, find(res.states.(state)));
        id = res.(fd)(trigIdx);
        t = res.t(trigIdx);
        [a, b]=xCorrPointProcess(t,t, .01, 101, 1);
        a(51)=0;
        plot(nexttile(s), b(51:end), a(51:end), 'Color', cols(r, :));
        acorrs.(states(s))(:, r) = b(51:end);
    end
end
for a = 1:3
    ax = nexttile(a);
    xlim([0,.5]);
end
%% Write acorr data to sourcedata file
fname = fullfile(S.codeRoot, "figure_scripts", "source_data", "ED11.xlsx");
varnames = "rat "+(1:numel(recs));
varnames = [varnames, "meanAcorr", "stdAcorr"];
varnames = vec(states+"_"+varnames');
data = [];
for s = states
    data = [data, acorrs.(s), mean(acorrs.(s), 2), std(acorrs.(s), [], 2)];
end
tbl = array2table(data, "VariableNames",varnames);
writetable(tbl, fname, 'Sheet', 'c', 'WriteMode', 'overwritesheet')
%% Plot distribution of interflicker angles and autocorrelation of id at ID peaks (Fig 5c)
% fig = S.figure([400, 250], "All alignement hists");
nshuff = 10;
figure
tl = tiledlayout(3,3, "TileSpacing","compact", "TileIndexing","rowmajor");
states = ["run", "rem", "sws"];
clear acorrs
for a = 1:3
    axs(a) = nexttile(a);
    set(gca, 'XLim', [-180,180]);
    yticks([]);
    xticks([-180, 0, 180])
    xline(0)
    title(upper(states(a)))
%     axis square
end

gve = S.gve.angular;
gv = rad2deg(S.gv.angular);
pb = ProgressBar;

clear dirstats
%
for r = 1:numel(recs)
    res = recs(r);
    disp(res.recName)
    for fd = "id_corr"       
        for s = 1:numel(states)        
            state = states(s);
            sc = gsmooth(res.idspk, 2);
            [pks, iloc] = findpeaks(sc);
            trigIdx = intersect(iloc, find(res.states.(state)));

            id = res.(fd)(trigIdx);
            t = res.t(trigIdx);

            tdiff = diff(t);
            id_diff = circ_diff(id);
            id_diff(tdiff>.25|tdiff<.02) = nan;
            
            % compute alternation
            diffsign = sign(id_diff);
            diffdiff = diff(diffsign);
            isalternation = abs(diffdiff)==2;
            npossible = sum(~isnan(diffdiff));
            flickering = sum(isalternation)./npossible;
            dirstats.(state).alternation(r, 1) = 100.*sum(isalternation)./npossible;
            [dirstats.(state).modes(r, :),concentration,significance] = computeModesInterangle(id_diff, 1);
            
            % Shuffled alternation
            egoid = id_diff;
            egoidnonan = egoid(~isnan(egoid));
            n = numel(egoidnonan);
            prcAltern = zeros(1000, 1);
            modes = zeros(1000, 2);
            disp("shuffling")
            tic
            for i = 1:nshuff
                  egoidrand = egoidnonan(randperm(n));
                  prcAltern(i) = computeAlternationPercent(egoidrand, 1);
                  modes(i, :) = computeModesInterangle(egoidrand, 1);
            end
            dirstats.(state).nhigher(r) = sum(prcAltern>flickering);
            dirstats.(state).maxshuff(r) = max(prcAltern);
            dirstats.(state).alternationShuffle(r) = mean(prcAltern);
            dirstats.(state).modesShuffle(r, :) = circ_mean(modes);
            toc
            
            %
            diffprev = circshift(id_diff, 1);
            h = histcounts(id_diff(diffprev>0), gve);
            h = gsmooth(h, 1);
            plot(nexttile(s), gv, h./sum(h), 'b');
            h = histcounts(id_diff(diffprev<0), gve);
            h = gsmooth(h, 1);
            plot(nexttile(s), gv, h./sum(h), 'r');
%             
            [lags, acorrs.(state)(:, r)] = circAlternationAcorrAdjacent(id_diff, 7);
        
        end
    end
    pb.update(r/numel(recs));
end

% plot acorrs
for a = 4:6
    axs(a) = nexttile(a);
    scatterEgodirAcorrs(acorrs = acorrs.(states(a-3))', plotHalf=0, ax=axs(a));ylim([-1,1])
    xlim([-5.5,5.5])
%     set(gca, "FontSize", 10)
    if a>4
        ylabel([])
    end
end

%% Print stats
% REM
modes = mean((dirstats.rem.modes.*[-1,1]), 2);
fprintf("\n REM alternation %.2f %.2f, Interangle %.2f %.2f", mean(dirstats.rem.alternation), std(dirstats.rem.alternation)./sqrt(numel(recs)), rad2deg(mean(modes)), rad2deg(std(modes)./sqrt(numel(recs))))
% REM shuffle
modes = mean((dirstats.rem.modesShuffle.*[-1,1]), 2);
fprintf("\n Shuffle REM alternation %.2f %.2f, Interangle %.2f %.2f\n", 100*mean(dirstats.rem.alternationShuffle), 100*std(dirstats.rem.alternationShuffle)./sqrt(numel(recs)), rad2deg(mean(modes)), rad2deg(std(modes)./sqrt(numel(recs))))

% SWS
modes = mean((dirstats.sws.modes.*[-1,1]), 2);
fprintf("\n SWS alternation %.2f %.2f, Interangle %.2f %.2f", mean(dirstats.sws.alternation), std(dirstats.sws.alternation)./sqrt(numel(recs)), rad2deg(mean(modes)), rad2deg(std(modes)./sqrt(numel(recs))))
% SWS shuffle
modes = mean((dirstats.sws.modesShuffle.*[-1,1]), 2); 
fprintf("\n Shuffle SWS alternation %.2f %.2f, Interangle %.2f %.2f\n", 100*mean(dirstats.sws.alternationShuffle), max(dirstats.sws.nhigher), rad2deg(mean(modes)), rad2deg(std(modes)./sqrt(numel(recs))))


%% Write acorr data to sourcedata file
fname = fullfile(S.codeRoot, "figure_scripts", "source_data", "fig5.xlsx");
varnames = "rat "+(1:numel(recs));
varnames = [varnames, "meanAcorr", "stdAcorr"];
varnames = vec(states+"_"+varnames');
data = [];
for s = states
    data = [data, acorrs.(s), mean(acorrs.(s), 2), std(acorrs.(s), [], 2)];
end
tbl = array2table(data, "VariableNames",varnames);
writetable(tbl, fname, 'Sheet', 'c', 'WriteMode', 'overwritesheet')

%% Alignment between sweeps and ID across states
i = 0;
nshuff = 1000;
clear stats hall statsshuff
pb = ProgressBar;
for r = 1:numel(recs)
    D = recs(r);
    chk = D.thetaChunks;
    nmods = numel(D.mods);
    sc = gsmooth(D.idspk, 2);
    disp(D.recName)
    
    sc = gsmooth(D.idspk, 2);
    
    for m = 1:nmods
        if~isempty(D.mods{m})
        i = i+1;
        for s = 1:numel(states)
            vt = D.states.(states(s));
            sweeps = D.mods{m}.sweeps;
            trigIdx = intersect(iloc, find(D.states.(state)));
            idx = vt([sweeps.iStart]);
            sweeps = sweeps(vt([sweeps.iStart]));

            
            vswp = [sweeps.nvalid]'>3 & [sweeps.straight]'>.3;
            sd = [sweeps.travelDirection]';
            icen = chk.icen;
            hd = D.id_corr(icen(idx));
            sd(~vswp) = nan;
            egosd=circ_dist(sd, hd);
            h = histcounts(egosd, gve);
            hall.(states(s))(i, :) = h;
            clear tmp
            [tmp.muoffset, tmp.rho, tmp.pval, muabsoffset, rez] = computeDirAlignment(sd,hd, 1);
            tmp.ray = rez.rayleigh;
            tmp.mvl = rez.mvl;
            tmp.length = mean([sweeps(vswp).length]);
            tmp.length_sd = std([sweeps(vswp).length]);
            tmp.straight = mean([sweeps([sweeps.nvalid]>3).straight]);
            tmp.prcswp = 100*sum(vswp./numel(vswp));
            stats.(states(s))(i) = tmp;
            
            valid = ~isnan(sd) & ~isnan(hd);
            hd=hd(valid); sd = sd(valid);
            clear tmpshuff tmp
            imax = numel(sd);
            pbs = ProgressBar();
            for ii =1:nshuff
                sdshuff = circshift(sd, randi(imax));
                [tmp.muoffset(ii), tmp.rho(ii), tmp.pval(ii)] = computeDirAlignment(sdshuff,hd, 1);
                pbs.update(ii/500);
            end
            tmp.muoffset = abs(tmp.muoffset);
            for fd = fieldnamesstr(tmp)
                tmpshuff.(fd) = mean(tmp.(fd), 'omitnan');
            end
            statsshuff.(states(s))(i) = tmpshuff;
        end
        end
    end
    disp('Done')
end
%% Print sweep stats
for s = states
    fprintf("\n%s\n", upper(s));
    fprintf("%s sweeps: %.2f %.2f\n", "%", mean([stats.(s).prcswp]), std([stats.(s).prcswp])./numel(recs))
    len_mu = 100*mean([stats.(s).length]);
    len_sem = 100*std([stats.(s).length])./numel(recs);
    fprintf("Sweep length: %.2f %.3f\n", len_mu, len_sem)

    len_mu = rad2deg(mean(abs([stats.(s).muoffset])));
    len_sem = rad2deg(std(abs([stats.(s).muoffset]))./numel(recs));
    fprintf("Angular offset: %.2f %.2f deg", len_mu, len_sem)

    pvalid = [stats.(s).ray]<.01;
    rhos = [stats.(s).mvl];
    fprintf("\nmvl %.2f, %.2f, p<.01 in %u/%u\n", mean(rhos(pvalid)), std(rhos(pvalid))./sum(pvalid), sum(pvalid), numel(pvalid));

    pvalid = [stats.(s).pval]<.01;
    rhos = [stats.(s).rho];
    fprintf("corr %.2f, %.2f, p<.01 in %u/%u\n", mean(rhos(pvalid)), std(rhos(pvalid))./sum(pvalid), sum(pvalid), numel(pvalid));

    pvalid = [statsshuff.(s).pval]<.01;
    rhos = [statsshuff.(s).rho]';
    fprintf("corr (shuffled) %.2f, %.2f, p<.01 in %u/%u\n", mean(rhos(:), 'omitnan'), std(rhos(:))./sum(pvalid), sum(pvalid), numel(pvalid));
    fprintf("Angular offset: %.2f %.2f deg, mvl %.2f, %.2f\n", len_mu, len_sem)
end
%%
% figure
cols = pink(36);
a = 7;
for s = states
    axs(a) = nexttile(a);
    cla
    h = hall.(s);
    h = gsmooth(h', 1)';
    h = h./sum(h, 2);
    for i = 1:size(h, 1)
        plot(gv, h(i, :), 'Color',cols(i, :))
    end
    set(gca, 'XLim', [-180,180]);
    yticks([]);
    xticks([-180, 0, 180])
    xline(0)
    a = a+1;
end
%% Write data to sourcedata file
fname = fullfile(S.codeRoot, "figure_scripts", "source_data", "fig5.xlsx");

i = 1;
clear varnames
for r = 1:numel(recs)
    D = recs(r);
    nmods = numel(D.mods);
    mm = 1;
    for m = 1:nmods
        if~isempty(D.mods{m})
            varnames(i) = sprintf("_rat%u_mod%u", r, mm);
            i = i+1;
            mm = mm+1;
        end
    end
end
varnames = vec(states+varnames');
data = [hall.run;hall.rem;hall.sws]';
tbl = array2table(data, "VariableNames",varnames);
writetable(tbl, fname, 'Sheet', 'e', 'WriteMode', 'overwritesheet')
%% %%%%%%%% EXAMPLE PLOTS %%%%%%%%%%%%
%% Plot 1000 sweeps for one example mod (roger mod 3)
clf
fig5.singleSweepsSleep(recs(1), 2)
%% Plot averaged sweeps for all mods
figure
clear cols;
fds = ["right", "left"];
cols.right = [1,0,0, .4];
cols.left = [0,0,1, .4];

tl = tiledlayout(1,3, "TileSpacing","compact", "TileIndexing","rowmajor");
states = ["run", "rem", "sws"];
clear acorrs
for a = 1:3
    axs(a) = nexttile(a);
    xline(0); yline(0);
    xlim([-.15,.15]);
    ylim([-.05,.25]);
    axis square off
    if a ==1
        plot([.15,.15, .05], [.15,.25, .25], 'k','LineWidth',2)
    end
    title(upper(states(a)))
end

nrecs = numel(recs);
pb = ProgressBar();
for r = 1:nrecs
    D = recs(r);
    dec = D;
    nmods = numel(dec.mods);
    for m = 1:nmods
        if isempty(dec.mods{m});continue;end
        gmod = dec.mods{m};
        for s = 1:3
            [tmp] = fig5.trigSweepsSleep(D,gmod, D.states.(states(s)));

            
            for fd = fds               
                plot(axs(s), tmp.(fd)(:, 1), tmp.(fd)(:, 2), Color=cols.(fd))
                scatter(axs(s), tmp.(fd)(end, 1), tmp.(fd)(end, 2), 20,cols.(fd)(1:3), 'MarkerFaceAlpha',.7)
            end
        end
    end
    pb.update(r/numel(recs))
end
%% Plot REM example
%plotRemExample
%% Plot SWS example
%plotSwsExample

