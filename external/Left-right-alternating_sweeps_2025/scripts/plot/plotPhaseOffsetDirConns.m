function [res, recs] = plotPhaseOffsetDirConns(recs, combs, doplot)
% Plots prefdir of pre and phase offset dir for given cell type pairings
ncombs = numel(combs);
nrecs = numel(recs);

for comb = combs
    predir = [];
    postdir = [];
    phasedist = [];
    predir = nan;
    postdir = nan;
    phasedist = nan;
    phaseoffset = [nan, nan];
    nrats = 0;
    for r = 1:nrecs
        if isfield(recs(r).tuning.gridphase, comb)
            conninds = recs(r).tuning.gridphase.(comb).linearInds;
            valid = find(recs(r).isconn(conninds));
%             valid = 1:numel(conninds);
            predir = [predir; recs(r).tuning.gridphase.(comb).predir(valid)];
            postdir = [postdir; recs(r).tuning.gridphase.(comb).phaseoffsetdir(valid)];
            shiftdist = recs(r).tuning.gridphase.(comb).phasedist(valid)*2.5; %2.5 is bin size? yes
            shiftdir = recs(r).tuning.gridphase.(comb).phaseoffsetdir(valid);
            % spacing = recs(r).tuning.gridphase.(comb).spacing(valid)*100;
%             phasedist = [phasedist; shiftdist./spacing];
            phasedist = [phasedist; shiftdist];
            if ~isempty(shiftdir)
                clear offset
                % [offset(:, 1), offset(:, 2)] = pol2cart(shiftdir, shiftdist./spacing);
                [offset(:, 1), offset(:, 2)] = pol2cart(shiftdir, shiftdist);
                phaseoffset = [phaseoffset; offset];
            end
            nrats = nrats+1;
            recs(r).tuning.gridphase.(comb).linearIndsValid = valid;
            disp(r)
        end
    end
    postdir = wrapToPi(postdir-pi);

    res.(comb).predir = double(predir(2:end));
    res.(comb).postdir = double(postdir(2:end));
    res.(comb).phasedist = phasedist(2:end);
    res.(comb).phaseoffset = phaseoffset(2:end, :);
    
    if numel(predir)>1
        [res.(comb).muoffset, res.(comb).rho, res.(comb).pval, res.(comb).muabsoffset, stats] = computeDirAlignment(double(predir),double(postdir));
        res.(comb).mvl = stats.mvl;
        res.(comb).ray_p = stats.rayleigh;
        res.(comb).std = rad2deg(stats.std);
        res.(comb).sem = rad2deg(stats.sem);
        res.(comb).muoffset = rad2deg(res.(comb).muoffset);
        v = ~isnan(predir) & ~isnan(postdir);
        res.(comb).offsets = circ_dist(predir(v), postdir(v));
        res.(comb).muabsoffset = rad2deg(res.(comb).muabsoffset);
    end

    v = ~isnan(predir) & ~isnan(postdir);
    [rho, pval] = circ_corrcc_no_rotation(predir(v), postdir(v));
    % [rho, pval] = circ_corrcc(predir(v), postdir(v));
    npairs = numel(predir);
    predir = predir+[0,2*pi, 2*pi,0];
    postdir = postdir+[0,0, 2*pi,2*pi];
    predir = predir(:);
    postdir = postdir(:);
    cols = repmat(mapcolors(phasedist(:), [1, 40], 'turbo'), 4, 1);
    if doplot
        nexttile
        xlim([-pi,3*pi]), ylim([-pi, 3*pi]);
        d = diagline('r');
        plot(d.XData+2*pi, d.YData, 'r');
        plot(d.XData, d.YData+2*pi, 'r');

        scatter(predir, postdir, 3,'k');
%         scatter(predir, postdir, 3,cols);
        
        title(comb)
        xlabel("Prefdir. Conj")
        ylabel("Grid phase offset dir")
        text (2*pi, 4.6*pi, ...
            sprintf("%.0f pairs from %.0frats\nr:%.2f, p:%.4f", ...
            npairs, nrats, rho, pval))
        xticks(-pi:2*pi:3*pi)
        yticks(-pi:2*pi:3*pi)
        xticklabels(0:360:720);
        yticklabels(0:360:720);
        ax = gca; ax.FontSize = 11;
        axis square
    end

end