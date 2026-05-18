function [res, recs] = plotDirTuningConns(recs, combs, doplot)
% Plots prefdir of pre and postsynaptic cell for given cell type pairings
ncombs = numel(combs);
nrecs = numel(recs);
for comb = combs
    predir = nan;
    postdir = nan;
    nrats = 0;
    for r = 1:nrecs
        if isfield(recs(r).tuning.prefdir, comb)
            conninds = recs(r).tuning.prefdir.(comb).linearInds;
            valid = find(recs(r).isconn(conninds));
%             valid = 1:numel(conninds);
            predir = [predir; recs(r).tuning.prefdir.(comb).predir(valid)];
            postdir = [postdir; recs(r).tuning.prefdir.(comb).postdir(valid)];
            recs(r).tuning.prefdir.(comb).linearIndsValid = valid;
            % disp(r)
            nrats = nrats+1;
        end
    end
    res.(comb).predir = predir(2:end);
    res.(comb).postdir = postdir(2:end);
    
    if numel(predir)>1
        [res.(comb).muoffset, res.(comb).rho, res.(comb).pval, res.(comb).muabsoffset, stats] = computeDirAlignment(predir,postdir);
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
    % [rho, pval] = circ_corrcc(predir, postdir);
    pfdoffset = circ_dist(predir, postdir);
    npairs = numel(predir);
    disp(npairs)
    predir = predir+[0,2*pi, 2*pi,0];
    postdir = postdir+[0,0, 2*pi,2*pi];
    predir = predir(:);
    postdir = postdir(:);
    if doplot
        nexttile
        xlim([-pi,3*pi]), ylim([-pi, 3*pi]);
        % scatter(predir, postdir, 3, 'k');
        d = diagline('r');
        plot(d.XData+2*pi, d.YData, 'r');
        plot(d.XData, d.YData+2*pi, 'r');

        scatter(predir, postdir, 3, 'k');
        title(comb)
        if ncombs>1
        ctypes = strsplit(comb, "_");
        xlabel("Prefdir. "+ctypes{1})
        ylabel("Prefdir. "+ctypes{2})
        else
        xlabel("Prefdir. Directional")
        ylabel("Prefdir. Conjunctive")
        end
        text (2*pi, 4.6*pi, ...
            sprintf("%.0f pairs from %.0frats\nr:%.3f, p:%.4f", ...
            npairs, nrats, rho, pval))
        xlim([-pi,3*pi]), ylim([-pi, 3*pi]);
        xticks(-pi:2*pi:3*pi)
        yticks(-pi:2*pi:3*pi)
        xticklabels(0:360:720);
        yticklabels(0:360:720);
        ax = gca; ax.FontSize = 11;
        axis square
    end

end
