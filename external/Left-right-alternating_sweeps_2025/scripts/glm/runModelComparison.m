function [L, scores, fitData] = runModelComparison(spikeCounts, dt, variablesD, interpD, fspaces)
% Run model-selection comparison for non-shift

shiftData = [];

for useHd = [0, 1]
    for usePos = [0, 1]
        for useId = [0, 1]
            
            betaTerms = {'intercept'};
            if useHd,  betaTerms{end+1} = 'hd'; end
            if useId,  betaTerms{end+1} = 'id'; end
            if usePos, betaTerms{end+1} = 'pos'; end
            
            mdlName = join(string(betaTerms), "_");
            
            fitDataTmp = fitPosShiftGlm(variablesD, interpD, fspaces, shiftData, spikeCounts, dt, ...
                "alphaTerms", {}, ...
                "betaTerms", betaTerms, ...
                "plot", false, ...
                "display", "off");
            
            fitDataTmp = rmfield(fitDataTmp, ["opts", "fitPosShiftGlmParams"]);
            
            fitData.(mdlName) = fitDataTmp;
            
        end
    end
end

% fds = fieldnames(fitData);
% for f = 1:numel(fds)
%     fd = fds{f};
%     L.(fd) = fitData.(fd).llh;
% end
% 
% L = structfun(@(s) s.llh

L.null      = fitData.intercept.llh;
L.hd        = fitData.intercept_hd.llh;
L.id        = fitData.intercept_id.llh;
L.hd_id     = fitData.intercept_hd_id.llh;
L.hd_id_pos = fitData.intercept_hd_id_pos.llh;

% HD modulation: HD vs. null model
scores.hd = L.null ./ L.hd;

% ID modulation: ID vs. null model
scores.id = L.null ./ L.id;

% Conjunctive HD+ID modulation: HD+ID vs. best of ID and HD
scores.hd_id = L.null ./ L.hd_id;

% Pos modulation: HD+ID+POS vs. HD+ID
scores.pos = L.null ./ L.hd_id_pos;

end