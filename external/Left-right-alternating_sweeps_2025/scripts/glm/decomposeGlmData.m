function [D, nanrows] = decomposeGlmData(fspaces, variables, vt, varargin)

if nargin < 3, vt = []; end

inp = inputParser();
inp.addParameter("plot", false);
inp.addParameter("checkRank", true);
inp.addParameter("nanAction", "error");
inp.addParameter("floatClass", "double");
inp.addParameter("useGpu", false);
inp.parse(varargin{:});
P = inp.Results;

% only use fields if we have both the data and the fspace
% fields = intersect(fieldnames(variables), fieldnames(fspaces));
fields = fieldnames(variables);

nFields = numel(fields);

if P.plot
    fig = gcf();
    sply = floor(sqrt(nFields));
    splx = ceil(nFields/sply);
end

for f = 1:nFields
    fd = fields{f};
    X = variables.(fd);
    X = cast(X, P.floatClass);
    if P.useGpu
        X = gpuArray(X);
    end
    
    % all variables should have the same number of rows
    if f==1
        X0 = X;
        if isempty(vt)
            vt = true(size(X, 1), 1);
        else
            assert(numel(vt) == size(X0,1));
        end
        nanrows = false(sum(vt), 1);
    else
        assert( ...
            size(X, 1) == size(X0, 1), ...
            'Number of rows in variables "%s" and "%s" do not match', ...
            fd, fields{1});
    end
    
    X = X(vt, :);
    Y = fspaces.(fd).evaluate(X);
    
    xnan = any(isnan(X), 2);
    nanrows = nanrows | xnan;
    if any(xnan)
        if P.nanAction == "error"
            error('NaN values were found in input data.');
        elseif P.nanAction == "warn"
            warning('NaN values were found in input data. These will be replaced with zeros in the output.');
        elseif P.nanAction ~= "none"
            error('Invalid value for parameter "nanAction"');
        end
        Y(xnan, :) = 0;
    end
    
    nrank = rank(Y);
    if P.checkRank && nrank < size(Y, 2)
        warning('glm:rankDeficient', 'Decomposition of variable "%s" is rank-deficient (ncols = %u, nrank = %u)', fd, size(Y,2), nrank);
    end
    
    if P.plot
        ax = subplot(sply, splx, f, 'parent', fig);
        imagesc(ax, corr(Y(~nanrows, :)));
        axis(ax, 'xy', 'equal', 'tight');
        xlabel(ax, 'Covariate #1');
        ylabel(ax, 'Covaraiate #2');
        title(fd);
        ax.CLim = [0, 1];
    end
    
    D.(fd).X = X;
    D.(fd).Y = Y;
    
end

end
