classdef LmtModel < LatentVariableModel

    properties

        invC

        Finit

        XboundPercentile = 0.1
        Fgrid
        Fgridv
        trimXToGrid = true
        XclampIter = 0;
        dt = 0.1;

        % Time-bin chunking
        ntseg = 10000
        % ntseg = [];
        nseg
        maxChunkNT
        binTimes                  % use this if bins are not uniformly spaced
        tgrid

        % Annealing parameters
        learnRate = 0.98
        annealPeriod = 0
        annealIter
        annealIterStart  % start annealing at this iteration (should be same as XclampIter)
        annealCycle = 0
        annealNCycles = 1

        % P-GPLVM HYPERPARAMS
        hparams = struct( ...
            "rhoff", 1, ...
            "rhoxx", 1e4, ...
            "lenxx", 0.015, ...
            "lamff", 0, ...
            "lenffR", 1, ...
            "sigma", 1);

        hparamRanges = struct( ...
            "rhoff", [], ...
            "rhoxx", [], ...
            "lenxx", [], ...
            "lamff", [], ...
            "sigma", [3, 1]);

    end

    properties (Dependent)
        FgridBounds
    end

    properties (Constant)
        MAX_CHUNK_NUMEL = 2e7
    end

    methods

        function self = LmtModel(varargin)
            self@LatentVariableModel(varargin{:});
            self.nf = 30;
        end

        function onInit(self)

            if self.isCircular
                self.enableXAlignment = false;
            end

            self.onInit@LatentVariableModel();

            if ~isempty(self.binTimes)
                assert(numel(self.binTimes)==self.nt, "Number of values in properites 'binTimes' must equal the number of time bins in Y and X");
            end

            % Define the max chunk size to cap the max. amount of memory
            % required
            self.maxChunkNT = self.MAX_CHUNK_NUMEL / max(self.nunits, self.nf.^self.ndims);
            self.initSegments();
            self.initF();

            % Annealing
            if isempty(self.annealIterStart)
                if self.enableXStep
                    self.annealIterStart = self.XclampIter;
                else
                    self.annealIterStart = 1;
                end
            end
            self.annealIter = 0;

            self.initHparams();
            self.checkHparams();
            self.Fgrid = [];
            self.Fgridv = [];
            self.X = self.updateFgrid(self.X);

            % Initialize invC (normally happens at the start of the
            % F-step, but we do it here too, just in case the F-step is
            % disabled)
            [~, self.invC] = self.calcXWeights();
        end

        function initSegments(self)

            % Either the number of segments or segment length can be set
            % manually.
            nseg = self.nseg;
            ntseg = self.ntseg;

            if isempty(nseg) && isempty(ntseg)
                % no user-specified segment length, calculate
                % automatically
                [~, nseg] = self.calcLargestSegmentSize();
                ntseg = self.nt / nseg;
            elseif isempty(nseg) && ~isempty(ntseg)
                % segment length manually set
                nseg = self.nt / ntseg;
            elseif ~isempty(nseg) && isempty(ntseg)
                % number of segments manually set
                ntseg = self.nt / nseg;
            end

            assert(rem(ntseg, 1)==0, "Segment length must be a whole number of samples");
            assert(rem(nseg, 1)==0, "Number of segments must be an integer");
           
             
            self.nseg = nseg;
            self.ntseg = ntseg;
            self.tgrid = [];

            for s = 1:nseg
                inds = self.segmentInds(s);
                if isempty(self.binTimes)
                    self.tgrid(:, s) = inds * self.dt;
                else
                    self.tgrid(:, s) = self.binTimes(inds);
                end
            end
            fprintf("Total %u data points, splitting into %u chunks of length %u\n", self.nt, nseg, ntseg);
        end

        function inds = segmentInds(self, iseg)
            n = self.ntseg;
            i0 = n * (iseg-1);
            inds = i0 + (1:n)';
        end

        function checkHparams(self)
            props = ["hparams", "hparamRanges"];
            for n = 1:2
                h = self.(props(n));
                if isfield(h, "lambdaff")
                    h = renameStructField(h, 'lambdaff', 'lamff');
                    self.(props(n)) = h;
                end

                names = fieldnames(h);
                valid = ismember(names, ["rhoff", "rhoxx", "lenxx", "lamff", "lenffR", "sigma", "lenff"]);

                if any(~valid)
                    idx = find(~valid, 1);
                    error("Invalid fieldname '%s' in %s", names{idx}, props(n));
                end
            end
        end

        function initHparams(self)
            hr = self.hparamRanges;
            hpnames = string(fieldnames(hr))';
            for name = hpnames
                r = hr.(name);
                if ~isempty(r)
                    self.hparams.(name) = r(1);
                end
            end
        end

        function startNewAnnealCycle(self)
            self.initHparams();
            self.annealIter = 0;
            self.annealCycle = self.annealCycle + 1;
        end

        function annealStep(self)

            self.annealIter = self.annealIter+1;
            ranges = self.hparamRanges;
            names = string(fieldnames(ranges))';

            for name = names
                r = ranges.(name);
                % r will be empty if the hyperparameter value is constant
                if ~isempty(r)
                    val = annealHparam(r, self.annealIter, self.learnRate);
                    self.hparams.(name) = val;
                end
            end
        end

        function onIterInit(self)
            % Runs before every optimization iteration

            iAnnealCyc = self.annealCycle;
            annealNCyc = self.annealNCycles;

            doAnneal = annealNCyc>0 && (self.iter >= self.annealIterStart);

            if doAnneal
                if isequal(self.annealPeriod, "auto")
                    newCycle = self.hparams.sigma <= self.sigmaFinal;
                else
                    newCycle = rem(self.annealIter, self.annealPeriod) == 0;
                end
                newCycle = newCycle & iAnnealCyc<annealNCyc;

                if newCycle
                    self.startNewAnnealCycle();
                end
                self.annealStep();
            end

        end

        function [L, Fnew] = onOptimizeF(self, F, inpR)
            % F-step subroutine to optimize tuning curves

            [wX, invC] = self.calcXWeights();
            self.invC = invC; % we use this again when optimizing X

            Y = self.castArr(self.Y);
            F = F(:);

            costFcn = @(prm) StateSpaceModelsofSpikeTrains_ref_rg_inp( ...
                prm, Y, invC, wX, inpR, self.hparams.lamff);

            [Fnew, L] = minFunc(costFcn, F, self.minFuncOptions);
            Fnew = reshape(Fnew, size(self.F));
            Fnew = Fnew-mean(Fnew); % enforce mean zero for each cell
        end

        function [L, X] = onOptimizeX(self, inpLogR)

            clampX = self.iter <= self.XclampIter;

            if clampX
                L = NaN;
                X = self.X;
                return;
            end

            Fk = self.invC*self.F; % F adjusted for gaussian covariance

            nd = self.ndims;
            nt = self.nt;
            nts = self.ntseg;
            ns = self.nseg;
            nu = self.nunits;

            inpLogRSegs = reshape(inpLogR, nts, ns, nu);
            Ysegs = reshape(self.Y, nts, ns, nu);

            XsegsNew = zeros(nts, ns, nd, "like", self.X);
            L = 0;

            h = self.hparams;
            covfcn = self.covfun();

            tgrid = self.tgrid;
            tgrid = tgrid - tgrid(1, :); % all times relative to segment start

            Xsegs = reshape(double(gather(self.X)), nts, ns, nd);

            parfor s = 1:ns
            % for s = 1:ns

                % Get the functions for mapping X <--> U.
                % This transformation applies the gaussian prior to X.
                [Bfun, BTfun] = prior_kernel_sp(h.rhoxx, h.lenxx, tgrid(:, s));
                Bfun  = permuteBfun(Bfun);
                BTfun = permuteBfun(BTfun);

                Xseg = squeeze(Xsegs(:, s, :)); % t, seg, dim
                Xseg3 = reshape(Xseg, 1, nts, nd); % (1, nt, nf)

                uu = Bfun(Xseg3,1); % X -> U
                if size(uu, 2) < nts
                    uu(:, (end+1):nts, :) = 0;
                end
                uu = reshape(uu, nts, nd);

                inpLogRSeg = squeeze(inpLogRSegs(:, s, :));

                Yseg = squeeze(Ysegs(:, s, :));
                FkY = -Fk*Yseg';

                lmlifun = @(u) logmargli_gplvm_se_block_ref_nogrid_rg_inp( ...
                    u, self.Fgrid, Fk, Yseg, Bfun, covfcn, nd, BTfun, 1, ...
                    inpLogRSeg, FkY, self.isCircular);

                [USegNew, Lseg] = minFunc(@(u) lmlifun(u), uu(:), self.XminFuncOptions);

                L = L + Lseg;

                USegNew = reshape(USegNew, 1, nts, nd);
                XsegNew = Bfun(USegNew, 0); % U -- > X
                XsegNew = reshape(XsegNew, nts, nd);
                XsegsNew(:, s, :) = XsegNew;
            end

            X = reshape(XsegsNew, nt, nd);
            X = self.updateFgrid(X); % N.B. if X is circular, it gets wrapped here
        end

        function yh = onPredictLogY(self, icol, X)
            if nargin < 2 || isempty(icol), icol = 1:self.nunits; end
            if nargin < 3 || isempty(X), X = self.X; end
            % predict log-firing rates from tuning curves
            wX = self.calcXWeights(X);
            n = self.nf ^ self.ndims;
            F = zeros(n, numel(icol), "like", self.F);
            vcol = icol~=0;
            isrc = icol(vcol);
            F(:, vcol) = self.F(:, isrc);
            yh = wX*F;
        end

        function [wX, invC, covfun] = calcXWeights(self, X)
            % Calculate the weighting of each X observation onto each bin
            % of the F grid.
            if nargin < 2 || isempty(X), X = self.X; end
            [covfun, invC] = self.covfun();
            Xg = self.Fgrid;
            wX = covfun(X, Xg) * invC; % original
        end

        function F = logTuningCurvesF(self)
            % Convert F into log-firing-rate "tuning curve"
            Xg = self.Fgrid;
            F = self.onPredictLogY([], Xg);
        end

        function str = onIterString(self)
            hprm = self.hparams;
            fds = fieldnames(hprm);
            for f = 1:numel(fds)
                strs(f) = sprintf("%6s=%.3g", fds{f}, hprm.(fds{f}));
            end
            str = sprintf("%-13s", strs);
            str = strip(str);

            if self.iter <= self.XclampIter
                str = str + sprintf(" <CLAMPED X>");
            end
        end

        function initF(self)
            gridsz = self.nf^self.ndims;
            F = 1e-4 * randn(gridsz, self.nunits, "like", self.Y);
            self.F = F;
            self.Finit = F;
        end

        function [ntseg, nseg] = calcLargestSegmentSize(self)
            % Calculate the largest allowed segment size
            divisors = [1, 2, 5, 10, 20, 50, 100];
            ntseg = self.nt;
            c = 1;
            while ntseg > self.maxChunkNT
                c = c+1;
                ntseg = self.nt/divisors(c);
            end
            nseg = divisors(c);
        end

        function X = updateFgrid(self, X)
            % TODO: some operations in this method aren't closely related
            % to the grid update, and may be better put somewhere else.

            if self.isCircular
                % If X is circular, the grid has a fixed range of [-pi, pi]
                % therefore it only needs computing once
                if isempty(self.Fgrid)
                    gv = linspace(-pi, pi, self.nf+1);
                    gv = edg2cen(gv);
                    gv = cast(gv, 'like', X);
                    gv = repmat({gv}, 1, self.ndims);
                    self.Fgridv = gv;
                    [gg{1:self.ndims}] = ndgrid(gv{:}); % N.B. this behaves differently from meshgrid!
                    gg = cellfun(@(x) x(:), gg, "uni", false);
                    self.Fgrid = cat(2, gg{:});
                end

            else
                [self.Fgrid, ~, self.Fgridv] = getGrid(X, self.nf, self.XboundPercentile);
            end

            % determine new value of hyperparameter lenff, as a function
            gbnd = self.FgridBounds;
            grng = diff(gbnd, [], 2);
            lenffR = self.hparams.lenffR;
            self.hparams.lenff = mean(grng)/lenffR;

            if self.trimXToGrid
                % could break this off into a new method
                if self.isCircular
                    X = wrapToPi(X);
                else
                    for d = 1:self.ndims
                        b = gbnd(d, :);
                        Xc = X(:, d);
                        Xc(Xc<b(1)) = b(1);
                        Xc(Xc>b(2)) = b(2);
                        X(:, d) = Xc;
                    end
                end
            end

        end

        function [fcn, invC, C] = covfun(self)
            % Latent manifold (spatial) covariance function
            % outputs:
            % fcn - covariance function, set by hyperparameters
            %   hparams.rhoff - variance
            %   hparams.lenff - length scale
            %
            % C, invC - covariance matrix and inverse, calculated for all
            % points on the latent manifold grid.
            %
            % N.B. output C (and hence invC) are adjusted by addition of
            % % "sigma", to the diagonal of C, which effectively adds
            % white noise

            h = self.hparams;

            if self.isCircular
                fcn = @(x1,x2) covarianceSECirc(h.lenff, h.rhoff, x1, x2);
            else
                fcn = @(x1,x2) covarianceSE(h.lenff, h.rhoff, x1, x2);
            end

            if nargout >= 2
                Fg = self.Fgrid;
                C = fcn(Fg, Fg);

                % add sigma to diagonal of covariance matrix (adds white
                % noise?)
                sigma = self.hparams.sigma; % should this be sigma^2?
                sdiag = C(1,1)*sigma*eye(size(C));
                C = C + sdiag;
                invC = pdinv(C);
            end

        end

        function cen = FcenterOfMass(self)
            import rg.spatial.grid.*
            R = self.logTuningCurvesF();
            R = exp(R);
            if self.isCircular
                [gg{1:self.ndims}] = ndgrid(self.Fgridv{:});
                for dim = 1:self.ndims
                    cen(:, dim) = circ_mean(gg{dim}(:), R);
                end
            else
                for dim = 1:self.ndims
                    cen(:, dim) = sum(self.Fgrid(:, dim) .* R) ./ sum(R);
                end
            end
        end

        function val = get.FgridBounds(self)
            if self.isCircular
                val = repmat([-pi, pi], self.ndims, 1);
            else
                bnd = cellfun(@(x) {x([1, end])}, self.Fgridv);
                val = cat(1, bnd{:});
            end
        end

    end

end