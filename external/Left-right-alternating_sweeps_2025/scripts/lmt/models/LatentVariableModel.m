classdef (Abstract) LatentVariableModel < TuningModel
% Abstract class adding latent-variable functionality to TuningModel
    properties
        ndims
        X                   % the latent variable (matrix, dimensions time x ndims)
        Xinit               % initial state of X%
        XL = NaN
        XstoreIter = false
        XAligned
        XAlignSmooth = 50   % smooth X over this number of samples for aligning to Xinit
        XTform              % transformation matrix for mapping X -> Xinit
        XAlignInds
        enableXStep = true      % enable optimization of X?
        enableXAlignment = true
        XminFuncOptions = defaultMinFuncOptions();
        isCircular = false  % true for special case where X is a 1-D circular variable
    end

    methods (Abstract)
        [L, X] = onOptimizeX(self, inpLogR)
    end

    methods
        function self = LatentVariableModel(varargin)
            self@TuningModel(varargin{:});
            if nargin
                self.Xinit = self.castArr(varargin{2});
                self.ndims = size(self.Xinit, 2);
            end
        end

        function onInit(self)
            self.X = self.Xinit;
            if self.isCircular
                self.X = wrapToPi(self.X);
            end
            self.XL = [];
            if self.enableXAlignment
                self.calcXTform();
            end
            self.XAligned = self.alignX();
        end

        function L = onStep(self, inpLogR)
            L = self.onStep@TuningModel(inpLogR);
            if self.enableXStep
                [L, self.X] = self.onOptimizeX(inpLogR);
                if self.enableXAlignment
                    self.calcXTform();
                end
                self.XAligned = self.alignX(self.X);
            end
        end


        function Xout = alignX(self, Xin, direction)
            % Transform coordinates between the initial and current X
            % reference frames.
            %
            % "forward" maps X -> Xinit
            % "reverse" maps Xinit -> X

            if nargin < 2 || isempty(Xin), Xin = self.X; end
            if nargin < 3 || isempty(direction), direction = "forward"; end
            direction = lower(direction);

            if ~self.enableXAlignment
                Xout = Xin;
                return;
            end

            if self.isCircular
                if self.nf == 1
                    % If X is circular, we just subtract the offset angle(s)
                    alpha = self.XTform;
                    if direction == "forward"
                        % X -> Xinit (remove the offset)
                        Xout = Xin - alpha;
                    elseif direction == "reverse"
                        % Xinit -> X (add the offset)
                        Xout = Xin + alpha;
                    end
                    Xout = wrapToPi(Xout);
                else
                    % Don't try to align other "circular" data (e.g. torus)
                    Xout = Xin;
                end
            else
                T = self.XTform;
                % In the fitted transformation, X is the "moving" points
                % and Xinit is the "fixed" points.
                if direction == "forward"
                    % Transform X -> Xinit
                    Xout = T.transformPointsForward(Xin);
                elseif direction == "reverse"
                    % Xinit -> X
                    Xout = T.transformPointsInverse(Xin);
                end
            end

        end

        function calcXTform(self)
            % smooth X and find linear mapping between X <-> Xinit
            sm = self.XAlignSmooth;
            X = self.X;
            Xinit = self.Xinit;

            iscirc = self.isCircular && self.nf==1;

            if sm
                sigma = self.XAlignSmooth;
                if iscirc
                    X = gsmoothcirc(X, sigma);
                else
                    X = gsmooth(X, sigma);
                end
            end
            inds = self.XAlignInds;
            if ~isempty(inds)
                X = X(inds, :);
                Xinit = Xinit(inds, :);
            end

            if iscirc
                % transformation T corresponds to the mean circ. difference
                % X and Xinit (i.e. a positive-valued T indicates that X is
                % positively rotated w.r.t. Xinit)
                a = circ_dist(X, Xinit);
                T = circ_mean(a);
            else
                T = fitgeotrans(X, Xinit, "nonreflectivesimilarity"); % in this setup, transforming *FORWARD* maps X->Xinit
            end
            self.XTform = T;
        end

    end

end