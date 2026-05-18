classdef DRFunctionSpace < BaseFunctionSpace & matlab.mixin.Copyable
    %DRFUNCTIONSPACE a dimension-reduced function space
    
    properties(Hidden)
        drPoints
        drAlgorithm
        fsImpl
    end
    
    properties(Dependent)
        plotFcn
        defaultPlotX
        defaultPlotAxesGrid
        explained
        nDimsOriginal
        dimInds
    end
    
    methods
        
        function self = DRFunctionSpace(FSpace, DR)
            % DRFUNCTIONSPACE constructor
            %
            % FR = DRFUNCTIONSPACE(F) creates an empty DRFunctionSpace
            % object FR, based on the FunctionSpace object F. FR will use
            % the same basis functions as F, but will reduce the
            % dimensionality of F values by applying PCA to the values of 
            % its basis functions.
            %
            % FR = DRFUNCTIONSPACE(F, METHOD) specifies the dimensionality
            % reduction algorithm to be used. METHOD can be "pca" (standard
            % PCA), "pcacorr" (PCA based on correlation matrix), or "svd".
            
            if nargin
                args = {FSpace.nDims};
            else
                args = {};
            end
            self@BaseFunctionSpace(args{:});
            if nargin 
                self.fsImpl = FSpace.copy();
                self.fsImpl.name = [FSpace.name '_impl'];
                if nargin < 2 || isempty(DR)
                    DR = PCA();
                else
                    % DR can be a string or a DRAlgorithm obj
                    if ischar(DR)
                        switch lower(DR)
                            case 'pca'
                                DR = PCA();
                            case 'svd'
                                DR = SVD();
                            case 'pcacorr'
                                DR = PCACorr();
                        end
                    end
                end
                self.drAlgorithm = DR;
            end
        end
        
        function D = runDR(self, X)
            % Run dim-reduction on data in the basis functions' domain
            self.drPoints = X;
            Y = self.fsImpl.evaluate(X);
            self.drAlgorithm.run(Y);
            FS = self.copy();
            D = DRFDecomp(FS, X);
        end
        
        function indsKeep = discard(self, nKeep, mode, mult)
            if nargin < 4, mult = []; end
            if nargin < 3, mode = ''; end
            indsKeep = self.drAlgorithm.discard(nKeep, mode, mult);
            self.discardDims(indsKeep);
        end
        
        function indsKeep = discardThresh(self, thresh, mode, mult)
            if nargin < 4, mult = []; end
            if nargin < 3, mode = ''; end
            indsKeep = self.drAlgorithm.discardThresh(thresh, mode, mult);
            self.discardDims(indsKeep);
        end
        
        function Y = transformFwd(self, X)
            % Map codomain -> DR-codomain
            Y = self.drAlgorithm.transformFwd(X);
            Y = self.scaleForward(Y);
        end
        
        function X = transformRev(self, Y)
            % Map DR-codomain -> codomain
            Y = self.scaleReverse(Y);
            X = self.drAlgorithm.transformRev(Y);
        end
        
        function val = get.explained(self)
            val = self.drAlgorithm.explained;
        end
        
        function val = get.nDimsOriginal(self)
            val = self.fsImpl.nDims;
        end
        
        function val = get.dimInds(self)
            val = self.drAlgorithm.dimIndsR;
        end
        
        function val = get.plotFcn(self)
            val = self.fsImpl.plotFcn;
        end
        
        function set.plotFcn(self, val)
            self.fsImpl.plotFcn = val;
        end
        
        function val = get.defaultPlotX(self)
            val = self.fsImpl.defaultPlotX;
        end
        
        function set.defaultPlotX(self, val)
            self.fsImpl.defaultPlotX = val;
        end
        
        function val = get.defaultPlotAxesGrid(self)
            val = self.fsImpl.defaultPlotAxesGrid;
        end
        
        function set.defaultPlotAxesGrid(self, val)
            self.fsImpl.defaultPlotAxesGrid = val;
        end
        
    end
    
    methods (Access = protected)
        function cpSelf = copyElement(self)
            cpSelf = copyElement@matlab.mixin.Copyable(self);
            cpSelf.fsImpl = self.fsImpl.copy();
            cpSelf.drAlgorithm = self.drAlgorithm.copy();
        end
        
        function val = nDimsImpl(self)
            val = self.drAlgorithm.nDimsR;
        end
        
        function val = fcnNDimsInImpl(self)
            val = self.fsImpl.fcnNDimsIn;
        end
        
        function val = fcnNDimsOutImpl(self)
            val = self.fsImpl.fcnNDimsOut;
        end
        
        function discardDims(self, indsKeep)
            self.nDims = numel(indsKeep);
            self.dimScales = self.dimScales(indsKeep);
        end
        
        function Y = evaluateImpl(self, X, p)
            % Function evaluation at a point in DR function space
            
            % Evaluate original fcn space at all-ones position. No scaling
            % is done here
            Y = self.fsImpl.evaluate(X);
            
            % Transform back into DR space and apply weights
            Y = p .* self.drAlgorithm.transformFwd(Y);
        end
        
        function y = evaluate1Impl(self, X, idx)
            % Evaluate one axis in DR function space
            
            % Evaluate original fcn space at all-ones position
            Y0 = self.fsImpl.evaluate(X);
            % Transform and apply weighting to reduced dim
            Y = self.drAlgorithm.transformFwd(Y0);
            y = Y(:, idx);
        end
        
        function D = decompImpl(self, X)
            D = DRFDecomp(self, X);
        end
        
    end
    
end