classdef BasisFunction < matlab.mixin.Copyable
    
    properties
        fcn function_handle
        nDimsIn
        nDimsOut
        plotFcn
    end
    
    properties (Dependent)
        name
    end

    properties (Hidden)
        nameImpl = ''
    end
    
    methods
        function self = BasisFunction(fcn, nDimsIn, nDimsOut, name)
            if nargin
                self.fcn = fcn;
                self.nDimsIn = nDimsIn;
                self.nDimsOut = nDimsOut;
                self.name = name;
            end
        end
        
        function Y = evaluate(self, X)
            Y = self.fcn(X);
        end
        
        function bFcnW = wrap(self, bFcn)
            % WRAP wraps one BasisFunction around another
            %
            % BFWRAPPED = BF1.WRAP(BF2) wraps the function represented by
            % BF1 around the function represented by a second BasisFunction 
            % object BF2. The output BasisFunction object BFWRAPPED
            % represents the wrapped function.
            assert(self.nDimsIn == bFcn.nDimsOut);
            bFcnW = self.copy();
            bFcnW.nameImpl = strrep(self.nameImpl, '<ARG>', bFcn.nameImpl);
            bFcnW.fcn = @(X) self.fcn(bFcn.fcn(X));
        end
        
        function str = toString(self)
            for n = 1:numel(self)
                str = sprintf('BasisFunction "%s", %u-D domain, %u-D range\n', ...
                    self.mathStr, self.nDimsIn, self.nDimsOut);
            end
            str(end) = '';
        end
        
        function set.name(self, val)
            val = char(val);
            self.nameImpl = [val '(<ARG>)'];
        end
        
        function val = get.name(self)
            val = strrep(self.nameImpl, '<ARG>', 'X');
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % OVERLOADED OPERATORS
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function bFcn = plus(self, arg)
            bFcn = arithmeticSub(self, arg, @plus);
        end
        
        function bFcn = minus(self, arg)
            bFcn = arithmeticSub(self, arg, @minus);
        end
        
        function bFcn = times(self, arg)
            bFcn = arithmeticSub(self, arg, @times);
        end
        
        function bFcn = rdivide(self, arg)
            bFcn = arithmeticSub(self, arg, @rdivide);
        end
        
        function bFcn = power(self, arg)
            bFcn = arithmeticSub(self, arg, @power);
        end
        
    end
    
    methods (Hidden)
        
        function bFcn = arithmeticSub(self, bFcn, op)
            sz1 = size(self);
            sz2 = size(bFcn);
            sz = max([sz1; sz2]);
            bFcn1 = repmat(self, sz ./ sz1);
            bFcn2 = repmat(bFcn, sz ./ sz2);
            nObj = prod(sz);
            
            switch func2str(op)
                case 'times'
                    opStr = '.*';
                case 'rdivide'
                    opStr = './';
                case 'plus'
                    opStr = '+';
                case 'minus'
                    opStr = '-';
            end
            
            bFcns = BasisFunction.empty();
            for n = 1:nObj
                bf1 = bFcn1(n);
                bf2 = bFcn2(n);
                bFcn = BasisFunction();
                bFcn.fcn = @(X) op(bf1.evaluate(X), bf2.evaluate(X));
                bFcn.nDimsIn = bf1.nDimsIn;
                bFcn.nDimsOut = bf2.nDimsOut;
                bFcn.nameImpl = sprintf( ...
                    '%s %s %s', ...
                    bf1.nameImpl, ...
                    opStr, ...
                    bf2.nameImpl);
                bFcns(n) = bFcn;
            end
            bFcns = reshape(bFcns, sz);
        end
        
        function str = mathStrImpl(self, arg)
            if self.hasWrappedFcn()
                bf = self.wrappedBasisFcn;
                str0 = bf.mathStr;
            else
                str0 = arg;
            end
            str = sprintf('%s( %s )', self.name, str0);
        end
    end
end