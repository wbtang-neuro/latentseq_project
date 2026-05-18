function coords = resizeobj(h,prop,anchor,scalechildren)
%RESIZEOBJ changes the dimensions of a 2D graphics object by a
%proportion of their original size
%
% INPUTS
% h      - object handle, or vector of handles
% prop   - proportion of original size to resize to.  Either a scalar or a
%          two-element vector [x y] specifying the proportions to shrink
%          the x and y dimension respectively.  E.g. prop = 1.2 will
%          enlarge the object by 20% in both dimensions, while prop = [0.5 1]
%          will shrink the x dimension by 50% while keeping the y dimension
%          the same
% anchor - defines the anchor point for resizing: 'center', 'left', 'right',
%          'top', or 'bottom'.
% scalechildren - determines whether the positions of children of the
%          resized object will be also be rescaled (1), or left in their
%          original poisitions (0)
%
% Adapted from resizeaxes

if nargin < 4 || isempty(scalechildren), scalechildren = 1; end
if nargin < 3 || isempty(anchor), anchor = 'center';  end
if nargin < 2 || isempty(prop), prop = 1; end

if length(prop) == 1,             prop = [prop prop]; end

if isempty(h), h = gca; end

nh = numel(h);

for a = 1:nh
    
    isfig = strcmp(get(h(a),'type'),'figure');
    
    pos = get(h(a),'position');
    szx = pos(3)*prop(1);
    szy = pos(4)*prop(2);
    dszx = (prop(1)-1)*pos(3);
    dszy = (prop(2)-1)*pos(4);
    
    % Calculate new coordinates for object and children (child coord
    % shifting is only necessary if children are to remain in same
    % poisition and left/bottom anchoring is in use - otherwise we don't
    % need to change them
    switch lower(anchor)
        case {'c','center'}
            pos = [pos(1)-dszx/2 pos(2)-dszy/2 szx szy];
            chshift = [dszx/2 dszy/2];
        case {'l','left'}
            pos = [pos(1) pos(2)-dszy/2 szx szy];
            chshift = [0 dszy/2];
        case {'r','right'}
            pos = [pos(1)-dszx pos(2)-dszy/2 szx szy];
            chshift = [dszx dszy/2];
        case {'t','top'}
            pos = [pos(1)-dszx/2 pos(2)-dszy szx szy];
            chshift = [dszx/2 dszy];
        case {'b','bottom'}
            pos = [pos(1)-dszx/2 pos(2) szx szy];
            chshift = [dszx/2 0];
        otherwise
            error('Invalid anchor')
    end
    
    % (For figures only) shift coords of children to keep them in place on
    % the figure relative to the specified anchor
    if isfig
        % Find all first-generation descendents
        ch = get(h(a),'children');
        chunit = get(ch,'units');
        
        if scalechildren
            set(ch,'units','normalized')
        else
            set(ch,'units','centimeters')
            chpos = get(ch,'position');
            
            % Compensate for shift in parent object's coordinates
            for aa = 1:length(ch)
                tmp = chpos{aa};
                tmp([1 2]) = tmp([1 2]) + chshift;
                set(ch(aa),'position',tmp)
            end
        end
    end
    
    % Resize figure
    set(h(a),'position',pos)
    
    % Restore original units for children
    if isfig
        for aa = 1:length(ch)
            set(ch(aa),'units',chunit{aa})
        end
    end
    
    if nargout
        coords(:,a) = pos;
    end
    
end

end

