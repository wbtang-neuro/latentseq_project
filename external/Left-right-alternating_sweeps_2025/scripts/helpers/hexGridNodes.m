function points = hexGridNodes(phase, nRings, spacing, orientation)
% HEXGRIDNODES - generate XY coords of grid nodes
%
% This is an adaptation of Tor's code

% z = re^{i*theta}; x = Re|z|, y = Im|z|

% Example run: z = TSgridNodes3 (15, 20, 30);

% nRings = number of rings around center field
% spacing --- is spacing.
% orientation, given in radians

rInd = 1:nRings;
nInd = 6*nRings; % Total number of nodes, increases by 6 per ring
phi = spacing;

z = [];

for step = 1:length(rInd)
    
    a = rInd(step); % Ring number = unit distance along x-axis
    if nInd == 0 || step ~= length(rInd)
        b = 0:6*a-1; % # nodes in ring
    else
        b = 0:nInd-1;
    end
    
    c = mod(b,a);
    theta = atan((sqrt(3).*c)./(2*a-c)) + pi*(b-c)/(3*a) + orientation; % RG 21/03/2015
    r = phi .*sqrt((a-c).^2 +a.*c);
    z = [z r.*exp(1i .*theta)];
    
end

x = real(z)+phase(1);
y = imag(z)+phase(2);

points = [phase; x' y'];

end

