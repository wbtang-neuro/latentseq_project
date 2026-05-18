function [xy, dt] = createLinearTrajectory(numTimeSteps)
% Create a linear locomotion trajectory (position units are arbitrary)
runSpeed = 1;
dt = 1/8;
dx = runSpeed * dt;
xpos = cumsum(ones(numTimeSteps, 1)*dx); xpos = xpos-mean(xpos);
ypos = zeros(size(xpos));
xy = [xpos, ypos];
end