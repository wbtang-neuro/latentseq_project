function [prob] = decodeBayes(spikeCounts,tuning,dt,nsmooth)
%THEREALBAYESIAN Summary of this function goes here
%   Detailed explanation goes here
% spk     % Spike count matrix (t x neurons)
% tuning  % Tuning curves (bins x neurons)
% dt      % Time bin size
% AZV 2023
tuning = tuning'+eps;
logtuning = log(tuning);

llh = spikeCounts * logtuning - dt.*sum(tuning, 'omitnan'); 
llh = gsmooth(llh, nsmooth);
prob = exp(double(llh));
prob = prob./sum(prob, 2, 'omitnan');
end
