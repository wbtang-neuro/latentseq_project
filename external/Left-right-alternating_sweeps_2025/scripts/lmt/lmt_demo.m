% This script demonstrates how to fit the "composite"
% latent-manifolder-tuning (LMT) model to multi-neuron spike-train data

% It's highly recommended to run this on a machine with a GPU. Running on
% CPU is possible, but will be very slow.

clear
S = SweepsSettings;
% S.dataRoot_ = "~/Temp/sweeps/data"; % change this to the root folder where the data files are saved
S.dataRoot_ = '/Users/wt248/Downloads/Left-right-alternating_sweeps_2025/sample_data/';

%% Load the session dataset file

fn = fullfile(S.dataRoot, "28229_3.mat");
load(fn, "Dsession");

%% Generate the time-by-units matrix of spike counts

% units = Dsession.units.mec;
units = Dsession.units.hc;

for u = 1:numel(units)
    units(u).spikeCounts = spikeIndsToCounts(units(u).spikeInds, Dsession.nt);
end

%% Configure the composite model and run the fitting

% Enable the "useGpu" option if your machine has a NVIDIA GPU
compModel = createCompositeModel(Dsession, units, useGpu=false, plotInterval=1);
compModel.fit();

%% Export the fitted LMT data

% The 'lmtData' struct generated below has the same format as the contents
% of the 'lmt' field in 'Dsession'.

lmtData = struct();

for m = 1:compModel.nmodels
    mdl = compModel.models{m};
    if isa(mdl, 'LmtModel')
        lmtData.(mdl.name) = createLmtDataStruct(Dsession, units, mdl);
    end
end
