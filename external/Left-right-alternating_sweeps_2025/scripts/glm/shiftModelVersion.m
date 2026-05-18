function [ver, date] = shiftModelVersion()
% v2 2019-05-03
%   -exported units now thresholded for mean firing rate
%
% v3
%   -add option for multiple user-specified time ranges in prepareNpxData
%
% v4 2019-05-05
%   -add timestamp to fit files
%
% v5 2019-05-14
%   - npx export streamlined into single script
%   - automatic creation of slurm job scripts
%   - shift angle source is now recorded in fit files
%
% v6 2019-05-15
%   - single unit .mat files now contain: depth, amplitude, channel, ks-group
%   - non-master probe recordings are now precisely synced with tracking
%   - save data source paths in "rec" mat files
%
% v7 2019-05-21
%   - fix v6 bug causing failure to save ks2 unit labels
%   - add "checkHdAlignment" function
%   - prepareNpxData refactored as a function
%   - minor changes to fitting functions:
%       - fitOneUnit_2
%       - fitPosShiftGlm
%       - posShift
%
% v8 2019-05-29
%   - prepareNpxData bugfixes:
%       -LFP synchronization bug affecting "probe_2" recordings
%       -
%   - prepareNpxData now exports wideband LFP (glmData "lfp" field)
%   - change in behavior of decomposeGlmData (time bin selection arg added)
%
%
% v9 2019-06-07
%     - BD format changes
%         - User-defined time ranges are now saved
%         - Unit post-spike fspace objs removed (too bulky)
%
% v10 2019-08-28
%     - Format changes to support multi-probe datasets:
%         - Unit 'id' field no longer corresponds to kilosort ID
%     - Format changes to simplify comparisons across datasets:
%         - Universal basis data now lives at a fixed location: "N:\richarga\misc\to_abraham\misc\model_common_data"
%
% v11 2019-09-14
%     - Fix for HD alignment. N.B. datasets from previous versions may not
%       be correctly aligned
%     - Fix for ks label import
%
% v12 2019-10-02
%     - Add UMAP sweep-dir detection to prepareNpxData
%     - Theta segmentation and phase data now calculated from MUA
%     - Save info about corrections applied to tracking offsets
%     - Increase speed threshold from 0.025 to 0.05 m/s
%     - Added SD as a model variable and removed speed
%     - All angular variables now use the same basis set (6 PCs)
%
% 
%
% v13 2019-10-03
%    - Export more detailed data from SD UMAP detection
%
% v13 2019-10-07
%    - More efficient memory use during variable decomposition
%    - Fix bug in SD interpolation (affects pre-v13 'rec' structs) 
%
% v14 2020-10-07
%    - Testing use of GPU during data export and model fitting
%    - Rigid body error values now used for eliminating bad frames

ver = 14;
date = '2020-01-16';

end