This repository includes code designed to perform analyses, organized by figure numbers, on the demo data located in the "sample_data" subfolder, pertaining to the publication titled "Goal-directed hippocampal theta sweeps during memory-guided navigation" authored by Tang, Mei, Harvey, Carbajal-Leon, Netzer, Chang, Oliva, and Fernandez-Ruiz (Nature Neuroscience, 2026).

## Requirements
External Matlab packages (included in this repository):
- Left-right-alternating_sweeps_2025 (https://zenodo.org/records/14548054)
- LMT (https://github.com/waq1129/LMT)
- ln-model-of-mec-neurons (https://github.com/GiocomoLab/ln-model-of-mec-neurons)
- Widloski_model_14728054 (https://zenodo.org/records/14728054)
Note: the package was implemented on MATLAB 2023b. Backward compatibility with older versions is not guaranteed.
Operating system: Mac or Windows.

## Setup
Copy or clone the repository to your local machine. 
Add the Matlab subfolder to the Matlab search path.
Specify the path to the demo data and code by setting the properties "S.dataRoot_" and "S.codeRoot_" to the respective folders, either by modifying these properties in "SweepsSettings.m" or by configuring the settings accordingly.

## Main Scripts
### 1. 'demo_AlterThetaSweeps.m'
Demonstrate left-right alternating theta sweeps during random foraging.
### 2. 'demo_GoalThetaSweeps.m'
Demonstrate goal-directed theta sweeps during the cheeseboard task.
### 3. 'demo_LNmodel_GoalCell.m'
Demonstrate the identification of goal-direction cells using the LN model.
### 4. 'demo_ThetaModulation_GoalCells.m'
Demonstrate the theta modulation of (non-)goal-direction cell firing during goal theta sweeps.
### 5. 'model_thetaSequences_2D_simple.m'
Demonstrate the left-right alternating theta sweeps in the CAN model.
### 6. 'model_thetaSequences_2D_simple_Goal.m'
Demonstrate the goal theta sweeps in the CAN model.
### 7. 'demo_Replay_toGoal.m'
Demonstrate the replay sweeps to the goal during the cheeseboard task.
### 8. 'demo_CA1PFC_thetasweep_mismatch.m'
Demonstrate the decoding mismatch between CA1 and PFC during goal vs. lateral theta sweeps.
### 9. 'demo_LatentSweep.m'
Demonstrate theta sweeps on the latent maps.

## Support
Questions can be directed to wenbo.tang07@gmail.com
