This repository contains code for the paper "Left-right-alternating sweeps in entorhinal-hippocampal maps of space" by Vollan, Gardner, Moser & Moser (Nature, 2025). Code is intended to be compatible with the accompanying dataset deposited at EBRAINS (https://doi.org/10.25493/R5FR-EDG).

## Requirements
External Matlab packages (included in this repository):
- Chronux v2-12 (http://chronux.org/)
- CircStat2011f (https://github.com/circstat/circstat-matlab)
- matGeom (https://github.com/mattools/matGeom)
- minFunc (http://www.cs.ubc.ca/~schmidtm/Software/minFunc.html)
Matlab version: 2024a
Hardware: It is recommended to run the code on a system with at least 32GB of RAM. An NVIDIA GPU is highly recommended for some of the analyses.
Operating system: Windows or Mac.
 
## Setup
Copy or clone repository to your local machine. Run the script called "sweepsSetup.m" to add relevant directories to Matlab's search path and configure settings (returned as SweepsSettings object "S"). It is recommended to download the accompanying dataset (https://doi.org/10.25493/R5FR-EDG) and store it locally. The path to the dataset can be specified by setting the property "S.dataRoot_" to the dataset folder (either by modyfing the property "S.dataRoot_" in "SweepsSettings.m" or by setting the property).