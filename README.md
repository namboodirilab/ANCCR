# ANCCR
Here we provided three different types of codes in separate folders:
1. simulation:
  This folder contains matlab scripts that simulate different experiments with three models (compound serial component, microstimlus, and ANCCR models)
2. analysis:
  This folder contains matlab and python scripts that analyze experimental data. Every experimental data was converted to Neuroscience Without Border (NWB) format and posted on DANDI (https://dandiarchive.org/dandiset/000351).
  Please note that we initially used matlab to analyze data, but transitioned to python, so analysis codes for Experiment 1-6 are in matlab, but codes for Experiment 7 are in python. We could have rewriteen the codes in one language, but did not to keep the original scripts used to generate figures.
  Every python code directly imports NWB files from DANDI server and uses them for analysis. There is an unresolved issue in importing NWB files in matlab, so we used mat file which was regenerated from nwb file in python. 
  We are working to resolve the importing issue in coordination with nwb team (https://github.com/NeurodataWithoutBorders/helpdesk/discussions/36), and matlab codes will be updated to directly use NWB files as soon as possible. 
3. function:
  This folder contains matlab and python functions that are used in simulation and analysis scripts.
  
Lastly,
4. extension: This folder contains yaml files that are needed for NWB file extensions.

last update: 11/9/2022
