broadbandComputation

=========

This is a set of Matlab tools for simulating broadband signals in LFP time series and comparing various analysis choices.

In this repository you find one main mlx script that runs a single simulation and analysis:

* s_simulateBroadband.mlx

This script allows you to see the effects of various simulation and analysis parameters and visualizes each step involved.
Please see the contents of this script for more details on the simulation and analysis steps.

Other scripts in the main folder run systematic comparisons of parameters: 

* Scripts starting with 'q' compare various simulations. 
* Scripts starting with 'd' analyze real ECoG data.

The subfolders contain the following:
* support: a set of subfunctions used in simulation and analysis.
* sample_data: mat files with pre-saved response profiles. 

Note that the data analysis scripts make use of functions from the Winawerlab ECoG_utils repository:
https://github.com/WinawerLab/ECoG_utils
