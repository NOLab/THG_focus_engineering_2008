# THG_focus_engineering_2008
Third-harmonic generation microscopy with focus-engineered beams: a numerical Study


The same code is used for all figures (but since it is light, it is copied in every folder), the only difference is the script to run the simulation, and the script for the analysis.


The only change to the code was to clean up the scripts and fix the isotropic tensor (in function 'assembletout.m")


Any questions, please send an email to nicolas.olivier "at" polytechnique.edu


A better (but more limited for some aspects) version of this code is available on the repository https://github.com/LaboratoryOpticsBiosciences/FDTD_for_NL_microscopy


Post-processing: once the simulation is done running, the results can be quickly visualized using the script "Script_Post-processing.m" (that requires the 2 functions 'traitementdatafinal.m" and "decoupage.m" in the same folder)
