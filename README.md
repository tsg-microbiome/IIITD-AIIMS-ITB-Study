# IIITD-AIIMS-ITB-Study
The repository contains the data and the codes utilized for the investigations performed as part of a study profiling the gut microbiome and mycobiome alterations in Intestinal Tuberculosis and Crohns' Disease patients.

The analysis was performed in different stages, each focused on a specific objective. We have placed a detailed code-based pipeline for each stage as described below:

**Stage 1**

The codes corresponding to this part of the study are uploaded as *Codes_Stage1.r*. This part reads all the input data files placed in the *Data_Files* folder and some user defined functions placed in *function_library.R* file. The later part of the codes also performs preliminary checks into various alpha-diversity and investigation of prevalent microbes.

**Stage 2**

The codes corresponding to this part of the study are uploaded as *Codes_Stage2.r*. This part performs all the beta-diversity analyses (e.g. Principal Coordinates Analysis and PERMANOVA) on microbiome and mycobiome profiles using various distance measures. It uses the 'batch_dunns.R' that is already present in the *function_library.R* file. 
