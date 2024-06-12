This folder contains the code and data for the H3N2 antibody transmission project (G Hoy, contact ghoy@umich.edu)

Code folders: 
01 code_antibody_transmisison_working_binned: core model with binned antibody categories
02 code_antibody_transmission_working_additive: model with categorical high vs. low variable for each assay, summed (0, 1, 2, or 3)
03 code_antibody_transmission_working_nacontrast: model to test high vs. low NA in participants with detectable anti-HA immunity

Data files:
[All data files are deindentified and created from this SAS program: T:\Nica Projects\Family Cohort\Projects-Flu\2023_Hoy_AntibodyTransmission\Code\H3N2 dataset generation.sas]
flu_data_v2.txt - data file for core model (corresponds to 01 code_antibody_transmisison_working_binned)
flu_data_simulated.txt - data file for the simulation analyses of the core model
flu_data_v3_additive - data file for the additive model, corresponds to 02 code_antibody_transmission_working_additive
flu_data_v4_nacontrast.txt - data file for the NA contrast model, corresponds to 03 code_antibody_transmission_working_nacontrast
flu_data_v5_additive_na.txt - data file for the additive model splitting out those with high antibodies in 2+ assays with NA vs. those with high antibodies in 2+ assays without NA, corresponds to 04 code_antibody_transmission_working_additive_na