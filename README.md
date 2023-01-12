# CKISCREEN_analysis_2022
Codes to reproduce the analysis of the manuscript:
CKI are numbered from 01 to 58 and the map can be found in /data/KI_Annotations/

01_DMSO_Impulse_and_Clustering.R
- Code to performed time-resolved differential analysis and CRE clustering of H3K27ac Chipseq signal in control (DMSO) condition using ImpulseDE2, UMAP/Louvain.
- The code can be run using the singularity image: docker://fgualdr/envimpulse

02_Combine_ATLAS.R
- Code to combine ~200 TFs Chipseq (both signal and MCAS2 called peaks) together with H3K27ac LPS/IL4 time course and H3K27ac IFNB1 time course in control (unperturbed) conditions
- The code can be run using the singularity image: docker://fgualdr/envrgeneralg

03_Main_MFA_analysis_and_plots.R
- Code to perform the main analysis devoted to the integration of CKI effects employing MFA (Multiple Factor Analysis) from the FactomineR package.

Folders:
/data/
|--BMDM_CHIP_ATLAS : Contains .RData files relative to the pre-processed TFs Chipseq ATLAS (i.e. harmonized Peaks and signals per CREs)
|--cre_counts : Contains a sample description "Normalisation_Parameters.txt" and the normalized read counts per CREs "Count_normalized.txt" for the H3K27ac Chipseq of the Screen presented in the manuscript
|--cre_counts_ifnb1_stimulation : normalized read counts per CREs "Count_normalized.txt" for the H3K27ac Chipseq following IFNb1 stimulation (GSE56121).
|--cre_dmso_impulseDE2 : ImpulseDE2 results (relative to script 01_DMSO_Impulse_and_Clustering.R)
|--KI_Annotations : Designated and Kinobeads CATDS based annotation of each of the 58 CKI used
|--SMALE_GSE67357 : ImpulseDE2 results (stored in .RData files) of the re-processing of the GSE67357 dataset.
|--total_rna_seq_ckiscreen : Total RNAseq DEseq2 results for the CKI screen (2h LPS stimulation with DMSO pre-treatment or the individual CKI)
|--total_rna_seq_timecourse : Total RNAseq time-course ImpulseDE2 results

/Results/
|--ATLAS_Harmonization_feature_transformation: Results of the 02_Combine_ATLAS.R code divided into lps and il4 stimulation including the X_Transformed_Features.txt file and the individual plots showing the correlation/v-test between original and transformed features
|--DMSO_Stat_folder_ImpulseDE2: Results relative to 01_DMSO_Impulse_and_Clustering.R code including ImpulseDE2 results and Macro cluster identification
|--lps : MAin results relative to 03_Main_MFA_analysis_and_plots.R for the LPS condition
|--il4 : MAin results relative to 03_Main_MFA_analysis_and_plots.R for the IL4 condition
