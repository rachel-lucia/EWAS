Glyphosate Analysis Workflow
Folder: Glyphosate_Methylation_Data_Analysis in MEE/Data/Methylation Data

Main Analysis

These initial steps are the same for all analyses with all 400 samples
1.	3/2/20 Load raw IDATs; filter out failed and low-quality probes; noob and BMIQ normalization (both are within-array)
	R script: LoadFilter.R
	Console output from successful cluster run: LoadFilterOut.txt
	R workspace from successful cluster run: LoadFilter.RData
	Raw beta matrix with detection p > 0.05 masked: RawBetas.csv
	Normalized/filtered beta matrix: NormalizedBetas.csv
	RGChannelSet object with loaded IDATs: RawRGSet.RData
2.	Sample QC steps (results are in MEE/Data/Methylation Data/Data Validation)
	2a. Verified sex inferred from methylation array matches reported sex - all 400 samples were predicted female
		R script: CheckSex.R
		Output: PredictedSex.csv
	2b. Verified SNP identity (no accidental duplications)
		R script: CheckSNPIdentity.R
		Output: SNPSimilarityMatrix.csv (% identity for every possible sample pair)
		The only samples that had >~50% identity were intentional duplicates (blood/saliva)
	2c. No samples had greater than 1% failed probes (see console output: LoadFilterOut.txt; max 0.126%)
	2d. No samples had obvious problems on standard QC/control probes
	
These steps will differ based on phenotype of interest
3.	ComBat correction
	Use champ.SVD to visualize which technical variables to worry about
	Need to specify which phenotype variables to protect
	
	Missing: HEI
	Redundant: Column = Array = Sentrix_Position; Slide = Sentrix_ID = Row [sort of: rows can repeat]
	Remove: Study_ID
	
	Plate: Group of 96
	
	Example script: BatchCorrectGP.R (for correction) BatchVisualizeGP.R (for visualizing which vars are important and if correction was effective)
4.	Modeling: basic + full models
	In MasterSampleSheet: glyphosate and AMPA substituted with 1/sqrt(2) of the LOD if <LOD

	Subsampling analysis: handled by 2 shell scripts
	Governor.sh submits RunParallelScript.sh repeatedly with the specified parameters (1:20 in the attached version)
	RunParallelScript.sh runs the analysis script (Parallel_Base_Script_GP.R) with the parameter as input
	Parallel_Base_Script_GP.R runs the analysis dividing 1000 iterations into N chunks. Read the comments for more info
	There are more proper way to do this but this is my "parallelization for dummies" approach
	
	BTW, Parallel_Base_Script.R requires a csv with the IDs for the members of each subsample to be already present. Generate this with GetSubsampleMembersGP.R
	
5.	Index development
	Index_CrossVal_GP_Clean.R generates and tests a methylation index using lasso regression. See comments for details


Etc.
Various other useful scripts scattered through the "Glyphosate_Methylation_Data_Analysis" and "Density_ReAnalysis" folders: plotting, scripts for combining overlapping regions, downstream analyses, etc. Most are *somewhat* documented
	
	

