This readme.txt file was generated on 2020-04-04 by MARIANO ALVAREZ


GENERAL INFORMATION

1. Title of Dataset: Spartina alterniflora oil spill epiGBS data

2. Author Information
	A. Principal Investigator Contact Information
		Name: Christina Richards
		Institution: University of South Florida
		Address: 
			4202 East Fowler Avenue SCA 127
			NES 107 (shipping)
			Tampa, FL 33620

		Email: clr@usf.edu

	B. Alternate Contact Information
		Name: Mariano Alvarez
		Institution: Duke University
		Address: 
			130 Science Drive, Box 90388
			Durham, NC 27708
		Email: mfa20@duke.edu

3. Date of data collection (single date, range, approximate date): Summer 2012

4. Geographic location of data collection : Gulf Coast, USA

5. Information about funding sources that supported the collection of the data: NSF


SHARING/ACCESS INFORMATION

1. Licenses/restrictions placed on the data: N/A

2. Links to publications that cite or use the data: https://doi.org/10.1101/426569

3. Links to other publicly accessible locations of the data: N/A

4. Links/relationships to ancillary data sets: github.com/alvarezmf

5. Was data derived from another source? NO
	A. If yes, list source(s): 

6. Recommended citation for this dataset: https://doi.org/10.1101/426569


DATA & FILE OVERVIEW

1. File List: 
consensus_cluster.renamed.fa (reference derived from epiGBS pipeline)
coordinates.txt (GPS coordinates of sampling locations)
env_spartina2.csv (sample names, sampling location, and oil presence)
medsummclean.txt (raw gene expression values from Alvarez et al. 2018)
methylation.bed (methylation polymorphisms, as called from the epiGBS pipeline)
pctlnorm_amr_sig.txt (significantly differentially expressed genes from Alvarez et al. 2018)
snp.vcf.gz (nucleotide polymorphisms, as called from the epiGBS pipeline)
Analysis.R (R script for analysis of epiGBS data)


2. Relationship between files, if important: All files are used by the Analysis.R script

3. Additional related data collected that was not included in the current data package:\ N/A

4. Are there multiple versions of the dataset? no


METHODOLOGICAL INFORMATION

1. Description of methods used for collection/generation of data: 
https://doi.org/10.1101/426569

2. Methods for processing the data: 
All scripts to process raw data are contained in the Analysis.R file, which contains function corresponding to each step of the analysis.

3. Instrument- or software-specific information needed to interpret the data: 
N/A

4. Standards and calibration information, if appropriate: N/A

5. Environmental/experimental conditions: Crude oil exposure

6. Describe any quality-assurance procedures performed on the data: Analysis.R contains read filtering steps and quality control checks.

7. People involved with sample collection, processing, analysis and/or submission: 

Mariano Alvarez, Marta Robertson, Thomas van Gurp, Niels Wagemaker, Delphine Giraud, Malika L. Ainouche, Armel Salmon, Koen J. F. Verhoeven, Christina L. Richards
