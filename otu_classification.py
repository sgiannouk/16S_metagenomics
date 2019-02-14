#!/bin/bash

### Creator n' Maintainer Stavros Giannoukakos
### University of Granada

# To activate this environment, use
#
#     $ conda activate qiime2
#
# To deactivate an active environment, use
#
#     $ conda deactivate

#Version of the program
__version__ = "0.1.0"

import argparse
import subprocess
from Bio.Seq import Seq
from datetime import datetime
import shutil, fnmatch, glob, sys, os

# Configuration file needed for FastQ Screen
fastQscreen_config = "/home/stavros/playground/16S_metagenomics/subsidiary_files/fastq_screen.conf"

# Tracking time of analysis
start_time = datetime.now()

usage = "otu_classification [options] -i <input_directory/input_files>"
epilog = " -- January 2019 | Stavros Giannoukakos -- "
description = "DESCRIPTION\
           		\n-----------"

parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, usage=usage, description=description, epilog=epilog)
# Create required section in help
requiredArgs = parser.add_argument_group('required arguments')
# Input folder option
requiredArgs.add_argument('-i', '--input_dir', required=True, metavar='', 
						   help="Path of the input directory that contains the raw data.\nBoth forward and reverse reads are expected to be found\nin this directory.")
# Number of threads/CPUs to be used
parser.add_argument('-th', '--threads', dest='threads', default=50, metavar='', 
                	help="Number of threads to be used in the analysis")
# Number of threads/CPUs to be used
parser.add_argument('-fp', '--forwardPrimer', default="GTGCCAGCMGCCGCGGTAA", required=False, metavar='', 
                	help="Sequence of the forward primer")
# Number of threads/CPUs to be used
parser.add_argument('-rp', '--reversePrimer', default="GGACTACHVGGGTWTCTAAT", required=False, metavar='', 
                	help="Sequence of the reverse primer")
# Display the version of the pipeline 
parser.add_argument('-v', '--version', action='version', version='%(prog)s {0}'.format(__version__))
# Path of where the output folder should be located
parser.add_argument('-o', '--output_dir', metavar='', 
					help="Path of the directory that will host the analysis.\n(default: <current_directory>)")

# Get the options and return them
args = parser.parse_args()

""" All the necessary directory that will host the analysis are being created. These 
include 'preprocessed_files' that will host the filtered and quality controlled 
data and reports, where all reports from all software will be stored. """
if not os.path.exists(args.input_dir):
	sys.exit('The given input folder does not exist...')
else:
	inputDir = args.input_dir

# Main folder hosting the analysis
analysisDir = os.path.join(args.output_dir if args.output_dir else os.getcwd(), "otu_analysis")

# Main subdirectories
reportsDir = os.path.join(analysisDir, "reports")  # Reports directory
qiimeDir = os.path.join(analysisDir, "qiime_analysis")  # Directory hosting the main analysis 
preprocessedFiles = os.path.join(analysisDir, "preprocessed_files")  # Save processed .fastq files

# Secondary subdirectories
temp = os.path.join(preprocessedFiles, "temp")
filteredDir = os.path.join(preprocessedFiles, "filtered_data")

qiimeResults = os.path.join(qiimeDir, "denoising_results")

preprocessingReports = os.path.join(reportsDir, "preprocessing_reports")
analysisReports = os.path.join(reportsDir, "analysis_reports")

# Generation of the directories
for files in [temp, qiimeDir, filteredDir, preprocessingReports, analysisReports]:
	if not os.path.exists(files): os.makedirs(files)



def assess_input_data(input_directory):
	""" In this function the PE input data will be assesses for valid format and 
	for correct pairing. """
	input_files = []  # Output list that will contain the paired-input files
	for path, subdirs, files in os.walk(input_directory):
		for name in files:
			# Checking the format of the reads
			if not name.endswith((".fastq.gz", ".fq.gz")):  
				sys.exit('Unidentified input format in read: {0}'.format(name))
			# Verifying that all reads have their pairs
			elif not any(x in name.upper() for x in ["_R1_", "_R2_"]):
				sys.exit('Unidentified member of a pair in read: {0}'.format(name))  
			elif "_R1_" in name:  # Obtaining the paired-input files
				inR1 = os.path.join(os.path.abspath(path), name)
				inR2 = inR1.replace("_R1_","_R2_")
				assert (os.path.isfile(inR2)), 'Could not detect the pair of {0} ({1})'.format(inR1, inR2)
				input_files.append((inR1, inR2))
	return input_files

def mildQualityTrimming_primerRemoval(forwardRead, reverseRead, i, totNum):
	""" An initial very mild base quality trimming will be performed. In this step, we are trying to 
	discard very troublesome bases (whos quality is below Q15). That way we remove obvious trash and 
	trying to improve the chances of a proper merge. """
	forwardRead_output = os.path.join(filteredDir, os.path.basename(forwardRead).replace([x for x in [".fastq.gz", ".fq.gz"]\
						 if forwardRead.endswith(x)][0], ".fastq.gz"))
	reverseRead_output = os.path.join(filteredDir, os.path.basename(reverseRead).replace([x for x in [".fastq.gz", ".fq.gz"]\
						 if reverseRead.endswith(x)][0], ".fastq.gz"))
	
	print("{0}/{1}. Mild base quality filtering in {2}".format((i+1), totNum, os.path.basename(forwardRead.split("_")[-1])))
	cutadapt = ' '.join([
	"cutadapt",  # Call Cutadapt to preprocess the raw data
	"--cores", str(args.threads),  # Number of CPUs to use
	"--max-n", "0",  # Discard reads with 'N' bases. Dada2 cannot comprehend Ns. 
	"--trim-n",  # Trim N's on ends of reads
	"--no-indels",  # Not allowing indels in the alignments
	"-m", "20",  # Discard reads shorter than 20 bases
	"--quality-cutoff", "15",  # Trim low-quality bases from 3' end of each read (Q<15)
	"--overlap", str(len(args.forwardPrimer)-3), # Min overlap between read and adapter for an adapter to be found
	"--discard-untrimmed",  # Discard reads that do not contain a primer
	"--output", forwardRead_output,  # Export edited forward read to file
	"--paired-output", reverseRead_output,  # Export edited reverse read to file
	"-a", "{0}...{1}".format(args.forwardPrimer, Seq(args.reversePrimer).reverse_complement()),  # R1 linked adapter FWDPRIMER...RCREVPRIMER
	"-A", "{0}...{1}".format(args.reversePrimer, Seq(args.forwardPrimer).reverse_complement()),  # R2 linked adapter REVPRIMER...RCFWDPRIMER
	forwardRead,  # Input of the forward file
	reverseRead,  # Input of the reverse file
	"|", "tee", "--append", os.path.join(preprocessingReports, "cutadapt_mildQtrimNoPrimers_report.txt")])  # Output trimming report
	subprocess.run(cutadapt, shell=True)
	return 

def quality_control():
	""" Running FastQ  Screen software to identify possible  contaminations in our samples. 
	Additionally, use AfterQC to make a preliminary quality check of the processed PE reads. 
	Then MultiQC  will summarise the QC reports from  all samples into a summary report """
	# Obtaining the preprocessed reads
	mfiltered_data = ' '.join([f for f in glob.glob(os.path.join(filteredDir, "*R1*.fastq.gz"))])

	print("Checking random reads for possible contamination: in progress ..")
	fastQscreen = ' '.join([
	"fastq_screen",  # Call fastQ screen to check contamination in the processed data
	"--threads", str(args.threads),  # Number of threads to use
	"--outdir",  preprocessingReports,  # Directory in which the output files will be saved
	"--quiet",  # Suppress all progress reports on stderr and only report errors
	"--conf", fastQscreen_config,  # Location of the required configuration file
	mfiltered_data, mfiltered_data.replace("_R1_", "_R2_"),  # Input PE files
	"|", "tee", "--append", os.path.join(preprocessingReports, "fastQscreen_report.txt")])  # Output fastQ screen report
	subprocess.run(fastQscreen, shell=True)

	print("Quality Control reports for the following data are being generated: in progress ..")
	afterQC = ' '.join([
	"after.py",  # Call fastQC to quality control all processed data
	"--input_dir", filteredDir,  # Input directory to be process automatically
	"--good_output_folder", temp,  # Storing good reads
	"--bad_output_folder", temp,  # Storing bad reads
	"--report_output_folder", preprocessingReports,  # Create all output files in this specified output directory
	"|", "tee", "--append", os.path.join(preprocessingReports, "afterQC_report.txt")])  # Output fastQC report
	subprocess.run(afterQC, shell=True) 

	multiQC = " ".join([
	"multiqc",  # Call MultiQC
	"--quiet",  # Print only log warnings
	"--outdir", preprocessingReports,  # Create report in the FastQC reports directory
	"--filename", "summarised_report",  # Name of the output report 
	preprocessingReports,  # Directory where all FastQC and Cutadapt reports reside
	"|", "tee", "--append", os.path.join(preprocessingReports, "multiQC_report.txt")])  # Output multiQC report
	subprocess.run(multiQC, shell=True)

	os.system('mv {0}/*_report.txt {1}'.format(preprocessingReports, reportsDir))  # Moving all reports in the reports folder
	# os.system('rm -r {0}/*_report_data'.format(preprocessingReports))  # Removing MultiQC temporary folder
	# os.system("chmod 755 -R {0}".format(preprocessedFiles))	
	return

def otu_mainAnalysis():
	""" Importing and denoising the preprocessed PE reads """
	# print("Importing the preprocessed reads to the Qiime2 Artifact: in progress ..")
	importingSamplesToQiime2 =	' '.join([
	"qiime tools import",  # Run QIIME IMPORT to import data and create a new QIIME 2 Artifact
	"--type", "\'SampleData[PairedEndSequencesWithQuality]\'",  # The semantic type of the artifact that will be created upon importing
	"--input-format", "CasavaOneEightSingleLanePerSampleDirFmt",
	"--input-path", filteredDir,  # Path to the directory that should be imported
	"--output-path", os.path.join(qiimeDir, "input_data.qza"),  # Path where output artifact should be written
	"|", "tee", os.path.join(reportsDir, "qiime2_importingData_report.txt")])  # Output denoising report
	# subprocess.run(importingSamplesToQiime2, shell=True)


	""" Denoising is an attempt to correct reads with sequencing errors and then 
	remove chimeric sequences originating from different DNA templates."""
	# print("Denoising, dereplicating and filtering chimera sequences from the paired-end data: in progress ..")
	denoisingNmerging = ' '.join([
	"qiime dada2 denoise-paired",  # Call qiime dada2 to denoise the preprocessed data
	"--p-n-threads", str(args.threads),  # Number of threads to use
	"--verbose",  # Display verbose output to stdout and/or stderr during execution
	"--p-trunc-len-f", "0",  # No truncation will be performed cause we have already trimmed low quality ends
	"--p-trunc-len-r", "0",  # No truncation will be performed cause we have already trimmed low quality ends
	"--output-dir", qiimeResults,  # Output results to a directory
	"--i-demultiplexed-seqs", os.path.join(qiimeDir, "input_data.qza"),  # The paired-end demultiplexed sequences to be denoised
	"|", "tee", os.path.join(reportsDir, "dada2_denoising_report.txt")])  # Output denoising report
	# subprocess.run(denoisingNmerging, shell=True)

	""" At this stage, we have obtained the artifacts containing the feature table and corresponding feature sequences. 
	Now we will generate summary of the above features and proceed with visualisation of the data. """
	featureTableSummary = ' '.join([
	"qiime feature-table summarize",  # Calling qiime2 feature-table summarize function
	"--verbose",  # Display verbose output to stdout and/or stderr during execution
	"--i-table", os.path.join(qiimeResults, "table.qza"),  # The feature table to be summarized
	"--o-visualization", os.path.join(qiimeDir, "featureTable.qzv"),  # Output results to directory
	# "--m-sample-metadata-file", os.path.join(qiimeDir, "sample-metadata.tsv"),  # Metadata file
	"|", "tee", os.path.join(preprocessingReports, "qiime_featureTableSummary_report.txt")])  # Output featureTableSummary report
	# subprocess.run(featureTableSummary, shell=True)
	
	featureIdentifierToSeqMapping = ' '.join([
	"qiime feature-table tabulate-seqs",  # Calling qiime2 feature-table tabulate-seqs function
	"--verbose",  # Display verbose output to stdout and/or stderr during execution
	"--i-data", os.path.join(qiimeResults, "representative_sequences.qza"),  # The feature sequences to be tabulated
	"--o-visualization", os.path.join(qiimeResults, "representative_sequences.qzv"),  # Output file
	"|", "tee", os.path.join(preprocessingReports, "qiime_featureIdentifierToSeqMapping_report.txt")])  # Output featureIdentifierToSeqMapping report
	# subprocess.run(featureIdentifierToSeqMapping, shell=True)
	

	# print("Visualisation of the basic statistics regarding the denoising process: in progress ..")
	visualisationOfDenoisingStats = ' '.join([
	"qiime metadata tabulate",  # Calling qiime2 metadata tabulate function
  	"--m-input-file", os.path.join(qiimeResults, "denoising_stats.qza"),  # Metadata input file
  	"--o-visualization", os.path.join(qiimeResults, "denoising-stats.qzv"),  # Output visualization file 
  	"|", "tee", os.path.join(reportsDir, "dada2_visualisationOfDenoisingStats_report.txt")])  # Output denoising report
	# subprocess.run(visualisationOfDenoisingStats, shell=True)
	


	######################
	#### add visualisation and extraction of stats
	#####################

	""" Applying basic filters for how frequent a variant needs to be """
	filteringVariant = ' '.join([
	"qiime feature-table filter-features",  # Calling qiime2 feature-table filter-features function
	"--verbose",  # Display verbose output to stdout and/or stderr during execution
	"-i-table", os.path.join(qiimeResults, "table.qza"),  # Input feature table
	"-p-min-frequency", "10"  # Minimum observed frequency of a 
	"-o-filtered-table", os.path.join(qiimeResults, "feature_frequencyFiltered_table.qza"),  #
	"|", "tee", os.path.join(reportsDir, "qiime2_filteringVariant_report.txt")])  # Output denoising report
	# subprocess.run(filteringVariant, shell=True)

	""" Generate a heatmap representation of the filtered feature table """
	filteredFeaturesHeatmap = ' '.join([
	"qiime feature-table heatmap",  # Calling qiime2 feature-table heatmap function
	"--verbose",  # Display verbose output to stdout and/or stderr during execution
	"--p-color-scheme", "Paired",  # The color scheme of the heatmap  
	"--i-data", os.path.join(qiimeResults, "feature_frequencyFiltered_table.qza"),  # The input filtered feature table
	"--output-dir", qiimeResults,  # Output directory
	"|", "tee", os.path.join(preprocessingReports, "qiime_filteredFeaturesHeatmap_report.txt")])  # Output filteredFeaturesHeatmap report
	# subprocess.run(filteredFeaturesHeatmap, shell=True)
	return 

def main():

	pairedReads = assess_input_data(inputDir)
	# Obtaining the number of pair files
	totNum = len(pairedReads)

	"""
	## Preprocessing of the input data
	# Performing mild quality trimming and 
	# removal of all primers on both reads
	for i, read in enumerate(pairedReads):
		mildQualityTrimming_primerRemoval(read[0], read[1], i, totNum) 

	quality_control()  # Checking the quality of the merged reads
	"""

	## OTU analysis 
	otu_mainAnalysis()  # This function hosts the main otu analysis
	

	## Downstream analysis

if __name__ == "__main__": main()