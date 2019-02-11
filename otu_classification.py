#!/bin/bash

### Creator n' Maintainer Stavros Giannoukakos
### University of Granada

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

""" All the necessary folders that will host the analysis are being created. These 
include 'preprocessed_files' that will host the filtered and quality controlled 
data and reports, where all reports from all software will be stored. """
if not os.path.exists(args.input_dir):
	sys.exit('The given input folder does not exist...')
else:
	inputDir = args.input_dir

# Main folder hosting the analysis
analysisDir = os.path.join(args.output_dir if args.output_dir else os.getcwd(), "otu_analysis")

# Main subfolder 
preprocessedFiles = os.path.join(analysisDir, "preprocessed_files")  # Save processed fastq files
reportsDir = os.path.join(analysisDir, "reports")  # Reports folders
# Secondary subfolders
temp = os.path.join(preprocessedFiles, "temp")
qiimeDir = os.path.join(preprocessedFiles, "qiime_analysis")
filteredDir = os.path.join(preprocessedFiles, "filtered_data")
preprocessingReports = os.path.join(reportsDir, "preprocessing_reports")

# Generation of the folders
for files in [preprocessedFiles, preprocessingReports, temp, filteredDir]:
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
	forwardRead_output = os.path.join(temp, os.path.basename(forwardRead).replace([x for x in [".fastq.gz", ".fq.gz"]\
						 if forwardRead.endswith(x)][0], ".fastq.gz"))
	reverseRead_output = os.path.join(temp, os.path.basename(reverseRead).replace([x for x in [".fastq.gz", ".fq.gz"]\
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
	print(cutadapt)
	subprocess.run(cutadapt, shell=True)
	return 

def quality_control():
	""" Running fastQ  screen software to identify possible  contaminations in our samples. 
	Additionally, use afterQC to make a preliminary quality check of the processed PE reads. 
	Then MultiQC  will summarise the QC reports from  all samples into a summary report """
	# Obtaining the preprocessed reads
	mfiltered_data = ' '.join([f for f in glob.glob(os.path.join(temp, "*R1*.fastq.gz"))])

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

	print("Quality Control reports for the data are being generated: in progress ..")
	afterQC = ' '.join([
	"after.py",  # Call fastQC to quality control all processed data
	"--input_dir", temp,  # Input directory to be process automatically
	"--good_output_folder", filteredDir,  # Storing good reads
	"--bad_output_folder", filteredDir,  # Storing bad reads
	"--report_output_folder", preprocessingReports,  # Create all output files in this specified output directory
	"|", "tee", "--append", os.path.join(preprocessingReports, "afterQC_preQC_report.txt")])  # Output fastQC report
	subprocess.run(afterQC, shell=True) 

	multiQC = " ".join([
	"multiqc",  # Call MultiQC
	"--quiet",  # Print only log warnings
	"--outdir", preprocessingReports,  # Create report in the FastQC reports directory
	"--filename", "summarised_report",  # Name of the output report 
	preprocessingReports,  # Directory where all FastQC and Cutadapt reports reside
	"|", "tee", "--append", os.path.join(preprocessingReports, "multiQC_report.txt")])  # Output multiQC report
	subprocess.run(multiQC, shell=True)

	# os.system('rm {0}/*fastqc.zip'.format(preprocessingReports))  # Removing all 'fastqc.zip' temporary files
	# os.system("chmod 755 -R {0}".format(preprocessedFiles))	
	return

def denoidingAndMerning_reads():
	""" Importing and denoising the preprocessed PE reads """
	importingSamplesToQiime2 =	' '.join([
	"qiime tools import",  # Run QIIME IMPORT to import data and create a new QIIME 2 Artifact
	"--type", "\'SampleData[PairedEndSequencesWithQuality]\'",  # The semantic type of the artifact that will be created upon importing
	"--input-format", "CasavaOneEightSingleLanePerSampleDirFmt",
	"--input-path", temp,  # Path to the directory that should be imported
	"--output-path", os.path.join(qiimeDir, "mqualitrim_noprim_pe_data.qza")])  # Path where output artifact should be written
	subprocess.run(importingSamplesToQiime2, shell=True)

	denoising_n_merging = ' '.join([
	"qiime dada2 denoise-paired",  # Call qiime dada2 to denoise the preprocessed data
	"--p-n-threads", str(args.threads),  # Number of threads to use
	"--p-trunc-len-f", "0",  # No truncation will be performed cause we have already trimmed low quality ends
	"--p-trunc-len-r", "0",  # No truncation will be performed cause we have already trimmed low quality ends
	"--output-dir", os.path.join(qiimeDir, "denoised_results"),  # Output results to a directory
	"--i-demultiplexed-seqs", os.path.join(qiimeDir, "mqualitrim_noprim_pe_data.qza"),  # The paired-end demultiplexed sequences to be denoised
	"|", "tee", os.path.join(preprocessingReports, "dada2_denoising_report.txt")])  # Output denoising report
	# subprocess.run(denoising_n_merging, shell=True)

	# qiime feature-table filter-features 
	# –i-table trimmed/DADA2/table-dada2.qza 
	# –p-min-frequency 10 
	# –o-filtered-table trimmed/DADA2/dada2-feature-frequency-filtered-table.qza
	######################
	#### add visualisation and extraction of stats
	######################

	return

def quiime2():

	# To activate this environment, use
	#
	#     $ conda activate qiime2-2019.1
	#
	# To deactivate an active environment, use
	#
	#     $ conda deactivate
	return 

def main():

	pairedReads = assess_input_data(inputDir)
	# Obtaining the number of pair files
	totNum = len(pairedReads)

	## Preprocessing of the input data
	"""
	# Performing mild quality trimming and 
	# removal of all primers on both reads
	for i, read in enumerate(pairedReads):
		mildQualityTrimming_primerRemoval(read[0], read[1], i, totNum) 

	quality_control()  # Checking the quality of the merged reads
	"""
	denoidingAndMerning_reads()  # Merging the pair files

	
	
	## OTU analysis 

	## Downstream analysis

if __name__ == "__main__": main()