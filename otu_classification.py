#!/usr/bin/python3

### Creator n' Maintainer Stavros Giannoukakos
### University of Granada

#Version of the program
__version__ = "0.1.0"

import argparse
import subprocess
import shutil, sys, os
from datetime import datetime

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
parser.add_argument('-th', '--threads', dest='threads', default=1, metavar='', 
                	help="Number of threads to be used in the analysis")
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
reportsFolder = os.path.join(analysisDir, "reports")  # Reports folders
# Secondary subfolders
preprocessingReports = os.path.join(reportsFolder, "preprocessing_reports")
trimmingReports = os.path.join(preprocessingReports, "trimming_reports")
temp = os.path.join(preprocessedFiles, "temp")

# Generation of the folders
for files in [preprocessedFiles, preprocessingReports, trimmingReports, temp]:
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

def mild_quality_trimming(forwardRead, reverseRead, i, totNum):
	""" An initial very mild base quality trimming will be performed. In this step, we are trying to 
	discard very troublesome bases (whos quality is below Q18). That way we remove obvious trash and 
	trying to improve the chances of a proper merge. """
	print("{0}/{1}. Mild base quality filtering in {2}".format((i+1), totNum, os.path.basename(forwardRead.split("_")[0])))
	forwardRead_output = os.path.join(temp, os.path.basename(forwardRead).replace([x for x in [".fastq.gz", ".fq.gz"]\
						 if forwardRead.endswith(x)][0], "_mtrim.fq.gz"))
	reverseRead_output = os.path.join(temp, os.path.basename(reverseRead).replace([x for x in [".fastq.gz", ".fq.gz"]\
						 if reverseRead.endswith(x)][0], "_mtrim.fq.gz"))
	cutadapt = ' '.join([
	"cutadapt",  # Call Cutadapt to preprocess the raw data
	"--cores", str(args.threads),  # Number of CPU cores to use
	"--quality-cutoff", "18",  # Trim low-quality bases from 5' end of each read ()
	"--max-n", "0.20",  # Discard reads with more than 20% 'N' bases
	"--trim-n",  # Trim N's on ends of reads
	"--output", forwardRead_output,  # Export edited forward read to file
	"--paired-output", reverseRead_output,  # Export edited reverse read to file
	forwardRead,  # Input of the forward file
	reverseRead,  # Input of the reverse file
	"|", "tee", os.path.join(trimmingReports, "cutadapt_mildQtrim_report.txt")])  # Output trimming report
	subprocess.run(cutadapt, shell=True)
	return 

def pairEndMerge():


	return


def main():

	pairedReads = assess_input_data(inputDir)
	# Obtaining the number of pair files
	totNum = len(pairedReads)

	## Preprocessing of the input data
	
	# Performing mild quality trimming
	for i, read in enumerate(pairedReads):
		mild_quality_trimming(read[0], read[1], i, totNum) 

	pairEndMerge()  # Merging the pair files

	## OTU analysis 

	## Downstream analysis

if __name__ == "__main__": main()