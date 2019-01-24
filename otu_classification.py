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
# Number of threads/cores to be used
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
	inDir = args.input_dir

# Main folder hosting the analysis
analysisDir = os.path.join(os.getcwd(), "otu_analysis")

# Main subfolder 
preprocessedFiles = os.path.join(analysisDir, "preprocessed_files")  # Save processed fastq files
reportsFolder = os.path.join(analysisDir, "reports")  # Reports folders
# Secondary subfolders
preprocessingReports = os.path.join(reportsFolder, "preprocessing_reports")
alignmentReports = os.path.join(reportsFolder, "alignment_reports")
trimmingReports = os.path.join(preprocessingReports, "trimming_reports")
temp = os.path.join(preprocessingReports, "temp")

# Generation of the folders
for files in [preprocessedFiles, preprocessingReports, alignmentReports, trimmingReports, temp]:
	if not os.path.exists(files): os.makedirs(files)




def mild_quality_trimming()
	
	 -a ADAPT1 -A ADAPT2 [options] -o out1.fastq -p out2.fastq in1.fastq in2.fastq
	print("{0}/{1}. Preprocessing {2}".format((i+1), len(raw_data), os.path.basename(files)))
	cutadapt = ' '.join([
	"cutadapt",  # Call Cutadapt to preprocess the raw data
	"--output_dir", preprocessedFiles,  # Output directory of processed files
	files])  # Output directory of FastQC reports
	subprocess.run(trimGalore, shell=True) 
	return 












def main():

	check_input_reads()
	
	# Mild quality trimming
	mild_quality_trimming()





if __name__ == "__main__": main()