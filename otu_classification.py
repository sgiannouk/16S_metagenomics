#!/bin/bash

### Creator n' Maintainer Stavros Giannoukakos
### University of Granada

# To activate this environment, use
#
#     $ conda activate /opt/anaconda3/envs/qiime2
#
# To deactivate an active environment, use
#
#     $ conda deactivate

#Version of the program
__version__ = "0.1.1"

import argparse
from Bio import SeqIO
import subprocess, gzip
from Bio.Seq import Seq
from datetime import datetime
import shutil, fnmatch, glob, sys, os

inputDir = '/shared/projects/martyna_rrna_illumina_igtp/batch2/'
metadata_file = "/shared/projects/martyna_rrna_illumina_igtp/batch2/igtp_batch2_metadata.txt"

# Configuration file needed for FastQ Screen
fastQscreen_config = "/home/stavros/references/fastQscreen_references/fastq_screen.conf"
silva_99_classifier = "/home/stavros/playground/progs/16S_subsidiary_files/silva_132_99_v3v4.qza"

# Tracking time of analysis
start_time = datetime.now()

usage = "otu_classification [options] -i <input_directory/input_files>"
epilog = " -- January 2019 | Stavros Giannoukakos -- "
description = "DESCRIPTION"

parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, usage=usage, description=description, epilog=epilog)
# Create required section in help
requiredArgs = parser.add_argument_group('required arguments')
# Input folder option
requiredArgs.add_argument('-i', '--input_dir', required=False, metavar='', 
						   help="Path of the input directory that contains the raw data.\nBoth forward and reverse reads are expected to be found\nin this directory.")
# Type of analysis
parser.add_argument('-a', '--analysis_type', default='16S', choices=['16S', 'ITS'], metavar='16S/ITS',
                	help="Type of analysis to be performed. Bacteria analysis(16S) or Fungi analysis (ITS).")

# Number of threads/CPUs to be used
parser.add_argument('-th', '--threads', dest='threads', default=30, metavar='', 
                	help="Number of threads to be used in the analysis")
# Number of threads/CPUs to be used
parser.add_argument('-fp', '--forwardPrimer', default="CCTACGGGNGGCWGCAG", required=False, metavar='', 
                	help="Sequence of the forward primer")
# Number of threads/CPUs to be used
parser.add_argument('-rp', '--reversePrimer', default="GACTACHVGGGTATCTAATCC", required=False, metavar='', 
                	help="Sequence of the reverse primer")
# Metadata file
parser.add_argument('-m', '--metadata', required=False, metavar='', 
                	help="Metadata file containing several info conserning the\ninput data")
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
# if not os.path.exists(args.input_dir):
# 	sys.exit('The given input folder does not exist...')
# else:
# 	inputDir = args.input_dir

args.metadata = metadata_file
args.input_dir = inputDir

# Main folder hosting the analysis
analysisDir = os.path.join(args.output_dir if args.output_dir else os.getcwd(), "batch2_analysis")
# Main subdirectories
reportsDir = os.path.join(analysisDir, "reports")  # Reports directory
qiimeDir = os.path.join(analysisDir, "qiime2_analysis")  # Directory hosting the main analysis 
preprocessedFiles = os.path.join(analysisDir, "preprocessed_files")  # Save processed .fastq files
qiimeResults = os.path.join(qiimeDir, "denoising_analysis")
diversityAnalysis = os.path.join(qiimeDir, "diversity_analysis")
statisticalAnalysis = os.path.join(diversityAnalysis, "statistical_analysis")
taxinomicAnalysis = os.path.join(qiimeDir, "taxinomic_analysis")
preprocessingReports = os.path.join(reportsDir, "preprocessing_reports")



def assess_input_data(input_directory):
	""" In this function the PE input data will be assesses for valid format and 
	for correct pairing. """
	input_files = []  # Output list that will contain the paired-input files
	for path, subdirs, files in os.walk(input_directory):
		for name in files:
			# Verifying that all reads have their pairs
			if name.endswith((".fastq.gz", ".fq.gz")) and not any(x in name.upper() for x in ["_R1_", "_R2_"]):
				sys.exit('Unidentified member of a pair in read: {0}'.format(name))  
			elif name.endswith((".fastq.gz", ".fq.gz")) and "_R1_" in name:  # Obtaining the paired-input files
				inR1 = os.path.join(os.path.abspath(path), name)
				inR2 = inR1.replace("_R1_","_R2_")
				assert (os.path.isfile(inR2)), 'Could not detect the pair of {0} ({1})'.format(inR1, inR2)
				input_files.append((inR1, inR2))
	return input_files

def primer_removal(forwardRead, reverseRead, i, totNum):
	""" An initial very mild base quality trimming will be performed. In this step, we are trying to 
	discard very troublesome bases (whos quality is below Q20). That way we remove obvious trash and 
	trying to improve the chances of a proper merge. """
	[os.makedirs(files) for files in [preprocessedFiles, preprocessingReports] if not os.path.exists(files)]  # Generation of the directories

	forwardRead_output = os.path.join(preprocessedFiles, os.path.basename(forwardRead).replace([x for x in [".fastq.gz", ".fq.gz"]\
						 if forwardRead.endswith(x)][0], ".fastq.gz"))
	reverseRead_output = os.path.join(preprocessedFiles, os.path.basename(reverseRead).replace([x for x in [".fastq.gz", ".fq.gz"]\
						 if reverseRead.endswith(x)][0], ".fastq.gz"))
	# Obtaining the sample name
	sample_name = os.path.basename(forwardRead).split("_")[0]

	print("{0}/{1}. BBDUK - Quality filtering and primer trimming of {2}".format(i, totNum, sample_name))
	# Calculating the minimum and maximum acceptable length after primer and quality trimming
	minlength, maxlength = calculateMinMax(forwardRead)
	bbduk = ' '.join([
	"/home/stavros/anaconda3/bin/bbduk.sh",  # Call BBDuck (BBTools) to preprocess the raw data
	"threads={0}".format(str(args.threads)),  # Set number of threads to use
	"in={0}".format(forwardRead),  # Input of the forward file
	"in2={0}".format(reverseRead),  # Input of the reverse file
	"out={0}".format(forwardRead_output),  # Export edited forward read to file
	"out2={0}".format(reverseRead_output),  # Export edited reverse read to file
	"literal={0},{1}".format(args.forwardPrimer, args.reversePrimer),  # Providing the forward and reverse primers
	"ordered=t",  # Keeps the reads in the same order as we gave them to the software
	"minlength={0}".format(minlength), #220 Pairs (or reads) will be discarded if both are shorter than this after trimming
	"maxlength={0}".format(maxlength), # ~ 560 Pairs (or reads) will be discarded only if both are longer than this after trimming
	"k=10",  # Setting the kmer size we want to search for
	"mm=f",  # Looking for exact kmers and not mask the middle bases
	"mink=2",  # Specifies the smallest word size it will check against either edge of a read
	"rcomp=t",  # Look for reverse-complements of kmers in addition to forward kmers
	"copyundefined=t",  # Process non-AGCT IUPAC reference bases by making all possible unambiguous copies.
	"ktrim=l",  # Trim to the left of reads to remove bases matching reference kmers
	"trimq=20",  # Regions with average quality BELOW this will be trimmed
	"qtrim=r",  # Trim read ends (left end only) to remove bases with quality below trimq
	"2>>", os.path.join(preprocessingReports, "bbduk_report.log")])  # Output trimming report
	print(bbduk, file=open(os.path.join(preprocessingReports, "bbduk_report.log"),"a"))
	subprocess.run(bbduk, shell=True)
	return 

def quality_control():
	""" Running FastQ  Screen software to identify possible  contaminations in our samples. 
	Additionally, use AfterQC to make a preliminary quality check of the processed PE reads. 
	Then MultiQC  will summarise the QC reports from  all samples into a summary report """
	print("\n\t{0} QUALITY CONTROL OF THE PREPROCESSED SAMPLES".format(datetime.now().strftime("%d.%m.%Y %H:%M")))
	
	# Obtaining the preprocessed reads
	mfiltered_data = ' '.join([f for f in glob.glob(os.path.join(preprocessedFiles, "*_R1_*.fastq.gz"))])

	print("{0}  FastQScreen - Checking random reads for possible contamination: in progress ..".format(datetime.now().strftime("%d.%m.%Y %H:%M")))
	fastQscreen = ' '.join([
	"fastq_screen",  # Call fastQ screen to check contamination in the processed data
	"--threads", str(args.threads),  # Number of threads to use
	"--outdir",  preprocessingReports,  # Directory in which the output files will be saved
	"--quiet",  # Suppress all progress reports on stderr and only report errors
	"--conf", fastQscreen_config,  # Location of the required configuration file
	mfiltered_data, mfiltered_data.replace("_R1_", "_R2_"),  # Input PE files
	"2>>", os.path.join(preprocessingReports, "fastQscreen_report.log")])  # Output fastQ screen report
	subprocess.run(fastQscreen, shell=True)

	print("{0}  FastQC - Quality Control reports for the forward reads are being generated: in progress ..".format(datetime.now().strftime("%d.%m.%Y %H:%M")))
	fastQC_frw = ' '.join([
	"fastqc",  # Call fastQC to quality contol all processed data
	"--threads", str(args.threads),  # Number of threads to use
	"--quiet",  # Print only log warnings
	"--nogroup",  # Disable grouping of bases for reads >50bp
	"--outdir", preprocessingReports,  # Create all output files in this specified output directory
	mfiltered_data, mfiltered_data.replace("_R1_", "_R2_"),  # String containing all samples that are about to be checked
	"2>>", os.path.join(preprocessingReports, "fastQC_report.log")])  # Output fastQC report
	subprocess.run(fastQC_frw, shell=True)

	print("{0}  FastP - Quality Control reports for the input reads are being generated: in progress ..".format(datetime.now().strftime("%d.%m.%Y %H:%M")))
	for files in glob.glob(os.path.join(preprocessedFiles, "*_R1_*.fastq.gz")):
		fastP = ' '.join([
		"fastp",  # Call fastQC to quality control all processed data
		"--thread", str(args.threads),  # Number of threads to use
		"--in1", files,  # Input read1 file
		"--in2", files.replace("_R1_", "_R2_"),  # Input read2 file
		"--disable_adapter_trimming",  # Adapter trimming is disabled
		"--disable_quality_filtering",  # Quality filtering is disabled
		"--disable_length_filtering",  # Length filtering is disabled
		"--overrepresentation_analysis",  # Enable overrepresented sequence analysis
		"--html", os.path.join(preprocessingReports, "{0}_fastp.html".format(os.path.basename(files)[:-9])),  # Create ftml file in this specified output directory
		"--json", os.path.join(preprocessingReports, "{0}_fastp.json".format(os.path.basename(files)[:-9])),  # Create json output file in this specified output directory
		"2>>", os.path.join(preprocessingReports, "fastP_report.log")])  # Output fastP report
		subprocess.run(fastP, shell=True) 

	print("{0}  MultiQC - Summarised Quality Control report for all the input reads is being generated: in progress ..".format(datetime.now().strftime("%d.%m.%Y %H:%M")))
	multiQC = " ".join([
	"multiqc",  # Call MultiQC
	"--quiet",  # Print only log warnings
	"--outdir", preprocessingReports,  # Create report in the FastQC reports directory
	"--filename", "summarised_report",  # Name of the output report 
	preprocessingReports,  # Directory where all FastQC and Cutadapt reports reside
	"2>>", os.path.join(preprocessingReports, "multiQC_report.log")])  # Output multiQC report
	subprocess.run(multiQC, shell=True)

	os.system('mv {0}/*_report.log {1}'.format(preprocessingReports, reportsDir))  # Moving all reports in the reports folder
	os.system('mv {0}/summarised_report.html {1}'.format(preprocessingReports, reportsDir))  # Moving summary report in the reports folder
	# os.system('rm -r {0}/*_report_data'.format(preprocessingReports))  # Removing MultiQC temporary folder
	return

def qiime2_analysis():
	print("\n\t{0} QIIME2 ANALYSIS".format(datetime.now().strftime("%d.%m.%Y %H:%M")))

	if not os.path.exists(qiimeDir): os.makedirs(qiimeDir)  # Creating the "qiime2_analysis" directory
	""" Importing and denoising the preprocessed PE reads """
	print("{0} Importing the preprocessed reads to the Qiime2 Artifact: in progress ..".format(datetime.now().strftime("%d.%m.%Y %H:%M")))
	importingSamplesToQiime2 =	' '.join([
	"qiime tools import",  # Run QIIME IMPORT to import data and create a new QIIME 2 Artifact
	"--type", "\'SampleData[PairedEndSequencesWithQuality]\'",  # The semantic type of the artifact that will be created upon importing
	"--input-format", "CasavaOneEightSingleLanePerSampleDirFmt",
	"--input-path", preprocessedFiles,  # Path to the directory that should be imported
	"--output-path", os.path.join(qiimeDir, "input_data.qza"),  # Path where output artifact should be written
	"2>>", os.path.join(reportsDir, "qiime2_importingData_report.log")])  # Output denoising report
	subprocess.run(importingSamplesToQiime2, shell=True)

	print("{0} Generating interactive positional quality plots: in progress ..".format(datetime.now().strftime("%d.%m.%Y %H:%M")))
	importSamplesQC =	' '.join([
	"qiime demux summarize",  # Calling qiime demux summarize to quality of each sample
	"--quiet",  # Silence output if execution is successful
	"--i-data", os.path.join(qiimeDir, "input_data.qza"),  # Path where the input artifact is written
	"--o-visualization", os.path.join(qiimeDir, "inputData_QC.qzv"),  # Output reports
	"2>>", os.path.join(reportsDir, "qiime2_importSamplesQC_report.log")])  # Output importSamplesQC report
	subprocess.run(importSamplesQC, shell=True)
	export(os.path.join(qiimeDir, "inputData_QC.qzv"))

	""" Denoising is an attempt to correct reads with sequencing errors and then 
	remove chimeric sequences originating from different DNA templates. """
	print("{0} Denoising, dereplicating and filtering chimera sequences from the paired-end data: in progress ..".format(datetime.now().strftime("%d.%m.%Y %H:%M")))
	denoisingNmerging = ' '.join([
	"qiime dada2 denoise-paired",  # Call qiime dada2 to denoise the preprocessed data
	"--p-n-threads", str(args.threads),  # Number of threads to use
	"--quiet",  # Silence output if execution is successful
	"--p-trunc-len-f", "0",  # No truncation will be performed cause we have already trimmed low quality ends
	"--p-trunc-len-r", "0",  # No truncation will be performed cause we have already trimmed low quality ends
	"--output-dir", qiimeResults,  # Output results to a directory
	"--i-demultiplexed-seqs", os.path.join(qiimeDir, "input_data.qza"),  # The paired-end demultiplexed sequences to be denoised
	"2>>", os.path.join(reportsDir, "dada2_denoising_report.log")])  # Output denoising report
	subprocess.run(denoisingNmerging, shell=True)

	print("{0} Visualising the feature-table: in progress ..".format(datetime.now().strftime("%d.%m.%Y %H:%M")))
	featureTableSummary = ' '.join([
	"qiime feature-table summarize",  # Calling qiime2 feature-table summarize function
	"--m-sample-metadata-file", args.metadata,  # Metadata file
	"--quiet",  # Silence output if execution is successful
	"--i-table", os.path.join(qiimeResults, "table.qza"),  # The feature table to be summarized
	"--o-visualization", os.path.join(qiimeResults, "feature_table.qzv"),  # Output results to directory
	"2>>", os.path.join(preprocessingReports, "qiime2_featureTableSummary_report.log")])  # Output featureTableSummary report
	subprocess.run(featureTableSummary, shell=True)
	export(os.path.join(qiimeResults, "feature_table.qzv"))

	print("{0} Visualising the denoising stats: in progress ..".format(datetime.now().strftime("%d.%m.%Y %H:%M")))
	denoisingStats = ' '.join([
	"qiime metadata tabulate",  # Calling qiime2 metadata tabulate
	"--m-input-file", os.path.join(qiimeResults, "denoising_stats.qza"),  # Input stats file
	"--o-visualization", os.path.join(qiimeResults, "denoising_stats.qzv"),  # Output results to directory
	"2>>", os.path.join(preprocessingReports, "qiime2_denoisingVisualStats_report.log")])  # Output denoisingVisualStats report
	subprocess.run(denoisingStats, shell=True)
	export(os.path.join(qiimeResults, "denoising_stats.qzv"))

	""" Applying basic filters for how frequent a variant needs to be """
	print("{0} Applying basic filters for how frequent a variant needs to be: in progress ..".format(datetime.now().strftime("%d.%m.%Y %H:%M")))
	filteringVariant = ' '.join([
	"qiime feature-table filter-features",  # Calling qiime2 feature-table filter-features function
	"--quiet",  # Silence output if execution is successful
	"--i-table", os.path.join(qiimeResults, "table.qza"),  # Input feature table
	"--p-min-frequency", freqTheshold(os.path.join(qiimeResults, "feature_table.qzv")),  # Least frequency that a feature must have to be retained
	"--p-min-samples", "1",  # The minimum number of samples that a feature must be observed in to be retained
	"--o-filtered-table", os.path.join(qiimeResults, "feature_table_filtered.qza"),  # Output file
	"2>>", os.path.join(reportsDir, "qiime2_filteringVariant_report.log")])  # Output denoising report
	subprocess.run(filteringVariant, shell=True)


	##################### VISUALISATION AND DATA EXPORT ##########################################


	""" At this stage, we have obtained the artifacts containing the feature table and corresponding feature sequences. 
	Now we will generate summary of the above features and proceed with visualisation of the data. """
	# Generate a heatmap representation of the filtered feature table
	print("{0} Generating a heatmap representing the filtered feature table: in progress ..".format(datetime.now().strftime("%d.%m.%Y %H:%M")))
	filteredFeaturesHeatmap = ' '.join([
	"qiime feature-table heatmap",  # Calling qiime2 feature-table heatmap function
	"--quiet",  # Silence output if execution is successful
	"--p-color-scheme", "Paired",  # The color scheme of the heatmap
	"--p-cluster", "features",  # Perform the clusterring based on the features
	"--o-visualization", os.path.join(qiimeResults, "heatmap_filtered.qzv"),  # Output directory
	"--i-table", os.path.join(qiimeResults, "feature_table_filtered.qza"),  # The input filtered feature table
	"2>>", os.path.join(preprocessingReports, "qiime2_filteredFeaturesHeatmap_report.log")])  # Output filteredFeaturesHeatmap report
	subprocess.run(filteredFeaturesHeatmap, shell=True)
	export(os.path.join(qiimeResults, "heatmap_filtered.qzv"))

	# Obtaing the filtered sequences
	print("{0} Obtaing the filtered sequences: in progress ..".format(datetime.now().strftime("%d.%m.%Y %H:%M")))
	filteredFeaturesSeqs = ' '.join([
	"qiime feature-table filter-seqs",
    # "--m-sample-metadata-file", args.metadata,  # Metadata file
	"--i-data", os.path.join(qiimeResults, "representative_sequences.qza"),  # The sequences from which features should be filtered.
    "--i-table", os.path.join(qiimeResults, "feature_table_filtered.qza"),  # Input table containing feature ids used for id-based filtering
    "--o-filtered-data", os.path.join(qiimeResults, "representative_sequences_filtered.qza"),  # The output filtered sequences
	"2>>", os.path.join(preprocessingReports, "qiime2_filteredFeaturesSeqs_report.log")])  # Output filteredFeaturesSeqs report
	subprocess.run(filteredFeaturesSeqs, shell=True)

	# New summary of the filtered abundance table
	print("{0} Obtaing the filtered abundance table: in progress ..".format(datetime.now().strftime("%d.%m.%Y %H:%M")))
	filteredFeaturesSummary = ' '.join([
	"qiime feature-table summarize",  # Calling qiime2 feature-table summarize function
	"--m-sample-metadata-file", args.metadata,  # Metadata file
	"--i-table", os.path.join(qiimeResults, "feature_table_filtered.qza"),  # Input feature table to be summarized
	"--o-visualization", os.path.join(qiimeResults, "feature_table_filtered.qzv"),  # Output file
	"2>>", os.path.join(preprocessingReports, "qiime2_filteredFeaturesSummary_report.log")])  # Output filteredFeaturesSummary report
	subprocess.run(filteredFeaturesSummary, shell=True)
	export(os.path.join(qiimeResults, "feature_table_filtered.qzv"))
	return 

def phylogenetic_diversity_analysis():
	""" This pipeline will start by creating a sequence alignment using MAFFT,
  	after which any alignment columns that are phylogenetically uninformative
  	or  ambiguously aligned  will be removed  (masked). The resulting masked
  	alignment will be used to infer a phylogenetic tree and then subsequently
  	rooted at its  midpoint. Afterwards, a  collection of diversity metrics 
  	(both phylogenetic and non-phylogenetic) is being applied to the feature 
  	table. """
	print("\n\t{0} PHYLOGENETIC DIVERSITY ANALYSIS".format(datetime.now().strftime("%d.%m.%Y %H:%M")))
	
	if not os.path.exists(diversityAnalysis): os.makedirs(diversityAnalysis)  # Creating the directories which will host the analysis
	if not os.path.exists(statisticalAnalysis): os.makedirs(statisticalAnalysis)

	print("{0} Generating a rooted tree for phylogenetic diversity analysis: in progress ..".format(datetime.now().strftime("%d.%m.%Y %H:%M")))
	phylogeneticDiversityAnalysis =	' '.join([
	"qiime phylogeny align-to-tree-mafft-fasttree", # Calling qiime2 align-to-tree-mafft-fasttree function
	"--quiet",  # Silence output if execution is successful
	"--p-n-threads", str(args.threads),  # Number of threads to use
	"--i-sequences", os.path.join(qiimeResults, "representative_sequences_filtered.qza"),  # he sequences to be used for creating a phylogenetic tree
	"--o-alignment", os.path.join(diversityAnalysis, "aligned_representative_sequences.qza"),  # The aligned sequences
	"--o-masked-alignment", os.path.join(diversityAnalysis, "masked_aligned_representative_sequences.qza"),  # The masked alignment
	"--o-tree", os.path.join(diversityAnalysis, "unrooted_tree.qza"),  # The unrooted phylogenetic tree
	"--o-rooted-tree", os.path.join(diversityAnalysis, "rooted_tree.qza"),  # The rooted phylogenetic tree
	"2>>", os.path.join(reportsDir, "qiime2_phylogeneticDiversityAnalysis_report.log")])  # Output phylogeneticDiversityAnalysis report
	subprocess.run(phylogeneticDiversityAnalysis, shell=True)

	""" A key quality control step is to plot rarefaction curves for all 
	the samples to determine if performed sufficient sequencing """
	print("{0} Performing a rarefraction curves analysis: in progress ..".format(datetime.now().strftime("%d.%m.%Y %H:%M")))
	# max_depth, steps = maxDepthNpStepsThresholds()
	rarefactionCurvesAnalysis =	' '.join([
	"qiime diversity alpha-rarefaction",  # Calling qiime2 diversity alpha-rarefaction function
	"--quiet",  # Silence output if execution is successful
	# "--p-max-depth", "79900",  # The maximum rarefaction depth
	"--p-max-depth", "258000",  # The maximum rarefaction depth
	"--p-steps", "5000",  # The number of rarefaction depths to include between min_depth and max_depth
	"--m-metadata-file", args.metadata,  # Metadata file
	"--i-table", os.path.join(qiimeResults, "feature_table_filtered.qza"),  # Input filtered feature table
	"--i-phylogeny", os.path.join(diversityAnalysis, "rooted_tree.qza"),  #  Input phylogeny for phylogenetic metrics
	"--o-visualization", os.path.join(diversityAnalysis, "rarefaction_curves.qzv"),  # Output visualisation
	"2>>", os.path.join(reportsDir, "qiime2_rarefactionCurvesAnalysis_report.log")])  # Output rarefactionCurvesAnalysis report
	subprocess.run(rarefactionCurvesAnalysis, shell=True)
	export(os.path.join(diversityAnalysis, "rarefaction_curves.qzv"))

	""" Common alpha and beta-diversity metrics and ordination plots (such as PCoA plots for weighted UniFrac distances) 
	This command will also rarefy all samples to the sample sequencing depth before calculating these metrics 
	(X is a placeholder for the lowest reasonable sample depth; samples with depth below this cut-off will be excluded) """
	print("{0} Calculate multiple diversity metrics: in progress ..".format(datetime.now().strftime("%d.%m.%Y %H:%M")))
	diversityMetrics =	' '.join([
	"qiime diversity core-metrics-phylogenetic",  # Calling qiime 2diversity core-metrics-phylogenetic function
	"--m-metadata-file", args.metadata,  # Metadata file
	"--quiet",  # Silence output if execution is successful
	"--p-n-jobs", str(args.threads), # The number of CPUs to be used for the computation
	"--i-phylogeny", os.path.join(diversityAnalysis, "rooted_tree.qza"),  # The rooted phylogenetic tree
	"--i-table", os.path.join(qiimeResults, "feature_table_filtered.qza"),  # The input filtered feature table
	"--p-sampling-depth", "79900",  # The total frequency that each sample should be rarefied to prior to computing diversity metrics
	"--output-dir", os.path.join(diversityAnalysis, "core_metrics_results"),  # Output directory that will host the core metrics
	"2>>", os.path.join(reportsDir, "qiime2_diversityMetrics_report.log")])  # Output diversityMetrics report
	subprocess.run(diversityMetrics, shell=True)
	export(os.path.join(diversityAnalysis, "core_metrics_results"))



	""" Visually and statistically compare groups of alpha diversity values """
	print("{0} Creating Shannon diversity significance boxplot: in progress ..".format(datetime.now().strftime("%d.%m.%Y %H:%M")))
	shannonSignificance =	' '.join([
	"qiime diversity alpha-group-significance",  # Calling qiime diversity alpha-group-significance function
	"--m-metadata-file", args.metadata,  # Metadata file
	"--quiet",  # Silence output if execution is successful
	"--i-alpha-diversity", os.path.join(diversityAnalysis, "core_metrics_results/shannon_vector.qza"),  # Vector of alpha diversity values by sample
	"--o-visualization", os.path.join(statisticalAnalysis, "shannon_diversity-significance.qzv"),  # Output stats
	"2>>", os.path.join(reportsDir, "qiime2_shannonSignificance_report.log")])  # Output stats report
	subprocess.run(shannonSignificance, shell=True)
	export(os.path.join(statisticalAnalysis, "shannon_diversity-significance.qzv")) 

	print("{0} Creating Faith Phylogenetic Diversity (a measure of community richness) significance boxplot: in progress ..".format(datetime.now().strftime("%d.%m.%Y %H:%M")))
	faithSignificance =	' '.join([
	"qiime diversity alpha-group-significance",  # Calling qiime diversity alpha-group-significance function
	"--m-metadata-file", args.metadata,  # Metadata file
	"--quiet",  # Silence output if execution is successful
	"--i-alpha-diversity", os.path.join(diversityAnalysis, "core_metrics_results/faith_pd_vector.qza"),  # Vector of alpha diversity values by sample
	"--o-visualization", os.path.join(statisticalAnalysis, "faith_pd-significance.qzv"),  # Output stats
	"2>>", os.path.join(reportsDir, "qiime2_faithSignificance_report.log")])  # Output report
	subprocess.run(faithSignificance, shell=True)
	export(os.path.join(statisticalAnalysis, "faith_pd-significance.qzv")) 

	print("{0} Creating Evenness significance boxplot: in progress ..".format(datetime.now().strftime("%d.%m.%Y %H:%M")))
	evennessSignificance =	' '.join([
	"qiime diversity alpha-group-significance",  # Calling qiime diversity alpha-group-significance function
	"--m-metadata-file", args.metadata,  # Metadata file
	"--quiet",  # Silence output if execution is successful
	"--i-alpha-diversity", os.path.join(diversityAnalysis, "core_metrics_results/evenness_vector.qza"),  # Vector of alpha diversity values by sample
	"--o-visualization", os.path.join(statisticalAnalysis, "evenness-significance.qzv"),  # Output stats
	"2>>", os.path.join(reportsDir, "qiime2_evennessSignificance_report.log")])  # Output report
	subprocess.run(evennessSignificance, shell=True)
	export(os.path.join(statisticalAnalysis, "evenness-significance.qzv")) 


	""" Visually and statistically compare groups of beta diversity values """
	# Determine whether groups of samples are significantly different from one another using a permutation-based statistical test
	print("{0} Analysing sample composition in the context of Indication using PERMANOVA: in progress ..".format(datetime.now().strftime("%d.%m.%Y %H:%M")))
	permanova_indication =	' '.join([
	"qiime diversity beta-group-significance",  # Calling qiime diversity beta-group-significance function
	"--m-metadata-file", args.metadata,  # Metadata file
	"--p-pairwise",  # Perform pairwise tests between all pairs of groups in addition to the test across all groups
	"--m-metadata-column", "Indication",  # Categorical sample metadata column
	"--quiet",  # Silence output if execution is successful
	"--i-distance-matrix", os.path.join(diversityAnalysis, "core_metrics_results/unweighted_unifrac_distance_matrix.qza"),  # Vector of beta diversity values by sample
	"--o-visualization", os.path.join(statisticalAnalysis, "unweighted_unifrac_indication-significance.qzv"),  # Output stats
	"2>>", os.path.join(reportsDir, "qiime2_permanova_indication_report.log")])  # Output report
	subprocess.run(permanova_indication, shell=True)
	export(os.path.join(statisticalAnalysis, "unweighted_unifrac_indication-significance.qzv")) 

	print("{0} Analysing sample composition in the context of Indication_Code using PERMANOVA: in progress ..".format(datetime.now().strftime("%d.%m.%Y %H:%M")))
	permanova_indication_code =	' '.join([
	"qiime diversity beta-group-significance",  # Calling qiime diversity beta-group-significance function
	"--m-metadata-file", args.metadata,  # Metadata file
	"--p-pairwise",  # Perform pairwise tests between all pairs of groups in addition to the test across all groups
	"--m-metadata-column", "Indication_Code",  # Categorical sample metadata column
	"--quiet",  # Silence output if execution is successful
	"--i-distance-matrix", os.path.join(diversityAnalysis, "core_metrics_results/unweighted_unifrac_distance_matrix.qza"),  # Vector of beta diversity values by sample
	"--o-visualization", os.path.join(statisticalAnalysis, "unweighted_unifrac_indication_code-significance.qzv"),  # Output stats
	"2>>", os.path.join(reportsDir, "qiime2_permanova_indication_code_report.log")])  # Output report
	subprocess.run(permanova_indication_code, shell=True)
	export(os.path.join(statisticalAnalysis, "unweighted_unifrac_indication_code-significance.qzv")) 
	return 

def taxonomic_assignmnet():
	""" We will train the Naive Bayes classifier using SILVA (132) reference sequences 
	and classify the representative sequences from the input dataset """
	print("\n\t{0} TAXONOMIC ASSIGNMENT".format(datetime.now().strftime("%d.%m.%Y %H:%M")))

	# Importing SILVA reference taxonomy sequences
	if not os.path.exists(taxinomicAnalysis): os.makedirs(taxinomicAnalysis)  # Creating the directory which will host the analysis


	# print("Importing the reference sequences and the corresponding taxonomic classifications of SILVA123 database: in progress ..")
	# importSilvaReference = ' '.join([
	# "qiime tools import",  # Import function
	# "--type", "\'FeatureData[Sequence]\'",  # Type of imported data
 #  	"--input-path", silva_reference,  # Input SILVA 132 database
 #  	"--output-path", os.path.join(taxinomicAnalysis, "silva132_99_OTUs.qza"),  # Output file
 #  	"2>>", os.path.join(reportsDir, "qiime2_importSilvaReference_report.log")])  # Output importSilvaReference report
	# subprocess.run(importSilvaReference, shell=True)

	# # Importing SILVA reference taxonomy annotation
	# importSilvaRefTaxonomy = ' '.join([
	# "qiime tools import",  # Import function
	# "--type", "\'FeatureData[Taxonomy]\'",  # Type of imported data
	# "--input-format", "HeaderlessTSVTaxonomyFormat",  # Type of input file
 #  	"--input-path", silva_taxinomy,  # Input annotation file
 #  	"--output-path", os.path.join(taxinomicAnalysis, "silva132_99_OTU_taxonomy.qza"),  # Output artifact 
	# "2>>", os.path.join(reportsDir, "qiime2_importSilvaReference_report.log")])  # Output importSilvaRefTaxonomy report
	# subprocess.run(importSilvaRefTaxonomy, shell=True)

	""" It has been shown that taxonomic classification accuracy of 16S rRNA gene sequences 
	improves when a Naive Bayes classifier is trained on only the region of the target 
	sequences that was sequenced. Here we will extract the reference sequences. """
	# Extract sequencing-like reads from a reference database
	# print("Extract sequencing-like reads from the reference database: in progress ..")
	# extractRefReads = ' '.join([
	# "qiime feature-classifier extract-reads",
	# "--quiet",  # Silence output if execution is successful
	# "--p-trunc-len", "550",  # Minimum amplicon length
	# "--p-f-primer", args.forwardPrimer,  # Forward primer sequence
	# "--p-r-primer", args.reversePrimer,  # Reverse primer sequence
	# "--i-sequences", os.path.join(taxinomicAnalysis, "silva132_99_OTUs.qza"),  # Input reference seq artifact
	# "--o-reads", os.path.join(taxinomicAnalysis, "silva132_reference_sequences.qza"),  # Output ref sequencing-like reads
	# "2>>", os.path.join(reportsDir, "qiime2_importSilvaReference_report.log")])  # Output extractRefReads report
	# subprocess.run(extractRefReads, shell=True)

	""" We can now train a Naive Bayes classifier as follows, using 
	the reference reads and taxonomy that we just created """
	# print("Training the Naive Bayes classifier using the reference reads: in progress ..")
	# trainClassifier = ' '.join([
	# "qiime feature-classifier fit-classifier-naive-bayes",
	# "--quiet",  # Silence output if execution is successful
	# "--i-reference-reads", os.path.join(taxinomicAnalysis, "silva132_reference_sequences.qza"),
	# "--i-reference-taxonomy", os.path.join(taxinomicAnalysis, "silva132_99_OTU_taxonomy.qza"),
	# "--o-classifier", os.path.join(taxinomicAnalysis, "classifier.qza"), 
	# "2>>", os.path.join(reportsDir, "qiime2_trainClassifier_report.log")])  # Output trainClassifier report
	# subprocess.run(trainClassifier, shell=True)
	# export(os.path.join(taxinomicAnalysis, "classifier.qza"))

	""" Assign the taxonomy """
	print("Assign the taxonomy: in progress ..")
	assignTaxonomy = ' '.join([
	"qiime feature-classifier classify-sklearn",
	"--quiet",  # Silence output if execution is successful
	"--p-n-jobs", str(args.threads),  # Number of threads to use
	# "--i-classifier", os.path.join(taxinomicAnalysis, "classifier.qza"), 
	"--i-classifier", silva_99_classifier,  # Using Greengenes classifier
	"--i-reads", os.path.join(qiimeResults, "representative_sequences_filtered.qza"),  # The output filtered sequences
	"--o-classification", os.path.join(taxinomicAnalysis, "taxonomic_classification.qza"),
	"2>>", os.path.join(reportsDir, "qiime2_assignTaxonomy_report.log")])  # Output assignTaxonomy report
	subprocess.run(assignTaxonomy, shell=True)

	outputClassifications = ' '.join([
	"qiime metadata tabulate",
	"--quiet",  # Silence output if execution is successful
	"--m-input-file", os.path.join(taxinomicAnalysis, "taxonomic_classification.qza"),
	"--o-visualization", os.path.join(taxinomicAnalysis, "taxonomic_classification.qzv"),
	"2>>", os.path.join(reportsDir, "qiime2_assignTaxonomy_report.log")])  # Output outputClassifications report
	subprocess.run(outputClassifications, shell=True)
	export(os.path.join(taxinomicAnalysis, "taxonomic_classification.qzv"))

	print("Generating a barplot: in progress ..")
	barplotOfTaxonomy = ' '.join([
	"qiime taxa barplot", 
	"--quiet",  # Silence output if execution is successful
	"--i-table", os.path.join(qiimeResults, "feature_table_filtered.qza"),
	"--i-taxonomy", os.path.join(taxinomicAnalysis, "taxonomic_classification.qza"),
	"--m-metadata-file", args.metadata,  # Metadata file
	"--o-visualization", os.path.join(taxinomicAnalysis, "taxonomy_barplot.qzv"),
	"2>>", os.path.join(reportsDir, "qiime2_barplotOfTaxonomy_report.log")])  # Output barplotOfTaxonomy report
	subprocess.run(barplotOfTaxonomy, shell=True)
	export(os.path.join(taxinomicAnalysis, "taxonomy_barplot.qzv"))
	summarisation()
	return

def differential_abundance():

	# Filter out the samples with less than 2319 reads (in the case of the 18S rRNA data)
	# qiime feature-table filter-samples \
	# --i-table feature_table_samples_filtered.qza \
	# --p-min-frequency 2319 \
	# --o-filtered-table temp.qza

	# # Filter features that only appear in a single sample
	# qiime feature-table filter-features \
	# --i-table temp.qza \
	# --p-min-samples 2 \
	# --o-filtered-table temp2.qza

	# # Filter features by reads
	# qiime feature-table filter-features \
	# --i-table temp2.qza \
	# --p-min-frequency 10 \
	# --o-filtered-table temp3.qza

	# mkdir ANCOM

	# # Collapse table at genus level
	# qiime taxa collapse \
	# --i-table temp3.qza \
	# --i-taxonomy taxonomy/classified_rep_seqs.qza \
	# --p-level 6 \
	# --o-collapsed-table ANCOM/feature_table_for_ANCOM.qza

	# rm temp.qza
	# rm temp2.qza
	# rm temp3.qza

	# cd ANCOM

	# # Add pseudocount
	# qiime composition add-pseudocount \
	# --i-table feature_table_for_ANCOM.qza \
	# --o-composition-table ANCOM_ready_table.qza
	# # Run ANCOM on 18S data by cases/controls
	# qiime composition ancom \
	# --i-table ANCOM_ready_table.qza \
	# --m-metadata-file ../metadata_for_qiime2.txt \
	# --m-metadata-column group \
	# --o-visualization ancom_group_results.qzv \
	# --verbose

	# mkdir biplot

	# # Make the relative frequency table from the rarefied table
	# qiime feature-table relative-frequency \
	# --i-table rarefied_table.qza \
	# --o-relative-frequency-table biplot/rarefied_table_relative.qza

	# # Make the biplot for unweighted UniFrac
	# qiime diversity pcoa-biplot \
	# --i-pcoa unweighted_unifrac_pcoa_results.qza \
	# --i-features biplot/rarefied_table_relative.qza \
	# --o-biplot biplot/biplot_matrix_unweighted_unifrac.qza

	# cd biplot

	# # Turn this matrix into an emperor plot
	# qiime emperor biplot \
	# --i-biplot biplot_matrix_unweighted_unifrac.qza \
	# --m-sample-metadata-file ../../metadata_for_qiime2_with_18S.txt \
	# --m-feature-metadata-file ../../taxonomy/classified_rep_seqs.qza \
	# --o-visualization unweighted_unifrac_emperor_biplot.qzv
	return

def calculateMinMax(forwardRead):
	primers_averageLength = int((len(args.forwardPrimer) + len(args.reversePrimer))/2)

	read_length = 0
	with gzip.open(forwardRead, "rt") as handle:
	    for read in SeqIO.parse(handle, "fastq"):
	        read_length = len(read.seq)
	        break

	min_readLength = read_length - (primers_averageLength * 4)
	max_readLength = (read_length * 2) - (primers_averageLength * 2)
	return (min_readLength, max_readLength)

def freqTheshold(exportFile):
	""" Based on the summary we will calculate a cut-off for how frequent a variant needs to be for it to be retained. 
	We will remove all ASVs that have a frequency of less than 0.1% of the mean sample depth. This cut-off excludes 
	ASVs that are likely due to MiSeq bleed-through between runs (reported by Illumina to be 0.1% of reads). """
	if not os.path.exists(exportFile[:-4]):
		subprocess.run("qiime tools export --input-path {0} --output-path {1}".format(exportFile, exportFile[:-4]), shell=True)

	# Calculating the mean frequency and the min cut-off for excluding ASVs that are likely 
	# due to MiSeq bleed-through between runs (reported by Illumina to be 0.1% of reads)
	mean_freq = 0
	for path, subdir, folder in os.walk(exportFile[:-4]):
		for name in folder:
			file = os.path.join(path, name)
			if file.endswith("sample-frequency-detail.csv"):
				samples = 0
				with open(file) as fin:
					for line in fin:
						samples += 1
						mean_freq += float(line.split(",")[1])
				mean_freq = mean_freq/samples

	# Calculation of the cut-off point
	cut_off = str(int(mean_freq * 0.001))
	return cut_off

def maxDepthNpStepsThresholds():
	""" Calculating the maximum depth and step for the rarefraction experiment"""
	depth_file = os.path.join(qiimeResults, "feature_table_filtered/sample-frequency-detail.csv")
	if not os.path.exists(depth_file):
		sys.exit("The file \"{0}\" does NOT exist...").format(depth_file)

	depths = []
	with open(depth_file) as fin:
		for line in fin:
			depths.append(line.split(",")[1].strip())
	print(depths)

	return

def export(exportFile):
	if os.path.isfile(exportFile):
		subprocess.run("qiime tools export --input-path {0} --output-path {1}".format(exportFile, exportFile[:-4]), shell=True)
	else:
		for path, subdir, folder in os.walk(exportFile):
			for files in folder:
				file = os.path.join(path, files)
				subprocess.run("qiime tools export --input-path {0} --output-path {1}".format(file, file[:-4]), shell=True)
			
	return 

def summarisation():

	for path, subdir, folder in os.walk(analysisDir):
		for name in folder:
			file = os.path.join(path, name)
			if os.stat(file).st_size == 0:  
				print("Removing:\t", file)  
				os.remove(file)
	return 

def main():

	pairedReads = assess_input_data(args.input_dir)
	# Obtaining the number of pair files
	totNum = len(pairedReads)

	### Preprocessing of the input data
	# Performing quality trimming and removal of all primers on both reads
	print("\t{0} PREPROCESSING THE INPUT SAMPLES".format(datetime.now().strftime("%d.%m.%Y %H:%M")))
	for i, read in enumerate(pairedReads, 1):
		primer_removal(read[0], read[1], i, totNum) 

	quality_control()  # Checking the quality of the merged reads

	### Qiime2 analysis 
	qiime2_analysis()
	
	## Downstream analysis
	# phylogenetic_diversity_analysis()
	
	# taxonomic_assignmnet()

	# differential_abundance()

	

if __name__ == "__main__": main()