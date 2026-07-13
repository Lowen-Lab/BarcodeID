'''
barcoded amplicon processing pipeline - v1.2
###################################################################################################
#####################################    About this script    #####################################
###################################################################################################

This script is designed to take in amplicon sequences that have nucleotide barcodes at specific 
sites. There is no minimum or maximum number of barcodes, but this script currently only permits 
two possible nucleotides per site (e.g. a barcode at site 33 can be either A or T, but not C or G).
This script process the illumina sequence files to generate tables that summarize the frequency of
all barcodes detected in each sample, as well as reporting alpha and beta diversity indicies for
those samples.

Barcodes are identified in reads where the alleles you detect are: (1) from the list of 
possible nucleotides at barcode sites, (2) match all the expected nucleotide at all non-barcode site,
and (3) meet the base quality thresholds set by the user.

Features of this pipeline:
 -	It automatically detect and process all fastq files in a specified directory without needing to
 	explicitly list the filenames.
 -	Input sequences can be compressed or uncompressed .fastq files
 -	It has a parallel processing mode that takes advantace of multi-core processors.
 -	It has the option to consider base quality of barcode and non-barcode sites (a.k.a. backbone 
	sites) separately.
 -	A key feature is that you may interrupt the script running at any time, then resume , and will not attempt
	to redo steps that have already been completed, unless requested by the user.
		- Safeguards are in place to prevent incompletely generated output files from being saved.
		- Additional samples can be added to a folder at any time, and re-running BarcodeID will process
		  only the additional samples, before generating combined outputs
 -	All information about reads discarded because they contain unexpected bases are saved in separate
	tables to help identify positively selected mutants.
 -	The script has the ability to complete its analysis where the only required user input is the path to
	the directory where reads are located (default: "raw_data")

########################################    Requirements    #######################################
python3
bbtools (available in path)
Java (for running bbtools)
bash (for calling bbtools) (this is available by default on Mac and Linus operating systems)

For any questions, please contact the author of this script, Dave VanInsberghe, through the Lowen Lab
GitHub, or directly via dvanins@emory.edu.
If you find this script helpful, please cite the original publication to which it was published.

###################################################################################################
###################################   HOW TO RUN THIS SCRIPT    ###################################
###################################################################################################

In terminal, enter:
cd /location/of/script/
python barcodeID.py

The first time you run the script, it will pick a random sample and attempt to predict the barcode 
site locations and possible nucleotides. Answer the text prompts to confirm the prediction is 
correct, or enter the correct values. All remaining samples will then be processed.

###################################################################################################
#######################################    UPDATE NOTES    ########################################
###################################################################################################
03/16/2033 - v0.2 -	Added ability to process reads that have already been cleaned/merged
05/17/2023 - v0.3 -	Added auto-detection of barcode sites and expected amplicon length if not available
05/31/2023 - v0.4 -	Added auto-detection of forward and reverse read file extension and file pairs
					and steps to calculate alpha and beta diversity indicies
07/05/2023 - v0.5 - Added automatic unzip feature to check if new files need unzipping first
10/06/2023 - v0.6 - Added function to trim ends of reads to skip low quality ends
11/17/2023 - v1.0 - First publicly available version
04/25/2024 - v1.1 - Added generation of stacked bar plot figures and improved file suffix prediction
10/07/2024 - v1.2 - Added barcode quality correction and agglomeration (similar to DADA2 for 16S)
'''
###################################### Load required modules ######################################
import os
import sys
import math
import scipy
import numpy as np
import pandas as pd
import time
import random
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('Agg')

###################################### Establish variables #######################################
project_dir = "./"
output_dir = project_dir+"barcode_info/"
trim_dir = project_dir + 'trim/'
temp_dir = trim_dir+"temp/"
command_prefix = "bash "
###################################################################################################
##################################  USER DEFINED VARIABLES  #######################################
###################################################################################################
input_read_dir = project_dir+"raw_data/"

amplicon_info_infile_name = 'amplicon_info.txt'
titer_filename = 'sample_titers.txt' #this file is not necessary, but can be added to scale stacked box plots to log10 transformed viral titers. Tab separated file, one sample per line: sampleID \t titer \n

min_Q_score = 20
F_edge_ignore = 0 #these values are used to exclude the edge of reads if their quality is low and needlessly discarding excessive data due to sequencing quality issues
R_edge_ignore = 0 #these values are used to exclude the edge of reads if their quality is low and needlessly discarding excessive data due to sequencing quality issues
force_trim_pre_merge = 0 #use this value to force trimming the first x nucleotides from input reads, which is helpful when the beginning of reads are low quality
truncate_length_post_merge = 0 #use this value to truncate reads post-merging if amplicons are not always equal lengths - use only as a last resort to make processing complete, or for troubleshooting when no barcodes are detected in some samples

overwrite_existing_files = False #when False, pipeline will avoid re-running samples that have already been processed 
generate_stacked_bar_plots = True #when True, stacked bar plots will be generated from observed barcodes frequencies

reads_already_screened_for_quality = False # 'True' or 'False' --- If reads were already merged and screened for average quality and length (such as during a demultiplexing step), the script will not attempt to filter or merge the input reads
raw_read_input_type = "paired" # 'paired' or 'single' --- When 'single' the pipeline will look for barcodes in one file, when 'paired' the pipeline will look for forward and reverse reads. Single core use has not been extensively tested - recommend using "paired" and setting "parallel_max_cpu" to 1 if you only wish to use one core at a time
forward_read_suffix = "_R1.fastq" #if empty, the script will attempt to predict these values, but it only works under specific assumptions
reverse_read_suffix = "_R2.fastq" #this value is ignored if 'reads_already_screened_for_quality' is 'True' or 'raw_read_input_type' is 'single' - NOTE: this variable cannot be removed without causing missing value errors

parallel_process = True #'True' or 'False' --- when True, the script uses the joblib package to run multiple samples simultaneously. Currently not stable on windows
parallel_max_cpu = 150 # when parallel process true, this sets the max CPUs the script is allowed to use. If you want to use all available, set to 0. If you want to use all but one, set to -1 (or -2 to use all but 2 CPUs)

global verbose
verbose = True #when this value is 'True', the script will print status updates

num_subsamples = 5 #number of times barcode counts should be subsampled to calculate diversity indicies
min_barcodes_per_sample = 4000 #non-unique barcode counts that pass processing to be included in final output tables and statistics
############################################ FUNCTIONS ############################################
def run_command(command,mode="quiet"):
# def run_command(command,mode="verbose"):
	run_status = False
	if mode == "verbose":
		return_code = os.system(command)
	elif mode == "quiet":
		return_code = os.system(command+" >/dev/null 2>&1")
	# subprocess.check_output(command,stderr=subprocess.STDOUT)
	if return_code == 0:
		run_status = True
	return run_status

def unzip_gz(filename,path_to_file):
	command = "gzip -d "+path_to_file+filename
	out = run_command(command)
	if out == False:
		sys.exit()
	return filename

def unzip_tar_gz(filename,path_to_file):
	command = "tar -xvzf "+path_to_file+filename
	out = run_command(command)
	if out == False:
		sys.exit()
	return filename

def zip_gz(filename,path_to_file):
	command = "gzip -9 "+path_to_file+filename
	out = run_command(command)
	if out == False:
		sys.exit()
	return filename

def zip_tar_gz(filename,path_to_file):
	command = "tar -cvzf "+path_to_file+filename+".tar.gz "+path_to_file+filename
	out = run_command(command)
	if out == False:
		sys.exit()
	return filename

def load_barcode_info(amplicon_info_infile_name):
	expected_amplicon_seq = ''
	barcode_sites = {}
	amplicon_info_infile = open(amplicon_info_infile_name,"r")
	for line in amplicon_info_infile:
		line = line.strip()
		if len(line) >0:
			if line[0] != "#":
				line = line.split("\t")
				if line[0] == "amplicon_seq":
					expected_amplicon_seq = line[1]
				else:
					site = int(line[0])
					barcode_sites[site] = line[1].split(",")
	amplicon_info_infile.close()
	return expected_amplicon_seq,barcode_sites


def check_if_input_zipped(input_dir,parallel_process,num_cores):
	files_found = False
	#Unzip any files that need unzipping
	tar_gz_inputs = [f for f in os.listdir(input_dir) if f.endswith(".tar.gz")]
	if len(tar_gz_inputs) > 0:  #unzip .tar.gz files
		if parallel_process == True:
			processed_list = Parallel(n_jobs=num_cores)(delayed(unzip_tar_gz)(i,input_dir) for i in tar_gz_inputs)
		else:
			for targzfile in tar_gz_inputs:
				unzip_tar_gz(gzfile,input_dir)
	gz_inputs = [f for f in os.listdir(input_dir) if f.endswith(".gz")]
	if len(gz_inputs) >0:  #unzip .gz files
		if parallel_process == True:
			processed_list = Parallel(n_jobs=num_cores)(delayed(unzip_gz)(i,input_dir) for i in gz_inputs)
		else:
			for gzfile in gz_inputs:
				unzip_gz(gzfile,input_dir)
	fq_inputs = [f for f in os.listdir(input_dir) if f.endswith(".fastq") or f.endswith(".fq")]
	if len(fq_inputs) >0:
		files_found = True
	return files_found


def predict_paired_file_extensions(input_dir):
	exten_filelist = [f for f in os.listdir(input_dir) if f.endswith(".fastq") or f.endswith(".fq")]
	shortest_filename = exten_filelist[0]
	for f in range(0,len(exten_filelist)):
		filename = exten_filelist[f]
		if len(filename) < len(shortest_filename):
			shortest_filename = filename
	string_dict = {}
	for i in range(-1,(-1*len(shortest_filename)),-1):
		string_dict[i] = {}
		for f in range(0,len(exten_filelist)):
			filename = exten_filelist[f]
			string = filename[i:]
			try:
				string_dict[i][string] += 1
			except:
				string_dict[i][string] = 1
	first_dichotomous = 0
	last_dichotomous = 0
	for i in range(-1,(-1*len(shortest_filename)),-1):
		if len(string_dict[i]) == 2:
			if first_dichotomous == 0:
				first_dichotomous = i
			else:
				last_dichotomous = i
	exten_count = {}
	for f in range(0,len(exten_filelist)):
		filename = exten_filelist[f]
		try:
			exten_count[filename[last_dichotomous:]] += 1
		except:
			exten_count[filename[last_dichotomous:]] = 1
	exten_pair = {1:'',2:''}
	for exten in exten_count:
		try:
			exten_pair[int(exten[first_dichotomous])] = exten
		except:
			pass
	successful_prediction = False
	if len([f for f in os.listdir(input_dir) if f.endswith(exten_pair[1])]) == len([f for f in os.listdir(input_dir) if f.endswith(exten_pair[2])]):
		if len(exten_filelist)/2 == len([f for f in os.listdir(input_dir) if f.endswith(exten_pair[1])]):
			successful_prediction = True
			return exten_pair[1],exten_pair[2]
	if successful_prediction == False:
		sys.exit('Unable to predict forward and reverse read suffix values. Enter manually in "user defined variables" section in BarcodeID script to proceed.\nExiting.')

def predict_single_file_extension(input_dir):
	exten_filelist = [f for f in os.listdir(input_dir) if f.endswith(".fastq") or f.endswith(".fq")]
	shortest_filename = exten_filelist[0]
	for f in range(0,len(exten_filelist)):
		filename = exten_filelist[f]
		if len(filename) < len(shortest_filename):
			shortest_filename = filename
	string_dict = {}
	for i in range(-1,(-1*len(shortest_filename)),-1):
		string_dict[i] = {}
		for f in range(0,len(exten_filelist)):
			filename = exten_filelist[f]
			string = filename[i:]
			try:
				string_dict[i][string] += 1
			except:
				string_dict[i][string] = 1
	last_monomorphic = 0
	searching_for_last_monomorphic = True
	for i in range(-1,(-1*len(exten_filelist[0])),-1):
		if len(string_dict[i]) == 1 and searching_for_last_monomorphic == True:
			last_monomorphic = i
		elif len(string_dict[i]) > 1 and searching_for_last_monomorphic == True:
			searching_for_last_monomorphic = False
	exten_count = {}
	for f in range(0,len(exten_filelist)):
		filename = exten_filelist[f]
		try:
			exten_count[filename[last_monomorphic:]] += 1
		except:
			exten_count[filename[last_monomorphic:]] = 1
	successful_prediction = False
	if exten_count[filename[last_monomorphic:]] == len(exten_filelist):
		successful_prediction = True
		return filename[last_monomorphic:]
	if successful_prediction == False:
		sys.exit("Unable to predict read suffix value. Enter manually to proceed.\nExiting.")

def detect_barcode_content(seq_dict,amplicon_length):
	amplicon_expect = ''
	barcode_sites = {}
	nts = ['A','T','C','G']
	prop_lineout = 'loc\tA\tT\tC\tG\n'
	minor_freq_list = []
	minor_freq_loc_list = []
	for col in range(0,amplicon_length):
		prop_lineout += str(col)
		nt_count_dict = {'A':0,'T':0,'C':0,'G':0}
		for seq in seq_dict:
			if len(seq) == amplicon_length:
				try:
					cur_nt = seq[col]
					nt_count_dict[cur_nt] += 1
				except:
					if seq[col] != "N":
						sys.exit('Unexpected character found: "'+cur_nt+'" at site '+str(col))
		nt_list = []
		for n in range(0,len(nts)):
			nt = nts[n]
			nt_counted = nt_count_dict[nt]
			nt_prop = float(nt_counted)/float(len(seq_dict))
			prop_lineout += '\t'+str(nt_prop)
			temp_tup = (nt_count_dict[nt],nt)
			nt_list.append(temp_tup)
		nt_list = sorted(nt_list, reverse=True)
		major_allele = nt_list[0][1]
		minor_allele = nt_list[1][1]
		major_freq = float(nt_list[0][0])/float(len(seq_dict))
		minor_freq = float(nt_list[1][0])/float(len(seq_dict))
		minor_freq_list.append(minor_freq)
		loc_tup = (minor_freq,col,(major_allele,minor_allele))
		minor_freq_loc_list.append(loc_tup)
		prop_lineout += '\n'
	outfile = open(project_dir+"barcode_predict_info.txt","w")
	outfile.write(prop_lineout)
	outfile.close()
	
	amplicon_expect = "N"*amplicon_length
	minor_freq_loc_list = sorted(minor_freq_loc_list,reverse=True)
	med_freq = np.median(minor_freq_list)
	for i in range(0,len(minor_freq_loc_list)):
		local_freq = minor_freq_loc_list[i][0]
		loc = minor_freq_loc_list[i][1]
		allele_tup = minor_freq_loc_list[i][2]
		if (local_freq/med_freq)>=10 and local_freq >= 0.03:
			barcode_sites[loc] = [allele_tup[0],allele_tup[1]]
			# amplicon_expect += "N"
		else:
			amplicon_expect = amplicon_expect[0:loc]+allele_tup[0]+amplicon_expect[loc+1:amplicon_length]
			# amplicon_expect[loc] = allele_tup[0]
	return amplicon_expect, barcode_sites


def subset_fastq(filename,reads_to_calc_median_len=1000,unique_seqs_to_pull=100):	
	read_len_list = []
	trim_seq_subset_dict = {}
	num_added = 0
	num_seq_passed = 0
	line_counter = -1
	avg_read_len = 0
	infile = open(filename,"r")
	for line in infile:
		line = line.strip()
		line_counter += 1
		if line_counter == 0:
			header = line
			writing = False
		elif line_counter == 1:
			num_seq_passed += 1
			if num_seq_passed <= 5000:
				read_len_list.append(len(line))
				if num_seq_passed == 5000:
					avg_read_len = int(np.median(read_len_list))
					# print(avg_read_len)
			elif len(line) == avg_read_len:
				try:
					trim_seq_subset_dict[line]
				except:
					trim_seq_subset_dict[line] = ''
					num_added += 1
					if num_added >= unique_seqs_to_pull:
						break
		elif line_counter == 3:
			line_counter = -1
	infile.close()
	return trim_seq_subset_dict,avg_read_len


def filter_fastq(filename_in,filename_out,length_min,length_max):
	line_counter = -1
	infile = open(filename_in,"r")
	outfile = open(filename_out,"w")
	for line in infile:
		line = line.strip()
		line_counter += 1
		if line_counter == 0:
			header = line
			writing = False
		elif line_counter == 1:
			if len(line) >= length_min and len(line) <= length_max:
				outfile.write(header+"\n"+line+"\n")
				writing = True
			else:
				writing = False
		else:
			if line_counter == 3:
				line_counter = -1
			if writing == True:
				outfile.write(line+"\n")
	infile.close()
	outfile.close()


def calc_expect_read_len(filename_in):
	len_list = []
	line_counter = -1
	infile = open(filename_in,"r")
	for line in infile:
		line = line.strip()
		line_counter += 1
		if line_counter == 0:
			header = line
			writing = False
		elif line_counter == 1:
			len_list.append(len(line))
		elif line_counter == 3:
			line_counter = -1
			if len(len_list) >= 1000:
				break
	infile.close()
	try:
		return int(np.nanmedian(len_list))
	except:
		return 0


def check_adapter(filename_in):
	infile = open(filename_in,'r')
	adapter_pass = False
	for line in infile:
		line = line.strip()
		if len(line) >0:
			if line[0] == ">":
				header = line
			else:
				seq = line
				if seq != "N":
					adapter_pass = True
	return adapter_pass


def fastq_to_fasta(filename_in,filename_out):
	line_counter = -1
	infile = open(filename_in,"r")
	outfile = open(filename_out,"w")
	for line in infile:
		line = line.strip()
		line_counter += 1
		if line_counter == 0:
			outfile.write(line.split(" ")[0]+"\n")
		elif line_counter == 1:
			outfile.write(line+"\n")
		elif line_counter == 3:
			line_counter = -1
	infile.close()
	outfile.close()
	return filename_out

def pull_qual_info(input_seq_filename,max_read_count=1e6,min_qual_count=20):
	flank_ignore = 0
	infile = open(input_seq_filename,"r")
	line_counter = -1
	read_count = 0
	qual_sub_dict = {}
	loc_sub_dict = {}
	sub_per_read_dict = {}
	keep_reading = True
	for line in infile:
		line = line.strip()
		line_counter += 1
		if line_counter == 0:
			head = line
			read_count += 1
			if read_count == max_read_count:
				keep_reading = False
		elif line_counter == 1:
			nt_string = line
			nt_string = nt_string[0:len(nt_string)-flank_ignore]
		elif line_counter == 3:
			if keep_reading == True:
				line_counter = -1
			phreds = line
			phreds = phreds[0:len(phreds)-flank_ignore]
			read_sub_info = []
			mismatch_count = 0
			lowQ_count = 0
			if len(nt_string)==len(expected_amplicon_seq):
				for loc in range(0,len(nt_string)):
					backbone_nt = ''
					expected_nts = {}
					site_type = ''
					try:
						expected_nts = barcode_sites[loc]
						site_type = "barcode"
					except:
						expected_nts[expected_amplicon_seq[loc]] = ''
						backbone_nt = expected_amplicon_seq[loc]
						site_type = "backbone"
					if site_type == "backbone":
						nt = nt_string[loc]
						qu = ord(phreds[loc])-33
						sub_tup = (backbone_nt,nt,qu,loc)
						if qu >= min_qual_count:
							read_sub_info.append(sub_tup)
							if nt != backbone_nt:
								mismatch_count += 1
						else:
							lowQ_count += 1
				try:
					sub_per_read_dict[mismatch_count] += 1
				except:
					sub_per_read_dict[mismatch_count] = 1
				if mismatch_count <= 2 and lowQ_count == 0:
					for sub_tup in read_sub_info:
						backbone_nt = sub_tup[0]
						sub_nt = sub_tup[1]
						qu = sub_tup[2]
						loc = sub_tup[3]
						try:
							qual_sub_dict[backbone_nt][sub_nt][qu] += 1
						except:
							try:
								qual_sub_dict[backbone_nt][sub_nt][qu] = 1
							except:
								try:
									qual_sub_dict[backbone_nt][sub_nt] = {}
									qual_sub_dict[backbone_nt][sub_nt][qu] = 1
								except:
									qual_sub_dict[backbone_nt] = {}
									qual_sub_dict[backbone_nt][sub_nt] = {}
									qual_sub_dict[backbone_nt][sub_nt][qu] = 1
	output_tup = (input_seq_filename.split(".f")[0],qual_sub_dict,sub_per_read_dict)
	return output_tup

def read_in_barcodes(barcode_filename_in,mismatch_filename_in,nonbarcode_filename_in):
	barcode_count_dict = {}
	barcode_list = []
	infile = open(barcode_filename_in)
	for line in infile:
		line = line.strip().split("\t")
		if line[0] != "barcode":
			barcode_string = line[0]
			count = int(line[1])
			barcode_count_dict[barcode_string] = count
			barcode_list.append(barcode_string)
	infile.close()

	barcode_list = list(set(barcode_list))
	mismatch_by_site = {}
	infile = open(mismatch_filename_in)
	for line in infile:
		line = line.strip().split("\t")
		if line[0] != "Site_number":
			site = int(line[0])
			prop = float(line[1])
			mismatch_by_site[site] = prop
	infile.close()
	return barcode_count_dict, barcode_list,mismatch_by_site


def shannon_alpha(freq_array):
	H = 0
	S_obs = 0
	for freq in freq_array:
		if freq >0.0:
			H_i = -1*freq*np.log(freq)
			H += H_i
			S_obs += 1
	if S_obs >1:
		H_max = np.log(S_obs)
		E = H/H_max
		H = round(H,3)
		E = float(round_to_n_sig_figs(E,3))
	else:
		H,E = 0,0
	return H,E

def simpson_alpha(freq_array):
	p_sum = 0
	S = 0
	for freq in freq_array:
		if freq >0.0:
			S += 1
			p_sum += freq**2
	H = p_sum
	return float(round_to_n_sig_figs(H,3))

def bray_curtis_dissimilarity_log(array1,array2):
	C,S1,S2 = 0,0,0
	for i in range(0,len(array1)):
		val1 = array1[i]
		val2 = array2[i]
		if val1 > 0.0 and val2 > 0.0:
			minval = min(val1,val2)

			C += 1/(-1*np.log(min(val1,val2)))
		if val1 > 0.0:
			S1 += 1/(-1*np.log(val1))
		if val2 > 0.0:
			S2 += 1/(-1*np.log(val2))
	B = 1-((2*C)/(S1+S2))
	return float(round_to_n_sig_figs(B,3))

def bray_curtis_dissimilarity(array1,array2):
	C,S1,S2 = 0,0,0
	for i in range(0,len(array1)): #assumes both arrays are same size
		val1 = array1[i]
		val2 = array2[i]
		if val1 > 0.0 and val2 > 0.0:
			C += min(val1,val2)
		if val1 > 0.0:
			S1 += val1
		if val2 > 0.0:
			S2 += val2
	B = 1-((2*C)/(S1+S2))
	return float(round_to_n_sig_figs(B,3))

def jaccard_dissimilarity(array1,array2,min_count=1.0):
	AB,A,B,C = 0,0,0,0
	for i in range(0,len(array1)):
		val1 = array1[i]
		val2 = array2[i]
		if val1 >= min_count or val2 >= min_count:
			C += 1
		if val1 >= min_count and val2 >= min_count:
			AB += 1
	if C == 0:
		j = 1.
	else:
		J = 1. - (AB/C)
	return float(round_to_n_sig_figs(J,3))

def chao_1_richness(count_array,min_count=1):
	S_obs = 0
	S_single = 0
	S_double = 0
	for count in count_array:
		if count >= min_count:
			S_obs += 1
			if count == 1:
				S_single += 1
			elif count == 2:
				S_double += 1
	S_unobs = (S_single*(S_single-1))/(2*S_double+1)
	R = S_obs + S_unobs
	return round(R,1)

def true_richness(count_array,min_count=1):
	t_rich = 0
	for i in range(0,len(count_array)):
		count = count_array[i]
		if count >= min_count:
			t_rich += 1
	return t_rich

def is_number(n):
	try:
		float(n)
	except ValueError:
		return False
	return True

def round_to_n_sig_figs(val,num_sig_figs):
	if is_number(val):
		val = float(val)
		if val == 0.0:
			return '0.0'
		elif is_number(num_sig_figs):
			if val<0:
				multiplier = -1
			else:
				multiplier = 1
			val = val*multiplier
			num_sig_figs = float(num_sig_figs)
			if num_sig_figs.is_integer:
				num_sig_figs = int(num_sig_figs)
				if num_sig_figs == 0:
					num_sig_figs = 1
				sci_val = "{:.10e}".format(val)
				split_sci_val = sci_val.split("e")
				if len(split_sci_val) == 2:
					rounded_base_number = round(float(split_sci_val[0]),num_sig_figs-1)
					exponent = int(split_sci_val[1])
					if str(rounded_base_number) == '10.0':
						val_out = str(round(float(val),num_sig_figs))
					elif exponent == 0:
						val_out = str(rounded_base_number) + ((num_sig_figs)-1)*'0'
					elif exponent < 0:
						exponent*=-1
						val_out = '0.' + (exponent-1)*'0' + str(rounded_base_number).replace(".","")
						if exponent >3:
							val_out = str(float(val_out))
						elif len(val_out) >7:
							val_out = str(float(val_out))
					elif exponent > 0:
						val_out = str(rounded_base_number) +'e'+ (str(exponent))
					else:
						sys.exit("Unexpected error while rounding: "+str(val))
					if multiplier == -1:
						val_out = '-'+val_out
					return val_out
			else:
				sys.exit("Non-integer value for 'num_sig_figs' provided: "+str(num_sig_figs))
		else:
			sys.exit("Unable to round: '"+str(val) + "' to: '"+str(num_sig_figs)+"' decimals")
	else:
		sys.exit("Unable to round: '"+str(val) + "' to: '"+str(num_sig_figs)+"' decimals")

def rsquared(x, y):
	slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(x, y)
	return r_value**2

def subsample_counts(count_array,barcodeID_array,num_to_pick):
	num_to_pick = int(num_to_pick)
	rand_list = []
	for i in range(0,len(count_array)):
		count = count_array[i]
		barcodeID = barcodeID_array[i]
		if count >1:
			for j in range(0,count):
				rand_list.append(barcodeID)
	random.shuffle(rand_list)
	subset_list = rand_list[0:num_to_pick]
	count_dict_out = {}
	freq_dict_out = {}
	for i in range(0,len(subset_list)):
		barcodeID = subset_list[i]
		try:
			count_dict_out[barcodeID] += 1
		except:
			count_dict_out[barcodeID] = 1
	local_barcodeID_list = list(count_dict_out.keys())
	float_num_to_pick = float(num_to_pick)
	for i in range(0,len(local_barcodeID_list)):
		barcodeID = local_barcodeID_list[i]
		count = count_dict_out[barcodeID]
		barocde_freq = float(count)/float_num_to_pick
		freq_dict_out[barcodeID] = barocde_freq
	return count_dict_out,freq_dict_out

def subsample_alpha_func(input_tup,full_barcode_list,subsample_iterations=3):
	sampleID = input_tup[0]
	counts = input_tup[1]
	sample_total = np.sum(counts)
	pick_set_list = [1.e3,2.e3,4.e3,8.e3,1e4,2.e4,3e4,4.e4,5e4,6e4,7e4,8e4,9e4,1e5,2e5,4e5,8e5,1e6,2e6]
	diversity_dict = {}
	for s in range(0,len(pick_set_list)):
		subsample_size = pick_set_list[s]
		if sample_total >= subsample_size:
			diversity_dict[subsample_size] = {}
			for iter_num in range(0,subsample_iterations):
				rand_count_dict,rand_freq_dict = subsample_counts(counts,full_barcode_list,subsample_size)
				rand_count = []
				rand_freq = []
				local_barcode_list = []
				for b in range(0,len(full_barcode_list)):
					barcodeID = full_barcode_list[b]
					try:
						count = rand_count_dict[barcodeID]
					except:
						count = 0
					try:
						freq = rand_freq_dict[barcodeID]
					except:
						freq = 0
					if count > 0:
						local_barcode_list.append(barcodeID)
						rand_count.append(count)
						rand_freq.append(freq)
				local_shannon,local_even = shannon_alpha(rand_freq)
				local_simpson = simpson_alpha(rand_freq)
				local_chao = chao_1_richness(rand_count)
				try:
					diversity_dict[subsample_size]['shannon'].append(local_shannon)
					diversity_dict[subsample_size]['simpson'].append(local_simpson)
					diversity_dict[subsample_size]['even'].append(local_even)
					diversity_dict[subsample_size]['chao'].append(local_chao)
					diversity_dict[subsample_size]['rich'].append(int(len(local_barcode_list)))
				except:
					diversity_dict[subsample_size]['shannon'] = [local_shannon]
					diversity_dict[subsample_size]['simpson'] = [local_simpson]
					diversity_dict[subsample_size]['even'] = [local_even]
					diversity_dict[subsample_size]['chao'] = [local_chao]
					diversity_dict[subsample_size]['rich'] = [int(len(local_barcode_list))]
	subsample_count_list = list(diversity_dict.keys())
	subsample_count_list = sorted(subsample_count_list)
	median_diversity_dict = {}
	for s in range(0,len(subsample_count_list)):
		subsample_size = subsample_count_list[s]
		shannon_med = np.median(diversity_dict[subsample_size]['shannon'])
		simpson_med = np.median(diversity_dict[subsample_size]['simpson'])
		even_med = np.median(diversity_dict[subsample_size]['even'])
		chao_med = np.median(diversity_dict[subsample_size]['chao'])
		rich_med = np.median(diversity_dict[subsample_size]['rich'])
		try:
			median_diversity_dict['shannon'].append(shannon_med)
			median_diversity_dict['simpson'].append(simpson_med)
			median_diversity_dict['even'].append(even_med)
			median_diversity_dict['chao'].append(chao_med)
			median_diversity_dict['rich'].append(rich_med)
		except:
			median_diversity_dict['shannon'] = [shannon_med]
			median_diversity_dict['simpson'] = [simpson_med]
			median_diversity_dict['even'] = [even_med]
			median_diversity_dict['chao'] = [chao_med]
			median_diversity_dict['rich'] = [rich_med]
	fit_div_dict = {}
	if median_diversity_dict != {}:
		if len(median_diversity_dict['shannon']) >= 3:
			fit_div_dict['shannon'] = float(round_to_n_sig_figs(fit_alpha(median_diversity_dict['shannon'],subsample_count_list),4))
			fit_div_dict['simpson'] = float(round_to_n_sig_figs(fit_alpha(median_diversity_dict['simpson'],subsample_count_list),4))
			fit_div_dict['even'] = float(round_to_n_sig_figs(fit_alpha(median_diversity_dict['even'],subsample_count_list),4))
			fit_div_dict['chao'] = float(round_to_n_sig_figs(np.exp(fit_alpha(np.log(median_diversity_dict['chao']),subsample_count_list)),4))
			fit_div_dict['rich'] = float(round_to_n_sig_figs(np.exp(fit_alpha(np.log(median_diversity_dict['rich']),subsample_count_list)),4))
	out_tup = (sampleID,diversity_dict,fit_div_dict)
	return out_tup

def fit_alpha(alpha_vals,sub_num_list,poly_num=2,subsample_fit_level=1e4):
	sub_num_list = np.log(sub_num_list)
	poly_coeff = np.polyfit(sub_num_list,alpha_vals,2)
	poly_func = np.poly1d(poly_coeff)
	subsample_fit = poly_func(np.log(subsample_fit_level))
	return subsample_fit
	# min_rsquare = 0.95
	# poly_rsq = rsquared(alpha_vals,poly_func(sub_num_list))
	# if poly_rsq < min_rsquare:
	# 	return 0
	# subsample_fit = float(round_to_n_sig_figs(subsample_fit,4))
	# poly_rsq = float(round_to_n_sig_figs(poly_rsq,4))

def subsample_beta_func(input_tup,full_barcode_list,subsample_iterations=3,subsample_fit_level=1e4):
	sample_tup_1 = input_tup[0]
	sample_tup_2 = input_tup[1]
	sample_1,counts_1 = sample_tup_1[0],sample_tup_1[1]
	sample_2,counts_2 = sample_tup_2[0],sample_tup_2[1]
	sample_total_1 = np.sum(counts_1)
	sample_total_2 = np.sum(counts_2)
	
	# pick_set_list = [5e2,1e3,2e3,4e3,8e3,1.6e4,3.2e4,6.4e4,1.28e5]
	# pick_set_list = [4e3,6e3,8e3,1e4,2e4,3e4,4e4,5e4,6e4,8e4]
	pick_set_list = [8e3,1e4,2e4,4e4]

	dissimilarty_dict = {}
	for s in range(0,len(pick_set_list)):
		subsample_size = pick_set_list[s]
		if sample_total_1 >= subsample_size and sample_total_2 >= subsample_size:
			dissimilarty_dict[subsample_size] = {}
			for iter_num in range(0,subsample_iterations):
				rand_count_dict_1,rand_freq_dict_1 = subsample_counts(counts_1,full_barcode_list,subsample_size)
				rand_count_dict_2,rand_freq_dict_2 = subsample_counts(counts_2,full_barcode_list,subsample_size)
				rand_count_1 = []
				rand_count_2 = []
				rand_freq_1 = []
				rand_freq_2 = []
				local_barcode_list = []
				for b in range(0,len(full_barcode_list)):
					barcodeID = full_barcode_list[b]
					try:
						count1 = rand_count_dict_1[barcodeID]
					except:
						count1 = 0
					try:
						count2 = rand_count_dict_2[barcodeID]
					except:
						count2 = 0
					try:
						freq1 = rand_freq_dict_1[barcodeID]
					except:
						freq1 = 0
					try:
						freq2 = rand_freq_dict_2[barcodeID]
					except:
						freq2 = 0
					if count1 >0 or count2 >0:
						local_barcode_list.append(barcodeID)
						rand_count_1.append(count1)
						rand_count_2.append(count2)
						rand_freq_1.append(freq1)
						rand_freq_2.append(freq2)
				local_bray_dist = bray_curtis_dissimilarity(rand_freq_1,rand_freq_2)
				local_jaccard_dist = jaccard_dissimilarity(rand_count_1,rand_count_2)
				try:
					dissimilarty_dict[subsample_size]['bray'].append(local_bray_dist)
					dissimilarty_dict[subsample_size]['jaccard'].append(local_jaccard_dist)
				except:
					dissimilarty_dict[subsample_size]['bray'] = [local_bray_dist]
					dissimilarty_dict[subsample_size]['jaccard'] = [local_jaccard_dist]
	subsample_count_list = list(dissimilarty_dict.keys())
	subsample_count_list = sorted(subsample_count_list)
	median_dissimilarty_dict = {}
	for s in range(0,len(subsample_count_list)):
		subsample_size = subsample_count_list[s]
		bray_med = np.median(dissimilarty_dict[subsample_size]['bray'])
		jaccard_med = np.median(dissimilarty_dict[subsample_size]['jaccard'])
		try:
			median_dissimilarty_dict['bray'].append(bray_med)
			median_dissimilarty_dict['jaccard'].append(jaccard_med)
		except:
			median_dissimilarty_dict['bray'] = [bray_med]
			median_dissimilarty_dict['jaccard'] = [jaccard_med]
	fit_diss_dict = {}
	try:
		num_med = len(median_dissimilarty_dict['bray'])
	except:
		num_med = 0
	if num_med >=3:
		bray_list = median_dissimilarty_dict['bray']#[max(0,len(median_dissimilarty_dict['bray'])-8):len(median_dissimilarty_dict['bray'])]
		bray_count_list = subsample_count_list#[max(0,len(median_dissimilarty_dict['bray'])-8):len(median_dissimilarty_dict['bray'])]
		jaccard_list = median_dissimilarty_dict['jaccard']#[max(0,len(median_dissimilarty_dict['jaccard'])-8):len(median_dissimilarty_dict['jaccard'])]
		jaccard_count_list = subsample_count_list
		fit_diss_dict['bray'] = float(round_to_n_sig_figs(fit_alpha(bray_list,bray_count_list,1,subsample_fit_level),4))
		fit_diss_dict['jaccard'] = float(round_to_n_sig_figs(fit_alpha(jaccard_list,jaccard_count_list,1,subsample_fit_level),4))
	out_tup = (sample_1,sample_2,dissimilarty_dict,fit_diss_dict)
	return out_tup


def fit_qual_model(prop_sub_infile_name,count_sub_infile_name,fit_prop_sub_outfile_name,min_qual_val_to_fit=20,poly_num=1):
	fit_sub_dict = {}
	min_qual_obs = 999
	max_qual_obs = -1
	num_col = 2
	num_row = 2
	col = -1 #column
	row = -1 #row
	# fig_num = 0
	focal_nt_obs_dict = {}
	nt_count_dict = {}
	count_sub_infile = open(count_sub_infile_name,"r")
	first_line = True
	for line in count_sub_infile:
		line = line.rstrip('\n').split("\t")
		if first_line == True:
			qual_list = line
			first_line = False
		else:
			base_nt = line[0]
			sub_nt = line[1]
			for a in range(2,len(line)):
				qual_val = int(qual_list[a])
				count_val = int(line[a])
				try:
					focal_nt_obs_dict[base_nt+str(qual_val)] += count_val
				except:
					focal_nt_obs_dict[base_nt+str(qual_val)] = count_val

				try:
					nt_count_dict[base_nt+sub_nt+str(qual_val)] += count_val
				except:
					nt_count_dict[base_nt+sub_nt+str(qual_val)] = count_val
	count_sub_infile.close()

	bases = ['A','T','C','G']
	qual_val_list = {}

	for p in range(0,len(bases)):
		focal_nt = bases[p]
		col+=1
		if col%num_col==0:
			row+=1
			col = 0
		for h in range(0,len(bases)):
			focal_sub_nt = bases[h]
			prop_sub_infile = open(prop_sub_infile_name,"r")
			first_line = True
			for line in prop_sub_infile:
				line = line.rstrip('\n').split("\t")
				if first_line == True:
					qual_list = line
					first_line = False
				else:
					base_nt = line[0]
					sub_nt = line[1]
					if base_nt== focal_nt and sub_nt == focal_sub_nt:
						for a in range(2,len(line)):
							qual_val = int(qual_list[a])
							prop_val = float(round_to_n_sig_figs(float(line[a]),4))
							focal_nt_count = focal_nt_obs_dict[base_nt+str(qual_val)]
							sub_nt_count = nt_count_dict[base_nt+sub_nt+str(qual_val)]
							if prop_val > 0. and qual_val >= min_qual_val_to_fit and focal_nt_count >= 1e4 and sub_nt_count>=10:
								fit_sub_dict[focal_nt+focal_sub_nt+'-'+str(qual_val)] = prop_val
								qual_val_list[qual_val] = ''
			prop_sub_infile.close()
	qual_val_list = list(qual_val_list.keys())
	qual_val_list = sorted(qual_val_list)

	for qual_val in range(min_qual_val_to_fit,max(41,max(qual_val_list)+1)):
		for p in range(0,len(bases)):
			focal_nt = bases[p]
			for h in range(0,len(bases)):
				focal_sub_nt = bases[h]
				try:
					prop = fit_sub_dict[focal_nt+focal_sub_nt+'-'+str(qual_val)]
				except:
					fit_sub_dict[focal_nt+focal_sub_nt+'-'+str(qual_val)] = 0.0
			# prop_total = 0.
				# prop_total += prop
			# if prop_total > 0.:
			# 	remainder = 1.0 - prop_total
			# 	fit_sub_dict[focal_nt+focal_nt+'-'+str(qual_val)] = remainder
		# qual_val = qual_val_list[q]
	outlines = '\t'
	for qual_val in range(min_qual_val_to_fit,max(41,max(qual_val_list)+1)):
		outlines += '\t'+str(qual_val)
	outlines += '\n'
	for p in range(0,len(bases)):
		focal_nt = bases[p]
		prop_total = 0.
		for h in range(0,len(bases)):
			focal_sub_nt = bases[h]
			outlines += focal_nt+'\t'+focal_sub_nt
			for qual_val in range(min_qual_val_to_fit,max(41,max(qual_val_list)+1)):
			# for q in range(0,len(qual_val_list)):
			# 	qual_val = qual_val_list[q]
				try:
					prop = fit_sub_dict[focal_nt+focal_sub_nt+'-'+str(qual_val)]
				except:
					prop = 0.0
				outlines += '\t'+str(prop)
			outlines += '\n'
	fit_sub_outfile = open(fit_prop_sub_outfile_name,'w')
	fit_sub_outfile.write(outlines)
	fit_sub_outfile.close()
	return fit_sub_dict


def load_fit_qual(prop_infile_name):
	sub_dict = {}
	sub_infile = open(prop_infile_name,"r")
	first_line = True
	for line in sub_infile:
		line = line.rstrip('\n').split("\t")
		if first_line == True:
			qual_list = line
			first_line = False
		else:
			base_nt = line[0]
			sub_nt = line[1]
			for a in range(2,len(line)):
				qual_val = int(qual_list[a])
				prop_val = float(line[a])
				sub_dict[base_nt+sub_nt+'-'+str(qual_val)] = prop_val
	sub_infile.close()
	return sub_dict


def barcode_from_fullseq(barcode_dict,seq_in):
	barcode = ''
	site_list = []
	for barcode_site in barcode_dict:
		site_list.append(barcode_site)
	site_list = sorted(site_list)
	for s in range(0,len(site_list)):
		barcode_site = site_list[s]
		barcode += seq_in[barcode_site]
	return barcode

def find_all_possible_barcode_combinations(barcode_sites):
	full_barcode_list = []
	barcode_site_list = sorted(list(set(list(barcode_sites.keys()))))
	for b in range(0,len(barcode_site_list)):
		site = barcode_site_list[b]
		full_barcode_list = list(set(full_barcode_list))
		nt_list = barcode_sites[site]
		building_barcodes = []
		if full_barcode_list == []:
			for num in range(0,len(nt_list)): 
				building_barcodes.append(nt_list[num])
		else:
			for bar_num in range(0,len(full_barcode_list)):
				for nt_num in range(0,len(nt_list)):
					growing_nt_string = full_barcode_list[bar_num]+nt_list[nt_num]
					building_barcodes.append(growing_nt_string)
		full_barcode_list = building_barcodes
	full_barcode_list = sorted(list(set(full_barcode_list)))
	return full_barcode_list


def screen_barcode_site_qual(seq_in,qual_in,seq_expect,barcode_dict,F_edge_ignore,R_edge_ignore):
	barcode_mismatches = 0
	backbone_mismatches = 0
	# high_qual_mismatches = 0
	mismatch_sites = []
	expected_barcode_nt_list = []
	barcode_string = ''
	qual_list = []
	for site in range(0+F_edge_ignore,len(seq_expect)-R_edge_ignore):
		try:
			expected_barcode_nt_list = barcode_dict[site]
			barcode_site = True
		except:
			expected_barcode_nt_list = []
			barcode_site = False
		site_nt = seq_in[site]
		site_qual = ord(qual_in[site])-33
		if barcode_site == False:
			if site_nt != seq_expect[site]:
				backbone_mismatches += 1
				mismatch_sites.append(site)
		elif barcode_site == True:
			barcode_string += site_nt
			qual_list.append(site_qual)
			if site_nt != expected_barcode_nt_list[0] and site_nt != expected_barcode_nt_list[1]:
				barcode_mismatches += 1
				mismatch_sites.append(site)
	qual_list = tuple(qual_list)
	out_tup = (barcode_string,qual_list)
	# total_mismatches = barcode_mismatches+backbone_mismatches
	mismatch_sites = list(set(mismatch_sites))
	return out_tup,barcode_mismatches,backbone_mismatches,mismatch_sites#,high_qual_mismatches

# def read_in_barcode_qual_data(filepath,minQ_val=15):
# 	output_list = []
# 	infile = open(filepath)
# 	for line in infile:
# 		line = line.strip().split('\t')
# 		B = line[0]
# 		Qlist = line[1].split(',')
# 		Q = []
# 		for q in range(0,len(Qlist)):
# 			Q.append(int(Qlist[q]))
# 		Q = tuple(Q)
# 		min_Q = min(Q)
# 		max_Q = max(Q)
# 		Qsum = np.sum(Q)
# 		C = int(line[2])
# 		if min_Q >= minQ_val:
# 			out_tup = (C,min_Q,Qsum,B,Q)
# 			output_list.append(out_tup)
# 	return output_list

def cal_hamming_dist(seq1,seq2):
	mismatch_count = 0
	if len(seq1)!=len(seq2):
		sys.exit("Cannot calculate hamming distance between strings with different lengths: "+seq1+' '+seq2)
	for i in range(0,len(seq1)):
		if seq1[i] != seq2[i]:
			mismatch_count+=1
	return mismatch_count

def abundance_pval(intup_i, intup_j): #i is the focal barcode, j is the higher abundance barcode
	barcode_i, qual_tup_i, count_i = intup_i[3],intup_i[4],intup_i[0]
	barcode_j, count_j = intup_j[0],intup_j[1]
	p_list = []
	for l in range(0,len(barcode_i)):
		nt_i = barcode_i[l]
		nt_j = barcode_j[l]
		q_i = qual_tup_i[l]
		p = sub_prob_dict[nt_j+nt_i+'-'+str(q_i)]
		p_list.append(p)
	y_ij = np.prod(p_list)
	val1 = scipy.stats.poisson.pmf(k=0,mu=y_ij*count_j)
	val2 = scipy.stats.poisson.cdf(k=count_i,mu=y_ij*count_j)
	if val1 == 1.0:
		pA = 0.0
	else:
		pA = (1.0/(1.0-val1)) * (1.0-val2)
	if str(pA) == 'nan':
		pA = 1.0
	else:
		pA = float(round_to_n_sig_figs(pA,4))
	return pA

def divisive_partition(input_barcode_list, full_barcode_dict, pval_thresh=1e-40,max_hamming_dist=3):
	new_seed_counter = 0
	ranked_barcode_list = sorted(input_barcode_list,reverse=True)
	partition_dict = {}
	partition_count_dict = {}
	consolodated_dict = {}
	invalid_barcode_list = {}
	debug_lines = ''
	first_seed_barcode = ''
	active_ranked_barcode_list = []

	#initiate first seed barcode
	for i in range(0,len(ranked_barcode_list)):
		query_tup = ranked_barcode_list[i]
		query_C = query_tup[0]
		query_Qsum = query_tup[2]
		query_B = query_tup[3]
		query_Q = query_tup[4]
		valid_barcodeID = False
		try:
			full_barcode_dict[query_B]
			valid_barcodeID = True
		except:
			pass
		if valid_barcodeID == True:
			if partition_dict == {} and query_C > 1:
				partition_dict[query_B] = [query_tup]
				partition_count_dict[query_B] = query_C
				consolodated_dict[query_tup] = 'first-seed'
				first_seed_barcode = query_B
				new_seed_counter += 1
			else:
				if query_B == first_seed_barcode:
					partition_dict[query_B].append(query_tup)
					partition_count_dict[query_B] += query_C
					consolodated_dict[query_tup] = 'first-seed-collapse'
				elif query_C >= 1:
					active_ranked_barcode_list.append(query_tup)
	active_ranked_barcode_list = sorted(active_ranked_barcode_list,reverse=True)

	#Find subsequent valid barcodes and initiate seeds
	if len(active_ranked_barcode_list)==0:
		return partition_count_dict
	still_iterating = True
	iteration_count = 0
	while still_iterating == True:
		iteration_count +=1
		seed_barcode_list = list(partition_dict.keys())
		pA_rank_list = []
		new_seed_B = ''
		if len(active_ranked_barcode_list) >0:
			minimum_pA_threshold = pval_thresh/(float(len(active_ranked_barcode_list)))
			for i in range(0,len(active_ranked_barcode_list)):
				query_tup = active_ranked_barcode_list[i]
				query_C = query_tup[0]
				query_Qsum = query_tup[2]
				query_B = query_tup[3]
				query_Q = query_tup[4]
				valid_barcodeID = False
				try:
					full_barcode_dict[query_B]
					valid_barcodeID = True
				except:
					pass
				consolodated_barcodeID = False
				try:
					consolodated_dict[query_tup]
					consolodated_barcodeID= True
				except:
					pass
				if valid_barcodeID == True and query_C > 1 and consolodated_barcodeID == False:
					query_pA_list = []
					for b in range(0,len(seed_barcode_list)):
						seed_B = seed_barcode_list[b]
						seed_C = partition_count_dict[seed_B]
						seed_tup = (seed_B,seed_C)
						if seed_B != query_B:
							hdist = cal_hamming_dist(query_B,seed_B)
							if seed_C > query_C and hdist <= max_hamming_dist:
								pA = abundance_pval(query_tup, seed_tup)
							else:
								pA = 0.0
							dist_val = 1.0/hdist
							count_val = 1.0/query_C
							temp_tup = (pA,count_val,dist_val,seed_B)
							query_pA_list.append(temp_tup)
					query_pA_list = sorted(query_pA_list,reverse=True)
					min_tup = query_pA_list[0]
					min_pA = min_tup[0]
					min_seed_B = min_tup[3]
					min_hdist = cal_hamming_dist(query_B,min_seed_B)
					min_seed_C = partition_count_dict[min_seed_B]
					temp_tup = (min_pA,1.0/query_C,1.0/min_hdist,1.0/min_seed_C,min_seed_B,query_tup)
					pA_rank_list.append(temp_tup)
		if len(pA_rank_list) > 0:
			minimum_pA_threshold = pval_thresh/float(len(pA_rank_list))
			pA_rank_list = sorted(pA_rank_list,reverse=False)
			max_query_tup = pA_rank_list[0][-1]
			max_query_pA = pA_rank_list[0][0]
			max_query_B = max_query_tup[3]
			max_query_C = max_query_tup[0]
			max_query_Q = max_query_tup[4]
			valid_query_barcodeID = False
			try:
				full_barcode_dict[max_query_B]
				valid_query_barcodeID = True
			except:
				pass
			new_active_ranked_barcode_list = []
			if max_query_C > 1 and max_query_pA<=minimum_pA_threshold and valid_query_barcodeID == True:
				partition_dict[max_query_B] = [max_query_tup]
				partition_count_dict[max_query_B] = max_query_C
				consolodated_dict[max_query_tup] = 'new-seed'
				new_seed_B = max_query_B
				new_seed_counter += 1
				for b in range(0,len(active_ranked_barcode_list)):
					query_tup = active_ranked_barcode_list[b]
					query_C = query_tup[0]
					# query_Qsum = query_tup[2]
					query_B = query_tup[3]
					query_Q = query_tup[4]
					
					consolodated_query_barcodeID = False
					try:
						consolodated_dict[query_tup]
						consolodated_query_barcodeID = True
					except:
						pass
					if query_B == new_seed_B and query_tup != max_query_tup and consolodated_query_barcodeID == False:
						partition_dict[query_B].append(query_tup)
						partition_count_dict[query_B] += query_C
						consolodated_dict[query_tup] = 'new-seed-collapse'
					elif query_B != new_seed_B:
						new_active_ranked_barcode_list.append(query_tup)
			else:
				new_active_ranked_barcode_list = active_ranked_barcode_list
			new_active_ranked_barcode_list = list(set(new_active_ranked_barcode_list))
			if len(new_active_ranked_barcode_list) == len(active_ranked_barcode_list):
				still_iterating = False
			active_ranked_barcode_list = new_active_ranked_barcode_list
			active_ranked_barcode_list = sorted(active_ranked_barcode_list,reverse=True)
		else:
			active_ranked_barcode_list = []
			still_iterating = False
	for i in range(0,len(ranked_barcode_list)):
		query_tup = ranked_barcode_list[i]
		query_B = query_tup[3]
		if len(query_B) == len(query_B.replace('N','')):	
			try:
				consolodated_dict[query_tup]
			except:
				active_ranked_barcode_list.append(query_tup)
	active_ranked_barcode_list = list(set(active_ranked_barcode_list))
	active_ranked_barcode_list = sorted(active_ranked_barcode_list,reverse=True)

	#agglomerate remaining valid barcodes
	seed_barcode_list = list(partition_dict.keys())
	try:
		minimum_pA_threshold = pval_thresh/(float(len(active_ranked_barcode_list)))
	except:
		minimum_pA_threshold = 1.
	for i in range(0,len(active_ranked_barcode_list)):
		query_tup = active_ranked_barcode_list[i]
		query_C = query_tup[0]
		query_Qsum = query_tup[2]
		query_B = query_tup[3]
		query_Q = query_tup[4]
		try:
			full_barcode_dict[query_B]
		except:
			pass
		consolodated_barcodeID = False
		try:
			consolodated_dict[query_tup]
			consolodated_barcodeID = True
		except:
			pass
		valid_barcodeID = False
		try:
			full_barcode_dict[query_B]
			valid_barcodeID = True
		except:
			pass
		if consolodated_barcodeID == False:
			query_pA_list = []
			for b in range(0,len(seed_barcode_list)):
				seed_B = seed_barcode_list[b]
				seed_C = partition_count_dict[seed_B]
				seed_tup = (seed_B,seed_C)#partition_dict[seed_B][0]
				hdist = cal_hamming_dist(query_B,seed_B)
				if seed_B != query_B:
					if query_C < seed_C and hdist <= max_hamming_dist:
						pA = abundance_pval(query_tup, seed_tup)
					else:
						pA = 0.0
					temp_tup = (pA,1.0/hdist,0,seed_B)
					query_pA_list.append(temp_tup)
			if len(query_pA_list) >= 1:
				query_pA_list = sorted(query_pA_list,reverse=True)
				min_tup = query_pA_list[0]
				min_pA = min_tup[0]
				min_seed_B = min_tup[3]
				min_hdist = cal_hamming_dist(query_B,min_seed_B)
				min_seed_C = partition_count_dict[min_seed_B]
				minimum_pA_threshold = pval_thresh/(float(len(seed_barcode_list)))
				if min_pA>=minimum_pA_threshold and query_C > 1 and min_hdist <= max_hamming_dist:
					partition_dict[min_seed_B].append(query_tup)
					count_before = partition_count_dict[min_seed_B]
					partition_count_dict[min_seed_B] += query_C
					consolodated_dict[query_tup] = 'seed-collapse'
				elif valid_barcodeID == True:
					try:
						partition_dict[query_B].append(query_tup)
						count_before = partition_count_dict[query_B]
						partition_count_dict[query_B] += query_C
						consolodated_dict[query_tup] = 'remain-add'
					except:
						partition_dict[query_B] = [query_tup]
						partition_count_dict[query_B] = query_C
						consolodated_dict[query_tup] = 'remain'
			elif valid_barcodeID == True:
				try:
					partition_dict[query_B].append(query_tup)
					count_before = partition_count_dict[query_B]
					partition_count_dict[query_B] += query_C
					consolodated_dict[query_tup] = 'last-add'
				except:
					partition_dict[query_B] = [query_tup]
					partition_count_dict[query_B] = query_C
					consolodated_dict[query_tup] = 'last'
	return partition_count_dict

def raw_read_processing_single(accession,command_prefix,raw_reads_forward,trim_dir,temp_dir,detect_barcode_content_first,amplicon_length,force_trim_pre_merge,truncate_length_post_merge):
	global verbose
	if verbose == True:
		update_print_line("Processing input reads",accession)
	trim_fastq_filename = temp_dir+accession+".trim.fq"
	clean_filename_out = trim_dir+accession+".fq"
	if detect_barcode_content_first == True or amplicon_length == 0:
		read_len_expected = calc_expect_read_len(raw_reads_forward)
	else:
		read_len_expected = amplicon_length
	command = command_prefix+'bbduk.sh in="'+raw_reads_forward+'" out="'+trim_fastq_filename+'" qin=33 maq=20 minlength='+str(read_len_expected)+' maxlength='+str(read_len_expected)+' overwrite=t'
	out = run_command(command+ ' overwrite=t threads=2 -Xmx1000m')
	if out == False:
		sys.exit()
	time.sleep(1)
	os.rename(trim_fastq_filename,clean_filename_out)
	return clean_filename_out

def truncate_merged_reads(filename_in,truncate_length):
	infile = open(filename_in,"r")
	outfile = open(filename_in.replace(".fq",".temp.fq"),"w")
	line_counter = -1
	seq_count = 0
	for line in infile:
		line = line.strip()
		line_counter += 1
		if line_counter == 0:
			header = line
			outfile.write(line+"\n")
		elif line_counter == 1: #base calls
			nucleotide_line = line[0:min(truncate_length,len(line))]
			outfile.write(nucleotide_line+"\n")
		elif line_counter == 2:
			outfile.write(line+"\n")
		elif line_counter == 3: #quality scores
			qual_line = line[0:min(truncate_length,len(line))]
			outfile.write(qual_line+"\n")
			line_counter = -1
	outfile.close()
	os.remove(filename_in)
	os.rename(filename_in.replace(".fq",".temp.fq"),filename_in)


def raw_read_processing_paired(accession,command_prefix,raw_reads_forward,raw_reads_reverse,trim_dir,temp_dir,detect_barcode_content_first,amplicon_length,force_trim_pre_merge,truncate_length_post_merge):
	global verbose
	if verbose == True:
		update_print_line("Processing paired input reads",accession)
	repair_fastq_filename = temp_dir+accession+".repair.fq"
	temp_adapter_filename = temp_dir+accession+".adapters.fa"
	pretrim_fastq_filename = temp_dir+accession+".merge_pretrim.fq"
	merged_trimmed_filename = temp_dir+accession+".merge.fq"
	clean_filename_out = trim_dir+accession+".fq"

	if detect_barcode_content_first == True or amplicon_length == 0:
		median_forward_read_len = calc_expect_read_len(raw_reads_forward)
		median_reverse_read_len = calc_expect_read_len(raw_reads_reverse)
		diff_median_len = np.absolute(median_forward_read_len-median_reverse_read_len)
		read_len_expected = max(median_forward_read_len,median_reverse_read_len)
	else:
		read_len_expected = amplicon_length
		diff_median_len = 0
	command = command_prefix+'repair.sh in="'+raw_reads_forward+'" in2="'+raw_reads_reverse+'" out="'+repair_fastq_filename+'" fint=t repair=t threads=1 overwrite=t -da -Xmx1000m'
	out = run_command(command)
	if out == False:
		sys.exit("repair.sh failed: "+accession)
	time.sleep(1)

	len_diff_buffer = 50
	command = command_prefix+'bbmerge.sh in="'+repair_fastq_filename+'" outadapter="'+temp_adapter_filename+'" minlength='+str(read_len_expected-(len_diff_buffer+diff_median_len))+' maxlength='+str(read_len_expected+(len_diff_buffer+diff_median_len))+' overwrite=t' #forcetrimleft=9
	if force_trim_pre_merge >0:
		command += " forcetrimleft="+str(force_trim_pre_merge)#+" forcetrimright=130"
	out = run_command(command+" threads=1 -Xmx2000m")
	if out == False:
		sys.exit("bbmerge.sh adapter screen failed: "+accession)
	time.sleep(1)

	adapter_pass = check_adapter(temp_adapter_filename)
	if adapter_pass == True:
		command = command_prefix+'bbmerge.sh in="'+repair_fastq_filename+'" out="'+pretrim_fastq_filename+'" qin=33 mix=f adapters="'+temp_adapter_filename+'" minlength='+str(read_len_expected-(len_diff_buffer+diff_median_len))+' maxlength='+str(read_len_expected+(len_diff_buffer+diff_median_len)) #forcetrimleft=9 
		if force_trim_pre_merge >0:
			command += " forcetrimleft="+str(force_trim_pre_merge)#+" forcetrimright=130"
		out = run_command(command+" overwrite=t threads=1 -Xmx2000m")
		if out == False:
			sys.exit("bbmerge.sh merge with adapter trim failed: "+accession)
		time.sleep(1)
	else:
		command = command_prefix+'bbmerge.sh in="'+repair_fastq_filename+'" out="'+pretrim_fastq_filename+'" qin=33 mix=f minlength='+str(read_len_expected-(len_diff_buffer+diff_median_len))+' maxlength='+str(read_len_expected+(len_diff_buffer+diff_median_len))+' overwrite=t' #forcetrimleft=9 
		if force_trim_pre_merge >0:
			command += " forcetrimleft="+str(force_trim_pre_merge)#+" forcetrimright=130"
		out = run_command(command+" overwrite=t threads=1 -Xmx2000m")
		if out == False:
			sys.exit("bbmerge.sh merge failed")
		time.sleep(1)

	if amplicon_length == 0:
		median_merged_len = calc_expect_read_len(pretrim_fastq_filename)
	else:
		median_merged_len = amplicon_length
	command = command_prefix+'bbduk.sh in="'+pretrim_fastq_filename+'" out="'+merged_trimmed_filename+'" qin=33 maq=20 minlength='+str(median_merged_len)+' maxlength='+str(median_merged_len)+' tpe=f tbo=f'
	out = run_command(command+ ' overwrite=t threads=1 -Xmx1000m')
	if out == False:
		sys.exit("bbduk.sh failed: "+accession)
	time.sleep(1)

	if truncate_length_post_merge > 0:
		truncate_merged_reads(merged_trimmed_filename,truncate_length_post_merge)
	
	os.remove(repair_fastq_filename)
	os.remove(temp_adapter_filename)
	os.remove(pretrim_fastq_filename)
	os.rename(merged_trimmed_filename,clean_filename_out)
	
	return clean_filename_out

def subsample_read_counts(barcode_qual_count_dict_in,subsample_count):
	rand_list = []
	for barcode_qual_tup in barcode_qual_count_dict_in:
		barcode_string = barcode_qual_tup[0]
		barcode_qual_list = barcode_qual_tup[1]
		count = barcode_qual_count_dict_in[barcode_qual_tup]
		for i in range(0,count):
			rand_list.append(barcode_qual_tup)
	random.shuffle(rand_list)
	subset_list = rand_list[0:int(subsample_count)]
	barcode_qual_count_dict_out = {}
	for i in range(0,len(subset_list)):
		barcode_qual_tup = subset_list[i]
		try:
			barcode_qual_count_dict_out[barcode_qual_tup] += 1
		except:
			barcode_qual_count_dict_out[barcode_qual_tup] = 1
	return barcode_qual_count_dict_out

def barcode_screen_trimmed_reads(accession,sequence_filename,expected_amplicon_seq,barcode_sites,min_Q_score,temp_dir,output_dir,F_edge_ignore,R_edge_ignore):
	reads_per_sample_ceiling = 1e5
	### Load in processed reads and collect sequence to quality info to screen for barcodes
	global verbose
	if verbose == True:
		update_print_line("Screening reads for barcodes",accession)
	barcode_qual_count_filename = trim_dir+accession+'.barcode_qual_counts.txt'
	# if os.path.isfile() == True and overwrite_existing_files == False:
	# 	barcode_qual_list = read_in_barcode_qual_data(barcode_qual_count_filename)
	# else:
	# mismatched_barcode_count_dict = {}
	# high_qual_nonbarcode_count_dict = {}
	fastq_infile = open(sequence_filename,"r")
	mismatch_by_site = {}
	barcode_qual_count_dict = {}
	uncorrected_barcode_count_dict = {}
	barcode_list = []
	line_counter = -1
	seq_count = 0
	pass_seq_count = 0
	for line in fastq_infile:
		line = line.strip()
		line_counter += 1
		if line_counter == 0:
			header = line
			writing = False
		elif line_counter == 1: #base calls
			nucleotide_line = line
		elif line_counter == 3: #quality scores
			qual_line = line
			if len(nucleotide_line) == len(expected_amplicon_seq):
				seq_count += 1
				if verbose == True:
					if seq_count%1e5==0 and seq_count > 1:
						update_print_line("Reads screened: "+str(seq_count),accession)
						# print('\tReads screened from: "'+accession+'":\t'+str(seq_count))
				# mismatch_count,read_mismatch_locs,high_qual_mismatches = screen_barcode_site_qual(nucleotide_line,qual_line,expected_amplicon_seq,barcode_sites,min_qual_barcode_sites,min_qual_backbone,F_edge_ignore,R_edge_ignore)#(seq_in,qual_in,seq_expect,barcode_dict,min_barcode_site_qual,min_other_site_qual)
				read_barcode_tup,barcode_mismatches,backbone_mismatches,mismatch_list = screen_barcode_site_qual(nucleotide_line,qual_line,expected_amplicon_seq,barcode_sites,F_edge_ignore,R_edge_ignore)#(seq_in,qual_in,seq_expect,barcode_dict,min_barcode_site_qual,min_other_site_qual)
				if min(read_barcode_tup[1]) >= min_Q_score:
					total_mismatches = barcode_mismatches+backbone_mismatches
					if backbone_mismatches == 0 and barcode_mismatches <= 1:
						try:
							barcode_qual_count_dict[read_barcode_tup] += 1
						except:
							barcode_qual_count_dict[read_barcode_tup] = 1
						pass_seq_count += 1
						if barcode_mismatches == 0:
							try:
								uncorrected_barcode_count_dict[read_barcode_tup[0]] += 1
							except:
								uncorrected_barcode_count_dict[read_barcode_tup[0]] = 1
					if len(mismatch_list)>0:
						for m in range(0,len(mismatch_list)):
							try:
								mismatch_by_site[mismatch_list[m]] += 1
							except:
								mismatch_by_site[mismatch_list[m]] = 1
			line_counter = -1
	# barcode_list = sorted(list(set(barcode_list)))
	fastq_infile.close()

	if pass_seq_count > reads_per_sample_ceiling:
		barcode_qual_count_dict = subsample_read_counts(barcode_qual_count_dict,reads_per_sample_ceiling)
		if verbose == True:
			update_print_line("Subsampled read counts to "+str(int(reads_per_sample_ceiling)),accession)

	barcode_qual_list = []
	outstring = ''
	for read_barcode_tup in barcode_qual_count_dict:
		barcode_string = read_barcode_tup[0]
		barcode_qual_tup = read_barcode_tup[1]
		count = barcode_qual_count_dict[read_barcode_tup]
		min_Q = min(barcode_qual_tup)
		max_Q = max(barcode_qual_tup)
		Qsum = np.sum(barcode_qual_tup)
		if min_Q >= min_Q_score:
			out_tup = (count,min_Q,Qsum,barcode_string,barcode_qual_tup)
			barcode_qual_list.append(out_tup)
		outstring += barcode_string+'\t'+str(barcode_qual_tup[0])
		for i in range(1,len(barcode_qual_tup)):
			outstring += ','+str(barcode_qual_tup[i])
		outstring += '\t'+str(count)+'\n'
	barcode_qual_count_file = open(barcode_qual_count_filename,'w')
	barcode_qual_count_file.write(outstring)
	barcode_qual_count_file.close()

	full_barcode_list = find_all_possible_barcode_combinations(barcode_sites)
	full_barcode_dict = {}
	for f in range(0,len(full_barcode_list)):
		barcode = full_barcode_list[f]
		full_barcode_dict[barcode] = None

	barcode_count_dict = divisive_partition(barcode_qual_list,full_barcode_dict)
	#### count unique barcodes
	raw_barcode_list = sorted(list(uncorrected_barcode_count_dict.keys()))
	raw_barcode_outfile = open(temp_dir+accession+".raw_barcode_count.txt","w")
	raw_barcode_outfile.write("barcode\tcount\n")
	for num in range(0,len(raw_barcode_list)):
		barcode_string = raw_barcode_list[num]
		raw_barcode_count = uncorrected_barcode_count_dict[barcode_string]
		raw_barcode_outfile.write(str(barcode_string)+"\t"+str(raw_barcode_count)+"\n")
	raw_barcode_outfile.close()

	barcode_list = sorted(list(barcode_count_dict.keys()))
	barcode_outfile = open(temp_dir+accession+".barcode_count.txt","w")
	barcode_outfile.write("barcode\tcount\n")
	for num in range(0,len(barcode_list)):
		barcode_string = barcode_list[num]
		barcode_count = barcode_count_dict[barcode_string]
		barcode_outfile.write(str(barcode_string)+"\t"+str(barcode_count)+"\n")
	barcode_outfile.close()

	### summarize mismatch count
	mismatch_prop_by_site = {}
	mismatch_by_site_outfile = open(temp_dir+accession+".mismatch_info.txt","w")
	mismatch_by_site_outfile.write("Site_number\tmismatches\n")
	for col in range(0,len(expected_amplicon_seq)):
		try:
			mismatch_count = mismatch_by_site[col]
		except:
			mismatch_count = 0
		if mismatch_count > 0:
			mismatch_prop = mismatch_count/seq_count
			mismatch_prop_by_site[col] = mismatch_prop
			mismatch_by_site_outfile.write(str(col)+"\t"+str(mismatch_prop)+"\n")
	mismatch_by_site_outfile.close()

	### write file with full read count
	read_count_outfile = open(temp_dir+accession+".total_read_count.txt","w")
	read_count_outfile.write(str(pass_seq_count)+'\n')
	read_count_outfile.close()

	# #### count high-quality non-barcode sites
	# nonbarcode_outfile = open(temp_dir+accession+".nonbarcode_count.txt","w")
	# nonbarcode_outfile.write("non-barcode_seq\tcount\n")
	# for nucl_string in high_qual_nonbarcode_count_dict:
	# 	nonbarcode_count = high_qual_nonbarcode_count_dict[nucl_string]
	# 	nonbarcode_outfile.write(str(nucl_string)+"\t"+str(nonbarcode_count)+"\n")
	# nonbarcode_outfile.close()
	# os.rename(temp_dir+accession+".nonbarcode_count.txt",output_dir+"samples/"+accession+".nonbarcode_count.txt")

	os.rename(temp_dir+accession+".mismatch_info.txt",output_dir+"samples/"+accession+".mismatch_info.txt")
	os.rename(temp_dir+accession+".barcode_count.txt",output_dir+"samples/"+accession+".barcode_count.txt")
	os.rename(temp_dir+accession+".raw_barcode_count.txt",output_dir+"samples/"+accession+".raw_barcode_count.txt")
	os.rename(temp_dir+accession+".total_read_count.txt",output_dir+"samples/"+accession+".total_read_count.txt")

	return barcode_count_dict,barcode_list,mismatch_prop_by_site#,high_qual_nonbarcode_count_dict


def clean_files_for_one_sample(accession,output_dir,trim_dir):
	sample_barcode_info_filename = output_dir+"samples/"+accession+".barcode_count.txt"
	barcode_mismatch_info_filename = output_dir+"samples/"+accession+".mismatch_info.txt"
	sample_nonbarcode_info_filename = output_dir+"samples/"+accession+".nonbarcode_count.txt"
	processed_read_filename = trim_dir+accession+".fq"
	
	#remove any files that exist
	if os.path.isfile(sample_barcode_info_filename) == True:
		os.remove(sample_barcode_info_filename)
	if os.path.isfile(barcode_mismatch_info_filename) == True:
		os.remove(barcode_mismatch_info_filename)
	if os.path.isfile(sample_nonbarcode_info_filename) == True:
		os.remove(sample_nonbarcode_info_filename)
	if os.path.isfile(processed_read_filename) == True:
		os.remove(processed_read_filename)

def update_print_line(string_in,accession):
	global longest_accession_length
	print(accession+" "*(longest_accession_length-len(accession))+" - "+string_in)

##########################      Core pipeline function      ##########################

def reads_to_barcodes_one_sample(accession,command_prefix,forward_read_suffix,reverse_read_suffix,input_read_dir,trim_dir,temp_dir,output_dir,expected_amplicon_seq,barcode_sites,min_Q_score,overwrite_existing_files,raw_read_input_type,reads_already_screened_for_quality,detect_barcode_content_first,F_edge_ignore,R_edge_ignore,force_trim_pre_merge,truncate_length_post_merge,trim_all_files_for_qual_learning=False):
	# update_print_line("Sample processing start",accession)
	processed_something = False
	raw_reads_forward = input_read_dir+accession+forward_read_suffix
	raw_reads_reverse = input_read_dir+accession+reverse_read_suffix

	amplicon_length = len(expected_amplicon_seq)
	sample_barcode_info_filename = output_dir+"samples/"+accession+".barcode_count.txt"
	barcode_mismatch_info_filename = output_dir+"samples/"+accession+".mismatch_info.txt"
	sample_nonbarcode_info_filename = output_dir+"samples/"+accession+".nonbarcode_count.txt"
	processed_read_filename = trim_dir+accession+".fq"
	#clear the temporary files folder from any previously interrupted runs
	temp_files_list = [f for f in os.listdir(temp_dir)]
	for files in temp_files_list:
		if accession+"." in files:
			os.remove(temp_dir+files)
	#remove any files that already exist if 'overwrite_existing_files' is set to 'True'
	if overwrite_existing_files == True:
		if os.path.isfile(sample_barcode_info_filename) == True:
			os.remove(sample_barcode_info_filename)
		if os.path.isfile(barcode_mismatch_info_filename) == True:
			os.remove(barcode_mismatch_info_filename)
		if os.path.isfile(sample_nonbarcode_info_filename) == True:
			os.remove(sample_nonbarcode_info_filename)
		if os.path.isfile(processed_read_filename) == True:
			os.remove(processed_read_filename)
	if reads_already_screened_for_quality == False:
		if os.path.isfile(processed_read_filename) == False:
			if raw_read_input_type == "paired":
				processed_read_filename = raw_read_processing_paired(accession,command_prefix,raw_reads_forward,raw_reads_reverse,trim_dir,temp_dir,detect_barcode_content_first,amplicon_length,force_trim_pre_merge,truncate_length_post_merge)
				processed_something = True
			elif raw_read_input_type == "single":
				processed_read_filename = raw_read_processing_single(accession,command_prefix,raw_reads_forward,trim_dir,temp_dir,detect_barcode_content_first,amplicon_length,force_trim_pre_merge,truncate_length_post_merge)
				processed_something = True
	elif reads_already_screened_for_quality == True:
		processed_read_filename = raw_reads_forward
	if trim_all_files_for_qual_learning == True:
		qual_summary_tup = pull_qual_info(processed_read_filename)
		return qual_summary_tup
	else:
		if os.path.isfile(sample_barcode_info_filename) == False:
			barcode_count_dict,barcode_list,mismatch_by_site = barcode_screen_trimmed_reads(accession,processed_read_filename,expected_amplicon_seq,barcode_sites,min_Q_score,temp_dir,output_dir,F_edge_ignore,R_edge_ignore) #,high_qual_nonbarcode_dict
			processed_something = True
		else:
			barcode_count_dict,barcode_list,mismatch_by_site = read_in_barcodes(sample_barcode_info_filename,barcode_mismatch_info_filename,sample_nonbarcode_info_filename) #,high_qual_nonbarcode_dict
		out_tup = (accession,barcode_count_dict,barcode_list,mismatch_by_site)#,high_qual_nonbarcode_dict)
		if processed_something == True:
			update_print_line("Completed",accession)
		return out_tup


#####################################   MAIN   ######################################

#Create output and temporary working directories if they do not exist already
if not os.path.exists(trim_dir):
	os.makedirs(trim_dir)
if not os.path.exists(output_dir):
	os.makedirs(output_dir)
if not os.path.exists(output_dir+"samples/"):
	os.makedirs(output_dir+"samples/")
if not os.path.exists(temp_dir):
	os.makedirs(temp_dir)

#Set up parallel processing if enabled
if parallel_process == True:
	import multiprocessing
	from joblib import Parallel, delayed
	if parallel_max_cpu == 0:
		num_cores = multiprocessing.cpu_count()
	elif parallel_max_cpu > 0:
		num_cores = min(parallel_max_cpu,multiprocessing.cpu_count())
	elif parallel_max_cpu < 0:
		num_cores = multiprocessing.cpu_count()-parallel_max_cpu
else:
	num_cores = 1

### Check if the files in the input directory need to be unzipped
files_found = check_if_input_zipped(input_read_dir,parallel_process,num_cores)
if files_found == False:
	sys.exit("No input files with '.fastq', '.fq', '.gz', or '.tar.gz' extensions found. Exiting.")

### Predict forward and reverse read suffix if left blank
continue_running = False
if forward_read_suffix == '':
	suffix_info_file_path = temp_dir+"file_extensions.txt"
	try:
		file_suffix_infile = open(suffix_info_file_path)
		for line in file_suffix_infile:
			line = line.strip().split("\t")
			if line[0] == "forward":
				forward_read_suffix = line[1]
			elif line[0] == "reverse":
				try:
					reverse_read_suffix = line[1]
				except:
					reverse_read_suffix = ''
			continue_running = True
	except:
		if raw_read_input_type == "paired" and reads_already_screened_for_quality == False:
			forward_read_suffix,reverse_read_suffix = predict_paired_file_extensions(input_read_dir)
			print('Predicted read suffix\n\tForward: "'+str(forward_read_suffix)+'"\n\tReverse: "'+str(reverse_read_suffix)+'"')
			print("Is this correct?")
			decision = input("(yes or no)\n")
			if decision == "Y" or decision == "y" or decision == "Yes" or decision == "yes" or decision == "YES":
				continue_running = True
			else:
				forward_read_suffix = input('Enter forward read suffix (ex "_L001_R1_001.fastq") or "exit" to cancel:\n')
				if forward_read_suffix == "exit":
					sys.exit()
				reverse_read_suffix	= input('Enter Reverse read suffix (ex "_L001_R2_001.fastq") or "exit" to cancel:\n')
				if reverse_read_suffix == "exit":
					sys.exit()
				if forward_read_suffix != "exit" and reverse_read_suffix != "exit":
					continue_running = True
			if continue_running == True:
				suffix_info_file = open(suffix_info_file_path,"w")
				suffix_info_file.write("forward\t"+forward_read_suffix+"\nreverse\t"+reverse_read_suffix+"\n")
				suffix_info_file.close()
		elif raw_read_input_type == "single" or reads_already_screened_for_quality == True:
			forward_read_suffix = predict_single_file_extension(input_read_dir)
			print('Predicted read suffix: "'+str(forward_read_suffix)+'"')
			print("Is this correct?")
			decision = input("(yes or no)\n")
			if decision == "Y" or decision == "y" or decision == "Yes" or decision == "yes" or decision == "YES":
				continue_running = True
			else:
				forward_read_suffix = input('Enter read suffix (ex "_L001_001.fastq") or "exit" to cancel:\n')
				if forward_read_suffix == "exit":
					sys.exit()
				else:
					continue_running = True
			if continue_running == True:
				suffix_info_file = open(suffix_info_file_path,"w")
				suffix_info_file.write("forward\t"+forward_read_suffix+"\nreverse\t"+reverse_read_suffix+"\n")
				suffix_info_file.close()

if raw_read_input_type == "paired":
	if forward_read_suffix != "" and reverse_read_suffix != "":
		continue_running = True
	else:
		continue_running = False
elif raw_read_input_type == "single":
	if forward_read_suffix != "":
		continue_running = True
	else:
		continue_running = False

if continue_running == False:
	sys.exit("No read suffix supplied, please re-run and enter information when prompted")


### Create list of samples to process
file_list = [f for f in os.listdir(input_read_dir) if f.endswith(forward_read_suffix)]
# if len(file_list) == 0:
# 	file_list = [f for f in os.listdir(input_read_dir) if f.endswith(forward_read_suffix+".gz")]
# 	if len(file_list) == 0:
# 		file_list = [f for f in os.listdir(input_read_dir) if f.endswith(forward_read_suffix+".tar.gz")]

accession_list = []
for files in file_list:
	accession = files.split(forward_read_suffix)[0]
	accession_list.append(accession)
accession_list = list(set(accession_list))
	# if accession != 'KH50_02':

# accession_list = []

sorted_accession_list = sorted(accession_list)
print_line_dict = {}
longest_accession_length = 0
for a in range(0,len(sorted_accession_list)):
	accession = sorted_accession_list[a]
	print_line_dict[accession] = a
	if len(accession)>longest_accession_length:
		longest_accession_length = len(accession)

### load barcode info. If unavailable, pick a random file and attempt to predict barcode info
detect_barcode_content_first = True
amplicon_info_infile_path = project_dir+amplicon_info_infile_name
if os.path.isfile(amplicon_info_infile_path):
	expected_amplicon_seq,barcode_sites = load_barcode_info(project_dir+amplicon_info_infile_path)
	if expected_amplicon_seq != "" and barcode_sites != {}:
		detect_barcode_content_first = False
else:
	print("\n\nCould not detect barcode info input file.")
if detect_barcode_content_first == True:
	'''	This section will attempt to predict the barcode info for your experiment if it is not already given.
	It will run through the entire pipeline with this one sample, then you need to check the output to determine
	if it correctly identified the barcode sites. If you agree with the prediction, remove ".auto" from
	the "barcode.info.auto.txt" file that was generated, then re-run the script to process all remaining samples.'''
	overwrite_existing_files = True
	accession_entered = input('Enter sample name for predicting barcode info (or "random" to pick random sample):\n')
	temp_cont = False
	if accession_entered == "random":
		accession = accession_list[random.randint(0,len(accession_list))]
	elif accession_entered in sorted_accession_list:
		accession = accession_entered
	else:
		sys.exit("sample ID not valid.")

	print('Attempting to predict barcode using sample: "'+accession+'"\n')
	num_expected_barcode_sites = input('How many barcode sites should be present in these samples? (enter number below)\n')
	num_expected_barcode_sites = num_expected_barcode_sites.strip()
	if is_number(num_expected_barcode_sites)== True:
		num_expected_barcode_sites = int(num_expected_barcode_sites)
	else:
		sys.exit("Value entered is not a number")
	raw_reads_forward = input_read_dir+accession+forward_read_suffix
	if raw_read_input_type == "paired":
		raw_reads_reverse = input_read_dir+accession+reverse_read_suffix
		processed_read_filename = raw_read_processing_paired(accession,command_prefix,raw_reads_forward,raw_reads_reverse,trim_dir,temp_dir,detect_barcode_content_first,0,force_trim_pre_merge,truncate_length_post_merge)
	elif raw_read_input_type == "single":
		processed_read_filename = raw_read_processing_single(accession,command_prefix,raw_reads_forward,trim_dir,temp_dir,detect_barcode_content_first,0,force_trim_pre_merge,truncate_length_post_merge)
	temp_seq_dict,temp_read_length = subset_fastq(processed_read_filename,10000,500000)
	expected_amplicon_seq,barcode_sites = detect_barcode_content(temp_seq_dict,temp_read_length)
	continue_running = False
	if len(barcode_sites) == num_expected_barcode_sites:
		print('\n\nThe correct number of expected barcode sites were automatically detected.\n')
		continue_running = False
		print("amplicon sequence:\n\t"+expected_amplicon_seq)
		temp_barcode_site_list = list(barcode_sites.keys())
		temp_barcode_site_list = sorted(temp_barcode_site_list)
		for site_num in range(0,len(temp_barcode_site_list)):
			site = temp_barcode_site_list[site_num]
			print("\t"+str(site)+"\t"+barcode_sites[site][0]+" or "+barcode_sites[site][1])
		print('\nAre these sites and possible nucleotides correct?')
		decision = input("(yes or no)\n")
		if decision == "Y" or decision == "y" or decision == "Yes" or decision == "yes" or decision == "YES":
			continue_running = True
			outfile = open(project_dir+"amplicon_info.txt","w")
		else:
			continue_running = False
			clean_files_for_one_sample(accession,output_dir,trim_dir)
			outfile = open(project_dir+"amplicon_info.auto.txt","w")
	else:
		clean_files_for_one_sample(accession,output_dir,trim_dir)
		print('\n\nThe incorrect number of expected barcode sites were predicted.\nExpected '+str(num_expected_barcode_sites)+", found: "+str(len(barcode_sites)))
		print("\namplicon sequence:\n\t"+expected_amplicon_seq)
		temp_barcode_site_list = list(barcode_sites.keys())
		temp_barcode_site_list = sorted(temp_barcode_site_list)
		for site_num in range(0,len(temp_barcode_site_list)):
			site = temp_barcode_site_list[site_num]
			print("\t"+str(site)+"\t"+barcode_sites[site][0]+" or "+barcode_sites[site][1])
		print('\nThe predicted barcode sites will be written to the file: "amplicon_info.auto.txt"\nEither try re-running the script to re-predict barcode content using a different random sample, or edit this file manually, then rename it to "amplicon_info.txt" before re-running this script to proceed.')
		outfile = open(project_dir+"amplicon_info.auto.txt","w")
	outfile.write("amplicon_seq\t"+expected_amplicon_seq+"\n")
	temp_barcode_site_list = list(barcode_sites.keys())
	temp_barcode_site_list = sorted(temp_barcode_site_list)
	for site_num in range(0,len(temp_barcode_site_list)):
		site = temp_barcode_site_list[site_num]
		outfile.write(str(site)+"\t"+barcode_sites[site][0]+","+barcode_sites[site][1]+"\n")
	outfile.close()

	if expected_amplicon_seq == "" or barcode_sites == {} or continue_running == False:
		clean_files_for_one_sample(accession,output_dir,trim_dir)
		sys.exit("\n\nExiting.")

### Trim read for all samples to summarize substitution frequency given quality
lane_prop_sub_qual_filename = output_dir+"sub_qual_info.prop.txt"
lane_count_sub_qual_filename = output_dir+"sub_qual_info.count.txt"
lane_prop_fit_qual_filename = output_dir+"sub_qual_fit.prop.txt"
if os.path.isfile(lane_prop_sub_qual_filename) == False:
	if parallel_max_cpu == 0:
		pre_trim_numcore = int(multiprocessing.cpu_count()/5)
	elif parallel_max_cpu > 0:
		pre_trim_numcore = min(parallel_max_cpu,int(multiprocessing.cpu_count()/5))
	elif parallel_max_cpu < 0:
		pre_trim_numcore = int(multiprocessing.cpu_count()/5)-parallel_max_cpu
	processed_list = []
	processed_list = Parallel(n_jobs=pre_trim_numcore)(delayed(reads_to_barcodes_one_sample)(accession,command_prefix,forward_read_suffix,reverse_read_suffix,input_read_dir,trim_dir,temp_dir,output_dir,expected_amplicon_seq,barcode_sites,min_Q_score,overwrite_existing_files,raw_read_input_type,reads_already_screened_for_quality,detect_barcode_content_first,F_edge_ignore,R_edge_ignore,force_trim_pre_merge,truncate_length_post_merge,True) for accession in accession_list)
	full_qual_sub_dict = {}
	qual_vals = {}
	for output_tup in processed_list:
		sampleID = output_tup[0]
		local_qual_sub_dict = output_tup[1]
		# sub_count_dict = output_tup[2]
		for backbone_nt in local_qual_sub_dict:
			for sub_nt in local_qual_sub_dict[backbone_nt]:
				for qual in local_qual_sub_dict[backbone_nt][sub_nt]:
					try:
						qual_vals[qual]
					except:
						qual_vals[qual] = ''
					try:
						full_qual_sub_dict[backbone_nt][sub_nt][qual] +=  local_qual_sub_dict[backbone_nt][sub_nt][qual]
					except:
						try:
							full_qual_sub_dict[backbone_nt][sub_nt][qual] =  local_qual_sub_dict[backbone_nt][sub_nt][qual]
						except:
							try:
								full_qual_sub_dict[backbone_nt][sub_nt] = {}
								full_qual_sub_dict[backbone_nt][sub_nt][qual] = local_qual_sub_dict[backbone_nt][sub_nt][qual]
							except:
								full_qual_sub_dict[backbone_nt] = {}
								full_qual_sub_dict[backbone_nt][sub_nt] = {}
								full_qual_sub_dict[backbone_nt][sub_nt][qual] = local_qual_sub_dict[backbone_nt][sub_nt][qual]
	qual_bins = sorted(list(qual_vals.keys()))
	count_output_string = '\t'
	prop_output_string = '\t'
	for q in range(0,len(qual_bins)):
		qual = qual_bins[q]
		count_output_string += '\t'+str(qual)
		prop_output_string += '\t'+str(qual)
	count_output_string += '\n'
	prop_output_string += '\n'
	bases = ['A','T','C','G']
	for b in range(0,len(bases)):
		full_count_total_dict = {}
		backbone_nt = bases[b]
		for n in range(0,len(bases)):
			sub_nt = bases[n]
			count_output_string += backbone_nt+'\t'+sub_nt
			for q in range(0,len(qual_bins)):
				qual = qual_bins[q]
				try:
					f_count = full_qual_sub_dict[backbone_nt][sub_nt][qual]
				except:
					f_count = 0
				count_output_string += '\t'+str(f_count)
				try:
					full_count_total_dict[qual] += f_count
				except:
					full_count_total_dict[qual] = f_count
			count_output_string += '\n'
		for n in range(0,len(bases)):
			sub_nt = bases[n]
			prop_output_string += backbone_nt+'\t'+sub_nt
			for q in range(0,len(qual_bins)):
				qual = qual_bins[q]
				try:
					f_count = full_qual_sub_dict[backbone_nt][sub_nt][qual]
				except:
					f_count = 0
				f_total = float(full_count_total_dict[qual])
				if f_total >0.0:
					f_prop = float(f_count)/f_total
				else:
					f_prop = 0.0
				prop_output_string += '\t'+str(f_prop)
			prop_output_string += '\n'

	outfile = open(lane_count_sub_qual_filename,"w")
	outfile.write(count_output_string)
	outfile.close()
	outfile = open(lane_prop_sub_qual_filename,"w")
	outfile.write(prop_output_string)
	outfile.close()
if os.path.isfile(lane_prop_fit_qual_filename) == False:
	sub_prob_dict = fit_qual_model(lane_prop_sub_qual_filename,lane_count_sub_qual_filename,lane_prop_fit_qual_filename,min_Q_score)
else:
	sub_prob_dict = load_fit_qual(lane_prop_fit_qual_filename)


### Fully process all samples
processed_list = []
if parallel_process == True:
	processed_list = Parallel(n_jobs=num_cores)(delayed(reads_to_barcodes_one_sample)(accession,command_prefix,forward_read_suffix,reverse_read_suffix,input_read_dir,trim_dir,temp_dir,output_dir,expected_amplicon_seq,barcode_sites,min_Q_score,overwrite_existing_files,raw_read_input_type,reads_already_screened_for_quality,detect_barcode_content_first,F_edge_ignore,R_edge_ignore,force_trim_pre_merge,truncate_length_post_merge) for accession in accession_list)
elif parallel_process == False:
	for accession in accession_list:
		barcode_tup = reads_to_barcodes_one_sample(accession,command_prefix,forward_read_suffix,reverse_read_suffix,input_read_dir,trim_dir,temp_dir,output_dir,expected_amplicon_seq,barcode_sites,min_Q_score,overwrite_existing_files,raw_read_input_type,reads_already_screened_for_quality,detect_barcode_content_first,F_edge_ignore,R_edge_ignore,force_trim_pre_merge,truncate_length_post_merge)
		processed_list.append(barcode_tup)


##### create list of all possible barcodes
bases = ['A','T','C','G']
temp_full_barcode_list = find_all_possible_barcode_combinations(barcode_sites)
full_barcode_list = []

barcodes_to_exclude = {}
# barcodes_to_exclude = {'TCTATATGCG':'','TGGACTGCGG':''}

for b in range(0,len(temp_full_barcode_list)):
	barcodeID = temp_full_barcode_list[b]
	skip = False
	try:
		barcodes_to_exclude[barcodeID]
		skip = True
	except:
		skip = False
	if skip == False:
		full_barcode_list.append(barcodeID)

##### merge all barcode info once all samples finish running
full_barcode_count_dict = {}
full_mismatch_dict = {}
full_nonbarcode_dict = {}
full_read_count = {}
nosingleton_read_count = {}

# full_nonbarcode_seq_list = []
for accession_output in processed_list:
	accession = accession_output[0]
	full_barcode_count_dict[accession] = accession_output[1]
	accession_barcode_list = accession_output[2]
	full_mismatch_dict[accession] = accession_output[3]

	accession_barcode_count = 0
	accession_barcode_count_nosingle = 0
	for barcodeID in full_barcode_count_dict[accession]:
		skip = False
		try:
			barcodes_to_exclude[barcodeID]
			skip = True
		except:
			skip = False
		if skip == False:
			local_count = full_barcode_count_dict[accession][barcodeID]
			if local_count > 0:
				accession_barcode_count += local_count
				if local_count > 1:
					accession_barcode_count_nosingle += local_count
	full_read_count[accession] = accession_barcode_count
	nosingleton_read_count[accession] = accession_barcode_count_nosingle
accession_list_passcount = []
for a in range(0,len(accession_list)):
	accession = accession_list[a]
	try:
		local_count = full_read_count[accession]
	except:
		local_count = 0
	if local_count>=min_barcodes_per_sample:
		accession_list_passcount.append(accession)
accession_list = accession_list_passcount
del accession_list_passcount

### Write final summary output tables
print("Writing barcode output summary tables")
accession_list = sorted(accession_list)
barcode_count_dict = {}
barcode_freq_dict = {}
barcode_count_dict_nosingleton = {}
barcode_freq_dict_nosingleton = {}
table_outfile = open(output_dir+"all_barcode_count.table.txt","w")
freq_table_outfile = open(output_dir+"all_barcode_freq.table.txt","w")
for a in range(0,len(accession_list)):
	accession = accession_list[a]
	table_outfile.write("\t"+accession)
	freq_table_outfile.write("\t"+accession)
table_outfile.write("\n")
freq_table_outfile.write("\n")
for b in range(0,len(full_barcode_list)):
	barcodeID = full_barcode_list[b]
	table_outfile.write(barcodeID)
	freq_table_outfile.write(barcodeID)
	for a in range(0,len(accession_list)):
		accession = accession_list[a]
		read_count = full_read_count[accession]
		read_count_nosingleton = nosingleton_read_count[accession]
		try:
			barcode_count = full_barcode_count_dict[accession][barcodeID]
			barcode_freq = (float(barcode_count)/float(read_count))
			barcode_freq_nosingleton = (float(barcode_count)/float(read_count_nosingleton))
		except:
			barcode_count = 0
			barcode_freq = 0
			barcode_freq_nosingleton = 0
		table_outfile.write("\t"+str(barcode_count))
		freq_table_outfile.write("\t"+str(barcode_freq))
		try:
			barcode_freq_dict[accession].append(barcode_freq)
			barcode_count_dict[accession].append(barcode_count)
		except:
			barcode_freq_dict[accession] = [barcode_freq]
			barcode_count_dict[accession] = [barcode_count]
		if read_count > 1:
			try:
				barcode_freq_dict_nosingleton[accession].append(barcode_freq_nosingleton)
				barcode_count_dict_nosingleton[accession].append(barcode_count)
			except:
				barcode_freq_dict_nosingleton[accession] = [barcode_freq_nosingleton]
				barcode_count_dict_nosingleton[accession] = [barcode_count]
	table_outfile.write("\n")
	freq_table_outfile.write("\n")
table_outfile.close()
freq_table_outfile.close()

max_mismatch_dict = {}
mismatch_table_outfile = open(output_dir+"all_mismatch_info.table.txt","w")
for a in range(0,len(accession_list)):
	accession = accession_list[a]
	mismatch_table_outfile.write("\t"+accession)
mismatch_table_outfile.write("\n")
for site in range(0,len(expected_amplicon_seq)):
	mismatch_table_outfile.write(str(site))
	for a in range(0,len(accession_list)):
		max_mismatch = 0
		accession = accession_list[a]
		try:
			mismatch_prop = round(float(full_mismatch_dict[accession][site]),5)
		except:
			mismatch_prop = 0.0
		mismatch_table_outfile.write("\t"+str(mismatch_prop))
		if mismatch_prop > max_mismatch:
			max_mismatch = mismatch_prop
		max_mismatch_dict[accession] = max_mismatch
	mismatch_table_outfile.write("\n")
mismatch_table_outfile.close()


print("Calculating alpha-diversity indicies")
all_process_inputs = []
for i in range(0,len(accession_list)):
	sample1 = accession_list[i]
	count1 = barcode_count_dict[sample1]
	process_tup = (sample1,count1)#,full_barcode_list)
	all_process_inputs.append(process_tup)
processed_list = Parallel(n_jobs=num_cores)(delayed(subsample_alpha_func)(process_tup,full_barcode_list,num_subsamples) for process_tup in all_process_inputs)

sub_size_list = []
alpha_subsample_infodict = {}
alpha_fit_dict = {}
for tup in processed_list:
	s1 = tup[0]
	alpha_subsample_infodict[s1] = tup[1]
	alpha_fit_dict[s1] = tup[2]
	local_sub_size_list = list(tup[1].keys())
	if len(local_sub_size_list) > len(sub_size_list):
		sub_size_list = sorted(local_sub_size_list)

dist_function_name_list = ['shannon','simpson','even','rich','chao']
for d in range(0,len(dist_function_name_list)):
	dist_function_name = dist_function_name_list[d]
	string_out = ''
	for i in range(0,len(sub_size_list)):
		subsample_size = sub_size_list[i]
		string_out += '\t'+str(subsample_size)
	string_out += '\n'
	for p in range(0,len(accession_list)):
		sampleID = accession_list[p]
		string_out += sampleID
		for i in range(0,len(sub_size_list)):
			subsample_size = sub_size_list[i]
			try:
				local_dist_list = alpha_subsample_infodict[sampleID][subsample_size][dist_function_name]
				m = np.median(local_dist_list)
			except:
				m = 'nan'
			string_out += '\t'+str(m)
		string_out += '\n'
	outfile = open(output_dir+"subsample_regress."+dist_function_name+".n"+str(num_subsamples)+".txt",'w')
	outfile.write(string_out)
	outfile.close()

alpha_div_outfile = open(output_dir+"all_barcode.alpha_diversity.txt","w")
alpha_div_outfile.write("sampleID\tshannon_div\tsimpson_div\tchao_richness\tshannon_evenness\tobs_richness\n")
for a in range(0,len(accession_list)):
	sampleID = accession_list[a]
	try:
		sample_fit_dict = alpha_fit_dict[sampleID]
	except:
		sample_fit_dict = {}
	if sample_fit_dict != {}:
		H_shannon = sample_fit_dict['shannon']
		H_simpson = sample_fit_dict['simpson']
		E_shannon = sample_fit_dict['even']
		R_chao = sample_fit_dict['chao']
		R_rich = sample_fit_dict['rich']
	else:
		H_shannon,E_shannon = shannon_alpha(barcode_freq_dict_nosingleton[sampleID])
		H_simpson = simpson_alpha(barcode_freq_dict_nosingleton[sampleID])
		R_chao = chao_1_richness(barcode_count_dict_nosingleton[sampleID],2) #value of 2 to not count singletons
		R_rich = true_richness(barcode_count_dict_nosingleton[sampleID],2) #value of 2 to not count singletons
	alpha_div_outfile.write(sampleID+"\t"+str(H_shannon)+"\t"+str(H_simpson)+"\t"+str(R_chao)+"\t"+str(E_shannon)+"\t"+str(R_rich)+"\n")
alpha_div_outfile.close()


print("Calculating beta-diversity dissimilarity indicies")
all_process_inputs = []
for i in range(0,len(accession_list)):
	sample1 = accession_list[i]
	for j in range(i+1,len(accession_list)):
		sample2 = accession_list[j]
		count1 = barcode_count_dict[sample1]
		count2 = barcode_count_dict[sample2]
		process_tup = ((sample1,count1),(sample2,count2))
		all_process_inputs.append(process_tup)
processed_list = Parallel(n_jobs=num_cores)(delayed(subsample_beta_func)(process_tup,full_barcode_list,num_subsamples) for process_tup in all_process_inputs)

sub_size_list = []
beta_subsample_infodict = {}
beta_fit_dict = {}
for tup in processed_list:
	s1 = tup[0]
	s2 = tup[1]
	beta_subsample_infodict[s1+s2] = tup[2]
	beta_fit_dict[s1+s2] = tup[3]
	local_sub_size_list = list(tup[2].keys())
	if len(local_sub_size_list) > len(sub_size_list):
		sub_size_list = sorted(local_sub_size_list)

dist_function_name_list = ['bray','jaccard']
for d in range(0,len(dist_function_name_list)):
	dist_function_name = dist_function_name_list[d]
	string_out = '\t'
	for i in range(0,len(sub_size_list)):
		subsample_size = sub_size_list[i]
		string_out += '\t'+str(subsample_size)
	string_out += '\n'
	for p in range(0,len(accession_list)):
		sample1 = accession_list[p]
		for q in range(p,len(accession_list)):
			sample2 = accession_list[q]
			string_out += sample1+'\t'+sample2
			for i in range(0,len(sub_size_list)):
				subsample_size = sub_size_list[i]
				try:
					local_dist_list = beta_subsample_infodict[sample1+sample2][subsample_size][dist_function_name]
				except:
					try:
						local_dist_list = beta_subsample_infodict[sample2+sample1][subsample_size][dist_function_name]
					except:
						local_dist_list = []
				if len(local_dist_list)>0:
					m = np.median(local_dist_list)
				else:
					m = 'nan'
				string_out += '\t'+str(m)
			string_out += '\n'
	outfile = open(output_dir+"subsample_regress."+dist_function_name+".n"+str(num_subsamples)+".txt",'w')
	outfile.write(string_out)
	outfile.close()



	bray_outfile = open(output_dir+"all_barcode.bray-curtis_dissimilarity.txt","w")
	for j in range(0,len(accession_list)):
		bray_outfile.write("\t"+accession_list[j])
	bray_outfile.write("\n")
	for i in range(0,len(accession_list)):
		bray_outfile.write(accession_list[i])
		for j in range(0,len(accession_list)):
			accession1 = accession_list[i]
			accession2 = accession_list[j]
			try:
				dist = beta_fit_dict[accession1+accession2]['bray']
			except:
				try:
					dist = beta_fit_dict[accession2+accession1]['bray']
				except:
					try:
						dist = bray_curtis_dissimilarity(barcode_freq_dict_nosingleton[accession1],barcode_freq_dict_nosingleton[accession2])
					except:
						dist = 'nan'
			bray_outfile.write("\t"+str(dist))
		bray_outfile.write("\n")
	bray_outfile.close()

	jaccard_outfile = open(output_dir+"all_barcode.jaccard_dissimilarity.txt","w")
	for j in range(0,len(accession_list)):
		jaccard_outfile.write("\t"+accession_list[j])
	jaccard_outfile.write("\n")
	for i in range(0,len(accession_list)):
		jaccard_outfile.write(accession_list[i])
		for j in range(0,len(accession_list)):
			accession1 = accession_list[i]
			accession2 = accession_list[j]
			try:
				dist = beta_fit_dict[accession1+accession2]['jaccard']
			except:
				try:
					dist = beta_fit_dict[accession2+accession1]['jaccard']
				except:
					try:
						dist = jaccard_dissimilarity(barcode_count_dict_nosingleton[accession1],barcode_count_dict_nosingleton[accession2],2) #value of 2 to exclude singletons
					except:
						dist = 'nan'
			jaccard_outfile.write("\t"+str(dist))
		jaccard_outfile.write("\n")
	jaccard_outfile.close()


####################################### Plot stacked bar plots #######################################
if generate_stacked_bar_plots == True:
	print("Creating stacked bar plots")
	
	#make color palette
	colorblind_friendly_palette = ['#E8ECFB', '#D9CCE3', '#D1BBD7', '#CAACCB', '#BA8DB4', '#AE76A3', '#AA6F9E', '#994F88', '#882E72', '#1965B0', '#437DBF', '#5289C7', '#6195CF', '#7BAFDE', '#4EB265', '#90C987', '#CAE0AB', '#F7F056', '#F7CB45', '#F6C141', '#F4A736', '#F1932D', '#EE8026', '#E8601C', '#E65518', '#DC050C', '#A5170E', '#72190E', '#42150A']
	full_barcode_color_list = []
	full_barcode_color_dict = {}
	for num in range(0,len(full_barcode_list)):
		barcode = full_barcode_list[num]
		color_index = num % len(colorblind_friendly_palette)
		color_hex_value = colorblind_friendly_palette[color_index]
		full_barcode_color_list.append(color_hex_value)
		full_barcode_color_dict[barcode] = color_hex_value

	#read in sample titers, if present
	titer_dict = {}
	if os.path.isfile(titer_filename) == True:
		titer_list = []
		titer_infile = open(titer_filename,"r")
		for line in titer_infile:
			line = line.strip()
			if len(line)>0:
				if line[0] !='#':
					line = line.split("\t")
					if is_number(line[1]) == True:
						accession = line[0]
						titer = float(line[1])
						titer_dict[accession] = titer
						titer_list.append(titer)
		titer_infile.close()
		titer_list = np.array(titer_list)
		max_obs_titer = np.max(titer_list)
		min_obs_titer = np.min(np.where(titer_list==0,999,titer_list))
		if max_obs_titer/min_obs_titer > 100:
			print("Sample do not appear to be log transformed. Log10 transforming values now.")
			log_titer_dict = {}
			for accession in titer_dict:
				log_titer_dict[accession] = -1*np.log10(titer_dict[accession])
			titer_dict = log_titer_dict
			del log_titer_dict
			del titer_list

	#make the plot
	bars_per_row = len(accession_list)
	width_per_bar = 0.4
	num_col = 1
	num_row = 1#math.ceil(len(accession_list)/bars_per_row)


	fig, axs = plt.subplots(num_row, num_col,sharex=True, sharey=True)##num_row, 
	fig.set_size_inches(num_col*(width_per_bar*bars_per_row), num_row*7)


	plot_freq_array = []
	max_mismatch_array = []
	titer_not_found = []
	df_key = ['Sample']
	df_key.extend(full_barcode_list)
	max_titer = 0
	for a in range(0,len(accession_list)):
		accession = accession_list[a]
		try:
			max_mismatch = max_mismatch_dict[accession]
		except:
			max_mismatch = 0
		max_mismatch_array.append(max_mismatch)
		try:
			accession_titer = titer_dict[accession]
		except:
			accession_titer = 1.0
			titer_not_found.append(accession)
		if accession_titer > max_titer:
			max_titer = accession_titer
		barcode_freq_array = barcode_freq_dict[accession]
		### Make scaled bar plots with barcodes scaled to sample titers
		titer_adjusted_freq_list = [accession]
		for i in range(0,len(full_barcode_list)):
			try:
				titer_adjusted_freq = barcode_freq_array[i]*accession_titer
			except:
				titer_adjusted_freq = 0.0
			titer_adjusted_freq_list.append(titer_adjusted_freq)
		plot_freq_array.append(titer_adjusted_freq_list)

	if max_titer == 1.0:
		y_axis_label = 'Proportion of viral population'
		print("\tStacked bar plots will be scaled to 1.0 - add file with viral titers to scale plots to viral titer instead.")
	else:
		y_axis_label = "Viral titer"
		if len(titer_not_found)>0:
			print("\nUnable to pair the following samples with a titer:")
			for t in range(0,len(titer_not_found)):
				print('\t'+titer_not_found[t])
	
	barcode_df = pd.DataFrame(plot_freq_array,columns=df_key)
	print('\t...plot rendering starting - this can take a long time to complete with a lot of samples - please do not interrupt')
	bar_plot = barcode_df.plot(x='Sample', kind='bar', stacked=True,legend=False,color=full_barcode_color_list,ax=axs,width=0.9, ylim=(0, max_titer*1.05),ylabel=y_axis_label)
	if max_titer == 1.0:
		axs.set_yticks([0.0,0.25,0.5,0.75,1.0])
	mismatch_axis = axs.twinx()
	mismatch_axis.plot(max_mismatch_array,color='black',marker='o', markersize=3, mfc='grey',linestyle='None')#linewidth=0.5, 
	mismatch_axis.set_ylim(0.0,0.55)
	mismatch_axis.set_yticks([0.0,0.1,0.2,0.3,0.4,0.5])
	mismatch_axis.set_ylabel("Max proportion of reads discarded")
	plt.tight_layout()
	if max_titer == 1.0:
		plt.savefig(output_dir+"all_barcode.stacked_bars.pdf")
	else:
		plt.savefig(output_dir+"all_barcode.stacked_bars.log10_freq.pdf")
	# plt.show()


print("Processing Complete")
