'''
barcoded amplicon processing pipeline - v1.0
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
python barcodeID_v1.0.py

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
'''
###################################### Load required modules ######################################
import os
import sys
import numpy as np
import time

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
min_qual_backbone = 25
min_qual_barcode_sites = 30
amplicon_info_infile_name = 'amplicon_info.txt'

F_edge_ignore = 0 #these values are used to exclude the edge of reads if their quality is low and needlessly discarding excessive data due to sequencing quality issues
R_edge_ignore = 0 #these values are used to exclude the edge of reads if their quality is low and needlessly discarding excessive data due to sequencing quality issues
force_trim_pre_merge = 0 #use this value to force trimming the first x nucleotides from input reads, which is helpful when the beginning of reads are low quality
truncate_length_post_merge = 0 #use this value to truncate reads post-merging if amplicons are not always equal lengths - use only as a last resort to make processing complete, or for troubleshooting when no barcodes are detected in some samples

overwrite_existing_files = False #when False, pipeline will avoid re-running samples that have already been processed 

reads_already_screened_for_quality = False # 'True' or 'False' --- If reads were already merged and screened for average quality and length (such as during a demultiplexing step), the script will not attempt to filter or merge the input reads
raw_read_input_type = "paired" # 'paired' or 'single' --- When 'single' the pipeline will look for barcodes in one file, when 'paired' the pipeline will look for forward and reverse reads. Single core use has not been extensively tested - recommend using "paired" and setting "parallel_max_cpu" to 1 if you only wish to use one core at a time
forward_read_suffix = "" #if empty, the script will attempt to predict these values, but it only works under specific assumptions
reverse_read_suffix = "" #this value is ignored if 'reads_already_screened_for_quality' is 'True' or 'raw_read_input_type' is 'single' - NOTE: this variable cannot be removed without causing missing value errors

parallel_process = True #'True' or 'False' --- when True, the script uses the joblib package to run multiple samples simultaneously. Currently not stable on windows
parallel_max_cpu = -1 # when parallel process true, this sets the max CPUs the script is allowed to use. If you want to use all available, set to 0. If you want to use all but one, set to -1 (or -2 to use all but 2 CPUs)

global verbose
verbose = True #when this value is 'True', the script will print status updates

############################################ FUNCTIONS ############################################
def run_command(command):
	run_status = False
	# return_code = os.system(command)  #unhash this line and hash out the line below if you wish to see all updated output from command line applications for troubleshooting
	return_code = os.system(command+" >/dev/null 2>&1")
	if return_code == 0: #if return_code is something other than zero, python was not able to successfully run your analysis
		run_status = True
	return run_status


def unzip_gz(filename,path_to_file):
	command = "gzip -d "+path_to_file+filename
	out = run_command(command)
	if out == False:
		sys.exit()
	return filename

def unzip_tar_gz(filename,path_to_file):
	command = "tar -xvf "+path_to_file+filename
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

def zip_tar_gz(filename,path_to_file,zipped_file_extension):
	command = "tar -cvzf "+path_to_file+filename.split(zipped_file_extension)[0]+".tar.gz "+path_to_file+filename
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
	gz_inputs = [f for f in os.listdir(input_dir) if f.endswith(".gz")]
	if len(gz_inputs) >0:  #unzip .gz files
		if parallel_process == True:
			processed_list = Parallel(n_jobs=num_cores)(delayed(unzip_gz)(i,input_dir) for i in gz_inputs)
		else:
			for gzfile in gz_inputs:
				unzip_gz(gzfile,input_dir)
	tar_gz_inputs = [f for f in os.listdir(input_dir) if f.endswith(".tar.gz")]
	if len(tar_gz_inputs) > 0:  #unzip .tar.gz files
		if parallel_process == True:
			processed_list = Parallel(n_jobs=num_cores)(delayed(unzip_tar_gz)(i,input_dir) for i in tar_gz_inputs)
		else:
			for targzfile in tar_gz_inputs:
				unzip_tar_gz(gzfile,input_dir)
	fq_inputs = [f for f in os.listdir(input_dir) if f.endswith(".fastq") or f.endswith(".fq")]
	if len(fq_inputs) >0:
		files_found = True
	return files_found


def predict_paired_file_extensions(input_dir):
	exten_filelist = [f for f in os.listdir(input_dir) if f.endswith(".fastq") or f.endswith(".fq")]
	char_dict = {}
	for i in range(-1,(-1*len(exten_filelist[0])),-1):
		char_dict[i] = {}
		for f in range(0,len(exten_filelist)):
			filename = exten_filelist[f]
			try:
				char = filename[i]
			except:
				char = "na"
			try:
				char_dict[i][char] += 1
			except:
				char_dict[i][char] = 1
	first_dichotomous = 0
	dichotomous_sites = []
	last_monomorphic = 0
	searching_for_last_monomorphic = False
	for i in range(-1,(-1*len(exten_filelist[0])),-1):
		if len(char_dict[i]) == 2:
			dichotomous_sites.append(i)
			if first_dichotomous == 0:
				first_dichotomous = i
			searching_for_last_monomorphic = True
		elif len(char_dict[i]) >1 and searching_for_last_monomorphic == True:
			searching_for_last_monomorphic = False
		elif len(char_dict[i]) == 1 and searching_for_last_monomorphic == True:
			last_monomorphic = i
	exten_count = {}
	for f in range(0,len(exten_filelist)):
		filename = exten_filelist[f]
		try:
			exten_count[filename[last_monomorphic:]] += 1
		except:
			exten_count[filename[last_monomorphic:]] = 1
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
	char_dict = {}
	for i in range(-1,(-1*len(exten_filelist[0])),-1):
		char_dict[i] = {}
		for f in range(0,len(exten_filelist)):
			filename = exten_filelist[f]
			try:
				char = filename[i]
			except:
				char = "na"
			try:
				char_dict[i][char] += 1
			except:
				char_dict[i][char] = 1
	last_monomorphic = 0
	searching_for_last_monomorphic = True
	for i in range(-1,(-1*len(exten_filelist[0])),-1):
		if len(char_dict[i]) == 1 and searching_for_last_monomorphic == True:
			last_monomorphic = i
		elif len(char_dict[i]) > 1 and searching_for_last_monomorphic == True:
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
	for col in range(0,amplicon_length):
		nt_count = {'A':0,'T':0,'C':0,'G':0}
		for seq in seq_dict:
			if len(seq) == amplicon_length:
				try:
					cur_nt = seq[col]
					nt_count[cur_nt] += 1
				except:
					if seq[col] != "N":
						skip_to_print_bottom()
						sys.exit('Unexpected character found: "'+cur_nt+'" at site '+str(col))
		nt_list = []
		for nt in nt_count:
			temp_tup = (nt_count[nt],nt)
			nt_list.append(temp_tup)
		nt_list = sorted(nt_list, reverse=True)
		if float(nt_list[0][0])/float(len(seq_dict)) >= 0.7:
			amplicon_expect += nt_list[0][1]
		else:
			amplicon_expect += "N"
			barcode_sites[col] = [nt_list[0][1],nt_list[1][1]]
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


def read_in_barcodes(barcode_filename_in,mismatch_filename_in,nonbarcode_filename_in):
	global verbose
	if verbose == True:
		update_print_line("Reading in barcodes",accession)
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

	high_qual_nonbarcode_dict = {}
	infile = open(nonbarcode_filename_in)
	for line in infile:
		line = line.strip().split("\t")
		if line[0] != "non-barcode_seq":
			seq = line[0]
			count = int(line[1])
			high_qual_nonbarcode_dict[seq] = count
	infile.close()
	if verbose == True:
		update_print_line("Finished reading barcodes",accession)
	return barcode_count_dict, barcode_list,mismatch_by_site,high_qual_nonbarcode_dict


def shannon_alpha(freq_array):
	H = 0
	S_obs = 0
	for freq in freq_array:
		if freq >0.0:
			H_i = -1*freq*np.log(freq)
			H += H_i
			S_obs += 1
	H_max = np.log(S_obs)
	E = H/H_max
	H = round(H,3)
	E = float(round_to_n_sig_figs(E,3))
	return H,E

def simpson_alpha(freq_array):
	p_sum = 0
	S = 0
	for freq in freq_array:
		if freq >0.0:
			S += 1
			p_sum += freq**2
	H = p_sum
	return round_to_n_sig_figs(H,3)

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
	return round_to_n_sig_figs(B,3)

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
	return round_to_n_sig_figs(B,3)

def jaccard_dissimilarity(array1,array2,min_freq=0.0):
	AB,A,B,C = 0,0,0,0
	for i in range(0,len(array1)):
		val1 = array1[i]
		val2 = array2[i]
		if val1 > min_freq or val2 > min_freq:
			C += 1
		if val1 > min_freq and val2 > min_freq:
			AB += 1
	J = 1 - (AB/C)
	return round_to_n_sig_figs(J,3)

def chao_1_richness(count_array):
	S_obs = 0
	S_single = 0
	S_double = 0
	for count in count_array:
		if count >0:
			S_obs += 1
			if count == 1:
				S_single += 1
			elif count == 2:
				S_double += 1
	S_unobs = (S_single*(S_single-1))/(2*S_double+1)
	R = S_obs + S_unobs
	return round(R,1)

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
		num_sig_figs = int(num_sig_figs)
		if num_sig_figs == 0:
			num_sig_figs = 1
		sci_val = "{:.10e}".format(val)
		split_sci_val = sci_val.split("e")
		if len(split_sci_val) == 2:
			rounded_base_number = round(float(split_sci_val[0]),num_sig_figs-1)
			exponent = int(split_sci_val[1])
			if exponent == 0:
				val_out = str(rounded_base_number) + ((num_sig_figs)-1)*'0'
			elif exponent < 0:
				exponent*=-1
				val_out = '0.' + (exponent-1)*'0' + str(rounded_base_number).replace(".","")
				val_out = str(float(val_out))
			elif exponent > 0:
				val_out = str(rounded_base_number) +'e'+ (str(exponent))
			return val_out
		else:
			return val
	else:
		return val

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


def screen_barcode_site_qual(seq_in,qual_in,seq_expect,barcode_dict,min_barcode_site_qual,min_other_site_qual,F_edge_ignore,R_edge_ignore):
	barcode_mismatches = 0
	backbone_mismatches = 0
	high_qual_mismatches = 0
	mismatch_sites = []
	expected_barcode_nt_list = []
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
			if site_qual < min_other_site_qual:
				backbone_mismatches += 1
				mismatch_sites.append(site)
			elif site_nt != seq_expect[site]:
				backbone_mismatches += 1
				mismatch_sites.append(site)
				if site_qual >= min_qual_backbone:
					high_qual_mismatches += 1
		elif barcode_site == True:
			if site_nt == expected_barcode_nt_list[0] or site_nt == expected_barcode_nt_list[1]:
				if site_qual < min_barcode_site_qual:
					barcode_mismatches += 1
					mismatch_sites.append(site)
			else:
				barcode_mismatches += 1
				mismatch_sites.append(site)
				if site_qual >= min_barcode_site_qual:
					high_qual_mismatches += 1

	total_mismatches = barcode_mismatches+backbone_mismatches
	mismatch_sites = list(set(mismatch_sites))
	return total_mismatches,mismatch_sites,high_qual_mismatches


def raw_read_processing_single(accession,command_prefix,raw_reads_forward,trim_dir,temp_dir,amplicon_length,force_trim_pre_merge,truncate_length_post_merge):
	global verbose
	if verbose == True:
		update_print_line("Processing input reads",accession)
	trim_fastq_filename = temp_dir+accession+".trim.fq"
	clean_filename_out = trim_dir+accession+".fq"
	command = command_prefix+'bbduk.sh in="'+raw_reads_forward+'" out="'+trim_fastq_filename+'" qin=33 maq=20 minlength='+str(amplicon_length)+' maxlength='+str(amplicon_length)+' overwrite=t'
	out = run_command(command)
	if out == False:
		sys.exit()
	time.sleep(5)
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
	command = command_prefix+'repair.sh in="'+raw_reads_forward+'" in2="'+raw_reads_reverse+'" out="'+repair_fastq_filename+'" fint=t repair=t overwrite=t -da -Xmx1000m'
	out = run_command(command)
	if out == False:
		sys.exit()
	time.sleep(1)

	len_diff_buffer = 50
	command = command_prefix+'bbmerge.sh in="'+repair_fastq_filename+'" outadapter="'+temp_adapter_filename+'" minlength='+str(read_len_expected-(len_diff_buffer+diff_median_len))+' maxlength='+str(read_len_expected+(len_diff_buffer+diff_median_len))+' overwrite=t' #forcetrimleft=9
	if force_trim_pre_merge >0:
		command += " forcetrimleft="+str(force_trim_pre_merge)#+" forcetrimright=130"
	out = run_command(command+" -Xmx2000m")

	if out == False:
		sys.exit()
	time.sleep(1)
	adapter_pass = check_adapter(temp_adapter_filename)
	if adapter_pass == True:
		command = command_prefix+'bbmerge.sh in="'+repair_fastq_filename+'" out="'+pretrim_fastq_filename+'" qin=33 mix=f adapters="'+temp_adapter_filename+'" minlength='+str(read_len_expected-(len_diff_buffer+diff_median_len))+' maxlength='+str(read_len_expected+(len_diff_buffer+diff_median_len)) #forcetrimleft=9 
		if force_trim_pre_merge >0:
			command += " forcetrimleft="+str(force_trim_pre_merge)#+" forcetrimright=130"
		out = run_command(command+" overwrite=t -Xmx2000m")
		if out == False:
			sys.exit()
		time.sleep(1)
	else:
		command = command_prefix+'bbmerge.sh in="'+repair_fastq_filename+'" out="'+pretrim_fastq_filename+'" qin=33 mix=f minlength='+str(read_len_expected-(len_diff_buffer+diff_median_len))+' maxlength='+str(read_len_expected+(len_diff_buffer+diff_median_len))+' overwrite=t' #forcetrimleft=9 
		if force_trim_pre_merge >0:
			command += " forcetrimleft="+str(force_trim_pre_merge)#+" forcetrimright=130"
		out = run_command(command+" overwrite=t -Xmx2000m")
		if out == False:
			sys.exit()
		time.sleep(1)

	if amplicon_length == 0:
		median_merged_len = calc_expect_read_len(pretrim_fastq_filename)
	else:
		median_merged_len = amplicon_length
	command = command_prefix+'bbduk.sh in="'+pretrim_fastq_filename+'" out="'+merged_trimmed_filename+'" qin=33 maq=20 minlength='+str(median_merged_len)+' maxlength='+str(median_merged_len)+' tpe=f tbo=f'
	out = run_command(command+ ' overwrite=t -Xmx1000m')
	if out == False:
		sys.exit()
	time.sleep(1)

	if truncate_length_post_merge > 0:
		truncate_merged_reads(merged_trimmed_filename,truncate_length_post_merge)
	
	os.remove(repair_fastq_filename)
	os.remove(temp_adapter_filename)
	os.remove(pretrim_fastq_filename)
	os.rename(merged_trimmed_filename,clean_filename_out)
	
	return clean_filename_out


def barcode_screen_trimmed_reads(accession,sequence_filename,expected_amplicon_seq,barcode_sites,min_qual_barcode_sites,min_qual_backbone,temp_dir,output_dir,F_edge_ignore,R_edge_ignore):
	### Load in processed reads and collect sequence to quality info to screen for barcodes
	global verbose
	if verbose == True:
		update_print_line("Screening reads for barcodes",accession)
	infile = open(sequence_filename,"r")
	mismatch_by_site = {}
	barcode_count_dict = {}
	high_qual_nonbarcode_count_dict = {}
	barcode_list = []
	line_counter = -1
	seq_count = 0
	for line in infile:
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
					if seq_count%100000==0 and seq_count > 1:
						update_print_line("Reads screened: "+str(seq_count),accession)
						# print('\tReads screened from: "'+accession+'":\t'+str(seq_count))
				mismatch_count,read_mismatch_locs,high_qual_mismatches = screen_barcode_site_qual(nucleotide_line,qual_line,expected_amplicon_seq,barcode_sites,min_qual_barcode_sites,min_qual_backbone,F_edge_ignore,R_edge_ignore)#(seq_in,qual_in,seq_expect,barcode_dict,min_barcode_site_qual,min_other_site_qual)
				if mismatch_count == 0:
					barcode_string = barcode_from_fullseq(barcode_sites,nucleotide_line)
					try:
						barcode_count_dict[barcode_string] += 1
					except:
						barcode_count_dict[barcode_string] = 1
						barcode_list.append(barcode_string)
				elif mismatch_count <= 3:
					for site_num in read_mismatch_locs:
						try:
							mismatch_by_site[site_num] += 1
						except:
							mismatch_by_site[site_num] = 1
					if high_qual_mismatches == mismatch_count:
						try:
							high_qual_nonbarcode_count_dict[nucleotide_line] += 1
						except:
							high_qual_nonbarcode_count_dict[nucleotide_line] = 1
			line_counter = -1
	barcode_list = sorted(list(set(barcode_list)))
	infile.close()

	#### summarize mismatch count
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

	#### count unique barcodes
	barcode_outfile = open(temp_dir+accession+".barcode_count.txt","w")
	barcode_outfile.write("barcode\tcount\n")
	for num in range(0,len(barcode_list)):
		barcode_string = barcode_list[num]
		barcode_count = barcode_count_dict[barcode_string]
		barcode_outfile.write(str(barcode_string)+"\t"+str(barcode_count)+"\n")
	barcode_outfile.close()

	#### count high-quality non-barcode sites
	nonbarcode_outfile = open(temp_dir+accession+".nonbarcode_count.txt","w")
	nonbarcode_outfile.write("non-barcode_seq\tcount\n")
	for nucl_string in high_qual_nonbarcode_count_dict:
		nonbarcode_count = high_qual_nonbarcode_count_dict[nucl_string]
		nonbarcode_outfile.write(str(nucl_string)+"\t"+str(nonbarcode_count)+"\n")
	nonbarcode_outfile.close()

	os.rename(temp_dir+accession+".mismatch_info.txt",output_dir+"samples/"+accession+".mismatch_info.txt")
	os.rename(temp_dir+accession+".barcode_count.txt",output_dir+"samples/"+accession+".barcode_count.txt")
	os.rename(temp_dir+accession+".nonbarcode_count.txt",output_dir+"samples/"+accession+".nonbarcode_count.txt")

	return barcode_count_dict,barcode_list,mismatch_prop_by_site,high_qual_nonbarcode_count_dict


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

def skip_to_print_bottom():
	return None


##########################      Core pipeline function      ##########################

def reads_to_barcodes_one_sample(accession,command_prefix,forward_read_suffix,reverse_read_suffix,input_read_dir,trim_dir,temp_dir,output_dir,expected_amplicon_seq,barcode_sites,min_qual_barcode_sites,min_qual_backbone,overwrite_existing_files,raw_read_input_type,reads_already_screened_for_quality,detect_barcode_content_first,F_edge_ignore,R_edge_ignore,force_trim_pre_merge,truncate_length_post_merge):
	
	update_print_line("Read processing start",accession)

	raw_reads_forward = input_read_dir+accession+forward_read_suffix
	raw_reads_reverse = input_read_dir+accession+reverse_read_suffix

	amplicon_length = len(expected_amplicon_seq)
	sample_barcode_info_filename = output_dir+"samples/"+accession+".barcode_count.txt"
	barcode_mismatch_info_filename = output_dir+"samples/"+accession+".mismatch_info.txt"
	sample_nonbarcode_info_filename = output_dir+"samples/"+accession+".nonbarcode_count.txt"
	processed_read_filename = trim_dir+accession+".fq"
	
	#clear the temporary files folder from previous interrupted runs
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
			elif raw_read_input_type == "single":
				processed_read_filename = raw_read_processing_single(accession,command_prefix,raw_reads_forward,raw_reads_reverse,trim_dir,temp_dir,detect_barcode_content_first,amplicon_length,force_trim_pre_merge,truncate_length_post_merge)
	elif reads_already_screened_for_quality == True:
		processed_read_filename = raw_reads_forward

	if os.path.isfile(sample_barcode_info_filename) == False:
		barcode_count_dict,barcode_list,mismatch_by_site,high_qual_nonbarcode_dict = barcode_screen_trimmed_reads(accession,processed_read_filename,expected_amplicon_seq,barcode_sites,min_qual_barcode_sites,min_qual_backbone,temp_dir,output_dir,F_edge_ignore,R_edge_ignore)
	else:
		barcode_count_dict,barcode_list,mismatch_by_site,high_qual_nonbarcode_dict = read_in_barcodes(sample_barcode_info_filename,barcode_mismatch_info_filename,sample_nonbarcode_info_filename)

	out_tup = (accession,barcode_count_dict,barcode_list,mismatch_by_site,high_qual_nonbarcode_dict)
	
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
				reverse_read_suffix = line[1]
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
accession_list = []
for files in file_list:
	accession = files.split(forward_read_suffix)[0]
	accession_list.append(accession)
accession_list = list(set(accession_list))

# accession_list = ['VR-022_dil1_NA1','VR-023_dil2_NA1','VR-024_dil3_NA1','VR-026_dil5_NA1','VR-027_dil6_NA1','VR-028_dil7_NA1','VR-029_dil1_HA1-NA1','VR-030_dil2_HA1-NA1','VR-031_dil3_HA1-NA1','VR-032_dil4_HA1-NA1','VR-033_dil5_HA1-NA1','VR-034_dil6_HA1-NA1','VR-035_dil7_HA1-NA1']
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
	print("\n\nCould not detect barcode info input file.\nA random sample will be selected to attempt to detect barcode sites and nucleotides automatically.\n")
if detect_barcode_content_first == True:
	'''	This section will attempt to predict the barcode info for your experiment if it is not already given.
	It will run through the entire pipeline with this one sample, then you need to check the output to determine
	if it correctly identified the barcode sites. If you agree with the prediction, remove ".auto" from
	the "barcode.info.auto.txt" file that was generated, then re-run the script to process all remaining samples.'''
	overwrite_existing_files = True
	accession_entered = input('Enter sample name for predicting barcode info (or "random" to pick random sample):\n')
	temp_cont = False
	if accession_entered == "random":
		accession_list = [accession_list[0]]	
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
	raw_reads_reverse = input_read_dir+accession+reverse_read_suffix
	processed_read_filename = raw_read_processing_paired(accession,command_prefix,raw_reads_forward,raw_reads_reverse,trim_dir,temp_dir,detect_barcode_content_first,0,force_trim_pre_merge,truncate_length_post_merge)
	temp_seq_dict,temp_read_length = subset_fastq(processed_read_filename,10000,50000)
	expected_amplicon_seq,barcode_sites = detect_barcode_content(temp_seq_dict,temp_read_length)
	continue_running = False
	if len(barcode_sites) == num_expected_barcode_sites:
		print('\n\nThe correct number of expected barcode sites were automatically detected.\n')
		continue_running = False
		print("amplicon sequence:\n\t"+expected_amplicon_seq)
		for site in barcode_sites:
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
		for site in barcode_sites:
			print("\t"+str(site)+"\t"+barcode_sites[site][0]+" or "+barcode_sites[site][1])
		print('\nThe predicted barcode sites will be written to the file: "amplicon_info.auto.txt"\nEither try re-running the script to re-predict barcode content using a different random sample, or edit this file manually, then rename it to "amplicon_info.txt" before re-running this script to proceed.')
		outfile = open(project_dir+"amplicon_info.auto.txt","w")
	outfile.write("amplicon_seq\t"+expected_amplicon_seq+"\n")
	for site in barcode_sites:
		outfile.write(str(site)+"\t"+barcode_sites[site][0]+","+barcode_sites[site][1]+"\n")
	outfile.close()

	if expected_amplicon_seq == "" or barcode_sites == {} or continue_running == False:
		skip_to_print_bottom()
		clean_files_for_one_sample(accession,output_dir,trim_dir)
		sys.exit("\n\nExiting.")

### Process all samples
processed_list = []
if parallel_process == True:
	processed_list = Parallel(n_jobs=num_cores, prefer="threads")(delayed(reads_to_barcodes_one_sample)(accession,command_prefix,forward_read_suffix,reverse_read_suffix,input_read_dir,trim_dir,temp_dir,output_dir,expected_amplicon_seq,barcode_sites,min_qual_barcode_sites,min_qual_backbone,overwrite_existing_files,raw_read_input_type,reads_already_screened_for_quality,detect_barcode_content_first,F_edge_ignore,R_edge_ignore,force_trim_pre_merge,truncate_length_post_merge) for accession in accession_list)
elif parallel_process == False:
	for accession in accession_list:
		barcode_tup = reads_to_barcodes_one_sample(accession,command_prefix,forward_read_suffix,reverse_read_suffix,input_read_dir,trim_dir,temp_dir,output_dir,expected_amplicon_seq,barcode_sites,min_qual_barcode_sites,min_qual_backbone,overwrite_existing_files,raw_read_input_type,reads_already_screened_for_quality,detect_barcode_content_first,F_edge_ignore,R_edge_ignore,force_trim_pre_merge,truncate_length_post_merge)
		processed_list.append(barcode_tup)


##### merge all barcode info once all samples finish running
full_barcode_count_dict = {}
full_barcode_list = []
full_mismatch_dict = {}
full_nonbarcode_dict = {}
full_read_count = {}
full_nonbarcode_seq_list = []
for accession_output in processed_list:
	accession = accession_output[0]
	full_barcode_count_dict[accession] = accession_output[1]
	accession_barcode_list = accession_output[2]
	full_mismatch_dict[accession] = accession_output[3]
	full_nonbarcode_dict[accession] = accession_output[4]

	full_barcode_list += accession_barcode_list
	full_barcode_list = list(set(full_barcode_list))
	accession_barcode_count = 0
	accession_nonbarcode_count = 0
	for barcodeID in full_barcode_count_dict[accession]:
		accession_barcode_count += full_barcode_count_dict[accession][barcodeID]
	for string in full_nonbarcode_dict[accession]:
		accession_nonbarcode_count += full_nonbarcode_dict[accession][string]
		full_nonbarcode_seq_list.append(string)
	full_read_count[accession] = (accession_barcode_count,accession_nonbarcode_count)
	full_nonbarcode_seq_list = list(set(full_nonbarcode_seq_list))

### Write final summary output tables
print("Writing barcode output summary tables")
accession_list = sorted(accession_list)
barcode_freq_dict = {}
barcode_count_dict = {}
table_outfile = open(output_dir+"all_barcode_count.table.txt","w")
freq_table_outfile = open(output_dir+"all_barcode_freq.table.txt","w")
full_barcode_list = sorted(full_barcode_list)
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
		read_count = full_read_count[accession][0]
		try:
			barcode_count = full_barcode_count_dict[accession][barcodeID]
			barcode_freq = (float(barcode_count)/float(read_count))
		except:
			barcode_count = 0
			barcode_freq = 0
		table_outfile.write("\t"+str(barcode_count))
		freq_table_outfile.write("\t"+str(barcode_freq))
		try:
			barcode_freq_dict[accession].append(barcode_freq)
			barcode_count_dict[accession].append(barcode_count)
		except:
			barcode_freq_dict[accession] = []
			barcode_freq_dict[accession].append(barcode_freq)
			barcode_count_dict[accession] = []
			barcode_count_dict[accession].append(barcode_count)
	table_outfile.write("\n")
	freq_table_outfile.write("\n")
table_outfile.close()
freq_table_outfile.close()

mismatch_table_outfile = open(output_dir+"all_mismatch_info.table.txt","w")
for a in range(0,len(accession_list)):
	accession = accession_list[a]
	mismatch_table_outfile.write("\t"+accession)
mismatch_table_outfile.write("\n")
for site in range(0,len(expected_amplicon_seq)):
	mismatch_table_outfile.write(str(site))
	for a in range(0,len(accession_list)):
		accession = accession_list[a]
		try:
			mismatch_prop = round(float(full_mismatch_dict[accession][site]),5)
		except:
			mismatch_prop = 0.0
		mismatch_table_outfile.write("\t"+str(mismatch_prop))
	mismatch_table_outfile.write("\n")
mismatch_table_outfile.close()


subID_dict = {}
for seq in full_nonbarcode_seq_list:
	if len(seq) == len(expected_amplicon_seq):
		for site in range(0,len(expected_amplicon_seq)):
			try:
				expected_barcode_nt_list = barcode_sites[site]
				barcode_site = True
			except:
				expected_barcode_nt_list = []
				barcode_site = False
			site_nt = seq[site]
			if barcode_site == False and site_nt != expected_amplicon_seq[site]:
				subID = str(site)+site_nt
				subID_dict[subID] = ''
			elif barcode_site == True:
				if site_nt != expected_barcode_nt_list[0] and site_nt != expected_barcode_nt_list[1]:
					subID = str(site)+site_nt
					subID_dict[subID] = ''


bases = ['A','T','C','G']
table_outfile = open(output_dir+"all_nonbarcode_freq.table.txt","w")
for a in range(0,len(accession_list)):
	accession = accession_list[a]
	table_outfile.write("\t"+accession)
table_outfile.write("\n")
for site in range(0,len(expected_amplicon_seq)):
	for sub_nt in bases:
		found = False
		subID = str(site)+sub_nt
		try:
			subID_dict[subID]
			found = True
		except:
			found = False
		if found == True:
			table_outfile.write(subID)
			for a in range(0,len(accession_list)):
				accession = accession_list[a]
				sub_count = 0
				for string in full_nonbarcode_dict[accession]:
					if string[site] == sub_nt:
						sub_count += full_nonbarcode_dict[accession][string]
				read_count = full_read_count[accession][0]+full_read_count[accession][1]
				if sub_count > 0:
					sub_freq = (float(sub_count)/float(read_count))
				else:
					sub_freq = 0
				table_outfile.write("\t"+str(sub_freq))
			table_outfile.write("\n")
table_outfile.close()

# if detect_barcode_content_first == False:
print("Calculating alpha-diversity indicies")
alpha_div_outfile = open(output_dir+"all_barcode.alpha_diversity.txt","w")
alpha_div_outfile.write("sampleID\tshannon_div\tsimpson_div\tchao_richness\tshannon_evenness\n")
for i in range(0,len(accession_list)):
	accession = accession_list[i]
	H_shannon,E_shannon = shannon_alpha(barcode_freq_dict[accession])
	H_simpson = simpson_alpha(barcode_freq_dict[accession])
	R_chao = chao_1_richness(barcode_count_dict[accession])
	alpha_div_outfile.write(accession+"\t"+str(H_shannon)+"\t"+str(H_simpson)+"\t"+str(R_chao)+"\t"+str(E_shannon)+"\n")
alpha_div_outfile.close()

print("Calculating Bray-Curtis dissimilarity matrix")
bray_outfile = open(output_dir+"all_barcode.bray-curtis_dissimilarity.txt","w")
bray_curtis_dict = {}
for j in range(0,len(accession_list)):
	bray_outfile.write("\t"+accession_list[j])
bray_outfile.write("\n")
for i in range(0,len(accession_list)):
	bray_outfile.write(accession_list[i])
	for j in range(0,len(accession_list)):
		if i == j:
			dist = 0
		else:
			accession1 = accession_list[i]
			accession2 = accession_list[j]
			try:
				dist = bray_curtis_dict[accession2][accession1]
			except:
				try:
					dist = bray_curtis_dissimilarity(barcode_freq_dict[accession1],barcode_freq_dict[accession2])
				except:
					dist = 'na'
				try:
					bray_curtis_dict[accession1][accession2] = dist
				except:
					bray_curtis_dict[accession1] = {}
					bray_curtis_dict[accession1][accession2] = dist
		bray_outfile.write("\t"+str(dist))
	bray_outfile.write("\n")
bray_outfile.close()

print("Calculating Jaccard dissimilarity matrix")
jaccard_outfile = open(output_dir+"all_barcode.jaccard_dissimilarity.txt","w")
jaccard_dict = {}
for j in range(0,len(accession_list)):
	jaccard_outfile.write("\t"+accession_list[j])
jaccard_outfile.write("\n")
for i in range(0,len(accession_list)):
	jaccard_outfile.write(accession_list[i])
	for j in range(0,len(accession_list)):
		if i == j:
			dist = 0
		else:
			accession1 = accession_list[i]
			accession2 = accession_list[j]
			try:
				dist = jaccard_dict[accession2][accession1]
			except:
				try:
					dist = jaccard_dissimilarity(barcode_freq_dict[accession1],barcode_freq_dict[accession2])
				except:
					dist = "na"
				try:
					jaccard_dict[accession1][accession2] = dist
				except:
					jaccard_dict[accession1] = {}
					jaccard_dict[accession1][accession2] = dist
		jaccard_outfile.write("\t"+str(dist))
	jaccard_outfile.write("\n")
jaccard_outfile.close()

print("Processing Complete")
