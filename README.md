# BarcodeID

## What this script does

This script is designed to take in amplicon sequences that have nucleotide barcodes at specific sites. There is no minimum or maximum number of barcodes, but this script currently only permits two possible nucleotides per site (e.g. a barcode at site 33 can be either A or T, but not C or G). This script process the illumina sequence files to generate tables that summarize the frequency of all barcodes detected in each sample, as well as reporting alpha and beta diversity indicies for those samples.

It is written to perform barcoded amplicon sequence analysis for the Lowen Lab's influenza research, but it should work just as well with any barcoded amplicon library with a similar design. 

# How to use this script

In it's simplest use case, the user can put their compressed sequence files in a folder named "raw_data" and immediately start the analysis by entering:
```
python barcodeID_v1.0.py
```
in the appropriate terminal and following the text prompts for required information.

* In this case, the script will automatically detect if sequence files are compressed (.tar, .tar.gz, or .gz), unzip them in parallel if necessary, then predict forward and reverse sequence pairs from any detectible sequence files (.fastq or .fq).
   * If the script is unable to predict the forward and reverse sequence file extensions (e.g. *"_R1.fastq"* and *"_R2.fastq"*), the script will ask the user to input them manually in the "User defined variables" section at the top of the script.
* BarcodeID will then predict the barcoded sites within the amplicon library after the user indicates the number of barcodes that they expect to find in their library, using either a random or user-supplied sample. If the user verifies the prediction is correct through a text prompt, BarcodeID will immediately perform analysis on all remaining samples.
   * In the event that BarcodeID is unable to detect the correct number of barcodes, or if the user indicates that they are incorrect based on their expectation, a file containing the predicted values will be generated called *"amplicon_info.auto.txt"*, which the user can use as a template to correct the error.


# How to install the dependencies of this script

The requirements to run this pipeline are:
1. Python3
   * We recommend using [miniconda](https://docs.conda.io/en/latest/miniconda.html) to install python and the required packages
2. [BBTools](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/)
   * Note that the bbtools available on anaconda is not currently a functional copy, so it needs to be downloaded from the JGI website linked above
3. Java (to run the BBTools)
4. Bash (also required to run BBTools)
   * This is available by default on Mac and Linus operating systems
   * For windows users, we recommend using [GitBash](https://git-scm.com/downloads). In our testing, WSL was not a reliable method to perform the analysis.
	  * To use GitBash, the user will need to edit their *".bashrc" file in *"C:/Users/username"* to GitBash to add the location of their conda profile, the location of their python3 installation, and add the location of their bbtools installation to the GitBash path
	     * An example of how you might edit your .bashrc file:
> . /c/ProgramData/miniconda3/etc/profile.d/conda.sh
> alias python="winpty /c/ProgramData/miniconda3/python.exe"
> PATH=$PATH:"/c/ProgramData/bbmap"

Required python packages:
1. numpy
2. multiprocessing
3. joblib

   * Using conda, these packages can be installed using the following conda commands:
      - conda install numpy
	  - conda install multiprocess
	  - conda install joblib

While *multiprocessing* and *joblib* are not necessary for single-core processing, they are highly recommend not only because parallel processing is substantially faster and can scale with the number of cores available, but the single core processing option has not been as rigorously tested as the parallel processing option.

# Interpreting the results

Output files are generated in the folder *"./barcode_info/"*
* After the pipeline has finished running, this folder will contain:
   1. *"all_barcode_count.table.txt"*
      * A table containing all of the raw barcode counts for all samples
   2. *"all_barcode_freq.table.txt"*
      * A barcode frequency table for all samples
   3. *"all_mismatch_info.table.txt"*
      * A table that shows the frequency at which reads were discarded for having an unexpected nucleotide while also passing the user specified quality threshold
	  * This table is a very useful tool for determining if reads are being discarded erroneously, or if there are non-barcode alleles sweeping your samples through positive selection
   4. *"all_nonbarcode_freq.table.txt"*
      * A table showing the frequency of unexpected alleles in your samples across all reads
   5. *all_barcode.alpha_diversity.txt*
      * This file contains Shannon diversity, Simpson diversity, Chao1 richness, and Shannon evenness values calculated for all samples
   6. *"all_barcode.bray-curtis_dissimilarity.txt"*
      * This table reports the bray-curtis dissimilarity between all samples (0 means identical, 1 means there is no overalp in barcodes, and intermediate values reflect differences in both barcode presence and barcode frequency.
   7. *"all_barcode.jaccard_dissimilarity.txt"*
      * This table reports the jaccard dissimilarity between all samples (0 means identical, 1 means there is no overalp in barcodes. This metric only considers overlap in barcode presence between samples, not their frequency).

The intermediate files generated for each sample are located in *"./barcode_info/samples/"*
