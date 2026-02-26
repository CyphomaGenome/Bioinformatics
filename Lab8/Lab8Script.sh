#!/usr/bin/env zsh

# change working directory to the root
cd /mnt/c/Bioinformatics

wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash ~/Miniconda3-latest-Linux-x86_64.sh
# I installed conda and the following command produces a list as desired.
conda list

wget https://raw.githubusercontent.com/faircloth-lab/phyluce/v1.7.3/distrib/phyluce-1.7.3-py36-Linux-conda.yml
conda env create -n phyluce-1.7.3 --file phyluce-1.7.3-py36-Linux-conda.yml
# I installed phyluce

# this next step sets up an x86 environment
CONDA_SUBDIR=osx-64 conda create -n phyluce-1.7.3-x86 python=3.6
# this line 'activated' the conda environment
conda activate phyluce-1.7.3-x86
# configuring the environment
conda config --env --set subdir linux-64
# installing from the yml file
conda env update --name phyluce-1.7.3-x86 --file phyluce-1.7.3-py36-Linux-conda.yml --prune
# if you are finished remove pound sign
conda deactivate

# check to see if you have the environment installed
conda env list 

# phyluce can then be activate on all future cases  
# using this line (you may need to update the version number):
conda activate phyluce-1.7.3-x86

# check that a few of the tools are working, and peruse the options:
samtools --help
raxml-ng --help
#both were working


# now, change your working directory to a folder for this project
cd /mnt/c/Bioinformatics/Bioinformatics
#I then realized that I was previously supposed to be in the home directory, but I don't believe this changed anything thankfully.

# create a project directory
mkdir Lab8

# change to that directory
cd Lab8

# download the .fastq data into a file named fastq.zip
wget -O fastq.zip https://ndownloader.figshare.com/articles/1284521/versions/2

# make a directory to hold the data
mkdir raw-fastq

# move the zip file into that directory
mv fastq.zip raw-fastq

# move into the directory we just created
cd raw-fastq

# unzip the fastq data
unzip fastq.zip

# if the unzip worked, you can delete the zip file
rm fastq.zip

# you should see 6 files in this directory now
ls -l
#there are 8 zip files but I am assuming that is ok

# count the read data
for i in *_R1_*.fastq.gz; do echo $i; gunzip -c $i | wc -l | awk '{print $1/4}'; done

# look at the first few lines of one of the read files
# I've given two options in case one doesn't work on your computer
zcat Alligator_mississippiensis_GGAGCTATGG_L001_R1_001.fastq.gz | head -n 20
gzip -cd Alligator_mississippiensis_GGAGCTATGG_L001_R1_001.fastq.gz | head -n 20
#both commands worked

# clean the read data
# first, move "up" one folder
cd ..

# make sure that the illumiprocessor.conf file is in your Lab8 folder
ls
#it is in the folder

# run illumiprocessor
# change the number of cores to what your computer has
# you can check how many you have with:
getconf _NPROCESSORS_ONLN
# first, look at the help page for illumiprocessor, especially the quality scores
# check the code above for how to access a help page

# then run illumiprocessor
illumiprocessor \
    --input raw-fastq/ \
    --output clean-fastq \
    --config illumiprocessor.conf \
    --cores 8

# look at the first few lines of one of the cleaned read files
# does it look different than the raw reads?
cd clean-fastq/alligator_mississippiensis/split-adapter-quality-trimmed/
# two ways to do this:
zcat alligator_mississippiensis-READ1.fastq.gz | head -n 20
gzip -cd alligator_mississippiensis-READ1.fastq.gz | head -n 20
#the reads are now trimmed, removing ends with low quality bases.

# move back up three levels to the Lab8 folder
cd ../../..

# make sure that you're in the Lab8 folder
pwd

# move to the directory holding our cleaned reads
cd clean-fastq/

# this code is for quality control
# it checks the lengths of all the read files in all your folders
# the file "clean_read_counts_header.csv" has the column header information for this output
for i in *;
do
    phyluce_assembly_get_fastq_lengths --input $i/split-adapter-quality-trimmed/ --csv;
done

# from the tutorial, copy and paste the assembly configuration data into a file 
# named assembly.conf and save it to your Lab8 folder

# move back up one folder to Lab8
cd ..

# run the assembly
# this will take a few minutes
phyluce_assembly_assemblo_spades \
    --conf Assembly.conf.txt \
    --output spades-assemblies \
    --cores 8

# if this doesn't work, try a different assembler
# https://phyluce.readthedocs.io/en/latest/daily-use/list-of-programs.html
#spades worked for me once the files were properly named

# check the quality of your assemblies
for i in spades-assemblies/contigs/*.fasta;
do
    phyluce_assembly_get_fasta_lengths --input $i --csv;
done

#I saved them to the csv file with this command
for i in spades-assemblies/contigs/*.fasta; 
do     
    phyluce_assembly_get_fasta_lengths --input $i --csv; 
done>>"clean_read_counts_header.csv"

#I then installed seqkit to find N50 and stored the statistics with more information
conda install -c bioconda seqkit
seqkit stats -a spades-assemblies/contigs/*.fasta -T | tr '\t' ',' >"clean_read_full_stats.csv"

#########################################################################
### THIS WAS THE LAST REQUIRED STEP OF THE LAB. THE REST IS OPTIONAL. ###
#########################################################################


# get the probe set
wget https://raw.githubusercontent.com/faircloth-lab/uce-probe-sets/master/uce-5k-probe-set/uce-5k-probes.fasta

# run the phyluce_assembly_match_contigs_to_probes program to find UCE loci
phyluce_assembly_match_contigs_to_probes \
    --contigs spades-assemblies/contigs \
    --probes uce-5k-probes.fasta \
    --output uce-search-results

# create an output directory for this taxon set - this just keeps
# things neat
mkdir -p taxon-sets/all

# create the data matrix configuration file
phyluce_assembly_get_match_counts \
    --locus-db uce-search-results/probe.matches.sqlite \
    --taxon-list-config taxon-set.conf.txt \
    --taxon-group 'all' \
    --incomplete-matrix \
    --output taxon-sets/all/all-taxa-incomplete.conf

# check out the contents. What is in this file?
head taxon-sets/all/all-taxa-incomplete.conf 
#contains organisms and loci

# change to the taxon-sets/all directory
cd taxon-sets/all

# make a log directory to hold our log files - this keeps things neat
mkdir log

# get FASTA data for taxa in our taxon set
phyluce_assembly_get_fastas_from_match_counts \
    --contigs ../../spades-assemblies/contigs \
    --locus-db ../../uce-search-results/probe.matches.sqlite \
    --match-count-output all-taxa-incomplete.conf \
    --output all-taxa-incomplete.fasta \
    --incomplete-matrix all-taxa-incomplete.incomplete \
    --log-path log

# explode the monolithic FASTA by taxon (you can also do by locus)
phyluce_assembly_explode_get_fastas_file \
    --input all-taxa-incomplete.fasta \
    --output exploded-fastas \
    --by-taxon

# get summary stats on the FASTAS
for i in exploded-fastas/*.fasta;
do
    phyluce_assembly_get_fasta_lengths --input $i --csv;
done

# make sure we are in the correct directory
pwd
#cd Lab8/taxon-sets/all

# align the UCE data
# make sure to adjust your cores
phyluce_align_seqcap_align \
    --input all-taxa-incomplete.fasta \
    --output mafft-nexus-edge-trimmed \
    --taxa 4 \
    --aligner mafft \
    --cores 8 \
    --incomplete-matrix \
    --log-path log

# you will get warnings. Those are loci for which you did not get data
# for some individuals

# last step is to open some of the resulting alignments in the text
# editor of your choice

# get the final concatenated alignment, if you want
phyluce_align_concatenate_alignments \
    --alignments mafft-nexus-edge-trimmed \
    --output mafft-nexus-edge-trimmed-raxml \
    --phylip \
    --log-path log

#I was able to succesfully finish the entire lab