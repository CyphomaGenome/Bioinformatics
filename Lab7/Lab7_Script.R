library(Biostrings)
library(seqinr)
library(pwalign)
library(msa)
library(rentrez)

#reading the fasta files into the document
sordidulus<-readDNAStringSet("Contopus_sordidulus_61800_ref.fasta") #converts the fasta file to a DNA strings set variable
names(sordidulus) <- "sordidulus" #makes name easier to read
sordidulus #print for visualization
virens<-readDNAStringSet("Contopus_virens_768_ref.fasta") #converts the fasta file to a DNA strings set variable
names(virens) <- "virens" #makes name easier to read
virens #print for visualization
individual_30620<-readDNAStringSet("Contopus_sp_30620_consensus_curated.fasta") #converts the fasta file to a DNA strings set variable
names(individual_30620) <- "individual_30620" #makes name easier to read
individual_30620 #print for visualization
individual_27325<-readDNAStringSet("Contopus_sp_27325_consensus_curated.fasta") #converts the fasta file to a DNA strings set variable
names(individual_27325) <- "individual_27325" #makes name easier to read
individual_27325 #print for visualization

#running an msa alignment
Alignment<-msaMuscle(c(sordidulus,virens,individual_27325,individual_30620), type="dna") #performs a multiple sequence alignment on the sequences using Muscle
Alignment

Alignment_seqinrformat<-msaConvert(Alignment, type="seqinr::alignment") #converts the multiple sequence alignment to an alignment within the seqinr package
Alignment_seqinrformat #print for visualization
d <- dist.alignment(Alignment_seqinrformat, "identity") #creates a similarity matrix of the sequences
d #print matrix to see percent identity scores

#The difference between 30620 and virens is 0.29836349 while the difference between 30620 and sordidulus is 0.08391814
#The difference between 27325 and virens is 0.29008400 while the difference between 27325 and sordidulus is 0.07753724
#This seems to suggest that both individual 30620 and individual 27325 are both of the species Contopus sordidulus