library(Biostrings)
library(seqinr)
library(pwalign)
library(msa)
library(rentrez)

Strings<-lapply(c("Assimine_sp.txt","Boreoelona_ussuriensis.txt","Haliotis_ovina.txt", "Haliotis_varia.txt", "Lacuna_carinifera.txt","Lepetodrilus_sp.txt", "Peasiella_habei.txt", "Pseudorimula_sp.txt", "Tylomelania neritiformis.txt", "Tylomelania_hannelorae.txt"), readDNAStringSet)
Strings
combined_sequences <- do.call(c, Strings)
combined_sequences
names(combined_sequences) <- c("Assimine","Boreoelona","Haliotis_o", "Haliotis_v", "Lacuna", "Lepetodril", "Peasiella", "Pseudorimu", "Tylomelan_n", "Tylomelan_h")
combined_sequences

library(DECIPHER)
codon_aligned <- AlignTranslation(combined_sequences)
codon_aligned
trimmed_dna <- subseq(codon_aligned, start = 1, end = 1545)
trimmed_dna
#BiocManager::install(c("IRanges", "XVector"), force=TRUE)
#remotes::install_github("andrew-hipp/RADami", upgrade = "never")
library(RADami)
write.DNAStringSet(x = trimmed_dna, format = "phylip", filename = "CytrochomeMSA.phy")


###################################################################################
STORAGE

#names(combined_sequences) <- c()
#combined_sequences
#writeXStringSet(combined_sequences, "botulinumcompletesequences.fasta")
#MSA_Alignment<-msaMuscle(combined_sequences)
#MSA_Alignment
#print(MSA_Alignment, show="complete")

#these are valid lines of code but not neccesary
phylip_format <- msaConvert(MSA_Alignment, type="phangorn::phyDat")
phylip_format
#next single line of code writes phylip file
write.phylip(MSA_Alignment, file="botulinum.phy")

#Then go to Bash
#cd /mnt/c/Bioinformatics/botulinum
#iqtree2 -s botulinum.phy
#botulinum.phy.iqtree has the information


#Commands used for iqtree installation
#wget https://github.com/iqtree/iqtree2/releases/download/v2.2.2.7/iqtree-2.2.2.7-Linux.tar.gz
#/mnt/c/Bioinformatics/iqtree-2.2.2.7-Linux/bin
#export PATH="$PATH:/mnt/c/Bioinformatics/iqtree-2.2.2.7-Linux/bin"

#NOTE: command for DNA with codons
#iqtree -s coding_gene.phy -st CODON 
#NOTE: command for DNA with codons using invertebrate mitochondrial DNA
#iqtree -s coding_gene.phy -st CODON5

#BiocManager::install("DECIPHER")
library(DECIPHER)
codon_aligned <- AlignTranslation(combined_sequences)
codon_aligned
trimmed_dna <- subseq(codon_aligned, start = 1, end = 4119)
trimmed_dna
#BiocManager::install(c("IRanges", "XVector"), force=TRUE)
#remotes::install_github("andrew-hipp/RADami", upgrade = "never")
library(RADami)
write.DNAStringSet(x = trimmed_dna, format = "phylip", filename = "botulinumcodons.phy")

#Then go to Bash
#cd /mnt/c/Bioinformatics/botulinum
#iqtree2 -s botulinumcodons.phy -st CODON
#botulinumcodons.phy.iqtree has the information

#next lines of code try to use better codon model
write.DNAStringSet(x = trimmed_dna, format = "phylip", filename = "botulinumcodonmodel.phy")
#iqtree2 -s botulinumcodonmodel.phy -m MFP+CODON -st CODON

Strings2<-lapply(c("botulinum.fasta","boNTG.txt","boNTE.txt","boNTF.txt","boNTB.txt", "boNTEn.fasta"), readDNAStringSet)
Strings2
combined_sequences2 <- do.call(c, Strings2)
combined_sequences2
names(combined_sequences2) <- c("E1", "E2", "E3", "E4", "E5", "E6", "E7", "C1", "A1", "A2", "B1", "B2", "F1", "B3", "B4", "G1", "G2", "E8", "F2", "B5", "En1")
combined_sequences2
writeXStringSet(combined_sequences2, "botulinumcompletesequenceswithancestral.fasta")
MSA_Alignment2<-msaMuscle(combined_sequences2)
MSA_Alignment2
print(MSA_Alignment2, show="complete")
codon_aligned2 <- AlignTranslation(combined_sequences2)
codon_aligned2
trimmed_dna2 <- subseq(codon_aligned2, start = 1, end = 4230)
trimmed_dna2
write.DNAStringSet(x = trimmed_dna2, format = "phylip", filename = "botulinumwithancestralcodons.phy")

#Then go to Bash
#cd /mnt/c/Bioinformatics/botulinum
#iqtree2 -s botulinumwithancestralcodons.phy -st CODON
#botulinumwithancestralcodons.phy.iqtree has the information


#additional commands
#/mnt/c/Bioinformatics/iqtree-1.6.12-Linux/iqtree-1.6.12-Linux/bin
#export PATH="$PATH:/mnt/c/Bioinformatics/iqtree-1.6.12-Linux/iqtree-1.6.12-Linux/bin"
#iqtree -s botulinumcodons.phy -st CODON
#still does not work

