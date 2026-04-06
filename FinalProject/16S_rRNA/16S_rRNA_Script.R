library(Biostrings)
library(seqinr)
library(pwalign)
library(msa)
library(rentrez)
library(ape)

fasta_files <- list.files(pattern = "\\.fasta$")
fasta_files
Strings<-lapply(fasta_files, readDNAStringSet)
Strings
combined_sequences <- do.call(c, Strings)
combined_sequences
names(combined_sequences)
names(combined_sequences) <- c("Aforia", "Antillophos", "Babyloni_j", "Babyloni_l", "Biplex", "Buccinum_a", "Buccinum_m", "Burnupena", "Cancell_ca", "Cancell_co", "Cancell_si", "Cantharus", "Colubrar_n", "Colubrar_o", "Cominella", "Comitas", "Conus", "Cronia", "Drupella", "Engina", "Ficus", "Fusinus", "Fusitriton", "Granulifus", "Hemifusus", "Japeuthria", "Kelletia_k", "Kelletia_l", "Latiromitr", "Littorina", "Melongena", "Metula", "Microvolut", "Mitra", "Mitrella", "Muricanthu", "Nassariu_f", "Nassariu_s", "Neobuccinu", "Neptunea_l", "Neptunea_p", "Nucella", "Oliva_m", "Oliva_s", "Olivella", "Paraeuthri", "Penion", "Peristerni", "Pisania_p", "Pisania_s", "Plesiotrit", "Pollia", "Polystira", "Pomacea", "Siphonalia", "Siratus", "Sylvanococ", "Thais_c", "Thais_h", "Thalessa", "Tonna", "Vexillum", "Volema", "Zeuxis")
combined_sequences
names(combined_sequences)
# remove "Sylvanococ"
combined_sequences <- combined_sequences[names(combined_sequences) != "Sylvanococ"]
# remove "Pomacea"
combined_sequences <- combined_sequences[names(combined_sequences) != "Pomacea"]
combined_sequences
names(combined_sequences)
MSA_Alignment<-msaMuscle(combined_sequences)
MSA_Alignment
print(MSA_Alignment, show="complete")
alignment_set <- as(MSA_Alignment, "DNAStringSet")
MSA_Subset<-subseq(alignment_set, start = 79, end = 1496)
MSA_Subset
#install.packages("ips")
library(ips)
alignment_bin <- as.DNAbin(as.matrix(MSA_Subset))
#On next step adjust multiplier for level of trimming (1=untrimmed, 0.5 means 90%=trim where 90% of sequences have gaps)
trimmed_alignment <- deleteGaps(alignment_bin, gap.max = nrow(alignment_bin) * 1)
trimmed_alignment
write.dna(trimmed_alignment, 
          file = "gastropod_16S_NoPomacea.phy", 
          format = "sequential", 
          nbcol = -1, 
          colsep = "")

#Navigate to proper folder
#cd /mnt/c/Bioinformatics/Bioinformatics/FinalProject/16S_rRNA

#Run raxml
#../../../raxml-ng_v2.0.0_linux_x86_64/raxml-ng --all --msa gastropod_16S_NoPomacea.phy --model GTR+G --tree pars{10} --bs-trees 200

#Once tree is finished
# 1. Load your tree
my_tree <- read.tree("gastropod_16S_NoPomacea.phy.raxml.support")

# 2. Create a translation vector
# The names must match exactly what is currently in the tree
old_names <- my_tree$tip.label
old_names
new_names <- c("Neptunea_purpurea", "Neptunea_lamellosa", "Buccinum_aniwanum", "Buccinum_morchianum", "Japeuthria_ferrea", 
               "Siphonalia_cassidariaeformis", "Kelletia_kelletii", "Kelletia_lischkei", "Penion_chathamensis", "Peristernia_sulcata", 
               "Neobuccinum_eatoni", "Pollia_tincta", "Engina_mendicaria", "Pisania_striata", "Pisania_pusio", 
               "Cantharus_multangulus", "Granulifusus_niponicus", "Fusinus_akitai", "Colubraria_nitidula", "Colubraria_obscura", 
               "Metula_amosi", "Nassarius_semisulcatus", "Zeuxis_siquijorensis", "Nassarius_festivus", "Microvoluta_sp", 
               "Vexillum_plicarium", "Latiromitra_sp", "Ficus_subintermedia", "Oliva_mustelina", "Oliva_spicata", 
               "Olivella_volutella", "Mitra_lens", "Polystira_picta", "Comitas_kaderlyi", "Aforia_magnifica", 
               "Conus_praecellens", "Tonna_luteostoma", "Biplex_perca", "Fusitriton_oregonensis", "Plesiotriton_sp", 
               "Littorina_brevicula", "Cancellaria_cooperi", "Cancellaria_cancellata", "Cancellaria_sinensis", "Babylonia_japonica", 
               "Babylonia_lutosa", "Mitrella_bicincta", "Thalessa_savignyi", "Thais_clavigera", "Thais_haemastoma", 
               "Muricanthus_radix", "Siratus_beauii", "Cronia_sp", "Drupella_cornus", "Nucella_lapillus", "Hemifusus_tuba", "Volema_myristica", "Melongena_patula", "Cominella_adspersa", "Antillophos_laevis", 
               "Paraeuthria_plumbea", "Burnupena_cincta")
new_names

# 3. Swap the names
my_tree$tip.label <- new_names

# 4. Save the renamed tree
write.tree(my_tree, "gastropod_final_named.tre")




###################################################################################
STORAGE


library(DECIPHER)
codon_aligned <- AlignTranslation(combined_sequences)
codon_aligned
trimmed_dna <- subseq(codon_aligned, start = 1, end = 1545)
trimmed_dna
#BiocManager::install(c("IRanges", "XVector"), force=TRUE)
#remotes::install_github("andrew-hipp/RADami", upgrade = "never")
library(RADami)
write.DNAStringSet(x = trimmed_dna, format = "phylip", filename = "CytrochomeMSA.phy")


#names(combined_sequences) <- c()
#combined_sequences
#writeXStringSet(combined_sequences, "botulinumcompletesequences.fasta")


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

