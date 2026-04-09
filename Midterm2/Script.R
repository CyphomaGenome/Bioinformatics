#Part 1: Phylogenetics

#Load packages (not all will be used)
library(Biostrings)
library(seqinr)
library(pwalign)
library(msa)
library(rentrez)
library(RADami)
library(ips)
library(DECIPHER)
library(UniprotR)
library(GenomicAlignments)
library(protti)
library(r3dmol)


###Convert to .phy format:
#Load in .fasta
Sequences<-readDNAStringSet("metazoa_alignment.5k.fasta")
Sequences #print
#Write .phy file
write.DNAStringSet(x = Sequences, format = "phylip", filename = "Sequences.phy")

#After I ran this, I realized that the last few steps are unnecessary as RAxML can handle .fasta files.

#Within Bash
#Navigate to proper folder
#cd /mnt/c/Bioinformatics/Bioinformatics/Midterm2

#I chose to use a maximum likelihood phylogenetic software (RAxML) instead of a Batesian one because the data set is fairly large.
#The analysis could be also done with a Batesian phylogenetic software, such as BEAST 2, but it would require much more computation time.
#This is a slight sacrifice of biological accuracy for a shorter computation time.

#Run raxml

#../../raxml-ng_v2.0.0_linux_x86_64/raxml-ng --all --msa Sequences.phy --model GTR+G --tree pars{10} --bs-trees 200

###Part 2: Protein

#read in the sequences as a DNA stringset
Part2Sequences<-readDNAStringSet("metazoa_alignment.gene.fasta")
Part2Sequences #print
#print names to find Homo_sapiens
names(Part2Sequences)
#Isolate Homo_sapiens
Homo_sapiens<-Part2Sequences[5]
Homo_sapiens #print
#These next few steps are a method for removing gaps that I figured out for the final project so I have adapted it here
#I would guess that this would also be possible without conversion to a DNAbin
alignment_bin <- as.DNAbin(as.matrix(Homo_sapiens))
trimmed_alignment <- deleteGaps(alignment_bin, gap.max = nrow(alignment_bin) * 1)
trimmed_alignment #print
#convert back to DNASringSet, apply is a loop function, with "1" making it go row by row and collapse = "" making the character strings collapse to one string without gaps
trimmed_stringset <- DNAStringSet(apply(as.character(trimmed_alignment), 1, paste, collapse = ""))
trimmed_stringset #print

#actually here is the code without the conversion (this uses the DECIPHER package)
trimmed_stringset2<-RemoveGaps(Homo_sapiens)
trimmed_stringset2 #print

identical(trimmed_stringset, trimmed_stringset2) #The two methods yield an identical result!

#convert to amino acids
Protein<-Biostrings::translate(trimmed_stringset)
Protein #print
class(Protein) #Check class
#Write .fasta file
writeXStringSet(Protein, "Protein.fasta", format="fasta")

#after putting the sequence in BLASTp, I obtained the following result
#DNA polymerase subunit gamma-1 [Homo sapiens]

#10. P54098 inside of UniProt
UniProtAccessionNumbers<-read.csv("UniProtAccessionNumbers.txt") #loads txt file into workspace
UniProtAccessionNumbers #print
class(UniProtAccessionNumbers) #check format (format was data frame)
AccessionString<-UniProtAccessionNumbers[, "AccessionNumbers"] #stores accession numbers in new variable
AccessionString #print
class(AccessionString) #check format (format was character which is correct)
#now I am realizing that I could have just written the accession number into GetProteinGOInfo
ProteinGOInfo<-GetProteinGOInfo(AccessionString) #gets gene ontology terms of UniProt accession numbers
ProteinGOInfo #print
#P54098
#biological.process: DNA-templated DNA replication [GO:0006261]
#molecular.function: DNA-directed DNA polymerase activity [GO:0003887]
#cellular.component: gamma DNA polymerase complex [GO:0005760]

class(ProteinGOInfo) #check format (format was data frame which is correct)
PlotGoInfo(ProteinGOInfo, directorypath="C:/Bioinformatics/Bioinformatics/Midterm2") #creates a graphical representation of the GO info (the graphic is in a file called "GOPlotInfo.jpeg")
#files are in GitHub
