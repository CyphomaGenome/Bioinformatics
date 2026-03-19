#Start by loading in all packages that may be used during the project.
#Not all of these will necessarily be used but it is helpful to have them available

#BiocManager::install("GenomicAlignments")
#install.packages("UniprotR")
#install.packages("protti")
#install.packages("r3dmol")
library(UniprotR)
library(GenomicAlignments)
library(protti)
library(r3dmol)
library(seqinr)
library(pwalign)
library(msa)
library(rentrez)
library(rBLAST)

#Unsuccessful attempt below

#Next, we take the amino acid data from lab 6 and load it into a variable.
AminoAcidSequence<-readAAStringSet("AminoAcidSequence.fasta") 
#The next few lines of code create a BLAST database and search for the lab 6 sequence.
#This is where I initially made a mistake, because this searches in NCBI, not UniProt, so the accession numbers will be incorrect.
Protein_Database<-blast(db="nr", remote=TRUE, type="blastp") #creates database
Protein_Results<-predict(Protein_Database,AminoAcidSequence) #searches database
Protein_Results #prints results
class(Protein_Results) #check to make sure class is data.frame
AccessionNumbers<-Protein_Results[, "sseqid"] #stores accession numbers in new variable
AccessionNumbers #prints accession numbers
class(AccessionNumbers) #check to make sure class is character
GetProteinGOInfo(AccessionNumbers) #gets gene ontology terms of UniProt accession numbers
#Result was many lines of "Bad request. The resource you requested doesn't exist or There is a problem with your input."
#This is when I realized that I searched in the wrong database, so I decided to try using UniProt ID mapping
cat(AccessionNumbers, file = "AccessionNumbers.txt") #creates text file of accession numbers
#I then used Uniprot ID mapping, however no results returned

#Attempt 2 (also somewhat unsuccessful)

#Since the previous method failed, I tried entering the amino acid sequence directly into UniProt
#The best results were around a 67% percent identity with E-values of zero. 
#The entries were all for the Cytochrome P450 protein in marine gastropods such as the rough periwinkle, though not any species in the same superfamily as Cyphoma gibbosum.
#I copied the top 5 accession numbers into a txt file in my project folder.
UniProtAccessionNumbers<-read.csv("UniProtAccessionNumbers.txt") #loads txt file into workspace
UniProtAccessionNumbers #print
class(UniProtAccessionNumbers) #check format (format was data frame)
AccessionString<-UniProtAccessionNumbers[, "AccessionNumbers"] #stores accession numbers in new variable
AccessionString #print
class(AccessionString) #check format (format was character which is correct)

ProteinGOInfo<-GetProteinGOInfo(AccessionString) #gets gene ontology terms of UniProt accession numbers
ProteinGOInfo #print
class(ProteinGOInfo) #check format (format was data frame which is correct)
#I then tried the following lines of code alongside many others with extensive error messages
PlotGoInfo(GOObj = ProteinGOInfo, directorypath="C:/Bioinformatics/Bioinformatics/Lab10")
PlotGOAll(GOObj = ProteinGOInfo, Top = 10, directorypath = getwd(), width = 8, height = 5)
GetPathology_Biotech(ProteinGOInfo) 
Get.diseases(ProteinGOInfo) 
#I believe that all of the N/A data entries were involved with the problem, but I wasn't able to resolve any of the errors.

#Attempt 3 (completes lab activity, but some data returns with little information or "NULL" response)

#At this point, I decided to use the provided accession numbers "P0A799" and "P08839" because my Cytochrome p450 sequences did not have enough available info
ProvidedAccessionNumbers<-c("P0A799", "P08839") #load accession numbers into a variable
ProvidedAccessionNumbers #print
class(ProvidedAccessionNumbers) #check format (format was character which is correct)

GOInfo<-GetProteinGOInfo(ProvidedAccessionNumbers) #gets gene ontology terms of UniProt accession numbers
GOInfo #print
class(GOInfo) #check format (format was data frame which is correct)
PlotGoInfo(GOInfo, directorypath="C:/Bioinformatics/Bioinformatics/Lab10") #creates a graphical representation of the GO info (the graphic is in a file called "GOPlotInfo.jpeg")
#this command worked now that there were not as many N/A data entries
PlotGOAll(GOObj = GOInfo, Top = 10, directorypath = getwd(), width = 8, height = 5) #handy visualization for publications, creates a file called "GO All.jpeg" (equipped with an empty space to destroy everything).

Pathology<-GetPathology_Biotech(ProvidedAccessionNumbers, directorypath = getwd()) #fetches pathology information
Pathology #print
class(Pathology)
Diseases<-Get.diseases(Pathology, directorypath = getwd()) #fetches disease information
Diseases #print
#this resulted in the response "NULL"

ProttiInfo<-fetch_uniprot(ProvidedAccessionNumbers) #access uniprot information
ProttiInfo #print
pdbs<-c(ProttiInfo$accession,ProttiInfo$xref_pdb) #find pdb accession numbers
pdbs #print pdb accession numbers
pdbinfo<-fetch_pdb(c("1EZA","1EZB", "1EZC", "1EZD", "1ZYM", "2EZA", "2EZB", "2EZC", "2HWG", "2KX9", "2L5H", "2MP0", "2N5T", "2XDF", "3EZA", "3EZB", "3EZE", "6V9K", "6VU0", "1ZMR")) 
pdbinfo #prints info about pdb accession numbers

threeDinfo<-fetch_alphafold_prediction(ProvidedAccessionNumbers) #fetch 3D info
threeDinfo #print 3D info

#I then isolated one of the models (independent of last few steps)
Model<-fetch_alphafold_prediction("P08839")
Model
class(Model) #shows variable is list, needs to be corrected
protein_table <- Model$P08839
protein_table
class(protein_table) #now shows tbl_df, "tbl", "data.frame"

#I could not get the following lines of code to work so I used the Alpha Fold Website to observe the structure of P08839 (Phosphoenolpyruvate-protein phosphotransferase)
r3dmol() %>%
  m_add_model(data = m_bio3d(protein_table)) %>%
  m_set_style(style = m_style_cartoon(color = "spectrum")) %>%
  m_zoom_to()