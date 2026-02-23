#Comments
###Answers
#####Questions

#Before starting the assignment, I loaded these packages (I already had most of them installed, but steps for installing them are also included here).
#BiocManager::install("Biostrings") #remove the pound sign if needed
#install.packages("seqinr") #remove the pound sign if needed
#install.packages("pwalign") #remove the pound sign if needed
#BiocManager::install("msa") #remove the pound sign if needed
#install.packages("rentrez") #remove the pound sign if needed
library(Biostrings)
library(seqinr)
library(pwalign)
library(msa)
library(rentrez)
#Not all of these packages will necessarily be used, but I want to have them all available if I need them.

#These next steps were my process of installing rBLAST and ensuring that it was functional.
#I recieved an error message when I tried to install rBLAST using install.packages("rBLAST")
#I went to https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ and downloaded the BLAST+ installer for windows and followed the installation process.
#I then refreshed R using .rs.restartR()
#install.packages("rBLAST") still did not work so I tried BiocManager::install("rBLAST") which did work.
library(rBLAST)
#I still had errors when running Sys.which("blastn"), so I added the bin file inside of the control panel on my computer, environmental variables, program variables.
#After restarting again (this time completely closing RStudio), Sys.which("blastn") returned a path.
#I am still not completely sure what problems were fixed or what fixed them, but now rBLAST works so I can run BLAST inside of R :)
#I believe the following two lines of code do the same thing
#install.packages("devtools")
#devtools::install_github("mhahsler/rBLAST")

#These next three lines of code were not used in this asignment because I would have had to register for an API code to use OMIM, so I just used the web database.
#To install romim to search the OMIM database, the following command worked for me.
#devtools::install_github("davetang/romim")
#library(romim)

#####1.	Import and align your DNA sequences

sequences<-readDNAStringSet("sequences.fasta") #converts the fasta file to a DNA strings set variable
sequences #print for visualization
Alignment<-msaMuscle(sequences, type="dna") #performs a multiple sequence alignment on the sequences using Muscle
Alignment #print for visualization

#####2.	Find a way to measure how “good” your alignment is. 
#####Is this a good alignment or a bad one? How can you tell?

Alignment_seqinrformat<-msaConvert(Alignment, type="seqinr::alignment") #converts the multiple sequence alignment to an alignment within the seqinr package
Alignment_seqinrformat #print for visualization
d <- dist.alignment(Alignment_seqinrformat, "identity") #creates a similarity matrix of the sequences
d #print matrix to see percent identity scores
###The alignment is good. The percent identity of 0 to 0.11 shows that the sequences are all homologous

#####3.	Calculate the consensus sequence for the alignment 

Alignment_consensus<-msaConsensusSequence(Alignment) #msaConsensus takes the consensus of the multiple sequence alignment.
Alignment_consensus #print for visualization

#####4.	Calculate the GC content for the overall alignment  

#The function GC requires the sequence to be a character vector
Alignment_List<-strsplit(Alignment_consensus, "") #this converts the sequence to a list
Alignment_Vector<-unlist(Alignment_List) #this converts the sequence to a character vector
GC(Alignment_Vector) #this calculates GC content of the alignment

#####5.	Check to see how different your samples are from one another. 
######Are any of them different from the rest? 
#####If so, what kinds of mutations do you observe in this individual (or individuals)?
###Homo_sapiens_6 is the most different of the sequences

print(Alignment, show="complete") #this allows us to visualize the sequences and look for mutations.
#Entering the sequence into a Clustal Omega web page can also give a better visualization of the mutations (https://www.ebi.ac.uk/jdispatcher/msa/clustalo?stype=dna).
###Homo_sapiens_6 has one deletion and a few substitutions

#####6.	You suspect that an individual (or individuals) in this population might have some mutations in this gene, but you don’t know what this gene might be. 
#####Compare your sequences to a database to figure out what the gene is. 
#####You may either A) Export your data, paste it into the relevant database search engine, and add your results to a comment line in R, or 
#####B) code those same steps with an R package. What is the gene? What is the accession number of the best match to your search?

#For this, I used rBLAST to compare sequences from the nt database to the consensus sequences 
consensus_stringset<- as(Alignment_consensus, "DNAStringSet") #the sequence must be converted to DNAStringSet format to run the blastn
consensus_stringset #print for visualization
Nucleotide_Database<-blast(db="nt", remote=TRUE, type="blastn") #Creating this database variable allows us to run blastn against the nt database.
Nucleotide_Results<-predict(Nucleotide_Database,consensus_stringset) #this runs the blastn command (may take several minutes)
head(Nucleotide_Results) #see closest matching results
#The best match not only has a percent identity score of 100, but it is also the exact same length.
#Running the consensus through BLAST on the web yielded very similar results with the same exact match.
#To do this, the consensus had to be converted to a fasta file using the following line of code.
writeXStringSet(consensus_stringset, "consensus.fasta")
###The sequence is part of the human gene for the hemoglobin beta subunit.
###The accession number of the best match is LC121775.1
###This corresponds to: Homo sapiens hbb gene for beta globin, partial cds, note: HbLimassol Cd8(AAG>AAC)

#I also used blastn with the individual that is the most different:
Individual6<-sequences[6] #this isolates sequence 6 (the most different sequence)
Individual6 #print for visualization
Individual6_NucleotideResults<-predict(Nucleotide_Database,Individual6) #this runs the blastn command (may take several minutes)
head(Individual6_NucleotideResults) #see closest matching results
writeXStringSet(Individual6, "Individual6_Nucleotide.fasta") #To run blastn on a web browser, you could also convert the sequence to a fast file like so.
###The accession number of the best match is AY356351.1
###This corresponds to: Homo sapiens mutant hemoglobin beta chain (HBB) gene, partial cds

#####7.	Find the individual that is the most different from the rest of the individuals in your dataset. 
#####Translate that sequence to protein. Write it to a .fasta file.

#I translated the sequence in all 6 reading frames (3 forward and 3 in reverse) to see if the sequence has multiple translated regions that use different reading frames.
#Translating protein using reading frame 1
Individual6_Translated1<-Biostrings::translate(Individual6) #this runs the blastp command (may take several minutes)
Individual6_Translated1 #print for visualization
writeXStringSet(Individual6_Translated1, "individual6protein_frame1.fasta")
#Translating protein using reading frame 2
Individual6_Translated2<-Biostrings::translate(subseq(Individual6, start=2)) #this runs the blastp command (may take several minutes)
Individual6_Translated2 #print for visualization
writeXStringSet(Individual6_Translated2, "individual6protein_frame2.fasta")
#Translating protein using reading frame 3
Individual6_Translated3<-Biostrings::translate(subseq(Individual6, start=3)) #this runs the blastp command (may take several minutes)
Individual6_Translated3 #print for visualization
writeXStringSet(Individual6_Translated3, "individual6protein_frame3.fasta")
#Translating protein using reverse reading frame 1
Individual6_TranslatedR1<-Biostrings::translate(subseq(reverseComplement(Individual6))) #this runs the blastp command (may take several minutes)
Individual6_TranslatedR1 #print for visualization
writeXStringSet(Individual6_TranslatedR1, "individual6protein_frameR1.fasta")
#Translating protein using reverse reading frame 2
Individual6_TranslatedR2<-Biostrings::translate(subseq(reverseComplement(Individual6) , start=2)) #this runs the blastp command (may take several minutes)
Individual6_TranslatedR2 #print for visualization
writeXStringSet(Individual6_TranslatedR2, "individual6protein_frameR2.fasta")
#Translating protein using reverse reading frame 3
Individual6_TranslatedR3<-Biostrings::translate(subseq(reverseComplement(Individual6) , start=3)) #this runs the blastp command (may take several minutes)
Individual6_TranslatedR3 #print for visualization
writeXStringSet(Individual6_TranslatedR3, "individual6protein_frameR3.fasta")

#####8.	Use a database to figure out what your protein matches to. 
#####Click on the record for the best match. What is the accession number of this entry?

#I used rBLAST to compare all of the reading frames to the nr database with blastp.
#It is important to note that by using blastx with the nucleotide sequence, the same results could be obtained with less effort.
Protein_Database<-blast(db="nr", remote=TRUE, type="blastp")
Protein_ResultsIndividual6_frame1<-predict(Protein_Database,Individual6_Translated1)
head(Protein_ResultsIndividual6_frame1)
#Best Match: KAI2558340.1   100.000 percent identity   length:30   E value:2.56e-09
Protein_ResultsIndividual6_frame2<-predict(Protein_Database,Individual6_Translated2)
head(Protein_ResultsIndividual6_frame2)
#Best Match: BAU68217.1    96.721 percent identity    length:61   E value:4.64e-33
Protein_ResultsIndividual6_frame3<-predict(Protein_Database,Individual6_Translated3)
head(Protein_ResultsIndividual6_frame3)
#Best Match: KAI4069682.1    100.000 percent identity    length:81   E value:1.28e-50
Protein_ResultsIndividual6_frameR1<-predict(Protein_Database,Individual6_TranslatedR1)
head(Protein_ResultsIndividual6_frameR1)
#Best Match: KAL4829972.1    77.027 percent identity    length:75   E value:1.27e-27
Protein_ResultsIndividual6_frameR2<-predict(Protein_Database,Individual6_TranslatedR2)
head(Protein_ResultsIndividual6_frameR2)
#Best Match: OWK17646.1    70.000 percent identity    length:50   E value:5.08e-07
Protein_ResultsIndividual6_frameR3<-predict(Protein_Database,Individual6_TranslatedR3)
head(Protein_ResultsIndividual6_frameR3)
#Best Match: No significant similarity found
###Best Match of all reading frames: KAI4069682.1    100.000 percent identity    length:81   E value:1.28e-50
###Protein: hemoglobin subunit beta [Homo sapiens] 
#Note: Entering the earlier created protein fasta files into the web browser for blastp yielded the same results. 
#The same results were also seen by entering the nucleotide sequence for Individual 6 into the web browser for blastx.

#Below is a method of obtaining the same results with blastx in R
blastx_Database<-blast(db="nr", remote=TRUE, type="blastx")
Individual6_blastxResults<-predict(blastx_Database,Individual6) #this runs the blastx command (may take several minutes)
head(Individual6_blastxResults) #see closest matching results
###This also found KAI4069682.1 to be the top result (actually the 3rd from top result, but top 2 are only around 90 percent identity and are not from Homo sapiens)

#####9.	Either using R or by searching in the database, what disease(s) is this gene associated with? 
#####Does the individual that you identified in question #7 have the disease?

#I tried to do this using the romim package, but I would have needed to register for an API code to use OMIM, so I just used the web database.
###This gene is very closely associated with Sickle Cell Anemia.
#To analyze the sequences to see if individual 6 has the point mutation that causes Sickle Cell Anemia, I created a fasta file with the sequences of the standard and mutant hemoglobin along with individuals 6 and 12.
#I created this file by pasting the sequences all into the same fasta file, but this could also be done with R.
#Based on information found on PubMed and OMIM, the important point mutation occurs on the 7th codon and swaps GAG for GTG (this is actually the 6th codon in the finished protein because the methionine is spliced off).
#I then ran a multiple sequence alignment to visualize the mutations.
#We will be able to tell if individual 6 has Sickle Cell based on whether it has that mutation.
Hemoglobin_sequences<-readDNAStringSet("MutantHemoglobinMSA.fasta") #converts the fasta file to a DNA strings set variable
Hemoglobin_sequences #print for visualization
Hemoglobin_alignment<-msaMuscle(Hemoglobin_sequences, type="dna") #performs a multiple sequence alignment on the sequences using Muscle
Hemoglobin_alignment #print for visualization
print(Hemoglobin_alignment, show="complete") #this allows us to look for the mutation
#Sure enough, the standard hemoglobin and Individual 12 have GAG at the 7th codon, while the mutant hemoglobin and Individual 6 have GTG at the 7th codon.
###This tells us that Individual 6 does indeed have Sickle Cell Anemia.




