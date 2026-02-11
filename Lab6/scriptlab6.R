library(biostrings)
library(seqinr)
library(pwalign)
library(msa)
seqs<-readDNAStringSet("Cyphoma_CytochromeP450.fasta")
names(seqs) <- c("Cyphoma2", "Cyphoma3", "Cyphoma9", "Cyphoma8", "Cyphoma7")
Align<-msaMuscle(seqs, type="dna")
#Align<-msa(seqs, method = "Muscle") ##this should also work
Align
print(Align, show="complete")
#next 2 steps not needed (just for fun)
aln1<-pwalign::pairwiseAlignment(seqs[1], seqs[2])
pwalign::pid(aln1)
#next few steps for number of gap spaces
Align_stringset <- as(Align, "DNAStringSet")
Align_stringset
letterFrequency(Align_stringset, letters = "-")
lapply(Align_stringset,length) #shows length of each sequence
ncol(Align) #this also shows length 
#(length can also be seen by printing the msa)
#Align_consensus1<-consensusString(Align) #this isolates consensus
Align_consensus<-msaConsensusSequence(Align) #this method also works
#note: when running identical(Align_consensus,Align_consensus1), response was FALSE
Align_consensus
Align_List<-strsplit(Align_consensus, "")
Align_Vector<-unlist(Align_List)
GC(Align_Vector)
Align_seqinrformat<-msaConvert(Align, type="seqinr::alignment")
d <- dist.alignment(Align_seqinrformat, "identity")
d
