###Create a comprehensive sequence db of Hepatozoon
library("parallel")
library("DECIPHER")
library("taxonomizr")
library("clstutils")

HepSeq <- Biostrings::readDNAStringSet("~/Phylo_Hepatozoon/data/ena.fasta")

## Clean sequences: 1) In an interval of lenght, 2) With gaps and 3) With bad characters. 
removeBadSeq <- function (DNAStringSet, min.length=0, max.length=0){
  DNAStringSet <- RemoveGaps(DNAStringSet)
  DNAStringSet <- DNAStringSet[width(DNAStringSet)>min.length]
  DNAStringSet <- DNAStringSet[width(DNAStringSet)<max.length]
  goodSeq <- dada2:::C_isACGT(as.character(DNAStringSet))
  DNAStringSet[goodSeq]
}

HepSeq <- removeBadSeq(HepSeq, 500, 700)

## rename with the accessions
names(HepSeq) <- gsub("ENA\\|.*?\\|(.*?\\.\\d+).*", "\\1", names(HepSeq))

## exclude duplicates
HepSeq <- HepSeq[!duplicated(names(HepSeq))]

tax.Id.seq <- accessionToTaxa(names(HepSeq), "/SAN/db/taxonomy/taxonomizr.sql")

## exclude Not Hepatozoon canis sequences
if(!exists("badTaxa")){
  source("~/Phylo_Hepatozoon/R/1_Metadata_sequences.R") ##   
}

badTaxaSeq <- names(HepSeq[tax.Id.seq%in%badTaxa])

HepSeqFinal<- HepSeq[!(names(HepSeq)%in%badTaxaSeq)]

writeXStringSet(HepSeqFinal, "~/Phylo_Hepatozoon/output/database/Hepatozoon_canis_final_Seq.fasta")

##Correct metadata generated before 
final_accession<- names(HepSeqFinal)

##Remove ".1" from the accession number to match with metadata 
final_accession<- gsub("\\..*", "", final_accession)
final_accession<- as.data.frame(final_accession)
colnames(final_accession)[1]<- "accession"

HepDBFinal<- plyr::join(final_accession, db.Hcanis, by= "accession")

write.csv(HepDBFinal, "~/Phylo_Hepatozoon/output/database/Metadata_H_canis_seq.csv")
