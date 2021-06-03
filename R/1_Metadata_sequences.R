###Create a comprehensive metadata of Hepatozoon

library("xml2")
library("rvest")
library("lifecycle", lib.loc="/usr/local/lib/R/site-library") 
library("dplyr")
library("Biostrings")
library("DECIPHER")
library("annotate")
library("taxonomizr")

##Database downloaded from ENA search
##https://www.ebi.ac.uk/ena/data/search?query=Hepatozoon+canis ##Accession date 05.06.2020

db <- read_xml("data/ena.xml") ##Read metadata extracted from ENA

# get all the entries
entries <- xml_children(db)

## Transform the nodeset into an R readable dataframe 
db.raw<- bind_rows(lapply(xml_attrs(entries), function(x) data.frame(as.list(x), stringsAsFactors=FALSE)))

##Get taxonomy 
##convert accession numbers to taxonomic IDs 
##SQLite database is already in Harriet
#taxIDs<- accessionToTaxa(as.character(db.raw$accession), sqlFile = "/SAN/db/taxonomy/taxonomizr.sql") #Not working :S
#accession_18S$Taxonomic_ID<- taxID_18S
##This method is not working... Let's do it in bash
##First save the accessions as a data frame
acces<- as.data.frame(db.raw$accession)
write.table(acces, file = "data/accession.txt", sep = "\t",
            row.names = F, col.names = F, quote = F)

###!/bin/bash
## cat ~/Phylo_Hepatozoon/data/accession.txt | epost -db nuccore | esummary -db nuccore | xtract -pattern DocumentSummary -element Caption,TaxId

tax.Id<- read.delim("data/taxID.txt", header = F, sep = "\t")
colnames(tax.Id)<- c("accession", "taxID") ### After taxID search just 1134 left... find out where are the 36 missing.

##get taxonomy using taxonomic IDs 
taxmy<- taxonomizr::getTaxonomy(as.character(tax.Id$taxID), sqlFile = "/SAN/db/taxonomy/taxonomizr.sql")
##Transform it to a database
taxmy<- as.data.frame(taxmy)
rownames(taxmy)<- c(1:nrow(taxmy))
 
###Merge tax.id and taxmy
taxmy <- merge(tax.Id, taxmy, by=0, all=TRUE)
###Join all the information
db.all<- plyr::join(db.raw, taxmy, by= "accession")
rm(acces, db.raw, taxmy, tax.Id)

###Let's clean the database
##Eliminate all the entries without complete information
db.compl<- db.all[complete.cases(db.all), ]

write.csv(db.compl, "output/database/Database_not_cleaned.csv", row.names = F)

##Summary of sequences
require(htmlTable)
seq.summary<- as.data.frame(summary(db.compl$species))
seq.summary<- htmlTable(seq.summary, 
                          header =  "Number of sequences",
                          caption="Summary of sequences by species in non-cleaned database")

write.table(seq.summary, "output/database/Summary_not_cleaned.html")

##Ok out of the total amount of sequences just 970 are assigned as H. canis
##Let's eliminate sequences that are not true H. canis and repeated 
db.compl%>%
  dplyr::filter(species == "Hepatozoon canis") %>%
  dplyr::distinct()-> db.Hcanis

##Get a vector with taxids not from H. canis 
db.compl%>%
  dplyr::filter(species != "Hepatozoon canis") %>%
  dplyr::distinct()-> db.NOT.Hcanis

badTaxa<- db.NOT.Hcanis$taxID

hist(as.numeric(db.Hcanis$sequenceLength), 
     xlim= c(min(as.numeric(db.Hcanis$sequenceLength)), 35000), 
     breaks= 150, main = "Frequency of sequences base on length", xlab = "Sequence length (bp)")

summary(as.numeric(db.Hcanis$sequenceLength))
