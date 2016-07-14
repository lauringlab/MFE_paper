#this analysis is for verifying mutations for the MFE paper
library(plyr)
#setwd("~/Desktop/PIBS/Lauring lab/2016-07-08_MFE_reviewer_comments/")

#verify The nucleotide mutation matches the clone ID

shawns_mutants<-read.table("shawn_master_database_mutant_and_genes.csv", sep=',', header=1)
table_1_mutants<-read.table("table_1_mutants.csv", sep=',', header=1)
table_3_mutants<-read.table("table_3_mutants.csv", sep=',', header=1)

#fix formatting
shawns_mutants$Mutant<-gsub("_","-",shawns_mutants$Mutant)

#make string for matching mutant and mutations betweent ables

shawns_mutants$matching<-paste(shawns_mutants$Mutant,shawns_mutants$gene.mutation,sep="_")
table_1_mutants$matching<-paste(table_1_mutants$Clon_.ID,table_1_mutants$Mutation,sep="_")
table_3_mutants$matching<-paste(table_3_mutants$Clone.ID,table_3_mutants$Mutation,sep="_")

#verify they match

sum(table_1_mutants$matching %in% shawns_mutants$matching)
#=38 which is the length of the # of mutants in table 1, so these are all correct

sum(table_3_mutants$matching %in% shawns_mutants$matching)
#sum is 90, meaning 1 does not match
table_3_mutants[!table_3_mutants$matching %in% shawns_mutants$matching,]
  # 74     NA-3    G758C NA-3_G758C in table 3
  #       NA-3  G752C   NA-3_G752C in my master database
#this one was found to be incorrect in the master database. the true mutation is NA-3    G758C NA-3_G758C 

