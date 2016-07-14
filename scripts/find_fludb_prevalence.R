
#This script finds the prevalence of mutaitons in our library that are in fludb

# FIND PREVALENCE OF MUTANTS IN FLUDB
#read in all of the tables from fludb
#add colnames
#delete extra column
#PB2
PB2_df<-read.table("PB2_fludb.txt",  sep = '\t')
PB2_df[,6]<-NULL
colnames(PB2_df)<-c("AA_position","consensus","score","details","number_of_seqs")

#PB1
PB1_df<-read.table("all_mutant_alignments/PB1_fludb.txt",  sep = '\t')
PB1_df[,7]<-NULL
PB1_df<-PB1_df[,-2]
colnames(PB1_df)<-c("AA_position","consensus","score","details","number_of_seqs")

#PB1_F2
PB1_F2_df<-read.table("all_mutant_alignments/PB1_F2_fludb.txt",  sep = '\t')
PB1_F2_df[,6]<-NULL
colnames(PB1_F2_df)<-c("AA_position","consensus","score","details","number_of_seqs")
#PA
PA_df<-read.table("all_mutant_alignments/PA_fludb.txt",  sep = '\t')
PA_df[,6]<-NULL
colnames(PA_df)<-c("AA_position","consensus","score","details","number_of_seqs")
#HA
HA_df<-read.table("all_mutant_alignments/HA_fludb.txt",  sep = '\t')
HA_df[,6]<-NULL
colnames(HA_df)<-c("AA_position","consensus","score","details","number_of_seqs")
#NP
NP_df<-read.table("all_mutant_alignments/NP_fludb.txt",  sep = '\t')
NP_df[,6]<-NULL
colnames(NP_df)<-c("AA_position","consensus","score","details","number_of_seqs")
#NA
NA_df<-read.table("all_mutant_alignments/NA_fludb.txt",  sep = '\t')
NA_df[,6]<-NULL
colnames(NA_df)<-c("AA_position","consensus","score","details","number_of_seqs")

#M1
M1_df<-read.table("all_mutant_alignments/M1_fludb.txt",  sep = '\t')
M1_df[,6]<-NULL
colnames(M1_df)<-c("AA_position","consensus","score","details","number_of_seqs")
#M2
M2_df<-read.table("all_mutant_alignments/M2_fludb.txt",  sep = '\t')
M2_df[,6]<-NULL
colnames(M2_df)<-c("AA_position","consensus","score","details","number_of_seqs")
#NS1
NS1_df<-read.table("all_mutant_alignments/NS1_fludb.txt",  sep = '\t')
NS1_df[,6]<-NULL
colnames(NS1_df)<-c("AA_position","consensus","score","details","number_of_seqs")
#NEP
NEP_df<-read.table("all_mutant_alignments/NS2_fludb.txt",  sep = '\t')
NEP_df[,6]<-NULL
colnames(NEP_df)<-c("AA_position","consensus","score","details","number_of_seqs")

#all AA substitution positions
PB2_pos = c(271, 511, 93, 138, 169, 285, 380, 408, 490, 500, 545, 609, 696 ) # STOP REVERSION mutation not included because not in alignment
PB1_pos<-c(76, 166, 175, 186, 192, 422, 433, 580, 721, 751, 217, 235, 367, 403, 415, 519)
PB1_F2_pos<-c(45)
PA_pos<-c(22, 23, 71, 72, 159, 285, 314, 339, 445, 554, 700, 80, 169, 434, 537, 553)
#wsn HA has a 1 base deletion in subject relative to fludb consensus at 147, so positions in fludb are our position+1 after 146.
HA_pos<-c( 8, 67 ,101, 112, 126, 171, 207, 224, 252, 257, 265, 296 ,304 ,326 ,341 ,343, 370, 372, 375, 382 ,400, 412, 421, 424, 439, 479, 495, 518, 520, 524 ,556)
NP_pos<-c(51, 127, 131, 137, 147, 372, 395, 480, 36, 159, 383, 487)
#NA has 16 base deletion in subject at 57, add 16 to our #s after 56 to get fludb numbering
NA_pos<-c(27,  30 , 47 , 53 , 77 ,121, 151, 160 ,161, 233, 241, 243, 263, 313, 333, 352 ,357 ,367, 387, 397, 399 ,448 ,457 ,462)
M1_pos<-c(161, 176, 50, 172, 212)
M2_pos<-c(50, 75)
NS1_pos<-c(9,67, 208) # NOTE THAT 208 EXISTS TWICE, BUT ONLY 1 LINE FOR IT, E208STOP, E208D, deleted E208D from excel file for merging
NEP_pos<-c(9,51,104,50) 

#subset all fludb df based on the AA positions in the position vectors

PB2_muts<-PB2_df[PB2_df$AA_position %in% PB2_pos,]
PB1_muts<-PB1_df[PB1_df$AA_position %in% PB1_pos,]
PB1_F2_muts<-PB1_F2_df[PB1_F2_df$AA_position %in% PB1_F2_pos,]
PA_muts<-PA_df[PA_df$AA_position %in% PA_pos,]
HA_muts<-HA_df[HA_df$AA_position %in% HA_pos,]
NP_muts<-NP_df[NP_df$AA_position %in% NP_pos,]
NA_muts<-NA_df[NA_df$AA_position %in% NA_pos,]
M1_muts<-M1_df[M1_df$AA_position %in% M1_pos,]
M2_muts<-M2_df[M2_df$AA_position %in% M2_pos,]
NS1_muts<-NS1_df[NS1_df$AA_position %in% NS1_pos,]
NEP_muts<-NEP_df[NEP_df$AA_position %in% NEP_pos,]

#make 1 big df
PB2_muts$segment<-"PB2"
PB1_muts$segment<-"PB1"
PB1_F2_muts$segment<-"PB1_F2"
PA_muts$segment<-"PA"
HA_muts$segment<-"HA"
NP_muts$segment<-"NP"
NA_muts$segment<-"NA_seg"
M1_muts$segment<-"M1"
M2_muts$segment<-"M2"
NS1_muts$segment<-"NS1"
NEP_muts$segment<-"NEP"

#merge dfs
all_muts_df<-rbind(PB2_muts,PB1_muts,PB1_F2_muts,PA_muts,HA_muts,NP_muts,NA_muts,M1_muts,M2_muts,NS1_muts,NEP_muts)
all_muts_df$score<-NULL
all_muts_df$AA_position<-gsub(" ", "", all_muts_df$AA_position)
all_muts_df$AA_position<-as.numeric(all_muts_df$AA_position)
all_muts_df<-all_muts_df[with(all_muts_df, order(segment,AA_position)), ] #sort by segment and position

# read in mutant list 
our_subs_df<-read.table("all_mutant_alignments/mutant_excel_df.csv", sep = ",", header=T)
# get position only column for sorting
our_subs_df$position<-gsub('[[:alpha:]]', '', our_subs_df$Amino.Acid.Change)
our_subs_df$position<-gsub(" ", "", our_subs_df$position)
our_subs_df$position<-as.numeric(our_subs_df$position)
our_subs_df<-our_subs_df[with(our_subs_df, order(Segment, position)), ] #sort by segment and position

#combine dataframes
mut_prevalence_df<-cbind(all_muts_df,our_subs_df)

#get variable for amino acid change
mut_prevalence_df$substitution<-gsub('[[:digit:]]', '', mut_prevalence_df$Amino.Acid.Change)
mut_prevalence_df$substitution<-gsub(" ","",mut_prevalence_df$substitution)
mut_prevalence_df$substitution<-substring(mut_prevalence_df$substitution,2)

#need lookup table to convert from single letter code to abbreviation from fludb
#read in codon vs single letter vs abbreviation table copied and pasted from wikipedia
codon_table<-read.table("all_mutant_alignments/codon_table.csv", sep = ",", header=T)
colnames(codon_table)<-c("amino_acid","substitution","abbrev")
mut_prevalence_df <-join(mut_prevalence_df,codon_table, by="substitution")

#remove redundant columns
mut_prevalence_df$Segment<-NULL #remove dup column
colnames(mut_prevalence_df)<-c("fludb_H1N1_AA_position","fludb_consensus","fludb_substitutions",
                               "number_of_sequences_with_position_in_fludb","segment","Amino_acid_change","Fitness",
                               "WSN33_mutant_backbone_position","substitution","Amino_acid","Amino_acid_abbreviation")
mut_prevalence_df$Amino_acid<-NULL
mut_prevalence_df$substitution<-NULL

#make appreviation a vector
mut_prevalence_df$Amino_acid_abbreviation<-as.vector(mut_prevalence_df$Amino_acid_abbreviation)


#get prevalence
# write to table and use python script to extract prevalence 
write.table(mut_prevalence_df,"all_mutant_alignments/2016-07-10_mut_prevalence_df.txt",sep="\t", quote=F, row.names =F)

#USE find_prevalence_column.py to find the substitutions that exist in fludb or prints none if they arent there.
#read in this column and cbind it to our df so we can calculate prevalence
prev_col<-read.table("all_mutant_alignments/prevalence_col.txt",sep='\n',header=F)

mut_prevalence_df<-cbind(mut_prevalence_df,prev_col)
mut_prevalence_df$occurence<-mut_prevalence_df$V1
mut_prevalence_df$V1<-NULL
#get only number of occurences in fludb
mut_prevalence_df$occurence<-gsub('[[:alpha:]]', '', mut_prevalence_df$occurence)
mut_prevalence_df$occurence<-gsub('=', '', mut_prevalence_df$occurence)
mut_prevalence_df$occurence<-gsub(',', '', mut_prevalence_df$occurence)
mut_prevalence_df$occurence<-gsub(' ', '', mut_prevalence_df$occurence)
mut_prevalence_df$occurence<-as.numeric(mut_prevalence_df$occurence)

#prevalence %
mut_prevalence_df$percent_prevalence_in_fludb<-(mut_prevalence_df$occurence / mut_prevalence_df$number_of_sequences_with_position_in_fludb)*100

write.table(mut_prevalence_df, "~/Desktop/2016-07-10_mutant_prevalence_df", quote=F, row.names =F, sep = '\t')


#correlation test
mut_prevalence_df$Fitness<-as.vector(mut_prevalence_df$Fitness)
mut_prevalence_df$Fitness[mut_prevalence_df$Fitness =="lethal"]<-0
mut_prevalence_df$Fitness<-as.numeric(mut_prevalence_df$Fitness)
cor.test(mut_prevalence_df$Fitness, mut_prevalence_df$percent_prevalence_in_fludb, method = "spearman", na.rm = T)
