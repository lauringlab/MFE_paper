
#USE find_prevalence_column.py to find the substitutions that exist in fludb or prints none if they arent there.
#read in this column and cbind it to our df so we can calculate prevalence
prev_col<-read.table("../data/prevalence_col.txt",sep='\n',header=F)
mut_prevalence_df<-read.table("../data/mut_prevalence_df_inprogess.txt",sep="\t",header = T)

mut_prevalence_df$occurence<-prev_col
mut_prevalence_df$V1<-NULL
#get only number of occurences in fludb
mut_prevalence_df$occurence<-gsub('[[:alpha:]]', '', mut_prevalence_df$occurence)
mut_prevalence_df$occurence<-gsub('=', '', mut_prevalence_df$occurence)
mut_prevalence_df$occurence<-gsub(',', '', mut_prevalence_df$occurence)
mut_prevalence_df$occurence<-gsub(' ', '', mut_prevalence_df$occurence)
mut_prevalence_df$occurence<-as.numeric(mut_prevalence_df$occurence)

#prevalence %
mut_prevalence_df$percent_prevalence_in_fludb<-(mut_prevalence_df$occurence / mut_prevalence_df$number_of_sequences_with_position_in_fludb)*100

write.csv(mut_prevalence_df, "../results/mutant_prevalence_df.csv", quote=F, row.names =F)


#correlation test
mut_prevalence_df$Fitness<-as.vector(mut_prevalence_df$Fitness)
mut_prevalence_df$Fitness[mut_prevalence_df$Fitness =="lethal"]<-0
mut_prevalence_df$Fitness<-as.numeric(mut_prevalence_df$Fitness)
#cor.test(mut_prevalence_df$Fitness, mut_prevalence_df$percent_prevalence_in_fludb, method = "spearman", na.rm = T)
