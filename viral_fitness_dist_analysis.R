# LOAD PACKAGES
library(moments)
library(nortest)
library(lattice)
library(qpcR)
library(plyr)

#this code is for comparing distributions of MFE for different viruses.
# written by Shawn Whitefield 2-2-2016 

setwd("~/Desktop/viral_fitness_dist_analysis/")


# READ IN FITNESS DATA FOR EACH VIRUS
  
  #polio
  polio_fitness_df<-read.csv("parsed_fitness_files/polio_fig_4_only.csv")
  polio_fitness_df$virus <-"polio"
  polio_fitness_df$W <-polio_fitness_df$Fitness

  #TEV
  TEV_fitness_df<-read.csv("parsed_fitness_files/TEV_parsing.csv")
  TEV_fitness_df$W<-sub("\\(.*$","",TEV_fitness_df$Relative_fitness) #add column for W
  TEV_fitness_df$W<-as.numeric(TEV_fitness_df$W)
  TEV_fitness_df$virus<-"TEV"
  
  #QB 
  #CHECK MATH
  QB_fitness_df<-read.csv("parsed_fitness_files/QB_parsed.csv")
  QB_fitness_df$virus <- "QB"
  #add column for W,= 1-log2(1-s) ?
  QB_fitness_df$one_minus_rf<-1-QB_fitness_df$relative_fitness
  QB_fitness_df$log_2_one_minus_rf<-log2(QB_fitness_df$one_minus_rf)
  QB_fitness_df$one_minus_log2_W<-1-QB_fitness_df$log_2_one_minus_rf
  QB_fitness_df$W<-QB_fitness_df$one_minus_log2_W

  #phiX174 
  phiX174_fitness_df<-read.csv("parsed_fitness_files/phi_x_174_parsed.csv")
  #add column for W, = 1-s
  phiX174_fitness_df$virus<-"phix174"
  phiX174_fitness_df$one_minus_rf<-1-phiX174_fitness_df$relative_fitness
  phiX174_fitness_df$log_2_one_minus_rf<-log2(phiX174_fitness_df$one_minus_rf)
  phiX174_fitness_df$one_minus_log2_W<-1-phiX174_fitness_df$log_2_one_minus_rf
  phiX174_fitness_df$W<-phiX174_fitness_df$one_minus_log2_W


  #f1
    #in progress parsing the data

  # flu
    #waiting on flu data
#-------------------------------------------------------------------------
# PLOT HISTOGRAM OF FITNESS FOR EACH VIRUS

  #shorter variable names
    phi<-phiX174_fitness_df$one_minus_log2_W
    QB<-QB_fitness_df$one_minus_log2_W
    TEV<-TEV_fitness_df$W
    polio<-polio_fitness_df$W

  #multiple histograms on same plot
    histogram(phi, col="red", xlim=c(0,max(c(phi,QB,TEV,polio))),freq=FALSE, xlab=list('PhiX174',cex=2), ylim=c(0,60), main=list(cex=3),breaks=15)
    histogram(QB, col="blue", xlim=c(0,max(c(phi,QB,TEV,polio))),freq=FALSE, xlab='Qβ',ylim=c(0,60),breaks=15)
    histogram(TEV, col="darkgreen",xlim=c(0,max(c(phi,QB,TEV,polio))),freq=FALSE, xlab='TEV',ylim=c(0,60),breaks=15)
    histogram(polio, col="grey",xlim=c(0,max(c(phi,QB,TEV,polio))),freq=FALSE, xlab='Polio',ylim=c(0,60),breaks=15)

  #make a dataframe of all of the Ws
  big_fitness_df<-c(phiX174_fitness_df$virus,QB_fitness_df$virus, polio_fitness_df$virus,TEV_fitness_df$virus)
  big_fitness_df<-as.data.frame(big_fitness_df)
  big_fitness_df$W<-c(phiX174_fitness_df$W,QB_fitness_df$W, polio_fitness_df$W,TEV_fitness_df$W)

  hist(big_fitness_df$W)

#-------------------------------------------------------------------------

# EXAMINE SKEWNESS
  skewness(phi)
  skewness(QB)
  skewness(TEV)
  skewness(polio)
  skewness(big_fitness_df$W) #all viruses together

# EXAMINE KURTOSIS
  kurtosis(phi)
  kurtosis(QB)
  kurtosis(TEV)
  kurtosis(polio)
  kurtosis(big_fitness_df$W) #all viruses together

#-------------------------------------------------------------------------

# PERFORM 2 SAMPLE KS TESTS 
  #perform all pairwise KS tests

#make dataframe of all fitness values with each virus as a column  
v <- list(phi,QB,TEV,polio)
n <- max(sapply(v, length))
all_fitness_df<-do.call(rbind, lapply(v, `[`, seq_len(n)))
rownames(all_fitness_df)<-c("phi","QB","TEV","polio")
all_fitness_df<-t(all_fitness_df)

#function to get all of the pairwise comparison pvalues in a martix
make_ks_results_table<-function(v1,v2,...){
  pvalue<-ks.test(v1,v2,...,)$p.value
  return(pvalue)
}

#use apply to get the pvalues

pairwise_ks_pvalue_mat<-apply(all_fitness_df, 2, function(x) apply(all_fitness_df, 2, function(y) make_ks_results_table(x, y))) 
write.xlsx(pairwise_ks_pvalue_mat,"./pairwise_ks_pvalue.xlsx")


#function to get all of the pairwise comparison the test statistic in a martix
make_ks_stat_table<-function(v1,v2,...){
  statistic<-ks.test(v1,v2,...,)$statistic
  return(statistic)
}

#use apply to get the test statistic, D
pairwise_ks_statistic_mat<-apply(all_fitness_df, 2, function(x) apply(all_fitness_df, 2, function(y) make_ks_stat_table(x, y))) 
write.xlsx(pairwise_ks_statistic_mat,"./pairwise_ks_statistic.xlsx")
# find the maximum test statistic ie biggest difference (most conservative)

max_ks_D<-max(pairwise_ks_statistic_mat)
                   
#randomly put all fitness values into 5 groups (or 6 when flu is here) of size n for each sample

make_ks_cdf<-function(all_W_list,phi,QB,TEV,polio){
  all_W_list<-as.numeric(all_fitness_df, na.rm=T) #make big list of all W values from which to sample
  phi_sample<-sample(all_W_list, size=length(phi), replace=T) #sample from the list of all Ws with replacement
  QB_sample<-sample(all_W_list, size=length(QB),replace=T)
  TEV_sample<-sample(all_W_list, size=length(TEV),replace=T)
  polio_sample<-sample(all_W_list, size=length(polio),replace=T)
}


#

  # find the largest KS distance again x 1000 or 2000 ish times
    # use replicate function

  # chect to see if the KS value for the original data is in the tail




