# MFE distribution analysis
Shawn Whitefield  
April 28, 2016  

##This is an R Markdown document for the MFE distribution analysis.

##### load in packages for analysis

```r
# LOAD PACKAGES
library(moments)
library(nortest)
library(lattice)
library(qpcR)
```

```
## Loading required package: MASS
```

```
## Loading required package: minpack.lm
```

```
## Loading required package: rgl
```

```
## Warning in rgl.init(initValue, onlyNULL): RGL: unable to open X11 display
```

```
## Warning: 'rgl_init' failed, running with rgl.useNULL = TRUE
```

```
## Loading required package: robustbase
```

```
## Loading required package: Matrix
```

```r
library(plyr)
library(ggplot2)
```

##### read in the parsed data files as dataframes


```r
setwd("~/Desktop/viral_fitness_dist_analysis/")

# READ IN FITNESS DATA FOR EACH VIRUS
  #note!!
  #Flu is in terms of fitness
  # S = fitness -1 , or Fitness = S + 1 
  
  #polio
    #already represented as fitness
  polio_fitness_df<-read.csv("parsed_fitness_files/polio_fig_4_only.csv")
  polio_fitness_df$virus <-"polio"
  polio_fitness_df$fitness <-polio_fitness_df$Fitness

  #these viruses are in therms of S so need to add 1
  #TEV
  TEV_fitness_df<-read.csv("parsed_fitness_files/TEV_parsing.csv")
  TEV_fitness_df$s<-sub("\\(.*$","",TEV_fitness_df$Relative_fitness) #add column for fitness
  TEV_fitness_df$fitness<-(as.numeric(TEV_fitness_df$s))+1
  TEV_fitness_df$virus<-"TEV"
  
  #QB 
  QB_fitness_df<-read.csv("parsed_fitness_files/QB_parsed.csv")
  QB_fitness_df$virus <- "QB"
  QB_fitness_df$fitness<-1+(QB_fitness_df$relative_fitness)

  #phiX174 
  phiX174_fitness_df<-read.csv("parsed_fitness_files/phi_x_174_parsed.csv")
  phiX174_fitness_df$virus<-"phix174"
  phiX174_fitness_df$fitness<-1+(phiX174_fitness_df$relative_fitness)

  # Flu 
    #there are 6 datasets here
    # +/- lethal fraction for:
      #1. All
      #2. HA + NA
      #3  External segments
  flu_fitness_df<- read.csv("parsed_fitness_files/flu.csv")
```

#####Examine Skew and Kurtosis in non-lethal fraction

```r
  #shorter variable names and remove lethal fractions
    phiX<-phiX174_fitness_df$fitness
    names(phiX)<-c(rep("phiX",length(phiX))) #name the vector
      phiX.CDF<-phiX[phiX!=0]
      
    QB<-QB_fitness_df$fitness
    names(QB)<-c(rep("QB",length(QB)))
      QB.CDF<-QB[QB!=0]

    TEV<-TEV_fitness_df$fitness
    names(TEV)<-c(rep("TEV",length(TEV)))
      TEV.CDF<-TEV[TEV!=0]

    polio<-polio_fitness_df$fitness
    names(polio)<-c(rep("polio",length(polio)))
      polio.CDF<-polio[polio!=0]

    # Without lethal fractions
    flu_total.CDF<-flu_fitness_df$Total.CDF #all
    flu_total.CDF<-subset(flu_total.CDF, flu_total.CDF > -10) # remove NAs by subsetting w/small #
    names(flu_total.CDF)<-c(rep("flu_total_CDF",length(flu_total.CDF)))

    flu_random.CDF<-flu_fitness_df$Random.CDF #this is all w/o lethal
    flu_random.CDF<-subset(flu_random.CDF, flu_random.CDF > -10)
    names(flu_random.CDF)<-c(rep("flu_random.CDF",length(flu_random.CDF)))

    flu_surface.CDF<-flu_fitness_df$Surface.CDF #HA + NA
    flu_surface.CDF<-subset(flu_surface.CDF, flu_surface.CDF > -10)
    names(flu_surface.CDF)<-c(rep("flu_surface_CDF",length(flu_surface.CDF)))

    flu_internal.CDF<-flu_fitness_df$Internal.CDF #internal
    flu_internal.CDF<-subset(flu_internal.CDF, flu_internal.CDF > -10)
    names(flu_internal.CDF)<-c(rep("flu_internal.CDF",length(flu_internal.CDF)))


# EXAMINE SKEWNESS
  skewness(phiX.CDF)
```

```
## [1] -1.873701
```

```r
  skewness(QB.CDF)
```

```
## [1] -1.110005
```

```r
  skewness(TEV.CDF)
```

```
## [1] 1.815592
```

```r
  skewness(polio.CDF)
```

```
## [1] -0.3176839
```

```r
  skewness(flu_total.CDF)
```

```
## [1] -1.170211
```

```r
  skewness(flu_random.CDF)
```

```
## [1] -1.120054
```

```r
  skewness(flu_surface.CDF)
```

```
## [1] -0.8869894
```

```r
  skewness(flu_internal.CDF)
```

```
## [1] -0.9867137
```

```r
# EXAMINE KURTOSIS
  kurtosis(phiX.CDF)
```

```
## [1] 6.321138
```

```r
  kurtosis(QB.CDF)
```

```
## [1] 3.011725
```

```r
  kurtosis(TEV.CDF)
```

```
## [1] 6.51863
```

```r
  kurtosis(polio.CDF)
```

```
## [1] 3.461289
```

```r
  kurtosis(flu_total.CDF)
```

```
## [1] 3.901619
```

```r
  kurtosis(flu_random.CDF)
```

```
## [1] 3.479086
```

```r
  kurtosis(flu_surface.CDF)
```

```
## [1] 3.803445
```

```r
  kurtosis(flu_internal.CDF)
```

```
## [1] 2.980413
```

```r
# MEAN
  mean(phiX.CDF)
```

```
## [1] 0.8739444
```

```r
  mean(QB.CDF)
```

```
## [1] 0.8972333
```

```r
  mean(TEV.CDF)
```

```
## [1] 1.513197
```

```r
  mean(polio.CDF)
```

```
## [1] 0.7467637
```

```r
  mean(flu_total.CDF)
```

```
## [1] 0.8268889
```

```r
  mean(flu_random.CDF)
```

```
## [1] 0.7979688
```

```r
  mean(flu_surface.CDF)
```

```
## [1] 0.8821951
```

```r
  mean(flu_internal.CDF)
```

```
## [1] 0.7806122
```

##### Histograms


```r
#make df with all of the CDF fitness

Big_fitness_df<-as.data.frame(c(phiX.CDF, QB.CDF,TEV.CDF, polio.CDF, flu_total.CDF, flu_random.CDF, flu_surface.CDF, flu_internal.CDF))
colnames(Big_fitness_df)<-"fitness"
Big_fitness_df$virus<-names(c(phiX.CDF, QB.CDF,TEV.CDF, polio.CDF, flu_total.CDF, flu_random.CDF, flu_surface.CDF, flu_internal.CDF))

#double check that these are only the non-lethals
sum(Big_fitness_df$fitness > -1)
```

```
## [1] 6920
```

```r
length(Big_fitness_df$fitness)
```

```
## [1] 6920
```

```r
  #yup! they are

#MAKE HISTOGRAM that we don't need
all_hist<-(ggplot(Big_fitness_df, aes(fitness, ..density.., colour = virus)) +
  geom_freqpoly(binwidth = .1)+
  theme_classic())
plot(all_hist)
```

![](./MFE_distribution_analysis_files/figure-html/unnamed-chunk-5-1.png)<!-- -->

```r
# Flu only histogram that we also don't need
flu_fitness<-as.data.frame(c(flu_total.CDF,flu_random.CDF,flu_surface.CDF,flu_internal.CDF))
colnames(flu_fitness)<-"fitness"
flu_fitness$virus<-names(c(flu_total.CDF,flu_random.CDF,flu_surface.CDF,flu_internal.CDF))

flu_hist<-(ggplot(flu_fitness, aes(fitness, ..density.., colour = virus)) +
  geom_freqpoly(binwidth = .1)+
  theme_classic())
plot(flu_hist)
```

![](./MFE_distribution_analysis_files/figure-html/unnamed-chunk-5-2.png)<!-- -->


