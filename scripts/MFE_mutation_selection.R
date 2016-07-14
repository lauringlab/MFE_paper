#This package contains a function to randomly select mutations to make within a nucleic acid sequence and then designs quickchange mutagenesis primers to create these mutations in an input gene. 

#Created Dec 14, 2012 by Matthew Pauly 

#This function "PrimerDesign" requires an input file of a gene sequence in a .fasta format, 
#and the number of random mutations that you would like to design primers for.  
#Options include defining maximum length of the primers,
#minimum GC content of the primers, and minimum melting temperature of the primers.
#Defaults for these are 45 bases, 40%, and 78 degrees C, respectively. (based on quickchange protocol recommendations)    

#Output of the function is a .csv file which can be opened in Excel or an equivalent program.  

#If no primer can be found matching the design criteria (Tm, GC content, etc. than you will get the output
#"No matching primer for this region" for that particular mutation.  



PrimerDesign<-function(GeneFile, numberOfMutations, maxPrimerLength=45, minGCcontent=40, minMeltTemp=78){

Gene1<-read.csv(GeneFile,sep=";", header=FALSE)	#reads the unzipped file

Gene2<-Gene1[-1,]      #removes header     

Gene3<-paste(Gene2)	#pastes the read-in gene into a list



mutSites<-round(runif(numberOfMutations,1,nchar(Gene3)))   #determines random sites to be mutated

Gene4<-strsplit(Gene3, "")  	#Splits gene into individual elements

same<-function(x){x}		#a function that does not change the input, useful for sapply


Gene5<-sapply(Gene4,same)	#Makes split genome into vectors

mutBases<-c("A","T","C","G")		#defines possible bases that can be mutated to


GeneMut<-Gene5				#creates genome in which mutations will be made


mutations<-matrix(rep(0,length(mutSites)),nrow=1)  #Table to hold the mutated bases

for (i in 1:length(mutSites)){		#Defines mutations in gene and disallows same base mutations (these would not change the gene)
	mut<-sample(mutBases,1)
	while(mut==Gene5[mutSites[i]]){
		mut<-sample(mutBases,1)}
	mutations[i]<-mut
	}

#Translating ORFs in protein sequences

trna<-function(x){
  y<-gsub("ATG","M",x)
  y<-gsub("TTT|TTC","F",y)
  y<-gsub("TTA|TTG|CTT|CTC|CTA|CTG","L",y)
  y<-gsub("ATT|ATC|ATA","I",y)
  y<-gsub("GTT|GTC|GTA|GTG","V",y)
  y<-gsub("TCT|TCC|TCA|TCG|AGT|AGC","S",y)
  y<-gsub("CCT|CCC|CCA|CCG","P",y)
  y<-gsub("ACT|ACC|ACA|ACG","T",y)
  y<-gsub("GCT|GCC|GCA|GCG","A",y)
  y<-gsub("TAT|TAC","Y",y)
  y<-gsub("TAA|TAG|TGA","stop",y)
  y<-gsub("CAT|CAC","H",y)
  y<-gsub("CAA|CAG","Q",y)
  y<-gsub("AAT|AAC","N",y)
  y<-gsub("AAA|AAG","K",y)
  y<-gsub("GAT|GAC","D",y)
  y<-gsub("GAA|GAG","E",y)
  y<-gsub("TGT|TGC","C",y)
  y<-gsub("TGG","W",y)
  y<-gsub("CGT|CGC|CGA|CGG|AGA|AGG","R",y)
  y<-gsub("GGT|GGC|GGA|GGG","G",y)
  y
}





PrimerTable<-matrix(rep(0,(length(mutSites)*9)),ncol=9)  #Creates matrix to be populated with primer info
colnames(PrimerTable)<-c("Mutation number", "Gene mutation", "Protein mutation","Forward primer", "Reverse primer", "Percent GC content","Tm","Asymmetry around mutation", "Number of bases")

for(z in 1:length(mutSites)){ 


#Finding Cs and Gs

sidedLength<-round(((maxPrimerLength-1)/2))

locCG<-gregexpr("C|G", Gene3)   #Locations of Cs and Gs in the gene

if((mutSites[z]-sidedLength)>0){primerStartSite<-(mutSites[z]-sidedLength)}
else{(primerStartSite<-1)}

if((mutSites[z]-10)>0){intermediate<-(mutSites[z]-10)}
else{intermediate<-2}

locCGF<-gregexpr("C|G", paste(Gene5[primerStartSite:intermediate], collapse=""))   #Locations of Cs and Gs in front of the mutation by 12-20 bases

locCGB<-gregexpr("C|G", paste(Gene5[(mutSites[z]+10):(mutSites[z]+sidedLength)], collapse=""))   #Locations of Cs and Gs in back of the mutation by 12-20 bases

locCGF1<-sapply(locCGF, same)  #Makes vectors of above values

locCGB1<-sapply(locCGB, same)  #Makes vectors of above values



#####Primer Design


GeneMut<-Gene5				#creates genome in which mutations will be made


GeneMut[mutSites[z]]<-mutations[z]  #Places the first mutation in


iprim<-matrix(rep(0, ((length(locCGB1)*length(locCGF1)))), ncol=1)		#initial primers

iprim2<-matrix(rep(0, ((length(locCGB1)*length(locCGF1)))*4), ncol=4)

					#For loop to create all primers possible with desired lengths
for(q in 1:length(locCGF1)){
	for(w in 1:length(locCGB1)){
		if((mutSites[z]-((sidedLength+1)-locCGF1[q]))>0){primStart<-(mutSites[z]-((sidedLength+1)-locCGF1[q]))}
		else{primStart<-1}
		iprim[(((q-1)*(length(locCGB1)))+w),1]<-paste(GeneMut[primStart:(mutSites[z]+(locCGB1[w]+9))], collapse="")
		iprim2[(((q-1)*(length(locCGB1)))+w),3]<-abs(((sidedLength+1)-locCGF1[q])-(locCGB1[w]+9))
		}
		}

for(e in 1:nrow(iprim)){
	iprim2[e,1]<-(length(sapply(gregexpr("C|G",iprim[e,1]),same))/nchar(iprim[e,1]))*100
	iprim2[e,2]<-(81.5+(0.41*(iprim2[e,1]))-(675/nchar(iprim[e,1]))-((1/nchar(iprim[e,1]))*100))
	iprim2[e,4]<-nchar(iprim[e,1])
	}




goodGC<-which(iprim2[,1]>=minGCcontent)			#Finds those primers with GC contents greater than 40

goodTm<-which(iprim2[,2]>=minMeltTemp) 	#Finds those primers with melting temps above 78

goodGCTm<-intersect(goodGC, goodTm)  #Finds primers with appropriate GC content and Tm

sortSymetry<-sort(iprim2[goodGCTm[1:length(goodGCTm)],3])  #Sorts the primers with good GC content and melt temp by symetry(numer of bases on either side of mutation)

goodSymetry<-which(iprim2[,3]==sortSymetry[1])  #selects primers that are the most symetrical, based on the most symetrical with good GC and Tm

goodGCTmSym<-intersect(goodGCTm, goodSymetry)  #Selects primers with good GC, Tm, and that are the most symetrical

sortLength<-sort(iprim2[goodGCTmSym[1:length(goodGCTmSym)],4])  #Sorts to find the shortest primer meeting the above criteria

goodLength<-which(iprim2[,4]==sortLength[1])  #Finds all primers with the shortest length determined above

bestPrimer<-intersect(goodGCTmSym, goodLength)

if(length(bestPrimer)>0) {

PrimerTable[z,4]<-iprim[bestPrimer]  #Prints forward primer to table

PrimerTable[z,6]<-iprim2[bestPrimer,1]  #Prints GC content of best primer to table

PrimerTable[z,7]<-iprim2[bestPrimer,2]  #Prints melting temp of best primer to table

PrimerTable[z,8]<-iprim2[bestPrimer,3]  #Prints asymmetry of best primer to table

PrimerTable[z,9]<-iprim2[bestPrimer,4]  #Prints length of best primer to table
}

else{
	PrimerTable[z,4]<-"No matching primer for this region"
	PrimerTable[z,5]<-"No matching primer for this region"
	PrimerTable[z,6]<-"-"
	PrimerTable[z,7]<-"-"
	PrimerTable[z,8]<-"-"
	PrimerTable[z,9]<-"-"	
	}

PrimerTable[z,2]<-sprintf("%s%s%s", Gene5[mutSites[z]], mutSites[z], GeneMut[mutSites[z]])  #prints out the gene mutation in standard base-location-base format

PrimerTable[z,1]<-z

rPrimer<-paste(rev(unlist(strsplit(iprim[bestPrimer],""))),collapse="") #Makes the reverse compliment of the primer
rPrimer<-gsub("A","t",rPrimer)
rPrimer<-gsub("T","a",rPrimer)
rPrimer<-gsub("C","g",rPrimer)
rPrimer<-gsub("G","c",rPrimer)
PrimerTable[z,5]<-revPrimer<-toupper(rPrimer)  			


#The next portion of the script is to get the ORFs to translate into protein sequence 
len1<-gregexpr("ATG(...){1,255}(TGA|TAA|TAG)",Gene3) 
length1<-attributes(gregexpr("ATG(...){1,255}(TGA|TAA|TAG)",Gene3)[[1]])
len2<-gregexpr("ATG(...){255}(...){1,255}(TGA|TAA|TAG)",Gene3)
length2<-attributes(gregexpr("ATG(...){255}(...){1,255}(TGA|TAA|TAG)",Gene3)[[1]])
len3<-gregexpr("ATG(...){255}(...){255}(...){1,255}(TGA|TAA|TAG)",Gene3)
length3<-attributes(gregexpr("ATG(...){255}(...){255}(...){1,255}(TGA|TAA|TAG)",Gene3)[[1]])
len4<-gregexpr("ATG(...){255}(...){255}(...){255}(...){1,255}(TGA|TAA|TAG)",Gene3)
length4<-attributes(gregexpr("ATG(...){255}(...){255}(...){255}(...){1,235}(TGA|TAA|TAG)",Gene3)[[1]])
 
#These next functions put all of the found start sites and lengths into vectors.  
 
 
len1a<-sapply(len1,same)
length1a<-sapply(length1[1],same)
len2a<-sapply(len2,same)
length2a<-sapply(length2[1],same)
len3a<-sapply(len3,same)
length3a<-sapply(length3[1],same)
len4a<-sapply(len4,same)
length4a<-sapply(length4[1],same)
 
start<-c(len1a,len2a,len3a,len4a)  #creates a vector for all of the start site data
len<-c(length1a,length2a,length3a,length4a)  #creates a vector containing all of the length data
 
ORF<-matrix(rep(0,(length(start)*2)),ncol=2)

ORF[,1]<-start
ORF[,2]<-len


rStart<-sort(ORF[which(ORF[,1]>0),1])  #Gets the base numbers of all ATG start sites in gene and sorts by position in gene

realStart<-which(ORF[,1]==rStart[1])  #Finds ORFs that have the earliest start site in the sequence

rLen<-sort(ORF[which(ORF[,1]==rStart[1]),2])  #sorts all of the Open reading frame lengths for the earliest start codon 

realLen<-which(ORF[,2]==rLen[1])	#Determines the ORF with the earlist start site and that is the shortest.

#This next line will create a character string with the genes ORF only
ORFgene<-paste(Gene5[ORF[realLen[1],1]:((ORF[realLen[1],1])+(ORF[realLen[1],2])-1)], collapse="")
mutORFgene<-paste(GeneMut[ORF[realLen[1],1]:((ORF[realLen[1],1])+(ORF[realLen[1],2])-1)], collapse="")


#Translating ORFs in protein sequences



#Next lines will break ORF into codons and then translate into a vector of aa sequence

  codons<-substring(ORFgene, seq(1,nchar(ORFgene)-1,3),seq(3,nchar(ORFgene),3))
  protein<-trna(codons)
  
  
mutcodons<-substring(mutORFgene, seq(1,nchar(mutORFgene)-1,3),seq(3,nchar(mutORFgene),3))
  mutprotein<-trna(mutcodons)
 
mutAAsite<-ceiling(((mutSites[z]-ORF[realLen,1]+1)/3))
 
if(mutSites[z]>ORF[realLen[1],2]){PrimerTable[z,3]<-"non-coding region"}
else{
PrimerTable[z,3]<-sprintf("%s%s%s", protein[mutAAsite], mutAAsite, mutprotein[mutAAsite])  #prints out the protein mutation in standard base-location-base format
}

}

#Saves the data frame as an CSV file which can be opened in Excel for easy viewing
print(write.table(PrimerTable, "MutagenicPrimers.csv", col.names = NA, sep = "," ))

return(PrimerTable)
}

