#! /usr/bin/env Rscript
#########################################################
######### LaSSO - Lariat Sequence Site Origin ###########
######### Author: Danny Asher Bitton UCL      ###########
#########	d.bitton@ucl.ac.uk	      ###########	
#########################################################

args <- commandArgs(TRUE); if("--help" %in% args){cat("

Name:

LaSSO - Lariat Sequence Site Origin 

Description:

An R script that creates a FASTA database containing all possible lariat signatures from a given set of introns.
It requires two text files, as follows:

	1) Sequence file containing each intron on a separate line (5'-3' direction; see also SampleSeqFile.txt)
	2) Header file containing the headers describing each intron on separate lines (same order as sequence file; see also SampleHeaderFile.txt)
	
	The Sequence file should be in the following format: 
	
	GTGCATCGTTTGCTAACTGTGTTTACTGCTCAAGATGCTCAG
	GTACGTTTCCTTCACACTGAAGCTGTTTTTTCTCTTTACTAATAGTTTATTGTAG
	GTATGTTTCCTTTATAGTGATGCTTTTTTTCTTTTTTTTTGCTAATAATTTAAAATAG
	GTAAGAACATTCAAGAAAGCGTTTGATATGAAGTCTGTTCCATTATTGCTCTTAGTGCCTGACAAACGTTGTTAGTCTGTTTTGAGTTTACAAGTTTGTTAAGCTAATAGCCATTCATTGATCTTCATTGTTACATGCTTTCCATGGAATAGTTGGAGGTGTTTACTTTTCAACTTGTTAGCACTCACTTGCTAACATTTTTAATAACAG
	
	The header file should be in the following format:(Intron_id;Chromsome;Strand;Start;Stop;5PrimeFlankingExon - 3PrimeFlankingExon )
	
	SPAC212.06c;I;1;18307;18348;SPAC212.06c.1:exon:1-SPAC212.06c.1:exon:2
	SPAC212.04c;I;1;22077;22131;SPAC212.04c.1:exon:1-SPAC212.04c.1:exon:2
	SPAC212.01c;I;1;29228;29285;SPAC212.01c.1:exon:1-SPAC212.01c.1:exon:2
	SPAC977.18;I;-1;31558;31767;SPAC977.18.1:exon:3-SPAC977.18.1:exon:2
	

LaSSO also requires the read length, branch-point letters and a database name to be specified
	    

Usage:

LaSSO --SeqFile=<SampleSeqFile.txt> --HeaderFile=<SampleHeaderFile.txt> --ReadLength=<51> --BPS=<A,T,G,C> --DatabaseName=<LariatSampleDatabase.fa>  
LaSSO --help

Arguments:

--SeqFile:

The input sequence file.

--HeaderFile:

The input header file.

--ReadLength:

The input read length.

--BPS:

The Branch Points to consider.

--DatabaseName:

The input database name.

--help:

Show the help file.

Details:

All arguments should be specified in the order shown above.

Publication:

Please cite: Bitton et al., 2014, Genome Research.
Author:

Danny Asher Bitton, Departemnt of Genetics Evolution and Environment, University College London 

Contact:

d.bitton@ucl.ac.uk 

");q()}


## Extract arguments
if(length(args) != 5L)
  stop("--SeqFile,--HeaderFile, --ReadLength --BPS and --DatabaseName. Please use --help for details.")

seqfile <- args[1]
headfile <- args[2]
readlength <- args[3]
BPS <- args[4]
databaseName <- args[5]

for(i in args)
  {
    if(substr(i, 1, 10) == "--SeqFile=")
        seqfile <- gsub("--SeqFile=", "", x = i)
    if(substr(i, 1, 13) == "--HeaderFile=")
        headfile <- gsub("--HeaderFile=", "", x = i)
    if(substr(i, 1, 13) == "--ReadLength=")
        readlength <- as.numeric(gsub("--ReadLength=", "", x = i))
    if(substr(i, 1, 6) == "--BPS=")
        BPS <- unlist(strsplit(gsub("--BPS=", "", x = i),","))
    if(substr(i, 1, 15) == "--DatabaseName=")
        databaseName  <- gsub("--DatabaseName=", "", x = i)
  }



#reads all introns ids file (1 id per line)
ids<-read.delim(headfile,stringsAsFactors=F,header=F)
int.spl<-ids[,1]

#read all intron sequence files (1 sequence per line)
intseq<-read.delim(seqfile ,stringsAsFactors=F,header=F)
int.seq<-intseq[,1]

#Remove database if already exist
file.dir<-list.files(path = ".", pattern = ".fa$")
	if(databaseName %in% file.dir){
		print("Database Name already exists, removing it")
		system(paste("rm",databaseName))
	}
	
####Write lariat database

for ( k in 1:length(int.seq)){
	
		for(i in 1:nchar(int.seq[k])){

			#setting the 5' and 3' parts of the larait sequence
		
			temp3prime<-int.seq[k]
			temp5prime<-substr(int.seq[k],0,(nchar(int.seq[k])-i))
			bp<-substr(int.seq[k],(nchar(int.seq[k])-i+1),(nchar(int.seq[k])-i+1))
			
			#Account for BP letters as specified
			if(bp %in% BPS ){
				#If intron length is greater that the readlength-1 take only the first readlength-1 bases (else use all intron sequence)  
								
				if(temp3prime >= (readlength-1) ) {
				
					temp3prime<-substr(int.seq[k],0,(readlength-1))
		
				}
			
				#If the sequence upstream of the branch site is greater that the readlength-1 take only the last readlength-1 bases (else use all upstream sequence)
			
				if(temp5prime >= (readlength-1) ) {
				
					temp5prime<-substr(temp5prime,(nchar(temp5prime)-(readlength-2)),nchar(temp5prime))
		
				}
			
				#concatentate all parts 

				lariseq<-paste(temp5prime,bp,temp3prime,sep="")

				#write into a  lariat database add unique id along with the exact branch point position
				write(paste(">",int.spl[k],"£id",bp,"@",i,"\n",lariseq,sep=""),databaseName,append=TRUE)
						
			}
		}
}

#####Lariat skiping

#get Gene IDs
gene.ids<-unlist(strsplit(int.spl,";"))[seq(1,length(unlist(strsplit(int.spl,";"))),6)]

# choose genes with more than two introns
forlariat<-table(gene.ids)
forlariat<-forlariat[forlariat >=2]

#names of all lariatskip genes
lariskip<-names(forlariat)

#make lists of ids and sequences for each gene
larilist<-list()
larilistids<-list()

for(i in 1:length(lariskip)){

		larilist<-c(larilist,list(int.seq[which(gene.ids==lariskip[i])]))
		larilistids<-c(larilistids,list(int.spl[which(gene.ids==lariskip[i])]))
}
names(larilist)<-lariskip
names(larilistids)<-lariskip



#add the larait skipping to database
the_bin<-lapply(seq_along(larilist),function(idx){
      

	for(i in 1:(length(larilist[[idx]])-1)){
	 
		 for(j in (i+1):length(larilist[[idx]])){

				intskip3primefull<-paste(larilist[[idx]][i],larilist[[idx]][j],sep="")
				intskip3prime<-paste(larilist[[idx]][i],larilist[[idx]][j],sep="")
				intids<-paste(names(larilistids[idx]),"int",i,j,sep="_")
				
				# check size according to read length as before 3' 
				if(intskip3prime >= (readlength-1) ) {
				
					intskip3prime<-substr(intskip3prime,0,(readlength-1))
		
				}
				
		
				for(l in 1:nchar(larilist[[idx]][j])){
				
					intskip5prime<-substr(intskip3primefull,0,(nchar(intskip3primefull)-l))
					bp<-substr(intskip3primefull,(nchar(intskip3primefull)-l+1),(nchar(intskip3primefull)-l+1))
					

					#Account for BP letters as specified
					if(bp %in% BPS){		

						# check size according to read length as before 5'
						if(intskip5prime >= (readlength-1) ) {
				
							intskip5prime<-substr(intskip5prime,(nchar(intskip5prime)-(readlength-2)),nchar(intskip5prime))
		
						}
					
						#concatenate sequences as before 
						larskip.seq<-paste(intskip5prime,bp,intskip3prime,sep="")
					
						#write into a  lariat database add unique id along with the exact branch point position within the 3' intron				
						write(paste(">",intids,"£id",bp,"@",l,"\n",larskip.seq,sep=""),databaseName,append=TRUE)
					}
			     }
		}		

	}
})



