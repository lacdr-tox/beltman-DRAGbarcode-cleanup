##############################################################################################
##############################################################################################
#
# File: tagcleanup.R
#
# Script to remove spurious barcodes in barcoding deepseq data generated from the
# genetic 'barcode mouse' construct as generated in the lab of Ton Schumacher (NKI, NL)
#
# Expected input: 1. data files containing the barcodes, tags and the frequency of the
# 								combination per sample within the directory 'rawdata'; two example input
# 								files are given in this directory (this is actually output from scripts
# 								by Arno Velds on the raw sequencing data; these file names end with
# 								'perfect-long-clipped-template2-e6.txt')
# 								2. a file containing the samples and sample information (such as sample index)
# 								that is also part	of the data file names
#
# Output: 				two types of output file are created; first, files containing the barcodes
# 								that remain after the cleaning procedure for each sample (file names start
# 								with 'cleaned'); second, files containing those barcodes splitted in
# 								V, D, J parts and deletions and additions (file names start with 'splitbcs')
#
#
# Author: Joost Beltman (LACDR, Leiden University, NL)
# Date: Jan 25, 2022
#
##############################################################################################
##############################################################################################


#load required packages
library("ggplot2")
library("plyr")
library(stringdist)
library(gplots)
library(RColorBrewer)

#threshold for read numbers per tag
freqthresh = 10

#to determine data paths
basedir = "../"
data.path.read = paste(basedir, "rawdata/", sep="")
data.path.write = paste(basedir, "output/", sep="")
data.path.R = paste(basedir, "R/", sep="")

#set file name to read samples and store sample ID names
samnames = read.table(paste(data.path.read, "samplenames.txt", sep=""), sep="\t", header=TRUE)
listsam.fn = paste('5290_', samnames$deepseqID, '_', samnames$samplecode, '_', samnames$index, '-perfect-long-clipped-template2-e6.txt', sep="")
listsam.sn = as.character(samnames$samplecode)


#store short and long sequence within constant part of J region
shortconst = "TAGCAAGCTC"
longconst = "TAGCAAGCTCGAGAGTAGACCTACTGGAATCAGACCGCCACCATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGACCTACGGCGTGCAGTGCTTCAGCCGCTACCCCGACCACATGAAGCAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGGACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTACAACAGCCACAACGTCTATATCATGGCCGACAAGCAGAAGAACGGCATCAAGGTGAACTTCAAGATCCGCCACAACATCGAGGACGGCAGCGTGCAGCTCGCCGACCACTACCAGCAGAACACCCCCATCGGCGACGGCCCCGTGCTGCTGCCCGACAACCACTACCTGAGCACCCAGTCCGCCCTGAGCAAAGACCCCAACGAGAAGCGCGATCACATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAG"


samnum = length(listsam.sn)
#loop over all samples if there is at least one sample
if (samnum >= 1) {
	for (samk in 1:samnum) {

		#read sample file and put in variable sam
		assign(listsam.sn[samk], read.table(paste(data.path.read, listsam.fn[samk], sep=""), sep="\t", header=FALSE))
		sam = get(listsam.sn[samk])
		names(sam) = c("var", "varconst", "tag", "freq")


		######################################
		#REMOVE BARCODES AND TAGS CONTAINING N
		######################################

		sam.rmN.bc = sam[!grepl("N", sam$var),]
		sam.rmN.tag = sam.rmN.bc[!grepl("N", sam.rmN.bc$tag),]


		#################################################
		# REMOVE BARCODES WITH VARIATION IN CONSTANT PART
		#################################################

		#select for short specific region early in constant part
		sam.firstselconst = sam.rmN.tag[grepl(shortconst, sam.rmN.tag$varconst),]

		#find position of first match
		matchpos = as.numeric(regexpr(shortconst, sam.firstselconst$varconst, fixed=T))

		#split in VDJ part and constant part
		VDJ = substr(sam.firstselconst$varconst, 1, matchpos-1)
		Jconst = substr(sam.firstselconst$varconst, matchpos, nchar(as.character(sam.firstselconst$varconst)))

		#find equal part of input sequence to compare to
		inputJconst = substr(rep(longconst, length(Jconst)), 1, nchar(Jconst))

		#find hamming distance between input and output sequence
		dist = stringdist(Jconst, inputJconst, method="hamming")

		#correct for erroneous splitting that occasionally occurs by the blast script of Arno Velds
		sam.firstselconst$var = VDJ

		#keep only sequences with no errors in the constant part
		sam.rmconstJ = sam.firstselconst[dist==0,]


		####################################################################
		# REMOVE TAGS BELOW THRESHOLD AND KEEP MOST FREQUENT BARCODE PER TAG
		####################################################################

		#first remove tags below threshold
		sam.rmtags = sam.rmconstJ[sam.rmconstJ$freq >= freqthresh,]

		#find all unique tags
		tags = unique(sam.rmtags$tag)

		tagsatend = as.data.frame(matrix(nrow=0, ncol=5))
		names(tagsatend) = c("var", "varconst", "tag", "freq")
		#loop over all tags to find the barcode with which the tag has the maximal read number
		for (j in 1:length(tags)) {

			tag = sam.rmtags[sam.rmtags$tag %in% tags[j],]
			topfreq = tag[tag$freq==max(tag$freq),]
			tagsatend = rbind(tagsatend, topfreq)
		}


		############################################################
		# FISH FOR LEFT-OVER BARCODES AND JOIN TAGS FOR SAME BARCODE
		############################################################

		left.withfish = sam[sam$var %in% tagsatend$var,]

		if (nrow(left.withfish) > 0) {

		  #count total reads per barcode
		  rds.fish = ddply(left.withfish, "var", summarize, rds = sum(freq))
	  	names(rds.fish)[2] = "rds.fish"

		  #count total number of molecules per barcode
	  	nummol.fish = ddply(left.withfish, "var", nrow)
  		names(nummol.fish)[2] = "nmol.fish"

  		halfsample = merge(rds.fish, nummol.fish, by="var")

  		} else {

  		  halfsample = as.data.frame(matrix(nrow=0, ncol=3))
        names(halfsample) = c("var", "rds.fish", "nmol.fish")

		}


		######################################################
		# REMOVE BARCODES WITH VARIATION IN VDJ CONSTANT PARTS
		######################################################

		#run script to split into V, D and J elements
		source(paste(data.path.R, "splitVDJ.R", sep=""))

		#print out the splitting of the barcodes that are considered consistent with the original construct
		write.table(add.del.ok, paste(data.path.write, "splitbcs_", listsam.sn[samk], ".txt", sep=""), sep="\t", quote=F, row.names=T)

		#find barcodes left after final cleaning step
		bcs.trueVDJ = rownames(add.del.ok)
		halfsam.trueVDJ = halfsample[halfsample$var %in% bcs.trueVDJ,]
		halfsam.falseVDJ = halfsample[!(halfsample$var %in% bcs.trueVDJ),]

		#print the barcodes left after cleaning along with their frequency
		halfsam.trueVDJ$var = as.character(halfsam.trueVDJ$var)
		write.table(halfsam.trueVDJ, paste(data.path.write, "cleaned_", listsam.sn[samk], ".txt", sep=""), sep="\t", quote=F, row.names=F)

	}
}
