#!/usr/bin/env Rscript

## automating the standard DADA2 ITS2 pipeline from:
## 		https://benjjneb.github.io/dada2/ITS_workflow.html
## Plus Trimming!!!
## Plus ITSx

## suggested command to run: 
##		conda activate cutadapt
##		./v20240723__standard.trimm_DADA2_ITS.R > v20240723__standard.trimm_DADA2_ITS.R.out 2>&1
##		conda deactivate



##########################################################################################################################################
######### SET THESE VARIABLES TO YOUR SPECIFICS ##########################################################################################
##########################################################################################################################################

# rm(list=ls())					                        # clear environment variables if needed

# for getting input fastq files
fqpath="~/fastq/clonostachys"	# folder containing fastq files
	## ^don't put trailing '/' for this path!
R1.pattern="_R1.fastq.gz"		                    # unique pattern to match your forward reads
R2.pattern="_R2.fastq.gz"		                    # unique pattern to match your forward reads

# for primer detection before/after cutadapt
FWD="GTGARTCATCGAATCTTTG"   	                  # fITS7
REV="TCCTCCGCTTATTGATATGC"  	                  # ITS4

# path to programs wrapped in this script
cutadapt="/home/kkyle/miniconda3/envs/cutadapt/bin/cutadapt"	# klassenlab1
Trimmomatic="java -jar /usr/local/bin/Trimmomatic-0.39"		# klassenlab1
ITSx="/usr/local/bin/ITSx"													      # klassenlab1

# argument to set threads for cmds that support it
n.threads=20

# parameters for Trimmomatic command
trimm.flags = paste("PE","-threads",n.threads)            # PE is paired mode
trimm.steps = paste("SLIDINGWINDOW:5:20"
				   )
trimm.string = "_Trimm0.39PE_SW5.20"		                  # for output file namings

# path to UNITE fungal db for taxa calling
unite.ref="~/ms_peptaibol/sh_general_release_dynamic_s_04.02.2020.fasta"	# klassenlab1

##########################################################################################################################################


## ----load libraries---------------------------------------------------------------------------------------------------------------------

# get latest version of dada2 if don't have it
#library("devtools")
#devtools::install_github("benjjneb/dada2", ref="v1.20") # change the ref argument to get other versions

# load libraries
library(dada2); packageVersion("dada2")
library(phyloseq); packageVersion("phyloseq")
library(ShortRead); packageVersion("ShortRead")
library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")
#library(dplyr)		# not necessary for DADA2 pipeline, but useful for dataframe manipulation

## ----get input files setup--------------------------------------------------------------------------------------------------------------

fnFs=sort(list.files(fqpath,pattern=R1.pattern,full.names=T))
fnRs=sort(list.files(fqpath,pattern=R2.pattern,full.names=T))


## ----get sample names-------------------------------------------------------------------------------------------------------------------

# Extract sample names, assuming filenames have format:
get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(fnFs, get.sample.name))
# make sure looks hunky dory
head(cbind(sample.names,fnFs,fnRs))

## ----Trimmomatic--------------------------------------------------------------------------------------------------------------

# output directory(s)
fqpath.trimm <- file.path("./trimmomatic")				# trimmed fastqs outdir
if(!dir.exists(fqpath.trimm)) dir.create(fqpath.trimm)
fqpath.trimm.out <- file.path(fqpath.trimm,"out")		# log, summary, and .out files outdir
if(!dir.exists(fqpath.trimm.out)) dir.create(fqpath.trimm.out)

# outfile names
fnFs.trimm.P <- file.path(fqpath.trimm, paste0(sample.names,trimm.string,"_R1__P.fastq.gz"))
fnFs.trimm.U <- file.path(fqpath.trimm, paste0(sample.names,trimm.string,"_R1__UnP.fastq.gz"))

fnRs.trimm.P <- file.path(fqpath.trimm, paste0(sample.names,trimm.string,"_R2__P.fastq.gz"))
fnRs.trimm.U <- file.path(fqpath.trimm, paste0(sample.names,trimm.string,"_R2__UnP.fastq.gz"))

trimlog <-  file.path(fqpath.trimm.out, paste0(sample.names,trimm.string,".trimlog"))
trimsumm <- file.path(fqpath.trimm.out, paste0(sample.names,trimm.string,".trimm.summary"))
trimm.out <-file.path(fqpath.trimm.out, paste0(sample.names,trimm.string,".stdoutstderr"))

# run Trimmomatic on each sample
for(i in seq_along(fnFs)) {
	system(paste(Trimmomatic,
					trimm.flags,						# PE/SE and threads	
					"-trimlog", trimlog[i],				# file name for trimlog output
					"-summary", trimsumm[i],			# file name for summary output
					fnFs[i], fnRs[i],					# input
  					fnFs.trimm.P[i], fnFs.trimm.U[i], 	# output F
  					fnRs.trimm.P[i], fnRs.trimm.U[i],	  # output R
  					trimm.steps,						# trimming options
  					"2>&1 >",							  # send stdout and stderr to same file
  					trimm.out[i]						# stdout/stderr output
	) )  
}

## ----identify/remove primers------------------------------------------------------------------------------------------------

allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}

primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}


## ---- get all primer orients and save to objects ------------------------------------------------------------------------------------------------

FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)


## ---- filter out reads w/ Ns (req for dada function.. and apparently primer searching..?) -----------------------------------------------------------------------

# Put N-filterd files in filtN/ subdirectory
fqpath.filtN <- file.path(fqpath.trimm, "filtN")
fnFs.filtN <- file.path(fqpath.filtN, gsub("\\.fastq\\.gz$","_filtN.fastq.gz",basename(fnFs.trimm.P))) 
fnRs.filtN <- file.path(fqpath.filtN, gsub("\\.fastq\\.gz$","_filtN.fastq.gz",basename(fnRs.trimm.P))) 

filterAndTrim(	fnFs.trimm.P, fnFs.filtN, 
				        fnRs.trimm.P, fnRs.filtN,
				        maxN = 0, 
				        multithread = TRUE,
				        compress=T
			        )


## ----check where primers are detected in reads------------------------------------------------------------------------------------------------

# may need to reverse complement and repeat if finding in the wrong places
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]),
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]),
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))


## ----remove primers (clipping) with cutadapt ------------------------------------------------------------------------------------------------

# Run shell commands from R
system2(cutadapt, args = "--version")
# setup cutadapt output dir and output filenames
fqpath.cut <- file.path(fqpath.filtN, "cutadapt")
if(!dir.exists(fqpath.cut)) dir.create(fqpath.cut)
fnFs.cut <- file.path(fqpath.cut, gsub("\\.fastq\\.gz$","_clip.fastq",basename(fnFs.filtN)))
fnRs.cut <- file.path(fqpath.cut, gsub("\\.fastq\\.gz$","_clip.fastq",basename(fnRs.filtN))) 
																  # ^note we dont want zipped after clipping

# setting reverse complements of primers
FWD.RC <- dada2:::rc(FWD)	
REV.RC <- dada2:::rc(REV)

# flags for cutadapt
R1.flags <- paste("-g", FWD, "-a", REV.RC) # Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC) # Trim REV and the reverse-complement of FWD off of R2 (reverse reads)

# Run Cutadapt
for(i in seq_along(fnFs)) {
  system2(cutadapt, 
  	args = c(   R1.flags, R2.flags,
  				      "-n", 2, # -n 2 required to remove FWD and REV from reads
                "-m", 1, # prevent zero-length seqs in output which mess up plotQualityProfile
                "-o", fnFs.cut[i], "-p", fnRs.cut[i], 	# output files
                fnFs.filtN[i], fnRs.filtN[i]			      # input files
                ,"-j 0" # auto-detect cores... apparently parallel not supported on python 2... :(
  	), 
  	stdout = paste0(fnFs.cut[i],".cutadapt.out"), 
  	stderr = paste0(fnFs.cut[i],".cutadapt.err")
  ) 
}


## ----check for no more primers detected after cutadapt------------------------------------------------------------------------------------

rbind(	FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
     	  FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]),
		    REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]),
       	REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))


## ----setting clipped filenames and get sample names-------------------------------------------------------------------------------------

# Forward and reverse fastq filenames have the format:
cutFs <- sort(list.files(fqpath.cut, pattern = ".*R1.*clip\\.fastq$", full.names = TRUE))
cutRs <- sort(list.files(fqpath.cut, pattern = ".*R2.*clip\\.fastq$", full.names = TRUE))


## ----forward clipped qual plots----------------------------------------------------------------

pdf("cutFs.qualProf.pdf")
plotQualityProfile(cutFs)
dev.off()
pdf("cutRs.qualProf.pdf")
plotQualityProfile(cutRs)
dev.off()


## ----filter and trim --------------------------------------------------------------------------------------------------------
															# end .gz because setting compress to TRUE in filterAndTrim
filtFs <- file.path(fqpath.cut, "filtered", gsub("\\.fastq$","_filt.fastq.gz",basename(cutFs)))	                                                                                         # compress=T in filterAndTrim here
names(filtFs) <- sample.names
filtRs <- file.path(fqpath.cut, "filtered", gsub("\\.fastq$","_filt.fastq.gz",basename(cutRs)))	                                                                                         # compress=T in filterAndTrim here
names(filtRs) <- sample.names

# written to file and commented out so dont have to redo a bunch, just read in filterAndTrim.out
out <- filterAndTrim(	cutFs, filtFs, 
						          cutRs, filtRs,
						          maxN = 0, 
						          maxEE = c(2,2), 
                      truncQ = 2, 
                      minLen = 50, 
                      rm.phix = TRUE, 
                      compress = TRUE, 
                      multithread = TRUE	# on windows, set multithread = FALSE
                    )  
write.table(out, "filterAndTrim.out")
#out = as.matrix(read.table("filterAndTrim.out"))
head(out) 


## ----learn error rates and plot---------------------------------------------------------------------------------------------

# setting randomize=T bc i dont feel like the first 10 samples in order necessarily represent all the data
errF <- learnErrors(filtFs, multithread = TRUE, randomize = T); saveRDS(errF,"errF.RDS")
errR <- learnErrors(filtRs, multithread = TRUE, randomize = T); saveRDS(errR,"errR.RDS")  

# can look at PDf if don't want to rerun errF/R
pdf("filt.plotErrors.pdf")
plotErrors(errF, nominalQ = TRUE)
plotErrors(errR, nominalQ = TRUE)
dev.off()


## ----derep identical reads--------------------------------------------------------------------------------------------------

# this step is not included in the v1.16 dada2 16S tutorial... see this: https://github.com/benjjneb/dada2/issues/1095
derepFs <- derepFastq(filtFs, verbose = FALSE); saveRDS(derepFs,"derepFs.RDS")
derepRs <- derepFastq(filtRs, verbose = FALSE); saveRDS(derepRs,"derepRs.RDS")


## ----sample inference (dada function!)----------------------------------------------------------------------------------------

dadaFs <- dada(derepFs, err = errF, multithread = TRUE, verbose = FALSE); saveRDS(dadaFs,"dadaFs.RDS")
dadaRs <- dada(derepRs, err = errR, multithread = TRUE, verbose = FALSE); saveRDS(dadaRs,"dadaRs.RDS")


## ----merge paired reads -----------------------------------------------------------------------------------------------------

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE); saveRDS(mergers,"mergers.RDS")


## ----construct sequence table-----------------------------------------------------------------------------------------------------------

seqtab <- makeSequenceTable(mergers)
rownames(seqtab) = sample.names; saveRDS(seqtab,"seqtab.RDS")
dim(seqtab)
table(nchar(getSequences(seqtab)))  	


## ----remove chimeras--------------------------------------------------------------------------------------------------------------------

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", 
									multithread=TRUE, verbose=FALSE
								   ); saveRDS(seqtab.nochim,"seqtab.nochim.RDS")
sum(seqtab.nochim)/sum(seqtab)			# 0.9669808
dim(seqtab.nochim)  					# ok, now 257 ASVs
table(nchar(getSequences(seqtab.nochim))) # mean 304.2 bp, med 300, min 50, max 415


## ----track reads through pipeline and plot---------------------------------------------------------------------------------------

getN <- function(x) sum(getUniques(x))

# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
if(nrow(out) == 1) {
  track <-cbind( out, 
                 getN(dadaFs), 
                 getN(dadaRs),
                 getN(mergers),
                 rowSums(seqtab.nochim)
  )
} else {
  track <-cbind(	out, 
                 sapply(dadaFs, getN), 
                 sapply(dadaRs, getN),
                 sapply(mergers, getN),
                 rowSums(seqtab.nochim)
  )
}
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names

# output csv and plot to pdf
write.csv(track,file = "track.csv")
pdf("track.barplot.pdf",width=12,height=6)
barplot(t(track),
		beside=T,legend.text=T,las=2,cex.names=1,
		col=c("black","red","blue","yellow","purple","green"),
		space=c(0.2,2),ylim=c(0,max(t(track))),
		args.legend = list(cex=0.6)
	)
dev.off()


## ----assign taxonomy--------------------------------------------------------------------------------------------------------

# Assign taxonomy 
taxa=assignTaxonomy(seqtab.nochim,unite.ref,multithread = T,tryRC = T); saveRDS(taxa,"taxa.RDS")
# inspect
taxa.print <- taxa  			# Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)  # looks like fungi?!


## ----phyloseq-----------------------------------------------------------------------------------------------------------------------------

ps = phyloseq( 
  otu_table(seqtab.nochim,taxa_are_rows = F), 
  tax_table(taxa)
  )

# refseq slot
sequences <- Biostrings::DNAStringSet(taxa_names(ps))
names(sequences) <- taxa_names(ps)
ps <- merge_phyloseq(ps, sequences)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))

# if you want to write these tables to read into excel
write.csv(seqtab.nochim, file = "ps.otutable.csv")
write.csv(taxa, file= "ps.taxtable.csv")
Biostrings::writeXStringSet(refseq(ps),"ps.ASV.fasta")

saveRDS(ps,"ps.RDS")

## ----plot-----------------------------------------------------------------------------------------------------------------------------

theme_set(theme_bw())
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))

# plot just top 20 ASVs
top20 = names(sort(taxa_sums(ps), decreasing = T))[1:20]
ps.prop.top20=prune_taxa(top20,ps.prop)

# color by genus
pdf("genus_stackedbar.pdf",width=12,height=6)
plot_bar(	ps.prop.top20, 
			fill="Genus", 
			title = "Top 20 taxa"
		) + theme(axis.text=element_text(size=14))
dev.off()

# color by family
pdf("fam_stackedbar.pdf",width=12,height=6)
plot_bar(	ps.prop.top20, 
			fill="Family", 
			title = "Top 20 taxa"
			) + theme(axis.text=element_text(size=14))
dev.off()


## ----save image-----------------------------------------------------------------------------------------------------------------------------

save.image(paste0(format(Sys.time(), "%b%d%Y_%X%Z"),".RData"))


## ----ITSx-----------------------------------------------------------------------------------------------------------------------------

system2(ITSx, 
		args=c("-i","ps.ASV.fasta",		# input fasta
				"-o","ITSx_psASV",		    # output prefix
				"--cpu",n.threads		      # multithread
			),
		stdout = "ITSx.out",
		stderr = "ITSx.err"
		)

