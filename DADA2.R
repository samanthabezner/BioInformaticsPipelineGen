#-------------------------------------------------------

# DADA2 Script, with Commentary
# Code adopted from Brendan Daisley, Sonja Drosdowech, and David Huyben
# Annotated by Samantha Bezner
# DADA2 is a bioinformatic pipeline made by 
# https://benjjneb.github.io/dada2/tutorial.html
# Further help can be found there!

#-------------------------------------------------------


#Prior to starting the DADA2 pipeline, make sure that RTools is downloaded
#And that RStudio & R are both up to date. Look at DADA2's download page to 
#Ensure what version of DADA2 is supported on your current version of R. 

#Install packages and load libraries



if (!requireNamespace("BiocManager", quietly = TRUE, force = TRUE))
  
  install.packages("BiocManager")

BiocManager::install("dada2", version = "3.22", force = TRUE)



library(dada2)
library(ggplot2)
library(stats)
library(dplyr)



#-------------------------------------------------------

# Before running

#-------------------------------------------------------

# 1) Demultiplex your samples if necessary using cutadapt



# 2) Save demultiplexed files in folder named "reads" and set working directory to parent directory:

#Set your working directory to the file where your reads are located and where outputs will go, too.
setwd()

dir.create("Outputs")

#-------------------------------------------------------

# Setup

#-------------------------------------------------------



#Paths to demultiplexed reads

reads <- "Reads"
reads



#Dump the R sessioninfo for later

#writeLines(capture.output(sessionInfo()), "03_dada2_bac/RsessionInfo_dada2.txt")

# Get the filenames with relative path

# sort to ensure same order of fwd/rev reads

fnFs <- sort(list.files(reads, pattern="_R1_001.fastq", full.names=TRUE))
fnFs
fnRs <- sort(list.files(reads, pattern="_R2_001.fastq", full.names=TRUE))
fnRs

# Note: fnFs is a variable name meaning forward read file paths. fnRs is a variable name reverse read file paths, sort means put things in alphabetical order (keeps the files matched forward and reverse). 
#( ( ) ) means do this first. The ones on the inside get done first like in math. pattern= means only show files with “_R1_001.fastq” (or whatever you wanted) in their file name. 
#There are a few other things in the folder that we don’t want for this purpose. full.names = TRUE means that dada2 needs to read the full file location (/Users/…..) not just the file name.

# Basically, this is saying: Look in the reads folder and find only files containing “_R1_001.fastq” and output the full file names, which is the output that goes into sort. 
#Sort command then sorts them alphabetically and the <- tells R to save them as “fnFs”

# Get sample names only (remove path, and everything after the first "-")
# Assuming filenames have format: SAMPLENAME-XXX.fastq

sample.names <- sapply(strsplit(basename(fnFs), "_S"), `[`, 1)
sample.names


#list the files (not required)

#list.files(taxpath)




message ("###### checking for duplicated sample names")

any(duplicated(sample.names)) #Should be false. 


#-------------------------------------------------------

# Check read quality

#-------------------------------------------------------

#This will pick a random subset of 4 samples to look at read quality
#note: you can also plot an aggregate of all fastqs instead using aggregate=TRUE

ids<-round(runif(12,1,length(sample.names)))
ids

pdf("Outputs/qualprofiles.pdf")

plotQualityProfile(fnFs[ids])

plotQualityProfile(fnRs[ids])

dev.off()

# Explanation of the plot: The grey scale heat map means the frequency of each quality score at each base position. The green line shows the mean quality score at each base position. 
#The red line is not important. The quality score is how confident the sequencer is that the base call is correct. A quality score above 30 is good (99.9% correct).


#NOTE: The forward reads will show a higher quality score (roughly 40) and drop
#off later in the sample. Reverse reads, due to degredation of sequences by the
#time that they are ran, means they will be lower quality and drop off faster.

#-------------------------------------------------------

# Filter reads based on QC

#-------------------------------------------------------

message ("###### Filtering reads based on QC")


# Make filenames for the filtered fastq files

filtFs <- paste0(reads, "/", sample.names, "_F-filt.fastq")
filtRs <- paste0(reads, "/", sample.names, "_R-filt.fastq")

filtFs
filtRs

# Note: Explanation of the command: FiltFs is the variable name to hold paths for filtered forward reads. FiltRs is the variable name to hold paths for filtered reverse reads. 
#file.path() builds a full file path in the “filtered” folder. paste0 means combine strings without separators to make the file names for each sample, 
#which will include the previously assigned sample names (for example, “sample1”) and add the parts in quotations.



# *************** On Mac set multithread=TRUE / On Windows set multithread=FALSE**********

# Multithread essentially speeds up the process by allowing multiple cores to process multiple files
# at once. This is possible on Linux and Mac computers, but not on Windows. 



#For all of the 

out<-filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                   
                   truncQ=20, # Enter quality score cut off (truncates after </=2 QC score)

                   trimLeft=c(21,19), #(FORWARD, REVERSE) Used to trim the primers from the left side. 
                   
                   maxN=0,  #DADA2 automatically sets to zero. 
                   
                   maxEE=c(2,2), #(FORWARD, REVERSE), expected errors. Points of lower quality "allowed" to happen before cut off
                   
                   rm.phix=TRUE, #Removes the genome of PhiX, a bacteriphage that may skew results
                   
                   compress=TRUE, verbose=TRUE, multithread=FALSE)



write.table(out, file = "Outputs/reads_filtered.txt", sep="\t", col.names=NA, quote=F)

# Note: Explanation of the command: Out <- means store a summary whats filtered in the variable out. FnFs and FnRs are the raw forward and reverse file paths. 
#FiltFs and FiltRs are where the filtered reads will be saved. Familiar users will notice no truncLen=c(F,R) command, and that is because of V3-V5's unique, 7 base pair of overlap. If you truncate even 1 bp on both sides, 
#Overlap lowers from, 7 bp to 5 bp. If your quality can handle it, you should truncate based on quality score. If quality descends below ~25, truncate the sequence and you should have enough average reads to still have sufficient
#Reads for downstream analysis. MaxN=0 means no N’s (unknown bases) allowed. MaxEE=c(2,2) means forward reads can have up to 2 expected errors, reverse reads up to 2 expected errors. 
#Any worse reads are thrown out (reads, not samples, and each sample can have hundreds of thousands of reads, each of which is one microbe). truncQ=2 means cut reads off at the first base with a quality score below 2. 
#rm.phix=TRUE means remove reads that match phix, a common control spike-in used in Illumina sequencing. compress=TRUE means save the files as .gz (zipped). 
#multithread=TRUE means run filtering on multiple CPU core, good for running on MAC or servers. Keep FALSE on Windows base. 
#head(out) means show the first few rows of the summary table saved as (out) at the beginning of the code here.

#-------------------------------------------------------

# Learn the error rates - SLOW !! ----30 mins for each *errF* and *errR* on this run

#-------------------------------------------------------



message ("###### Learning error rates - SLOW !!")



# *************** On Mac set multithread=TRUE / On Windows set multithread=FALSE**********

errF <- learnErrors(filtFs, multithread=FALSE, randomize=TRUE)
errR <- learnErrors(filtRs, multithread=FALSE, randomize=TRUE)

#randomize=TRUE #don't pick the first 1mil for the model, pick a random set
 


#Plot the error rates and CHECK THE FIT

# Do not proceed without a good fit (compare to red line)

pdf("Outputs/err.pdf")

plotErrors(errF, nominalQ=TRUE)

plotErrors(errR, nominalQ=TRUE)  

dev.off()

# Explanation of plot: X axis is the quality score of each base in a sequencing read (higher = better). Y axis is the frequency at which dada2 observed a base problem. 
#Each panel is a type of nucleotide substitution ( Thiamine to Adenine, for example). So, A2A is blank since you can’t substitute a nucleotide for itself. 
#The size of the dots shows the number of observations and the black line through it is a fitted error model. The red line is the same in all and is a reference point. 
#Here, error rates decrease as quality scores increase, which is good.



#message ("###### Saving your R session...")

#save.image("dada2.RData") #Insurance in case your script dies. Delete this later

message ("###### Dereplicating the reads")



# Dereplication combines all identical sequencing reads into "unique sequences" with a corresponding "abundance": the number of reads with that unique sequence. 
#Dereplication substantially reduces computation time by eliminating redundant comparisons. Basically, if you have 10,000 reads from one type of bacteria, this will bin those reads together for downstream analysis.




derepFs <- derepFastq(filtFs, verbose=TRUE)

derepRs <- derepFastq(filtRs, verbose=TRUE)

# Name the derep-class objects by the sample names

names(derepFs) <- sample.names

names(derepRs) <- sample.names



#message ("###### Saving your R session...")

#save.image("03_dada2_bac/dada2.RData")  #Insurance in case your script dies. Delete this later



#-------------------------------------------------------

# Sample inference, merge paired reads, remove chimeras

#-------------------------------------------------------

message ("###### Inferring the sequence variants in each sample - SLOW!!")



dadaFs <- dada(derepFs, err=errF, multithread=FALSE, pool=FALSE, verbose=TRUE)

dadaRs <- dada(derepRs, err=errR, multithread=FALSE, pool=FALSE, verbose=TRUE)

# Note: Command explanation: dada() infers amplicon sequence variants (which is the bin for identical reads) from the filtered reads. Err=errF means to use the error model that it learned for forward reads. 
#dadaFs <- means store the results under the variable named dadaFs.

# The DADA2 algorithm inferred number of true sequence variants from the unique sequences in the first sample. 
#There is much more to the dada-class return object than this (see help("dada-class") for some info), 
#including multiple diagnostics about the quality of each denoised sequence variant.



# overlap the ends of the forward and reverse reads

message ("###### merging the Fwd and Rev reads")

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, minOverlap=5, verbose=TRUE)

    #, justConcatenate=TRUE for V59
 
# Merge paired reads (and inspect first sample). Paired reads are the forward and reverse sequences that overlap and make one long read when put together.
#For V3-V5 region, because the maximum overlap of the two regions is 7bp, you need to manually set the minimum overlap otherwise it defaults to min. 20bp, which misses our window entirely. 


# make the sequence table, samples are by rows
seqtab <- makeSequenceTable(mergers)



message ("###### summarize the output by sequence length")
table(nchar(getSequences(seqtab)))

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", verbose=TRUE, multithread=FALSE) #minFoldParentOverAbundance=2
# Remove chimeras from the sample, which are the aggregated, erroneous reads that arise from PCR complications that create a unique, nonexistent sequence.


#First number is # of samples. Second is number of Unique sequences.  Should be less than 10,000 generally. IF diverse, may be above but unlikely. If super high, run next step.
dim(seqtab.nochim)

sum(seqtab.nochim)/sum(seqtab)


#If under 10,000, just get rid of singletons (>1); if over, filter up to 10. 
seqtab.nochim.filt <- as.data.frame(seqtab.nochim) %>% #Convert to dataframe from matrix
  sjmisc::rotate_df() %>%  #Trsansposing table such that sequences are rows
  filter(rowSums(.) > 5) %>% #filter out singletons, 10 is maximum conservatively 
  sjmisc::rotate_df() #TRANSPOSE IT BACK

dim(seqtab.nochim.filt) #(SAMPLE, ASV)



#seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% seq(420,440)]

#seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% seq(389,430)]



#let's write the table, just in case

#samples are rows

write.table(seqtab.nochim.filt, file="Outputs/temp_dada2_nochim.txt", sep="\t", col.names=NA, quote=F)

write.table((t(seqtab.nochim.filt)), file="Outputs/temp_dada2_nochim_transpose.txt", sep="\t", col.names=NA, quote=F)

# Or save the Rsession save.image("dada2.RData")

#message ("###### Saving your R session...")

#save.image("03_dada2_bac/dada2.RData")  #Insurance in case your script dies. Delete this later



#---



#-------------------------------------------------------

# Sanity check

#-------------------------------------------------------

message ("###### sanity check - how many reads made it: readsout.txt")

# Check how many reads made it through the pipeline

# This is good to report in your methods/results

getN <- function(x) sum(getUniques(x))

track <- cbind(out, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab.nochim), rowSums(seqtab.nochim.filt))

colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim", "nonchim_filtered_5")

rownames(track) <- sample.names

write.table(track, file="Outputs/readsout.txt", sep="\t", col.names=NA, quote=F)

track

#-------------------------------------------------------

# Assign taxonomy.

#-------------------------------------------------------


taxa <- assignTaxonomy(as.matrix(seqtab.nochim.filt), "C:/Users/huybendlab/Desktop/SB2025AxMicrobiome/silva_nr99_v138.2_toGenus_trainset.fa.gz", multithread=FALSE)
taxa <- addSpecies(taxa, "silva_v138.2_assignSpecies.fa.gz") #gets the taxonomy table down to species level. 

taxa.print <- taxa # Removing sequence row names for display only
rownames(taxa.print) <- NULL
head(taxa.print) 


#-------------------------------------------------------

# Export ASV + Taxa Table!

#-------------------------------------------------------
  
write.csv(seqtab.nochim.filt, file = "ASV.csv", row.names = TRUE, col.names = TRUE)
write.csv(taxa, file = "Taxa.csv", row.names = TRUE, col.names = TRUE)

  
