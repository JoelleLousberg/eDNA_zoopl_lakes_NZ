#Create a new folder for my fasta files on HDrive
rm(list= ls())
setwd("H:/Joelle.Lousberg/CAW-26-06_Joelle_Ch2_step3_fwh/CAW-26-06_fasta")
my_dirs <- list.files("H:/Joelle.Lousberg/CAW-26-06_Joelle_Ch2_step3_fwh/CAW-26-06_fasta", 
                      pattern = "CAW-26-06-",
                      recursive = TRUE, include.dirs = TRUE)

files <- sapply(my_dirs,list.files, full.names=TRUE)
new_dir <- "All_sequences"
dir.create(new_dir, recursive = TRUE)

for(file in files) {
  file.copy(file, new_dir)
}




##Run on HPC terminal 1





rm(list= ls())
#libraries
#if (!requireNamespace("BiocManager", quietly = TRUE))
 # install.packages("BiocManager")
#BiocManager::install("Biostrings")
#BiocManager::install("DECIPHER")
#BiocManager::install("dada2")
#BiocManager::install("phyloseq")


library(phyloseq)
#library(DECIPHER)
library(Rcpp)
library(dada2)
library(ShortRead)
library(Biostrings)
library(ggplot2)
library(stringr) # not strictly required but handy
library(readr)
###set.seed(106)
library(seqinr)
library(data.table)
library(plyr)
library(tidyverse)

###Choose your working directory
setwd("/srv/users/Joelle/CAW-26-06/")

path <- "/srv/users/Joelle/CAW-26-06/"

list.files(path)

#I only want certain samples from my list actually having the FwhF1/R1 primers,
#in this case samples 1-10, 31 & 32.

# Forward reads
#fnFs <- sort(list.files(
#  path,
#  pattern = "CAW-25-36-(00[1-9]|010|031|032).*R1_001.fastq.gz$",
#  full.names = TRUE
#))

# Reverse reads
#fnRs <- sort(list.files(
#  path,
#  pattern = "CAW-25-36-(00[1-9]|010|031|032).*R2_001.fastq.gz$",
#  full.names = TRUE
#))

###Instead here the old code that includes all files:
fnFs <- sort(list.files(path, pattern = "R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "R2_001.fastq.gz", full.names = TRUE))

FWD <- "YTCHACWAAYCAYAARGAYATYGG"
REV <- "ARTCARTTWCCRAAHCCHCC"

allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients

primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs[[1]]))



###Using cutadapt for trimming the primers off

cutadapt <- "/home/john/miniconda3/bin/cutadapt"

system2(cutadapt, args = "--version") 


###create path

path.cut <- file.path(path, "fwhF1R1/cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

# Trim FWD off of R1 (forward reads) - 
R1.flags <- paste0("-g", " ^", FWD) 
# Trim REV off of R2 (reverse reads)
R2.flags <- paste0("-G", " ^", REV) 
# Run Cutadapt
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c("-e 0.08 --discard-untrimmed --nextseq-trim=20", R1.flags, R2.flags,
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs[i], fnRs[i])) # input files
}

rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))

# Forward and reverse fastq filenames have the format:
cutFs <- sort(list.files(path.cut, pattern = "R1_001.fastq.gz", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "R2_001.fastq.gz", full.names = TRUE))

if(length(cutFs) == length(cutRs)) print("Forward and reverse files match. Go forth and explore")
if (length(cutFs) != length(cutRs)) stop("Forward and reverse files do not match. Better go back and have a check")

# Extract sample names, assuming filenames have format:
get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(cutFs, get.sample.name))
head(sample.names)

# plot the first 10 forward and reverse reads (or fewer if less than 10 exist)
n_show <- min(10, length(cutFs))

fwd_qual_plots <- plotQualityProfile(cutFs[1:n_show]) +
  scale_x_continuous(breaks = seq(0, 300, 10)) +
  scale_y_continuous(breaks = seq(0, 40, 2)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

rev_qual_plots <- plotQualityProfile(cutRs[1:n_show]) +
  scale_x_continuous(breaks = seq(0, 300, 10)) +
  scale_y_continuous(breaks = seq(0, 40, 2)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Show the plots
fwd_qual_plots
rev_qual_plots

jpeg(file="fwhF1R1/fwhF1R1.Quality.Plot.F.jpg",res=300, width=15, height=8, units="in")
fwd_qual_plots
dev.off()

jpeg(file="fwhF1R1/fwhF1R1.Quality.Plot.R.jpg",res=300, width=15, height=8, units="in")
rev_qual_plots
dev.off()



##Look at plots and decide truncation length __________________________________________________________


filtpathF <- file.path(path.cut, "filtered", basename(cutFs))
filtpathR <- file.path(path.cut, "filtered", basename(cutRs))

###fwhF1R1 is 178bp long, quality goes down at 210bp
# 2,4 for new run, default was 4,6 maxEE
out <- filterAndTrim(cutFs, filtpathF, cutRs, filtpathR,
                     truncLen=c(100,100), maxEE=c(2,4), truncQ=2, maxN=0, rm.phix=TRUE,
                     compress=TRUE, verbose=TRUE, multithread=TRUE)


out2 <- as.data.frame(out)
str(out2)
out2$perc <- (out2$reads.out/out2$reads.in)*100
out2

#check if rows = 0

data <- as.data.frame(out) 
rows_with_zero <- which(data$reads.in == 0 | data$reads.out == 0)
print(rows_with_zero)


sample.names <- sapply(strsplit(basename(filtpathF), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
sample.namesR <- sapply(strsplit(basename(filtpathR), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
if(identical(sample.names, sample.namesR)) {print("Files are still matching.....congratulations")
} else {stop("Forward and reverse files do not match.")}
names(filtpathF) <- sample.names
names(filtpathR) <- sample.namesR


set.seed(100) # set seed to ensure that randomized steps are replicatable

# Learn forward error rates

# Learn error rates

errF <- learnErrors(filtpathF, nbases=1e8,multithread=TRUE, verbose = TRUE)

errR <- learnErrors(filtpathR, nbases=1e8, multithread=TRUE, verbose = TRUE)

###plot errors

jpeg(file="fwhF1R1/fwhF1R1.Error.F.Plot.jpg",res=300, width=15, height=8, units="in")
plotErrors(errF, nominalQ=TRUE)
dev.off()

jpeg(file="fwhF1R1/fwhF1R1.Error.R.Plot.jpg",res=300, width=15, height=8, units="in")
plotErrors(errR, nominalQ=TRUE)
dev.off()

exists <- file.exists(filtpathF)
derepFs <- derepFastq(filtpathF[exists], verbose=TRUE)
derepRs <- derepFastq(filtpathR[exists], verbose=TRUE)

dadaF.pseudo <- dada(derepFs, err=errF, multithread=TRUE, pool="pseudo")
dadaR.pseudo <- dada(derepRs, err=errR, multithread=TRUE, pool="pseudo")

mergers <- mergePairs(dadaF.pseudo, derepFs, dadaR.pseudo, derepRs, maxMismatch = 1, minOverlap = 10, 
                      verbose=TRUE)

getN <- function(x) sum(getUniques(x))

out_fq1 = ShortRead::countFastq(fnFs) %>% dplyr::mutate(mean_width=nucleotides/records)%>%dplyr::select(records,mean_width)%>% rownames_to_column("sample_id") %>% setNames(c("sample_id","input","mean_width"))


out_1 <- out_fq1 %>%dplyr::select(sample_id,input)


out_2 <- out %>%
  as.data.frame() %>%
  rownames_to_column("sample_id") %>%
  setNames(c("sample_id", "demultiplexed", "filtered"))
head(out_2,100)


out_3 = data.frame(names(dadaF.pseudo),sapply(dadaF.pseudo, getN)) %>% setNames(c("sample_id","denoised"))
out_4 = data.frame(names(dadaF.pseudo),sapply(mergers, getN)) %>% setNames(c("sample_id","merged"))
#out_5 = data.frame(rownames(seqtab.nochim),rowSums(seqtab.nochim)) %>% setNames(c("sample_id","nonchim"))

out_1[,'sample_id'] = gsub("(.*)_S\\d+_L001.*","\\1", out_1[,'sample_id'] )
out_2[,'sample_id'] = gsub("(.*)_S\\d+_L001.*","\\1", out_2[,'sample_id'] )

out_final = dplyr::left_join(out_1,out_2, by="sample_id")
out_final = dplyr::left_join(out_final,out_3, by="sample_id")
out_final = dplyr::left_join(out_final,out_4, by="sample_id")
#out_final = dplyr::left_join(out_final,out_5, by="sample_id")
out_final = out_final %>% base::replace(is.na(.),0)

write.csv(out_final, "fwhF1R1/track.csv")

seqtab <- makeSequenceTable(mergers)

saveRDS(seqtab, "fwhF1R1/fwhF1R1.seqtab.rds")

#Read merge old and new files, chimera check them and merge them
seqtab1 <- readRDS("/srv/users/Joelle/CAW-25-36/fwhF1R1_CO1/fwhF1R1_CO1.seqtab.rds")
seqtab2 <- readRDS("fwhF1R1/fwhF1R1.seqtab.rds")


##checking for chimeras for each sequencing run separately before pooling them

#seqtab1 <- seqtab1[,nchar(colnames(seqtab1)) %in% seq(177, 178)]
#seqtab1.nochim <- removeBimeraDenovo(seqtab1, multithread=TRUE, verbose=TRUE)
#saveRDS(seqtab1.nochim, "fwhF1R1/fwhF1R1.seqtab1.nochim.rds")

#seqtab2 <- seqtab2[,nchar(colnames(seqtab2)) %in% seq(177, 178)]
#seqtab2.nochim <- removeBimeraDenovo(seqtab2, multithread=TRUE, verbose=TRUE)
#saveRDS(seqtab2.nochim, "fwhF1R1/fwhF1R1.seqtab2.nochim.rds")

merged.seqtab <- mergeSequenceTables(seqtab1, seqtab2)
seqtab <- merged.seqtab[,nchar(colnames(merged.seqtab)) %in% seq(177, 178)]
seqtab.nochim <- removeBimeraDenovo(merged.seqtab, multithread=TRUE, verbose=TRUE)
saveRDS(seqtab.nochim, "fwhF1R1/fwhF1R1.seqtab.nochim.rds")

#saveRDS(merged.seqtab, "fwhF1R1/fwhF1R1.seqtab.nochim_all.rds")

###Look at the length distribution of the sequences (histogram)

trimtable <- as.data.frame(table(nchar(getSequences(merged.seqtab))))
colnames(trimtable) <- c("Length.bp", "Frequency")
trimtable$Frequency <- as.numeric(trimtable$Frequency)
str(trimtable)

library(ggplot2)
g <- ggplot(trimtable, aes(x = Length.bp, y = Frequency)) +
  geom_bar(stat="identity")
jpeg(file="fwhF1R1/fwhF1R1.Chim.Dist.Plot.jpg",res=300, width=15, height=8, units="in")
g
dev.off()

###Do only Blast/Megan taxonomy. 

###turn rds into file that can be used for blast

seqtab.nochim <- readRDS("fwhF1R1/fwhF1R1.seqtab.nochim_all.rds")

asv_seqs <- colnames(seqtab.nochim)
head(asv_seqs)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

# making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
head(asv_fasta)
write.table(asv_fasta, "fwhF1R1/fwhF1R1_ASV_seq_IDs_and_sequences.tsv", sep="\t", quote=F, col.names=NA)
write(asv_fasta, "fwhF1R1/fwhF1R1.nochim_ASVs.fa")

# count table:
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "fwhF1R1/fwhF1R1.nochim_ASVs_counts.tsv", sep="\t", quote=F, col.names=NA)

##blast code - run in HPC terminal 1
#blastn -query fwhF1R1.nochim_ASVs.fa -db ~/../../srv/referencedb/COI/COI.NCBI.BOLD.blast.db -out fwhF1R1_CO1.chim_ASVs_blast_BOLD.txt -num_threads 20 -max_target_seqs 10 -perc_identity 70 -qcov_hsp_perc 80
#scp john@bio-tr:\srv\users\Joelle\CAW-26-06\fwhF1R1\* H:/Joelle.Lousberg/CAW-26-06_Joelle_Ch2_step3_fwh/fwh_nochim_together_2.4_errors





