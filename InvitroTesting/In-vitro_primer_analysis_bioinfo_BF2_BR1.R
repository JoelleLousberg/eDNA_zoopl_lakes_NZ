
library(phyloseq)
library(Rcpp)
library(dada2)
library(ShortRead)
library(Biostrings)
library(ggplot2)
library(stringr) # not strictly required but handy
library(readr)
library(seqinr)
library(data.table)
library(plyr)

###Choose your working directory
setwd("/srv/users/Joelle/CAW-25-36/")
path <- "/srv/users/Joelle/CAW-25-36/"

list.files(path)

#I only want certain samples from my list actually having the BF2/BR1 primers,
#in this case samples 21-30, 31 & 32.

# Forward reads
fnFs <- sort(list.files(
  path,
  pattern = "CAW-25-36-(02[1-9]|030|031|032).*R1_001.fastq.gz$",
  full.names = TRUE
))

# Reverse reads
fnRs <- sort(list.files(
  path,
  pattern = "CAW-25-36-(02[1-9]|030|031|032).*R2_001.fastq.gz$",
  full.names = TRUE
))

###Instead here the code that includes all files:
#fnFs <- sort(list.files(path, pattern = "R1_001.fastq.gz", full.names = TRUE))
#fnRs <- sort(list.files(path, pattern = "R2_001.fastq.gz", full.names = TRUE))

FWD <- "GCHCCDGAYATRGCHTTYCC"
REV <- "ARYATWGTRATDGCHCCHGC"

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

path.cut <- file.path(path, "BF2BR1_CO1/cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

# Trim FWD off of R1 (forward reads) - 
R1.flags <- paste0("-g", " ^", FWD) 
# Trim REV off of R2 (reverse reads)
R2.flags <- paste0("-G", " ^", REV) 
# Run Cutadapt
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c("-e 0.08 --discard-untrimmed", R1.flags, R2.flags,
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

jpeg(file="BF2BR1_CO1/BF2BR1_CO1.Quality.Plot.F.jpg",res=300, width=15, height=8, units="in")
fwd_qual_plots
dev.off()

jpeg(file="BF2BR1_CO1/BF2BR1_CO1.Quality.Plot.R.jpg",res=300, width=15, height=8, units="in")
rev_qual_plots
dev.off()



##______________________________________________________________________

filtpathF <- file.path(path.cut, "filtered", basename(cutFs))
filtpathR <- file.path(path.cut, "filtered", basename(cutRs))


out <- filterAndTrim(cutFs, filtpathF, cutRs, filtpathR,
                     truncLen=c(175,175), maxEE=c(4,6), truncQ=2, maxN=0, rm.phix=TRUE,
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


filtpathF2 <- filtpathF[-c(11, 12)] ##remove water blank - no sequences
filtpathR2 <- filtpathR[-c(11, 12)] ##remove water blank - no sequences

set.seed(100) # set seed to ensure that randomized steps are replicatable

# Learn forward error rates
library(magrittr)

loessErrfun_mod <- function (trans) {
  qq <- as.numeric(colnames(trans))
  est <- matrix(0, nrow = 0, ncol = length(qq))
  for (nti in c("A", "C", "G", "T")) {
    for (ntj in c("A", "C", "G", "T")) {
      if (nti != ntj) {
        errs <- trans[paste0(nti, "2", ntj), ]
        tot <- colSums(trans[paste0(nti, "2", c("A",
                                                "C", "G", "T")), ])
        rlogp <- log10((errs + 1)/tot)
        rlogp[is.infinite(rlogp)] <- NA
        df <- data.frame(q = qq, errs = errs, tot = tot,
                         rlogp = rlogp)
        mod.lo <- loess(rlogp ~ q, df, weights = log10(tot),span = 2)
        pred <- predict(mod.lo, qq)
        maxrli <- max(which(!is.na(pred)))
        minrli <- min(which(!is.na(pred)))
        pred[seq_along(pred) > maxrli] <- pred[[maxrli]]
        pred[seq_along(pred) < minrli] <- pred[[minrli]]
        est <- rbind(est, 10^pred)
      }
    }
  }
  MAX_ERROR_RATE <- 0.25
  MIN_ERROR_RATE <- 1e-07
  est[est > MAX_ERROR_RATE] <- MAX_ERROR_RATE
  est[est < MIN_ERROR_RATE] <- MIN_ERROR_RATE
  #Monotonicity
  # estorig <- est
  # est <- est %>%
  #   data.frame() %>%
  #   mutate_all(funs(case_when(. < X40 ~ X40,
  #          . >= X40 ~ .))) %>% as.matrix()
  #  rownames(est) <- rownames(estorig)
  #  colnames(est) <- colnames(estorig)
  #
  err <- rbind(1 - colSums(est[1:3, ]), est[1:3, ], est[4,
  ], 1 - colSums(est[4:6, ]), est[5:6, ], est[7:8, ], 1 -
    colSums(est[7:9, ]), est[9, ], est[10:12, ], 1 - colSums(est[10:12,
    ]))
  rownames(err) <- paste0(rep(c("A", "C", "G", "T"), each = 4),
                          "2", c("A", "C", "G", "T"))
  colnames(err) <- colnames(trans)
  return(err)
}
###
# Learn error rates

errF <- learnErrors(filtpathF2, nbases=1e8, errorEstimationFunction = loessErrfun_mod, multithread=TRUE, verbose = TRUE)

errR <- learnErrors(filtpathR2, nbases=1e8, errorEstimationFunction = loessErrfun_mod, multithread=TRUE, verbose = TRUE)

###plot errors

jpeg(file="BF2BR1_CO1/BF2BR1_CO1.Error.F.Plot.jpg",res=300, width=15, height=8, units="in")
plotErrors(errF, nominalQ=TRUE)
dev.off()

jpeg(file="BF2BR1_CO1/BF2BR1_CO1.Error.R.Plot.jpg",res=300, width=15, height=8, units="in")
plotErrors(errR, nominalQ=TRUE)
dev.off()

derepF <- derepFastq(filtpathF2, verbose=TRUE)
derepR <- derepFastq(filtpathR2, verbose=TRUE)

dadaF.pseudo <- dada(derepF, err=errF, multithread=TRUE, pool="pseudo")
dadaR.pseudo <- dada(derepR, err=errR, multithread=TRUE, pool="pseudo")

mergers <- mergePairs(dadaF.pseudo, derepF, dadaR.pseudo, derepR, maxMismatch = 1, minOverlap = 10, 
                      verbose=TRUE)

#head(mergers)

getN <- function(x) sum(getUniques(x))
track <- cbind(out[-c(11, 12), ], sapply(dadaF.pseudo, getN), sapply(dadaR.pseudo, getN), sapply(mergers, getN))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged")
rownames(track) <- sample.names[-c(11, 12)]

write.csv(track,"BF2BR1_CO1/BF2BR1_CO1.track.csv")

track


seqtab <- makeSequenceTable(mergers)

saveRDS(seqtab, "BF2BR1_CO1/BF2BR1_CO1.seqtab.rds")

#seqtab <- readRDS("/srv/users/Joelle/CAW-25-36/BF2BR1_CO1/BF2BR1_CO1.seqtab.rds")


#######################################################################################################################################################


####code ends

###Look at the length distribution of the sequences

trimtable <- as.data.frame(table(nchar(getSequences(seqtab))))
colnames(trimtable) <- c("Length.bp", "Frequency")
trimtable$Frequency <- as.numeric(trimtable$Frequency)
str(trimtable)

library(ggplot2)
g <- ggplot(trimtable, aes(x = Length.bp, y = Frequency)) +
  geom_bar(stat="identity")


jpeg(file="/srv/users/Joelle/CAW-25-36/BF2BR1_CO1/BF2BR1_CO1.Chim.Dist.Plot.jpg",res=300, width=30, height=8, units="in")
g
dev.off()


##checking for chimeras (sequences outside the expected size range) Look at table to figure out spread (min,max for amplicon size)

seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% seq(322, 322)]

seqtab.nochim <- removeBimeraDenovo(seqtab2, multithread=TRUE, verbose=TRUE)
saveRDS(seqtab.nochim, "/srv/users/Joelle/CAW-25-36/BF2BR1_CO1/BF2BR1_CO1.seqtab.nochim.rds")
#seqtab <- readRDS("BF2BR1_CO1/BF2BR1_CO1.seqtab.nonchim.rds")

###Produce a table which shows the number or reads at each stage
getN <- function(x) sum(getUniques(x))
track2 <- cbind(out[-c(11, 12), ],sapply(dadaF.pseudo, getN), sapply(dadaR.pseudo, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track2) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track2) <- sample.names[-c(11, 12)]

write.csv(track2, "/srv/users/Joelle/CAW-25-36/BF2BR1_CO1/BF2BR1_CO1.track2.csv")  

track2
###Do only Blast/Megan taxonomy. 

###turn rds into file that can be used for blast

seqtab.nochim <- readRDS("/srv/users/Joelle/CAW-25-36/BF2BR1_CO1/BF2BR1_CO1.seqtab.nochim.rds")

asv_seqs <- colnames(seqtab.nochim)
head(asv_seqs)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

# making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
head(asv_fasta)
write.table(asv_fasta, "/srv/users/Joelle/CAW-25-36/BF2BR1_CO1/BF2BR1_CO1_ASV_seq_IDs_and_sequences.tsv", sep="\t", quote=F, col.names=NA)
write(asv_fasta, "/srv/users/Joelle/CAW-25-36/BF2BR1_CO1/BF2BR1_CO1.nochim_ASVs.fa")

# count table:
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "/srv/users/Joelle/CAW-25-36/BF2BR1_CO1/BF2BR1_CO1.nochim_ASVs_counts.tsv", sep="\t", quote=F, col.names=NA)

##blast code
# blastn -query BF2BR1_CO1.nochim_ASVs.fa -db ~/../../srv/referencedb/COI/COI.NCBI.BOLD.blast.db -out BF2BR1_CO1.nochim_ASVs_blast_BOLD.txt -num_threads 20 -max_target_seqs 10 -perc_identity 70 -qcov_hsp_perc 80
