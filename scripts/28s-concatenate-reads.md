# Concatenate 28S pipeline


- [1. Load libraries](#1-load-libraries)
- [2. Define path to fastq files](#2-define-path-to-fastq-files)
- [3. Generate matched lists of the forward and reverse read files](#3-generate-matched-lists-of-the-forward-and-reverse-read-files)
- [4. Identify Primers](#4-identify-primers)
- [5. Verify the presence and orientation of the primers in the data](#5-verify-the-presence-and-orientation-of-the-primers-in-the-data)
  * [5.1 Create all orientantions of the primers](#51-create-all-orientantions-of-the-primers)
  * [5.1.1 Create all orientations of the input sequence](#511-create-all-orientations-of-the-input-sequence)
  * [5.2. Count the number of times the primers appear in the forward and reverse reads](#52-count-the-number-of-times-the-primers-appear-in-the-forward-and-reverse-reads)
- [5.3. Remove primers](#53-remove-primers)
    + [5.3.1 Call cutadapt in R-](#531-call-cutadapt-in-r-)
    + [5.3.2 Ceate output filenames for the cutadapt-ed files. Define the parameters of cutadapt](#532-ceate-output-filenames-for-the-cutadapt-ed-files-define-the-parameters-of-cutadapt)
    + [5.3.3 Cutadapt parameters](#533-cutadapt-parameters)
    + [5.3.4 Run cutadapt without filtering Ns](#534-run-cutadapt-without-filtering-ns)
    + [5.3.5 Sanity check. Look for primers in cutadapt-ed files (should be 0)](#535-sanity-check-look-for-primers-in-cutadapt-ed-files--should-be-0-)
- [6. DADA2 PIPELINE](#6-dada2-pipeline)
    + [1.0 Read in the names of the cutadapt-ed FASTQ](#10-read-in-the-names-of-the-cutadapt-ed-fastq)
    + [2.0 Inspect read quality profiles](#20-inspect-read-quality-profiles)
    + [3.0 Filter and trim](#30-filter-and-trim)
    + [4.0  Learn error rates](#40--learn-error-rates)
    + [5.0 Plot error rates](#50-plot-error-rates)
    + [6.0 Settings for trimming and merging](#60-settings-for-trimming-and-merging)
    + [7.0 Construct sequence (ASV) table](#70-construct-sequence--asv--table)
    + [8.0 Remove chimeras](#80-remove-chimeras)
    + [9.0 Last checks: track reads throughout the DADA2 pipeline](#90-last-checks--track-reads-throughout-the-dada2-pipeline)
    + [10. Write results to tsv file](#10-write-results-to-tsv-file)
- [7.0 Importing sequences to qiime2 for taxonomy assignment](#70-importing-sequences-to-qiime2-for-taxonomy-assignment)




## 1. Load libraries
```R
library(dada2)
library(ShortRead)
library(Biostrings)
library(biomformat)
library(qiime2R) 
```



## 2. Define path to fastq files

```R
path_seqs <- "/grp/animalbehaviour/microbiome/buzzard-microbiome/28s/seqs/"  ## CHANGE to the directory containing the fastq files.
list.files(path_seqs) 
```



## 3. Generate matched lists of the forward and reverse read files

```R
fnFs <- sort(list.files(path_seqs, pattern = "_R1_001.fastq.gz", full.names = TRUE)) # Just select forward read files
fnRs <- sort(list.files(path_seqs, pattern = "_R2_001.fastq.gz", full.names = TRUE)) # Just select reverse read files
```




## 4. Identify Primers
```R
FWD <- "GTAACTTCGGGAWAAGGATTGGCT" ## INPUT forward primer sequence
REV <- "AGAGTCAARCTCAACAGGGTCTT"  ## INPUT reverse primer sequence
```



## 5. Verify the presence and orientation of the primers in the data



### 5.1 Create all orientantions of the primers

```R
allorients <- function(primer) {
```

  ### 5.1.1 Create all orientations of the input sequence
```R
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allorients(FWD)
REV.orients <- allorients(REV)
FWD.orients
REV.orients
```

### 5.2. Count the number of times the primers appear in the forward and reverse reads

```R
#Considering all possible primer orientations. Identifying and counting the primers on one set of paired end FASTQ files is sufficient, #assuming all the files were created using the same library preparation, so we’ll just process the first sample.

#Counts number of reads in which the primer is found

primerHits <- function(primer, fn) {
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs[[20]]), # not filtering N at the moment 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs[[20]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs[[20]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs[[20]]))

#As expected, the FWD primer is found in the forward reads in its forward orientation, 
#and in some of the reverse reads in its reverse-complement orientation (due to read-through because 18S is short). 
#Similarly the REV primer is found with its expected orientations.

#Note: Orientation mixups are a common trip-up. If, for example, the REV primer is matching the Reverse reads 
#in its RevComp orientation, then replace REV with its reverse-complement orientation 

#(REV <- REV.orient[["RevComp"]]) before proceeding.
```



## 5.3. Remove primers



#### 5.3.1 Call cutadapt in R-

```R
cutadapt <- "/homes/hugoeira/.local/bin/cutadapt" # CHANGE to cutadapt path 
system2(cutadapt, args = "--version") # Run shell commands from R
```



#### 5.3.2 Ceate output filenames for the cutadapt-ed files. Define the parameters of cutadapt

```R
#Create output files

path.cut <- file.path(path_seqs, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs)) 
```



#### 5.3.3 Cutadapt parameters

```R
FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)

#Trim FWD and the reverse-complement of REV off of R1 (forward reads)

R1.flags <- paste("-g", FWD, "-a", REV.RC) 

#Trim REV and the reverse-complement of FWD off of R2 (reverse reads)

R2.flags <- paste("-G", REV, "-A", FWD.RC) 
```



#### 5.3.4 Run cutadapt without filtering Ns
```R
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2,"-m 50", "--match-read-wildcards", "--discard-untrimmed", # -n 2 required to remove FWD and REV from reads
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs[i], fnRs[i])) # input files
}
```




#### 5.3.5 Sanity check. Look for primers in cutadapt-ed files (should be 0)
```R
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[20]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[20]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[20]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[20]]))
```





## 6. DADA2 PIPELINE



#### 1.0 Read in the names of the cutadapt-ed FASTQ

```R
#Forward and reverse fastq filenames with same format

cutFs <- sort(list.files(path.cut, pattern = "_R1_001.fastq.gz", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "_R2_001.fastq.gz", full.names = TRUE))


#Extract sample names, assuming filenames have format

get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(cutFs, get.sample.name))
head(sample.names)
```



#### 2.0 Inspect read quality profiles

```R
plotQualityProfile(cutFs, aggregate = TRUE)
plotQualityProfile(cutRs, aggregate = TRUE) 
```



#### 3.0 Filter and trim

```R
Assigning the filenames for the output of the filtered reads to be stored as fastq.gz files

filtFs <- file.path(path.cut, "28S_filtered", basename(cutFs))
filtRs <- file.path(path.cut, "28S_filtered", basename(cutRs))

#Filter and trim parameters. 
#Used standard filtering parameters. Enforced a minLen here, to get rid of spurious very low-length sequences. 
#This was not needed in the 16S workflow because truncLen already served that purpose

filtered_out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs,truncLen =c(235,185), maxN = 0, maxEE = c(2, 2), 
                              truncQ = 2, minLen = 100, rm.phix = TRUE, compress = TRUE, multithread = TRUE)
head(filtered_out) 

# After filtering, check the quality of the filtered files to verify that we have solved the issues

plotQualityProfile(filtFs, aggregate = TRUE)
plotQualityProfile(filtRs[18:19])
```




#### 4.0  Learn error rates
```R
set.seed(12345)

path_filter_seqs <- "/grp/animalbehaviour/microbiome/buzzard-microbiome/28s/seqs/cutadapt/28S_filtered/"

filtFs <- sort(list.files(path_filter_seqs, pattern = "_R1_001.fastq.gz", full.names = TRUE)) # Just select forward read files
filtRs <- sort(list.files(path_filter_seqs, pattern = "_R2_001.fastq.gz", full.names = TRUE)) # Just select reverse read files

#need to give a new path to filtered seqs because some were excluded during
#filtering process. If not learn error gives an error

# Extract sample names again because new path was given

get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(filtFs, get.sample.name))
head(sample.names)

#Learn forward error rates
errF <- learnErrors(filtFs, nbases=1e8, multithread=TRUE)
#Learn reverse error rates
errR <- learnErrors(filtRs, nbases=1e8, multithread=TRUE)
```



#### 5.0 Plot error rates

```R
plotErrors(errF, nominalQ=TRUE)


dadaFs <- dada(filtFs, # the function will use the filtered files
               err=errF, # we pass on the calculated error rates
               multithread=TRUE # change if your computer does not allow multithreading
)
dadaRs <- dada(filtRs, # the function will use the filtered files
               err=errR, # we pass on the calculated error rates
               multithread=TRUE # change if your computer does not allow multithreading
)
```




#### 6.0 Settings for trimming and merging
```R
truncLen_r1 = 270
truncLen_r2 = 220
minOverlap <- 10

# merge denoised reads, keeping info about unmerged sequences in output

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, minOverlap = minOverlap, verbose = TRUE, returnRejects = TRUE)
concat <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, minOverlap = minOverlap, verbose = TRUE,  justConcatenate=TRUE)
concat$sequence <- gsub("NNNNNNNNNN", "", concat2$sequence)


#Inspect the merger data.frame from the first sample
names(concat) <- sample.names
head(concat[[1]])
```




#### 7.0 Construct sequence (ASV) table
```R
ASVtab <- makeSequenceTable(concat)
dim(ASVtab)

# Inspect distribution of sequence lengths

table(nchar(getSequences(ASVtab)))
```



#### 8.0 Remove chimeras

```R
ASVtab.nochim <- removeBimeraDenovo(ASVtab, method="consensus", multithread=TRUE, verbose=TRUE)

dim(ASVtab.nochim) # non chimera sequences

sum(ASVtab.nochim)/sum(ASVtab) #percentage of reads after chimera filtering

table(nchar(getSequences(ASVtab.nochim)))
```



#### 9.0 Last checks: track reads throughout the DADA2 pipeline

```R
filtered_out <- filtered_out[-c(2),] # use if some samples did not pass the filter. Remove those samples

#sample.names <- unname(sapply(filtFs, get.sample.name)) # need to get new sample names in case some samples were excluded

getN <- function(x) sum(getUniques(x))

track <- cbind(filtered_out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(concat, 

getN),rowSums(ASVtab.nochim),round(rowSums((ASVtab.nochim)/filtered_out[,1]*100, 2)))

colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim", "perc-reads-retained")

rownames(track) <- sample.names

head(track)
```



#### 10. Write results to tsv file

```R
write.table(track, "28SS-denoising-tracking.tsv", quote=FALSE, sep="\t", col.names=NA)
```




## 7.0 Importing sequences to qiime2 for taxonomy assignment

```R
# Giving our seq headers more manageable names (ASV_1, ASV_2...)

ASV_seqs <- colnames(ASVtab.nochim)
ASV_headers <- vector(dim(ASVtab.nochim)[2], mode="character")

for (i in 1:dim(ASVtab.nochim)[2]) {
  ASV_headers[i] <- paste(">ASV", i, sep="_")
}


# Making and writing out a fasta of our final ASV seqs

ASV_fasta <- c(rbind(ASV_headers, ASV_seqs))

write(ASV_fasta, "28S_ASVs_nochim.fasta")


# Export otu table as biom file

biom_table<- make_biom(ASVtab)

write_biom(biom_table, "table.biom")

```


