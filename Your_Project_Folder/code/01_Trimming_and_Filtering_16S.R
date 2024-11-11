#!/usr/bin/Rscript

# - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - #

cat("This script will take in .fastq files and: 
1) remove indeterminations (Ns) 
2) trim primers used for barcode amplification
3) filter and trim low quality nucleotide reads
\n")
cat("Working directory:", getwd(), "\n\n")

# - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - #

# Libraries ----
library(Biostrings) # fast manipulation of large biological sequences or sets of sequences
library(ShortRead) # sampling, iteration, and input of FASTQ files. filtering and trimming reads, 
                   # quality assessment report.
library(dplyr) # data management
library(dada2)

# - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - #

# Set general parameters ----
amplicon <- "16S"
cat(paste("Amplicon:", amplicon), "\n")
cat("Inputfiles' path: \n")
inpath <- paste("data/01_raw", amplicon, sep = "/"); print(inpath); cat("\n")
cat("Outputfiles' path: \n")
outpath <- paste("data/02_filtering", amplicon, sep = "/"); print(outpath); cat("\n")

# PRIMERS [CHANGE ME!]
cat("Forward and Reverse PRIMERS: \n")
FWD <- "CCTACGGGNGGCWGCAG"; print(FWD)
REV <- "GACTACHVGGGTATCTAATCC"; print(REV)
cat("\n")

# - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - #

# Input data ----

# .fastq files overview
print("We are going to analyse the following files:")
print(list.files(inpath, pattern = ".fastq.gz"))
cat("\n")
print(paste("Nº of files:", length(list.files(inpath, pattern = ".fastq.gz")))) ; cat("\n")

# Save forward and reverse reads file paths
fnFs <- sort(list.files(inpath, 
                        pattern = "_R1_001.fastq.gz", # [change me!] Forward filename pattern
                        full.names = TRUE)) 
fnRs <- sort(list.files(inpath, 
                        pattern = "_R2_001.fastq.gz", # [change me!] Reverse filename pattern
                        full.names = TRUE)) 

# - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - #

# Primer Trimming ----
  
  # - # - # - # - # - # - # - # - # - # 
  ## Find the primers in our reads ----
  # - # - # - # - # - # - # - # - # - # 

allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # Biostrings works w/ DNAString objects 
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna),
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}

FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
print("All possible orientations for forward and reverse primers:")
print(FWD.orients)
print(REV.orients)
cat("\n")

# We are now ready to count the number of times the primers appear in the 
# forward and reverse read, considering all possible primer orientations
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
cat("Scanning the samples for all possible orientations of the PCR primers:")
print(rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs[[1]]), 
            FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs[[1]]), 
            REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs[[1]]), 
            REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs[[1]]))
)
cat("\n")

  # - # - # - # - # - # - # - # - # - # 
  ## Primer removal w/ CUTADAPT ----
  # - # - # - # - # - # - # - # - # - #

# Specify CUTADAPT program location on computer 
# IF not installed, install in desired conda environment
cutadapt <- "/home/ecosystems/miniconda3/envs/metatax/bin/cutadapt"
print("Cutadapt version:")
system2(cutadapt, args = "--version")

# Specify CUTADAPT output directory
cutadapt.path <- paste(outpath, "cutadapt", sep = "/"); dir.create(cutadapt.path)
cat("The CUTADAPT output will go to the following directory:", cutadapt.path)
# File's paths
fnFs.cut <- file.path(cutadapt.path, basename(fnFs))
fnRs.cut <- file.path(cutadapt.path, basename(fnRs))
# Summary path
cutadapt_summary <-paste(cutadapt.path, "cutadapt_summary.txt", sep = "/") 
file.create(cutadapt_summary) # this creates an empty summary file

# create reverse complement sequences of primers to input as cutadapt arguments
FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
# Paste FWD and the reverse-complement of REV in R1 flags to remove 
R1.flags <- paste("-g", FWD, "-a", REV.RC) 
# Paste REV and the reverse-complement of FWD in R2 flags to remove
R2.flags <- paste("-G", REV, "-A", FWD.RC)

print("Running CUTADAPT (don't terminate script)")
# Run Cutadapt (runs outside of R (called by system2 command))
for(i in seq_along(fnFs)) {
  cmd_output <- system2(cutadapt, 
          args = c(R1.flags, # fwd primer & reverse-complement of rev primer
                   R2.flags, # rev primer & reverse-complement of fwd primer
                   "-n", 4, # -n number of primers to search for and trim in a read,
                            # before going on with next read
                   "-o", fnFs.cut[i],
                   "-p", fnRs.cut[i], # output file's PATH
                   fnFs[i],
                   fnRs[i], #input files
                   "--minimum-length", 10,
                   "--max-n", 0, # [IMPORTANT] Remove indeterminations (Ns)
                                 # Ns accumulate at reads' 5' and 3' ends, right 
                                 # where primers are located. This makes primer
                                 # trimming difficult to achieve
                   "--cores", 6), # nº of computing cores
          stdout = TRUE,
          stderr = TRUE)
  # Append command output to summary file
  cat(cmd_output, file = cutadapt_summary, append = TRUE, sep = "\n\n")
}

print("Check: CUTADAPT done") ; cat("\n")

  # - # - # - # - # - # - # - # - # - #
  ## Check remaining primers ----
  # - # - # - # - # - # - # - # - # - #

cat("Sanity check: look for presence of primers in our reads, that escaped the CUT: \n")
# Sanity check: we will count the presence of primers in the first cutadapt-ed sample:
print(rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
            FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]), 
            REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]), 
            REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]])) 
)

print("Yay! Data should be clean of primers!"); cat("\n")
cat("RECOMMENDED: Quality check with FASTQC + MULTIQC on the filtered reads under
    outputdirectory/cutadapt \n\n\n")

# - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - #

# Quality Trimming Filtering ----

#| Once the PRIMERS are removed, it is time to filter out low-quality nucleotide
#| reads from .fastq sequences. Deciding filtering parameters is tricky and 
#| requires finding a balance between obtainning high quality, trustable 
#| sequences and enough NUMBER and LENGTH of READS. The length of combined
#| forwards and reverse reads should sum up to the barcode's LENGTH 
#| (amplicon length - primers' lengths)

#| i.e.: Target region: 16S 341F-805R (length ~ 459). The sum of fwd & rev read
#| length should be ~ 460 MINUS the length of primers (~ 420) 

# Output directory for filtered .fastq files
filtered.path <- paste(outpath, "filterAndTrim", sep = "/"); dir.create(filtered.path)

filtFs <- file.path(filtered.path, basename(fnFs.cut))
filtRs <- file.path(filtered.path, basename(fnRs.cut))

  # - # - # - # - # - # - # 
  ## Check read QUALITY ---- 
  # - # - # - # - # - # - #

# Dada2 command for Read Quality Check
#plotQualityProfile(cutFs)
#plotQualityProfile(cutRs)

#| PERSONALY, these plots read poorly and it becomes increasingly difficult to 
#| assess read quality with increasing number of samples. To decide the filtering
#| parameters, use of bash CLI tools FASTQC + MULTIQC on input reads seems more
#| appropriate 

#| Go check cutadapted files QA: 
#| /data/microbial_community/QualityAnalysis/

#| Amplicon LENGTH: 341F-805R ==> 464 bps (+20bp overlap)
#| Fwd and Rev reads should be at least 238 bps long (mean between the two)
#| Fwd reads have better quality at the 3' end than Rev reads. We will trimm 
#| Rev reads more severely than Fwd ones.  

  # - # - # - # - # - #  
  ## Filter & Trimm ---- 
  # - # - # - # - # - #

# We’ll use standard filtering parameters: maxN=0 (DADA2 requires no Ns), 
# truncQ=2, rm.phix=TRUE and maxEE= ?. 
# The maxEE parameter sets the maximum number of “expected errors” allowed in a 
# read, which is a better filter than simply averaging quality scores.

cat("Filtering and trimming of low quality reads...")

out <- filterAndTrim(fnFs.cut, filtFs,
                     fnRs.cut, filtRs,
                     maxN = 0, # a MUST for DADA2: NO indeterminations (Ns)
                     maxEE = c(4, 6), 
                     truncQ = 2, 
                     minLen = 100, 
                     trimLeft = 5, # remove 10 first bps (low Quality in QA)
                     truncLen = c(267, 232), # [CHANGE ME according to QA
                                           # & length of amplicon region]
                     rm.phix = TRUE, 
                     compress = TRUE, 
                     multithread = 6, # bool or int with number of threads
                     matchIDs = TRUE,
                     verbose = TRUE
)

cat("Check: Quality filtering DONE")

print(out)

cat("Check final filtered reads QUALITY with FastQC + MultiQC.\n","Files in directory:", filtered.path)

# - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - #

print(sessionInfo())

rm(list = ls())
