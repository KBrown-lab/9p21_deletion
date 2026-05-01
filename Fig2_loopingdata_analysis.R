# ==============================================================================
# This script is used to prepare the HiChIP peak-to-peak files for bedtool intersect
# and extract out the loops that are in the 9p21 deletion region that connects to 
# p16, p14, CDKN2B, MTAP and DMRTA1
# Note:
# hiChIP fragment size is 2.5kb, output file used center coordinate -> +/-1.25kb from each side
# ==============================================================================

# Set up the working environment
library(ggplot2)
library(tidyverse)

# Read in the data
hifile <- read.delim("Melanocyte_HiChIP_noUV.interactions_FitHiC_Q0.01_WashU.bed",
                     header = FALSE) 
colnames(hifile) <- c("chr","start","end","loop_info","id","empty")

# Split the score column
hifile$scores <- sapply(hifile$loop_info, function(x) strsplit(x,",")[[1]][2])
hifile$scores <- as.numeric(hifile$scores)

# Subset to only keep chr9
hifile2 <- hifile[hifile$chr == "chr9",]

# Split the loop_info to get the coordinate for the loops
hifile2$n1 <- sapply(hifile2$loop_info, function(x) strsplit(x,"-")[[1]][1])
hifile2$n2 <- sapply(hifile2$loop_info, function(x) strsplit(x,"-")[[1]][2])
hifile2$target_start <- sapply(hifile2$n1, function(x) strsplit(x,":")[[1]][2])
hifile2$target_end <- sapply(hifile2$n2, function(x) strsplit(x,",")[[1]][1])
hifile2$n1 <- NULL
hifile2$n2 <- NULL

# Generate another file to use other end of loops as start/end
onion <- hifile2[,c(1,8,9)]
onion$loop_info <- paste0(hifile2$chr,":",hifile2$start,"-",hifile2$end,",",hifile2$scores)
onion$id <- hifile2$id
onion$empty <- "."
onion$scores <- hifile2$scores
colnames(onion) <- c("chr","start","end","loop_info","id","empty","scores")
onion$target_start <- hifile2$start
onion$target_end <- hifile2$end

# Write out files for bedtools intersect
write.table(onion, file = "hichIP_direction1_forbedtools.txt",
            row.names = FALSE,col.names = FALSE, quote = FALSE, sep = "\t")
write.table(hifile2, file = "hichIP_direction2_forbedtools.txt",
            row.names = FALSE,col.names = FALSE, quote = FALSE, sep = "\t")

# Use bedtools intersect to find overlapping between loop data and promoter data (biowulf)
# 2 steps bedtools intersect:
# 1. Intersect to extract loops overlapping with deletion region
# 2. Intersect to extract loops from step 1 that overlap gene promoters
# =============================================
# Extracting loops for p14, p16, CDKN2B, MTAP, DMRTA1 - Melanocyte promoters
# =============================================
# Read in the bedtools intersect output files
files <- list.files(path = "Results/Melanocytes", pattern="*.txt",
                    full.names = TRUE)

infile1 <- lapply(files, function(x) read_delim(x, col_names = FALSE))

# Remove 1 hichip file (no values)
infile1 <- infile1[-4]
infile1 <- infile1[-4]

# Add names to the dataframe
names(infile1) <- c("CaptureC","hichip1","hichip2","promoterC") #hichip1 : extend coordinates to +/-1.25kb, hichip2: extend coordinates to +/-2.5kb
# There is no difference in hichip1 and hichip2, remove one file to avoid writing out 2 identical files
infile1 <- infile1[-3]

# Remove 2 columns showing coordinate of other end of loops for CaptureC and hiChip, this was used for intersecting with promoters
infile1[[1]] <- infile1[[1]][,c(-2,-3)]
infile1[[2]] <- infile1[[2]][,c(-2,-3)]

# Rename the columns
new_names <- c("chr","start","end","loop_info","id","empty","loop_score","group",
               "chr.2","promoter_start","promoter_end",
               "gene_info","score","strand","bp_overlap")
infile1 <- lapply(infile1, setNames, new_names)

# Add in columns for geneid and transcriptid
infile1 <- lapply(infile1, function(x) x = cbind(x,
                                                 geneid = sapply(strsplit(x$gene_info,split="-") , "[[", 2),
                                                 transcriptid = sapply(strsplit(x$gene_info,split="-") , "[[", 4)))

# Remove quotes "" for the transcriptid column
infile1 <- lapply(infile1, function(x) x = cbind(x,
                                                 transcriptid2 = gsub("\"", "",x$transcriptid),
                                                 dataset = "Melanocyte Promoters"))

# Filter to keep only loops connect from deletion to p14, p16, CDKN2B, MTAP and DMRTA1
transcripts <- c("ENST00000304494.5", #p16
                 "ENST00000579755.1", #p14
                 "ENST00000428597", #CDKN2B-AS
                 "ENST00000276925.6", #CDKN2B p15
                 "ENST00000580900.1", #MTAP
                 "ENST00000325870.2", #DMRTA1
                 "ENST00000441769.2") #CDKN2A-AS (C9orf53)
outfile <- lapply(infile1, function(x) x=x[x$transcriptid2 %in% transcripts & 
                                             x$group == "Italian_deletion",])

# Write output as txt files
for (name in names(outfile)) {
  file_path <- paste0(name, "_Italiandeletion_loops_to_melanocytes_promoter.txt")
  write.table(outfile[[name]], file = file_path, 
              row.names = FALSE,col.names = TRUE, quote = FALSE, sep = "\t")
}


# Write output files for WashU epigenome browser
for (name in names(outfile)) {
  file_path <- paste0(name, "_Italiandeletion_loops_to_melanocytes_promoter_WashU.txt")
  write.table(outfile[[name]][,1:6], file = file_path, 
              row.names = FALSE,col.names = FALSE, quote = FALSE, sep = "\t")
}

# =============================================
# Extracting loops for p14, p16, CDKN2B, MTAP, DMRTA1 - Melanoma promoters
# =============================================
# Read in the bedtools intersect output files
files <- list.files(path = "Results/Melanoma/", pattern="*.txt",
                    full.names = TRUE)

infile1 <- lapply(files, function(x) read_delim(x, col_names = FALSE))

# Remove 1 hichip file (no values)
infile1 <- infile1[-3]

# Add names to the dataframe
names(infile1) <- c("CaptureC","hichip","promoterC")

# Remove 2 columns showing coordinate of other end of loops, this was used for intersecting with promoters
infile1[[1]] <- infile1[[1]][,c(-2,-3)]
infile1[[2]] <- infile1[[2]][,c(-2,-3)]

# Rename the columns
new_names <- c("chr","start","end","loop_info","id","empty","loop_score","group",
               "chr.2","promoter_start","promoter_end",
               "gene_info","score","strand","bp_overlap")
infile1 <- lapply(infile1, setNames, new_names)

# Add in columns for geneid and transcriptid
infile1 <- lapply(infile1, function(x) x = cbind(x,
                                                 geneid = sapply(strsplit(x$gene_info,split="-") , "[[", 2),
                                                 transcriptid = sapply(strsplit(x$gene_info,split="-") , "[[", 4)))

# Remove quotes "" for the transcriptid column
infile1 <- lapply(infile1, function(x) x = cbind(x,
                                                 transcriptid2 = gsub("\"", "",x$transcriptid),
                                                 dataset = "Melanoma Promoters"))

outfile <- lapply(infile1, function(x) x=x[x$transcriptid2 %in% transcripts & 
                                             x$group == "Italian_deletion",])

# Write output as txt files
for (name in names(outfile)) {
  file_path <- paste0(name, "_Italiandeletion_loops_to_promoter_melanoma.txt")
  write.table(outfile[[name]], file = file_path, 
              row.names = FALSE,col.names = TRUE, quote = FALSE, sep = "\t")
}

# Write output files for WashU epigenome browser
for (name in names(outfile)) {
  file_path <- paste0(name, "_Italiandeletion_loops_to_promoter__melanoma_WashU.txt")
  write.table(outfile[[name]][,1:6], file = file_path, 
              row.names = FALSE,col.names = FALSE, quote = FALSE, sep = "\t")
}
