# ==============================================================================
# Project: 9p21 deletion - Genoa (Italy)
# Description: Rscript for data analyses + generating Supplementary figures/Tables
# Author: Linh Bui-Raborn
# Finalized: 05/01/2026
# ==============================================================================

# ==============================================================================
# CTCF binding site analysis (Table S12, Figure S6)
# Data from: https://remap.univ-amu.fr/download_page#:~:text=hg38-,hg38,-CTCF
# Homo sapiens hg38
# ==============================================================================
# Read in the merged, non-redundant peak bed file
peakdata <- read_delim("CTCF_peaks/remap2022_nr_macs2_hg38_v1_0.bed", col_names = FALSE)
colnames(peakdata) <- c("chr","chromStart","chromEnd","name","score","strand",
                        "thickStart","thickEnd","itemRgb")

# Split the name column into TF and cell type
peakdata$TF_info <- sapply(strsplit(peakdata$name,split=":") , "[[", 1)
peakdata$Celltype <- sapply(strsplit(peakdata$name,split=":") , "[[", 2)

# Select for CTCF binding peaks only
peakdata2 <- peakdata[peakdata$TF_info == "CTCF",]
write.table(peakdata2, file = "remap2022_nr_hg38_ctcfpeaks.txt",
            row.names = FALSE, sep = "\t", quote = FALSE, col.names = FALSE)

# Filter out to keep only peaks in chr9
peakdata3 <- peakdata2[peakdata2$chr == "chr9",]

# Split the cell type column into multiple columns since some peaks are present in multiple cell types
peakdata3 <- peakdata3 %>% separate(Celltype, paste0('CT', c(1:20)), sep = ',', 
                                    remove = FALSE, extra = "merge")

# Count number of cell lines having CTCF peaks in deletion region
# Deletion coordinates: chr9:22,180,210 - 22,280,977
peakdata4 <- peakdata3[peakdata3$chromStart >= 22180210 & peakdata3$chromEnd <= 22280977,]

# Save output file
write.csv(peakdata4, 
          file = "11172025_ReMap2022_CTCFpeaks_nrpeaks_Italiandeletion.csv",
          row.names = FALSE)

# Prepare the peak file for HEK293 and H1-hESC 3D genome browser 
# biowulf: cat remap2022_nr_hg38_ctcfpeaks.txt | grep "hESC" > remap2022_nr_hg38_ctcf_hESC.bed
# biowulf: cat remap2022_nr_hg38_ctcfpeaks.txt | grep "HEK293" > remap2022_nr_hg38_ctcf_HEK293.bed
hek <- read_delim("CTCF_peaks/remap2022_nr_hg38_ctcf_HEK293.bed", col_names = FALSE)
hesc <- read_delim("CTCF_peaks/remap2022_nr_hg38_ctcf_hESC.bed", col_names = FALSE)
colnames(hek) <- c("chr","start","end","sample","score","blank","thick_start",
                   "thick_end","color_code","TF","cells")
colnames(hesc) <- c("chr","start","end","sample","score","blank","thick_start",
                    "thick_end","color_code","TF","cells")

hek$score2 <- hek$score/max(hek$score)*1000
hesc$score2 <- hesc$score/max(hesc$score)*1000

write.table(hek, file = "CTCF_peaks/remap2022_nr_hg38_ctcf_HEK293.scale.txt", 
            row.names = FALSE,col.names = FALSE, sep = "\t", quote = FALSE)
write.table(hesc, file = "CTCF_peaks/remap2022_nr_hg38_ctcf_HESC.scale.txt", 
            row.names = FALSE,col.names = FALSE, sep = "\t", quote = FALSE)

