# This script is to plot mean DP across 9p21 for the index sample carrying the 100 kb deletion

# Load the packages
library(ggplot)
library(dplyr)

# ==============================================================================
# Prior to running this script:
# Use samtools coverage to generate the file WGS_samtools_depth.txt using NIH biowulf cluster 
# samtools depth -r chr9:19900001-25600000 -q 20 -Q 20 $bamfile -o WGS_samtools_depth.txt
# ==============================================================================

# ==============================================================================
# Make coverage plots
# ==============================================================================
# Read in the samtools coverage file
data <- read_delim("WGS_samtools_depth.txt", 
                   delim = "\t", col_names = F)

# Set up colnames
colnames(data) <- c("SampleID","Chr","Position","DP")

# Calculate mean DP and coverage across 9p21
data$mean <- mean(data$DP)
data$depth <- data$DP/data$mean

# Plot the coverage
ggplot(data,aes(Position,depth))+
    geom_line(col="blue",size=0.3)+
    labs(x='\nChromosome position', y="Scaled coverage") 
 
# Adjust the plotting region to highlight the deletion
ggplot(data[data$Position >= 22117000 & data$Position <= 22442855,],
       aes(Position,depth))+
  geom_line(col="blue",size=0.3)+
  labs(x='\nChromosome position', y="Scaled coverage") 

