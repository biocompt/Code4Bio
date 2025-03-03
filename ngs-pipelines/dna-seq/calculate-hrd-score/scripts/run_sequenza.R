library(devtools)
library(sequenza)
library(copynumber)
library(data.table)

install_github('sztup/scarHRD')
library(scarHRD)

args <- commandArgs(trailingOnly = TRUE)
seqz_file <- args[1]
score_file <- args[2]

# Extract sample_id
sample_id <- basename(dirname(seqz_file))

# Files and directories 
out_dir <- dirname(seqz_file)
output_file <- paste0(out_dir, "/", sample_id, "_calculate_HRD.txt")

# Load and process sequenza files
Sys.setenv(VROOM_CONNECTION_SIZE = 1000000000)
source("utilities/scripts/preprocess.seqz.R")
seg <- preprocess.seqz(seqz_file)
write.table(seg, output_file, col.names = T, row.names = F, sep = "\t", quote = F)

# Run scarHRD
score <- scar_score(output_file, reference = "grch38", seqz = F, ploidy = T)
score_sample <- cbind(Sample_ID = sample_id, score)
write.table(score_sample, score_file, quote = F, row.names = F, col.names = T, sep = "\t")
