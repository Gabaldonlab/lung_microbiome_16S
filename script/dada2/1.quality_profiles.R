library(dada2); packageVersion("dada2")
library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")
library(ShortRead)

path <- "path/to/files"

list.files(path)

fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))

sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

sample.names.pre  <- sample.names[ grepl("-PRE", sample.names) ]
sample.names.post <- sample.names[ grepl("-POST", sample.names) ]
# get filenames for PRE and POST separately
fnFs.pre  <- fnFs[ grepl("-PRE", sample.names) ]
fnFs.post <- fnFs[ grepl("-POST", sample.names) ]
fnRs.pre  <- fnRs[ grepl("-PRE", sample.names) ]
fnRs.post <- fnRs[ grepl("-POST", sample.names) ]

png("quality_prof_F_26-55.png",
    width =
      465,
    height
    = 225, units='mm', res = 300)
plotQualityProfile(fnFs[25:55]) + scale_x_continuous(breaks = seq(0, 300, 25)) +
  scale_y_continuous(breaks = seq(0, 40, 5)) + theme_grey(base_size = 10)
dev.off()

png("quality_prof_R_26-55.png",
    width = 465,
    height = 225, units='mm', res = 300)
plotQualityProfile(fnRs[26:55]) + scale_x_continuous(breaks = seq(0, 300, 25)) +
  scale_y_continuous(breaks = seq(0, 40, 5)) + theme_grey(base_size = 10)
dev.off()
