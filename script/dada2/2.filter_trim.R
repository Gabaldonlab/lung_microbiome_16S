# load the needed libraries
library(dada2); packageVersion("dada2")
library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")
library(ShortRead)

path <- "path/to/files"
output_path <- "results/output/directory"

#list.files(path)
print("path ok")
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
print("fnFs ok")
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
print("fnRs ok")
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

sample.names.pre  <- sample.names[ grepl("-PRE", sample.names) ]
sample.names.post <- sample.names[ grepl("-POST", sample.names) ]
# get filenames for PRE and POST separately
fnFs.pre  <- fnFs[ grepl("-PRE", sample.names) ]
fnFs.post <- fnFs[ grepl("-POST", sample.names) ]
fnRs.pre  <- fnRs[ grepl("-PRE", sample.names) ]
fnRs.post <- fnRs[ grepl("-POST", sample.names) ]

filtFs <- file.path(path, "filtered3", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered3", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                     truncLen=c(270,270), trimLeft=c(10,10),
                     maxEE=c(2,2),
                     truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)
head(out)

mean(out[,2] / out[,1])

save.image(file.path(output_path, "filter_trim_output.RData"))
