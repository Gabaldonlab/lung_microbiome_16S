# load the needed libraries
library(dada2); packageVersion("dada2")            # sequence processing
library(phyloseq); packageVersion("phyloseq")      # data visualization and analysis
library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")        # data visualization and analysis
library(ShortRead)

path <- "path/to/files"

# list all the files in the filtered directory
all_filt_files <- list.files(sprintf("%s/filtered3", path), full.names = TRUE)

# Learn the errors
errF <- learnErrors(filtFs, multithread=TRUE)
print(paste0("errF:", errF))
errR <- learnErrors(filtRs, multithread=TRUE)
print(paste0("errR:", errR))

print("Errs learned")

png(sprintf("%s/results/error_rate_F_filt3.png", path), width = 465, height = 225, units='mm', res = 300)
plotErrors(errF, nominalQ = TRUE)
dev.off()

png(sprintf("%s/results/error_rate_R_filt3.png", path), width = 465, height =
      225, units='mm', res = 300)
plotErrors(errR, nominalQ = TRUE)
dev.off()

dada2:::checkConvergence(errF)
dada2:::checkConvergence(errR)

# Sample inference and merger of paired-end reads
mergers <- vector("list", length(sample.names))
names(mergers) <- sample.names

for(sam in sample.names) {
  derepF <- derepFastq(filtFs[ startsWith(basename(filtFs), sprintf("%s_F",sam)) ], verbose = T)
  ddF <- dada(derepF, err=errF, multithread=TRUE)
  
  derepR <- derepFastq(filtRs[ startsWith(basename(filtRs), sprintf("%s_R",sam)) ], verbose = T)
  ddR <- dada(derepR, err=errR, multithread=TRUE)
  
  merger <- mergePairs(ddF, derepF, ddR, derepR)
  mergers[[sam]] <- merger
}
rm(derepF); rm(derepR)

# Sample inference
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
head(mergers[[1]])

print(paste("mergers:", mergers))

## Construct sequence table
seqtab <- makeSequenceTable(mergers)
# save the sequence table
saveRDS(seqtab, sprintf("%s/results/seqtab.rds", path))
dim(seqtab)

save.image(sprintf("%s/results/seqtab_built.RData", path))

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

st.all <- seqtab

save.image(sprintf("%s/results/before_remove_chimeras_step.RData", path))

st.nochim <- removeBimeraDenovo(st.all,
                                method="consensus",
                                multithread=TRUE,
                                verbose=TRUE)
dim(st.all)
dim(st.nochim)
sum(st.nochim)/sum(st.all)

saveRDS(st.nochim, sprintf("%s/results/seqtab_nochim.rds", path))

save.image(sprintf("%s/results/chimeras_removed.RData", path))
print("chimeras have been removed")

## Taxonomy assignment
tax_spec <-
  assignTaxonomy(st.nochim,"silva_DB/silva_nr99_v138.1_wSpecies_train_set.fa.gz", # path to SILVA db
                 multithread = TRUE)

tax_spec[ grepl('\\.',tax_spec) ] <- gsub('\\.','_',tax_spec[grepl('\\.',tax_spec)])

saveRDS(tax_spec, sprintf("%s/results/tax_species_final.rds", path))

save.image(sprintf("%s/results/taxa_assigned.RData", path))
print("taxa have been assigned")

# Taxonomic assignment inspection
ts.classified <- matrix(0, nrow=6, ncol=3)
rownames(ts.classified) <- c("Phylum","Class","Order","Family","Genus","Species")
colnames(ts.classified) <- c("Classified","Unclassified","Percent_Classified")

for (tl in rownames(ts.classified)) {
  ts.classified[ tl, "Classified"] <- sum( ! is.na(tax_spec[,tl]))
  ts.classified[ tl, "Unclassified"] <- sum( is.na(tax_spec[,tl]))
  ts.classified[ tl, "Percent_Classified"] <- round(100 * sum( ! is.na(tax_spec[,tl])) /
                                                      nrow(tax_spec), 2)
}

ts.classified

## Get ASV table at the species level
tax.fixNames.spec <- matrix("", nrow = nrow(tax_spec), ncol = ncol(tax_spec))
rownames(tax.fixNames.spec) <- rownames(tax_spec)
colnames(tax.fixNames.spec) <- colnames(tax_spec)
for (tl in c("Kingdom","Phylum","Class","Order","Family","Genus","Species")) {
  tax.fixNames.spec[, tl] <- ifelse(is.na(tax_spec[,tl]), "unclassified", tax_spec[,tl])
}

# combine/adjust genus names when appropriate
tax.fixNames.spec[ tax.fixNames.spec[,"Genus"] %in% c("Clostridium_sensu_stricto_1"), "Genus"] <-
  "Clostridium"
tax.fixNames.spec[ tax.fixNames.spec[,"Genus"] %in% c("Escherichia/Shigella"), "Genus"] <-
  "Escherichia"
tax.fixNames.spec[ tax.fixNames.spec[,"Genus"] %in% c("Prevotella_9"), "Genus"] <- "Prevotella"
tax.fixNames.spec[ tax.fixNames.spec[,"Genus"] %in% c("Prevotella_7"), "Genus"] <- "Prevotella"
tax.fixNames.spec[ tax.fixNames.spec[,"Genus"] %in%
                     c("Ruminiclostridium_5","Ruminiclostridium_6","Ruminiclostridium_9"), "Genus"] <- "Ruminiclostridium"
tax.fixNames.spec[ tax.fixNames.spec[,"Genus"] %in% c("Ruminococcus_1","Ruminococcus_2"), "Genus"] <-
  "Ruminococcus"
tax.fixNames.spec[ tax.fixNames.spec[,"Genus"] %in% c("Treponema_2"), "Genus"] <- "Treponema"

# make species value scientific name by combining with genus name if not unclassified
tax.fixNames.spec[,"Species"] <- sapply(rownames(tax.fixNames.spec), function(ro) {
  if (tax.fixNames.spec[ro,"Species"] == "unclassified") {
    tax.fixNames.spec[ro,"Species"]
  } else {
    paste(tax.fixNames.spec[ro,"Genus"], tax.fixNames.spec[ro,"Species"], sep = " ")
  }
})

# add other levels to unclassified genera to have unique values
tax.fixNames.spec <- t(apply(tax.fixNames.spec, 1, function(x){
  all.levels.spec <- as.character(x)
  if ("unclassified" %in% all.levels.spec){
    new.row.spec <- sapply(2:ncol(tax.fixNames.spec), function(tl)
      ifelse(all.levels.spec[tl]=="unclassified",
             paste(c("unclassified", all.levels.spec[1:(tl-1)]),
                   collapse = '.'),
             all.levels.spec[tl]))
    new.row.spec <- c(all.levels.spec[1], new.row.spec)
  }
  else{
    all.levels.spec
  }
}
)
)

colnames(tax.fixNames.spec) <- colnames(tax_spec)
# get asv.spec table
asv.spec.table <- st.nochim
colnames(asv.spec.table) <- unname(tax.fixNames.spec[colnames(st.nochim) , "Species"])
non.bact.spec <- tax.fixNames.spec[ tax.fixNames.spec[,"Kingdom"] %in%
                                      c(NA,"Eukaryota","Archaea","unclassified"), "Species"]
asv.spec.table <- asv.spec.table[ , ! colnames(asv.spec.table) %in% non.bact.spec]


asv.spec.table <- t(asv.spec.table %*% sapply(unique(colnames(asv.spec.table)),"==",
                                              colnames(asv.spec.table)))
asv.spec.table.rel <- apply(asv.spec.table, 2, function(x) 100 * x/sum(x))
asv.spec.table.rel[is.nan(asv.spec.table.rel)] <- 0

write.csv(asv.spec.table, sprintf("%s/results/asv_spec_table.csv", path))
write.csv(asv.spec.table.rel, sprintf("%s/results/asv_spec_table_rel.csv", path))

# get tax table with rownames as species values
tax.spec.table <- unique(tax.fixNames.spec)
rownames(tax.spec.table) <- unname(tax.spec.table[,"Species"])
tax.spec.table <- tax.spec.table[ ! rownames(tax.spec.table) %in% non.bact.spec, ]
tax.spec.table <- tax.spec.table[ sort(rownames(tax.spec.table)), ]
write.csv(tax.spec.table, sprintf("%s/results/tax_spec_table.csv", path))

save.image(sprintf("%s/results/ASV_table_built.RData", path))
print("ASV table built")

# Update unclassified db
unclassified.db <-
  readRDS(sprintf("unclassified_identifier_db.rds"))

## UPDATE unclassified identifier database

tls <- c("Phylum","Class","Order","Family","Genus","Species")
for (tl in tls) {
  # new table subset to be added to unclassified identifier database
  table.subset <- tax.spec.table[ startsWith(tax.spec.table[ , tl], "unclassified"), 1:(match(tl,
                                                                                              tls)+1) ]
  # only update table if there are any unclassified tax at given level
  if ( length(table.subset) > 0 ) {
    if ( class(table.subset) == "character" ) {
      tmp_table <- rbind( unclassified.db[[ tl ]], t(as.matrix(table.subset)) )
    } else {
      tmp_table <- rbind( unclassified.db[[ tl ]], table.subset )
    }
    
    # reduce the expanded unclassified values at all levels
    tmp_table[ startsWith( tmp_table, "unclassified") ] <- "unclassified"
    
    # keep only those rows which are unique for the columns up to but excluding the current tax level
    # (which is how the unclassified values are being identified)
    if ( is.null( rownames(unique(tmp_table[ , 1:match(tl, tls)])) ) ) {
      # this case occurs for the Phylum level because it checks just unique values in Kingdom,
      # so instead of returning another table (as in the else below),
      # it returns an unnamed character vector, so cannot identify rows
      tmp_table <- unique(tmp_table)
    } else {
      tmp_table <- tmp_table[rownames(unique(tmp_table[ , 1:match(tl, tls)])),]
    }
    
    # update rownames to numbers,
    # will keep current order so that those that were already present will be labeled with the same number
    rownames(tmp_table) <- 1:nrow(tmp_table)
    # Use numbered rows plus letter of tax level to label the deepest level unclassified tax
    tletter <- substr(tl, 1, 1)
    tmp_table[ , tl] <- sapply(rownames(tmp_table), function(x) sprintf("unclassified.%s%s",tletter,x))
    rownames(tmp_table) <- tmp_table[ , tl]
    
    # finally, update the table for the given taxonomic level in the unclassified identifier database
    unclassified.db[[ tl ]] <- tmp_table
    
  }
}

# add the unclassified species to the database from Olfat, updated with the Aspergillosis unclassified
species

saveRDS(unclassified.db, sprintf("%s/results/unclassified_identifier_db_updated.rds",
                                 path))

# Update ASV and TAXA tables with new unclassified identifiers
t0 <- Sys.time()

asvs <- list()
asvs_rel <- list()
taxTables <- list()

# function for converting tax names to appropriate unclassified identifiers
get_unclass_IDs <- function(tax, level) {
  # get tax table up to the indicated level
  tmp_tax <- unique(tax.spec.table[ , 1:(match(level, tls)+1) ])
  
  uID <- sapply(tax,
                function(x) {
                  if (startsWith(tmp_tax[x, level], "unclassified")) {
                    # get vector of tax up to one level above current based on the extended
                    taxonomy <- tmp_tax[ x, 1:match(level, tls) ]
                    # remove upper tax levels from all unclassified names,
                    # in order to match to the unclassified identifier database in next step
                    taxonomy <- sapply(taxonomy, function(y) strsplit(y, '\\.')[[1]][1])
                    # get correct unclassified identifier by checking where taxonomy matches a row in the table for the given level
                    new.unclass <- rownames( unclassified.db[[ level ]] )[
                      sapply(rownames(unclassified.db[[ level ]]),
                             function(y)
                               sum(unclassified.db[[ level ]][y, 1:match(level, tls)] == taxonomy) == length(taxonomy) ) ]
                    return( new.unclass )
                  } else {
                    return( tmp_tax[x, level] )
                  }
                })
  if (level != "Species" & length( uID[ uID %in% uID[duplicated(uID)] ] ) > 0) {
    dups <- uID[uID %in% uID[duplicated(uID)]]
    lev.up.from.dups <- tmp_tax[ names(dups), ncol(tmp_tax)-1 ]
    uID[ uID %in% uID[duplicated(uID)] ] <- sprintf("%s.%s", dups, lev.up.from.dups)
  }
  return( unname(uID) )
}

# Start by including the Species tax and asv tables, update the rownames with unclassified identifiers
taxTables[["Species"]] <- tax.spec.table
rownames(taxTables[["Species"]]) <- get_unclass_IDs( rownames(taxTables[["Species"]]), "Species")
taxTables[["Species"]] <- unique(taxTables[["Species"]])
asvs[["Species"]] <- asv.spec.table
rownames(asvs[["Species"]]) <- get_unclass_IDs( rownames(asvs[["Species"]]), "Species")
asvs_rel[["Species"]] <- asv.spec.table.rel
rownames(asvs_rel[["Species"]]) <- get_unclass_IDs( rownames(asvs_rel[["Species"]]), "Species")

# then get asv_tables at each level by summing the values from same tax at that level within the species asv table
for (tl in c("Genus","Family","Order","Class","Phylum")) {
  print(tl)
  
  ### update tax table at given level first
  taxTables[[ tl ]] <- unique(tax.spec.table[ , 1:(match(tl, tls)+1) ])
  rownames(taxTables[[ tl ]]) <- get_unclass_IDs( rownames(taxTables[[ tl ]]), tl)
  
  ### get counts of all tax at given tl from Species counts table that are within each value of the given tl
  tax <- taxTables[[ tl ]][ , tl]
  asvs[[ tl ]] <- apply(asvs[["Species"]], 2,
                        function(x) sapply(tax,
                                           function(y) sum(x[ names( taxTables[["Species"]][,tl][ y ==
                                                                                                    taxTables[["Species"]][,tl] ] )]
                                           )
                        )
  )
  asvs_rel[[ tl ]] <- apply( asvs[[ tl ]], 2, function(x) 100 * x/sum(x))
  
}

# Make species names proper scientific names (by making them <Genus species>)

rownames(asvs[["Species"]]) <- unname(sapply(rownames(asvs[["Species"]]),
                                             function(x) {
                                               if (startsWith(x, "unclassified")) {
                                                 if (startsWith(taxTables[["Species"]][x,"Genus"],
                                                                "unclassified")) {
                                                   paste("unclassified", x, sep = ' ')
                                                 } else {
                                                   paste(taxTables[["Species"]][x,"Genus"], x, sep = '
')
                                                 }
                                               } else {
                                                 x
                                               }
                                               
                                             }))

rownames(asvs_rel[["Species"]]) <- unname(sapply(rownames(asvs_rel[["Species"]]),
                                                 function(x) {
                                                   if (startsWith(x, "unclassified")) {
                                                     if (startsWith(taxTables[["Species"]][x,"Genus"],
                                                                    "unclassified")) {
                                                       paste("unclassified", x, sep = ' ')
                                                     } else {
                                                       paste(taxTables[["Species"]][x,"Genus"], x, sep
                                                             = ' ')
                                                     }
                                                   } else {
                                                     x
                                                   }
                                                   
                                                 }))

rownames(taxTables[["Species"]]) <- unname(sapply(rownames(taxTables[["Species"]]),
                                                  function(x) {
                                                    if (startsWith(x, "unclassified")) {
                                                      if (startsWith(taxTables[["Species"]][x,"Genus"],
                                                                     "unclassified")) {
                                                        paste("unclassified", x, sep = ' ')
                                                      } else {
                                                        paste(taxTables[["Species"]][x,"Genus"], x, sep
                                                              = ' ')
                                                      }
                                                    } else {
                                                      x
                                                    }
                                                    
                                                  }))

# fix the problem with Family_XI, which is found in both the orders Bacillales and Clostridiales
for (tl in c("Species","Genus","Family")) {
  taxTables[[ tl ]][ , "Family" ] <- unname(sapply( 1:nrow(taxTables[[tl]]), function(x) {
    if (taxTables[[tl]][ x, "Family" ] == "Family_XI") {
      sprintf("%s.%s", taxTables[[tl]][ x, "Family" ], taxTables[[tl]][ x, "Order" ])
    } else {
      taxTables[[tl]][ x, "Family" ]
    }
  } ))
  
  if (tl == "Family" & "Family_XI" %in% rownames(taxTables[[ tl ]]))
    rownames(taxTables[[ tl ]]) <- unname(sapply( rownames(taxTables[[ tl ]]), function(x) {
      if (x == "Family_XI") {
        taxTables[[ tl ]][ x, "Family"]
      } else {
        x
      }
    } ))
}

# replace long unclassified names with correct unclassified IDs in each column at each level
for (tl in rev(tls)) {
  #  must go backwards starting with deepest levels because will rely on values within higher levels
  for (col_to_change in tls[ (match(tl, tls)):1 ]) {
    # for given column, check the rownames from the table of that tax level to see what the official label should be
    #   in the case that col_to_change is same as tl, will essentially just use that tables rownames
    taxTables[[ tl ]][ , col_to_change ] <- unname(sapply( unname(taxTables[[ tl ]][ , col_to_change
    ]), function(x)
      rownames(taxTables[[ col_to_change ]])[ taxTables[[ col_to_change ]][, col_to_change ] == x ] ))
  }
}

# save tax and asv table objects
saveRDS(taxTables, sprintf("%s/results/taxTables.rds", path))
saveRDS(asvs, sprintf("%s/results/asvs.rds", path))
saveRDS(asvs_rel, sprintf("%s/results/asvs_rel.rds", path))

print(Sys.time() - t0)

save.image(sprintf("%s/results/ASV_taxa_complete.RData", path))
print("ASV_taxa_tables_completed")
