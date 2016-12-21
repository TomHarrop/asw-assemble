#!/usr/bin/env Rscript

library(data.table)

# get input
command.args <- commandArgs(trailingOnly = TRUE)
# command.args <- c("-z", "data/ncbi/nucl_gb.accession2taxid.Rds",
#                   "-z", "data/ncbi/nodes.dmp.Rds",
#                   "-z", "data/ncbi/names.dmp.Rds")
parsed.args <- argparsR::ParseArguments(
    accepted.switches = list(
        `other.output` = "-z"),
    command.args)
nucl_gb_file <- grep("nucl_gb", parsed.args$other.output, value = TRUE)
names_file <- grep("names.dmp", parsed.args$other.output, value = TRUE)
nodes_file <- grep("nodes.dmp", parsed.args$other.output, value = TRUE)

outdir <- dirname(nucl_gb_file)

# download accession to taxid conversion table
rutils::GenerateMessage("Downloading accession2taxid")
nucl_gb.accession2taxid.url <-
    "ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz"
temp1 <- tempfile(fileext = ".gz")
download.file(nucl_gb.accession2taxid.url, temp1)

# read in with data.table
rutils::GenerateMessage("fread-ing accession2taxid")
nucl_gb.accession2taxid <- fread(paste("zcat", temp1))

# download taxdump (complete NCBI taxonomy)
rutils::GenerateMessage("Downloading taxdump")
taxdmp.url <- "ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdmp.zip"
temp2 <- tempfile(fileext = ".zip")
temp3 <- tempdir()
download.file(taxdmp.url, temp2)

# read nodes.dmp from inside taxdmp
rutils::GenerateMessage("fread-ing nodes.dmp")
nodes.dmp.file <- unzip(temp2, files = "nodes.dmp", exdir = temp3)
nodes.dmp.raw <- fread(nodes.dmp.file, sep = "|", quote = "\t", header = FALSE)

# tidy up nodes.dmp
rutils::GenerateMessage("Munging nodes.dmp")
nodes.dmp.raw[, V14 := NULL]
nodes.names <- c("tax_id", "parent_tax_id", "rank", "embl_code", "division_id",
                 "inherited_div_flag", "genetic_code_id", "inherited_GC",
                 "mitochondrial_genetic_code_id", "inherited_MGC_flag",
                 "GenBank_hidden_flag", "hidden_subtree_root_flag", "comments")
setnames(nodes.dmp.raw, names(nodes.dmp.raw), nodes.names)
nodes.dmp.raw[, tax_id := as.numeric(gsub("\t", "", tax_id, fixed = TRUE))]
setkey(nodes.dmp.raw, tax_id)

# read names.dmp from inside taxdmp
rutils::GenerateMessage("fread-ing names.dmp")
names.dmp.file <- unzip(temp2, files = "names.dmp", exdir = temp3)
names.dmp.raw <- fread(names.dmp.file, sep = "|", quote = "\t", header = FALSE)

# tidy up names.dmp
rutils::GenerateMessage("Munging names.dmp")
names.dmp.raw[, V5 := NULL]
names.names <- c("tax_id", "name_txt", "unique_name", "name_class")
setnames(names.dmp.raw, names(names.dmp.raw), names.names)
names.dmp.raw[, tax_id := as.numeric(gsub("\t", "", tax_id, fixed = TRUE))]
setkey(names.dmp.raw, tax_id)

# read merged.dmp from taxdmp
rutils::GenerateMessage("fread-ing merged.dmp")
merged.dmp.file <- unzip(temp2, files = "merged.dmp", exdir = temp3)
merged.dmp.raw <- fread(merged.dmp.file, sep = "|", quote = "\t", header = FALSE)

# tidy up merged.dmp
rutils::GenerateMessage("Munging merged.dmp")
merged.dmp.raw[, V3 := NULL]
merged.names <- c("old_tax_id", "new_tax_id")
setnames(merged.dmp.raw, names(merged.dmp.raw), merged.names)
merged.dmp.raw[, old_tax_id :=
                   as.numeric(gsub("\t", "", old_tax_id, fixed = TRUE))]
merged.dmp.raw[, new_tax_id := as.numeric(new_tax_id)]

# fix node names in accession2taxid
rutils::GenerateMessage("Using merged.dmp to fix taxid in accession2taxid")
nucl_gb.accession2taxid[, taxid := as.numeric(taxid)]
accession2taxid <- merge(nucl_gb.accession2taxid,
                         merged.dmp.raw,
                         by.x = "taxid",
                         by.y = "old_tax_id",
                         all.x = TRUE)
accession2taxid[, taxid.tmp := new_tax_id]
accession2taxid[is.na(taxid.tmp), taxid.tmp := taxid]
accession2taxid[, old_tax_id := taxid]
accession2taxid[, taxid := taxid.tmp]
accession2taxid[, taxid.tmp := NULL]
accession2taxid[taxid == old_tax_id, old_tax_id := NA]

# sort accession2taxid
rutils::GenerateMessage("Sorting accession2taxid")
setkey(accession2taxid, accession.version)

# save output
rutils::GenerateMessage("Writing output")
rutils::PrintF("outdir: %s\n", outdir)
if (!dir.exists(outdir)) {
    dir.create(outdir)
}

saveRDS(accession2taxid, nucl_gb_file)
saveRDS(nodes.dmp.raw, nodes_file)
saveRDS(names.dmp.raw, names_file)

# log metadata
rutils::GenerateMessage("Logging metadata")
log.file <- paste0(outdir, "/SessionInfo.txt")
s.inf <- rutils::GitSessionInfo()
writeLines(s.inf, log.file)

quit("no", 0)
