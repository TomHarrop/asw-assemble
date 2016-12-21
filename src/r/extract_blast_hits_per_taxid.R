#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)

#############
# arguments #
#############

# testing
# command.args <- c("-t", "output/blastqc/2125-01-06-1_mp_R1.table",
#                   "-t", "output/blastqc/2125-01-06-1_mp_R2.table",
#                   "-t", "output/blastqc/2125-01-06-1_pe_R1.table",
#                   "-t", "output/blastqc/2125-01-06-1_pe_R2.table",
#                   "-t", "output/blastqc/2125-01-06-1_se_R1.table",
#                   "-t", "output/blastqc/2125-01-06-1_se_R2.table",
#                   "-t", "output/blastqc/2125-01-06-1_unknown_R1.table",
#                   "-t", "output/blastqc/2125-01-06-1_unknown_R2.table",
#                   "-t", "output/blastqc/2125-01-11-1_R1.table",
#                   "-t", "output/blastqc/2125-01-11-1_R2.table",
#                   "-y", "data/ncbi/nucl_gb.accession2taxid.Rds",
#                   "-y", "data/ncbi/nodes.dmp.Rds",
#                   "-y", "data/ncbi/names.dmp.Rds",
#                   "-r", "output/blastqc/plots.pdf")

command.args <- commandArgs(trailingOnly = TRUE)

parsed.args <- argparsR::ParseArguments(
    accepted.switches = list(
        `input_table` = "-t",
        `other_input` = "-y",
        `output_pdf` = "-r"),
    command.args)

tax.names.file <- grep("names.dmp", parsed.args$other_input, value = TRUE)
tax.nodes.file <- grep("nodes.dmp", parsed.args$other_input, value = TRUE)
acc2taxid.file <- grep("accession2taxid", parsed.args$other_input,
                       value = TRUE)

#############
# functions #
#############

# extract the library name from the .table files
ExtractLibraryName <- function(x) {
    gsub(".*?([^_]+)_(R[[:digit:]]).table", "\\1_\\2", basename(x))
}

# extract the accession number from the blast hit_id
ExtractAccessionNumber <- function(x) {
    strsplit(x, split = "|", fixed = TRUE)[[1]][4]
}

# given a tax_id, find the "class" node
FindRootTaxid <- function(x) {
    # cat("search taxid ", x, "\n")
    my.taxid <- copy(x)
    while (!stop.nodes[, my.taxid %in% tax_id] & 
           my.taxid != 1) {
        # cat("root taxid", my.taxid, "\n")
        my.taxid <- tax.nodes[tax_id == my.taxid, as.numeric(parent_tax_id)]
    }
    # cat("parent taxid ", my.taxid, "\n")
    my.taxid
}

# given a root taxid, try to pick a good taxid name
PickTaxidName <- function(x) {
    possible.names <- tax.names[tax_id == x]
    if (possible.names[, "scientific name" %in% name_class]) {
        return(possible.names[name_class == "scientific name",
                                    name_txt][[1]])
    } else {
        return(possible.names[,name_txt][[1]])
    }
}

# given a taxid, return the node type
GetNodeType <- function(x) {
    tax.nodes[tax_id == x, rank]
}

########
# data #
########

# this takes a long time because the accession2taxid is large. That means it's
# best to run this script only once using all BLAST results simultaneously

# load taxonomy data
rutils::GenerateMessage("Loading data files. Usually takes a long time.")
tax.names <- readRDS(tax.names.file)
tax.nodes <- readRDS(tax.nodes.file)
acc2taxid <- readRDS(acc2taxid.file)

# nodes where we want to stop searching for parents
stop.ranks <- c("superkingdom", "kingdom", "phylum", "class")
stop.nodes <- tax.nodes[rank %in% stop.ranks]

########
# code #
########

# read parsed BLAST data (this is where we read input.files)
rutils::GenerateMessage("Loading parsed BLAST results")
parsed.blast.files <- parsed.args$input_table
names(parsed.blast.files) <- sapply(parsed.blast.files, ExtractLibraryName)
parsed.blast.list <- lapply(parsed.blast.files, fread)

# combine parsed BLAST data
rutils::GenerateMessage("Combining parsed BLAST results")
blast.data <- rbindlist(parsed.blast.list, idcol = "library_name")

# add accession number
rutils::GenerateMessage("Extracting accession numbers")
blast.data[, accession.version := ExtractAccessionNumber(hit_id), by = hit_id]

# get tax_id from acc2taxid by accession.version
rutils::GenerateMessage("Adding taxid")
setkey(blast.data, accession.version)
blast.data.with.acc <- acc2taxid[blast.data]

# get parent taxid
rutils::GenerateMessage("Searching nodes.dmp for root taxid (slow)")
root.taxids <- blast.data.with.acc[!is.na(taxid),
                                   .(root.taxid = FindRootTaxid(taxid)),
                                   by = taxid]

# add taxon names
rutils::GenerateMessage("Adding taxon names by taxid")
root.taxids[, root.taxid.name := PickTaxidName(root.taxid), by = root.taxid]
root.taxids[, taxid.name := PickTaxidName(taxid), by = taxid] # root nodes

# add parent taxid class
rutils::GenerateMessage("Adding taxon types")
root.taxids[, root.taxid.rank := GetNodeType(root.taxid), by = root.taxid]

# merge blast results with taxid results
rutils::GenerateMessage("Combining BLAST results with taxon data")
blast.data.with.taxinfo <- merge(blast.data.with.acc,
                                 root.taxids,
                                 by.x = "taxid",
                                 by.y = "taxid",
                                 all = TRUE)

# split library_name into library and read
rutils::GenerateMessage("Splitting library_name column")
blast.data.with.taxinfo[, c("library.name", "read") :=
                            tstrsplit(library_name, "_", fixed = TRUE)]

# pick the lowest evalue by group
rutils::GenerateMessage("Subsetting by evalue")
best.hit.per.query <- blast.data.with.taxinfo[
    blast.data.with.taxinfo[, .I[evalue == min(evalue)][1],
                            by = query_id][, V1]][evalue < 1]

# count the number of reads with blast hits for each library
rutils::GenerateMessage("BLAST hits with evalue < 1 per library:")
best.hit.per.query[,.(total.blast.hits = length(unique(query_id))),
                   by = .(library.name, read)]

# count the number of reads per "root" for each library
reads.per.root <- best.hit.per.query[!is.na(root.taxid),
                   .(root.hit = unique(root.taxid.name)),
                   by = .(library.name, read, query_id)][
                       ,
                       .(number.of.reads = length(unique(query_id))),
                       by = .(root.hit, library.name, read)]

# count the number of reads with blast hits AND parsed taxa for each library
rutils::GenerateMessage("BLAST hits with evalue < 1 and parsed root taxa:")
reads.per.root[, .(blast.hits.with.taxid = sum(number.of.reads)),
               by = .(library.name, read)]

#########
# plots #
#########

rutils::GenerateMessage("Generating plots")

# sum the hits < 5
other.roots <- reads.per.root[, sum(number.of.reads),
               by = root.hit][V1 < 10, as.character(root.hit)]
reads.per.root.summed <- rbind(reads.per.root[!root.hit %in% other.roots],
      reads.per.root[root.hit %in% other.roots,
                     .(root.hit = "other",
                       number.of.reads = sum(number.of.reads)),
                     by = .(library.name, read)])

# order bars
reads.per.root.summed[, total.hits := sum(number.of.reads), by = root.hit]
setkey(reads.per.root.summed, total.hits)
root.levels <- reads.per.root.summed[, rev(unique(root.hit))]
root.levels <- c(root.levels[root.levels != "other"], "other")
reads.per.root.summed[, root.hit := factor(root.hit, levels = root.levels)]
reads.per.root.summed[, total.hits := NULL]

gp <- ggplot(reads.per.root.summed, aes(x = root.hit, y = number.of.reads)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    geom_col() +
    facet_grid(library.name ~ read)

# write output
outdir <- dirname(parsed.args$output_pdf)
rutils::GenerateMessage("Writing output")
rutils::PrintF("outdir: %s\n", outdir)

ggsave(filename = parsed.args$output_pdf,
       device = cairo_pdf,
       plot = gp,
       width = 10, height = 7.5, units = "in")
saveRDS(blast.data.with.taxinfo,
        paste0(outdir, "/blast_data_with_taxinfo.Rds"))
saveRDS(best.hit.per.query,
        paste0(outdir, "/best_hit_per_query.Rds"))

# log metadata
rutils::GenerateMessage("Logging metadata")
log.file <- paste0(outdir, "/blast_plots.SessionInfo.txt")
s.inf <- rutils::GitSessionInfo()
writeLines(s.inf, log.file)

quit("no", 0)

        
        