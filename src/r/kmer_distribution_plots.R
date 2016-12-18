#!/usr/bin/env Rscript

rutils::GenerateMessage("Print kmer distribution plots")

library(data.table)
library(bit64)
library(ggplot2)
library(scales)

# file naming function
ExtractLibraryName <- function(x) {
    gsub(".*?([^_]+)_hist_[[:alnum:]]+.txt", "\\1", basename(x))
}

# for testing
# command.args <- c(
#     "--fq", "output/bbnorm/2125-01-06-1_mp.fastq.gz",
#     "--fq", "output/bbnorm/2125-01-06-1_pe.fastq.gz",
#     "--fq", "output/bbnorm/2125-01-06-1_se.fastq.gz",
#     "--fq", "output/bbnorm/2125-01-06-1_unknown.fastq.gz",
#     "--fq", "output/bbnorm/2125-01-11-1.fastq.gz",
#     "-r", "output/bbnorm/plots.pdf",
#     "-z", "output/bbnorm/plot_data.Rds")

# Parse arguments
parsed.args <- argparsR::ParseArguments(
    accepted.switches = list(
        "input_fq" = "--fq",
        "output_pdf" = "-r",
        "other_output" = "-z")
    # , command.args = command.args
)

# convert fastq files to hist
rutils::GenerateMessage("Loading input")

before.files <- gsub(".fastq.gz", "_hist_before.txt", parsed.args$input_fq)
rutils::PrintF("       Raw histogram file: %s\n", before.files)

after.files <- gsub(".fastq.gz", "_hist_after.txt", parsed.args$input_fq)
rutils::PrintF("Normalised histogram file: %s\n", after.files)

# test load all output data, need to parse with argparsR
# before.files <- list.files(path = "output/bbnorm",
#                            pattern = "hist_before",
#                            full.names = TRUE)
# after.files <- list.files(path = "output/bbnorm",
#                           pattern = "hist_after",
#                           full.names = TRUE)

# name the lists
names(before.files) <- paste0(
    sapply(before.files, ExtractLibraryName), ".Raw")
names(after.files) <- paste0(
    sapply(before.files, ExtractLibraryName), ".Normalised")

# read the data.tables
before.data.tables <- lapply(before.files, fread)
after.data.tables <- lapply(after.files, fread)

# make one big table
rutils::GenerateMessage("Combining data")
before.data <- rbindlist(before.data.tables, idcol = "lib.type")
after.data <- rbindlist(after.data.tables, idcol = "lib.type")
hist.data <- rbind(before.data, after.data)
hist.data[, c("lib", "type") := tstrsplit(lib.type, ".", fixed = TRUE)]
hist.data[, lib.type := NULL]

# order and rename data
hist.data[, type := factor(type, levels = c("Raw", "Normalised"))]
lib.order = c("2125-01-11-1" = "Paired end reads",
              "mp" = "Long mate pairs",
              "pe" = "Short pairs",
              "se" = "Singletons",
              "unknown" = "Unknown")
hist.data[, lib := factor(plyr::revalue(lib, lib.order),
                   levels = lib.order)]

# colour for plots
c.scale <- RColorBrewer::brewer.pal(9, "Set1")

# raw vs norm plot
rutils::GenerateMessage("Generating plot")
kmer_plot <- ggplot(hist.data, aes(x = `#Depth`, y = Unique_Kmers)) +
    theme_minimal() +
    facet_grid(lib ~ type) +
    geom_path(colour = c.scale[1]) +
    scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
    scale_x_continuous(trans = log_trans(base = 4)) +
                       # labels = trans_format(function(x) log(x, 4),
                       #                       math_format(4^.x))) +
    xlab("31-mer depth") + ylab("Number of unique 31-mers")

kmer_plot <- ggplot(hist.data, aes(x = `#Depth`, y = Unique_Kmers, colour = type)) +
    theme_minimal() +
    theme(legend.position = c(5/6, 1/4)) +
    facet_wrap(~ lib) +
    geom_path(alpha = 0.75) +
    scale_colour_brewer(palette = "Set1",
                        guide = guide_legend(title = NULL)) +
    scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
    scale_x_continuous(trans = log_trans(base = 4)) +
    # labels = trans_format(function(x) log(x, 4),
    #                       math_format(4^.x))) +
    xlab("31-mer depth") + ylab("Number of unique 31-mers")


# save output
rutils::GenerateMessage("Saving output")
ggsave(filename = parsed.args$output_pdf,
       plot = kmer_plot,
       width = 10, height = 7.5, units = "in")
saveRDS(hist.data, parsed.args$other_output)

# write logs
rutils::GenerateMessage("Writing logs")
out.dir <- dirname(parsed.args$other_output)
log.file <- paste0(out.dir, "/SessionInfo.kmer_distribution_plots.txt")
rutils::PrintF("Log file: %s\n", log.file)
s.inf <- rutils::GitSessionInfo()
writeLines(s.inf, log.file)

rutils::GenerateMessage("Done")
quit(save = "no", status = 0)
