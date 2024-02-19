# 2023 Benjamin J Perry
# MIT License
# Copyright (c) 2022 Benjamin J Perry
# Version: 1.0
# Maintainer: Benjamin J Perry
# Email: ben.perry@agresearch.co.nz
#
# Utility script to generate a proportional bar graph of reads from TaFFE QC.
# Intended use is with the TaFFE or GBS-TaFFE snakemake workflows.
#
# Inputs:
# 1. seqkit.report.raw.txt
# 2. seqkit.report.prinseq.txt
# 3. seqkit.report.bbduk.txt
# 4. seqkit.report.KDTrim.txt
# 5. seqkit.report.KDTRF.txt
# 6. seqkit.report.KDOvis.txt
# 7. seqkit.report.KDSILVA138.txt
# 8. seqkit.report.KDR.txt
# 9. keyfile - to link FID with sample
#
# Returns:
#   1. Reformatted matrix of merge reads counts for each step
#   2. A stacked bar plot of the read count proportions per library


### Import Libraries ###
suppressMessages(library("tidyverse"))
suppressMessages(library("optparse"))


### Parse Command Line Arguments ###
cli_args <- list( 
  make_option(c("-r", "--raw"), action="store", type = "character", default=NULL,
              help="path seqkit report: seqkit.report.raw.txt"),
  make_option(c("-b", "--bbduk"), action="store", type = "character", default=NULL,
              help="path seqkit report: seqkit.report.bbduk.txt"),
  make_option(c("-p", "--prinseq"), action="store", type = "character", default=NULL,
              help="path seqkit report: seqkit.report.prinseq.txt"),
  make_option(c("-t", "--kdtrim"), action="store", type = "character", default=NULL,
              help="path seqkit report: seqkit.report.KDTrim.txt"),
  make_option(c("-a", "--kdtrf"), action="store", type = "character", default=NULL,
              help="path seqkit report: seqkit.report.KDTRF.txt"),
  make_option(c("-o", "--kdhost"), action="store", type = "character", default=NULL,
              help="path seqkit report: seqkit.report.KDOvis.txt"),
  make_option(c("-s", "--kdsilva"), action="store", type = "character", default=NULL,
              help="path seqkit report: seqkit.report.KDSILVA138.txt"),
  make_option(c("-f", "--kdreads"), action="store", type = "character", default=NULL,
              help="path seqkit report: seqkit.report.KDR.txt"),
  make_option(c("-k", "--keyfile"), action="store", type = "character", default=NULL,
              help="keyfile containing FIDs and sampleID for all samples")
  )

parser <- OptionParser(
  usage = "%prog [--flags] path/to/seqkit.report.input.txt ",
  option_list = cli_args,
  formatter = IndentedHelpFormatter, 
  add_help_option = TRUE,
  prog = "make_read_proportions_plot.R")

args <- parse_args(parser)

print(args)


### Function Definitions ###
read_seqkit_report <- function(path = "") {
  # function to read and parse a seqkit stats report
  suppressMessages(library("tidyverse"))
  report <- as.tibble(read_table(file = path, col_names = TRUE))
  report$factid <- str_extract(report$file, "FID[0-9]+")
  return(report)
}

read_keyfile <- function(path = "") {
  # function to read and parse keyfile
  suppressMessages(library("tidyverse"))
  report <- as.tibble(read_tsv(file = path, col_names = TRUE))
  return(report)
}


### Main ###

raw <- read_seqkit_report(path = args$raw)
raw <- raw %>% dplyr::select(factid, raw_count = num_seqs)

bbduk <- read_seqkit_report(path = args$bbduk)
bbduk <- bbduk %>% dplyr::select(factid, bbduk_count = num_seqs)

prinseq <- read_seqkit_report(path = args$prinseq)
prinseq <- prinseq %>% dplyr::select(factid, prinseq_count = num_seqs)

trim <- read_seqkit_report(path = args$kdtrim)
trim <- trim %>% dplyr::select(factid, kdtrim_count = num_seqs)

trf <- read_seqkit_report(path = args$kdtrf)
trf <- trf %>% dplyr::select(factid, kdtrf_count = num_seqs)

host <- read_seqkit_report(path = args$kdhost)
host <- host %>% dplyr::select(factid, kdhost_count = num_seqs) 

silva <- read_seqkit_report(path = args$kdsilva)
silva <- silva %>% dplyr::select(factid, kdsilva_count = num_seqs)

kdr <- read_seqkit_report(path = args$kdreads)
kdr <- kdr %>% dplyr::select(factid, kdreads_count = num_seqs)  

keyfile <- read_keyfile(path = args$keyfile)

read_counts <- keyfile %>% dplyr::select(factid, animalid, uidtag, fullsamplename, sampleid = sample, species = genophyle_species, cohort = gbs_cohort, control)
read_counts$species_cohort <- paste(read_counts$species, read_counts$cohort, read_counts$control, sep = "-")

read_counts <- left_join(read_counts, raw, by = "factid")
read_counts <- left_join(read_counts, bbduk, by = "factid")
read_counts <- left_join(read_counts, prinseq, by = "factid")
read_counts <- left_join(read_counts, trim, by = "factid")
read_counts <- left_join(read_counts, trf, by = "factid")

read_counts$bbduk_removed <- read_counts$raw_count - read_counts$bbduk_count
read_counts$prinseq_removed <- read_counts$bbduk_count - read_counts$prinseq_count
read_counts$kdtrim_removed <- read_counts$prinseq_count - read_counts$kdtrim_count
read_counts$kdtrf_removed <- read_counts$kdtrim_count - read_counts$kdtrf_count

read_counts <- left_join(read_counts, silva, by = "factid")
read_counts <- left_join(read_counts, kdr, by = "factid")
read_counts <- left_join(read_counts, host, by = "factid")

read_counts$host_proportion <- read_counts$kdhost_count / read_counts$raw_count
read_counts$KDR_proportion <- read_counts$kdreads_count / read_counts$raw_count

write_csv(read_counts,"results/merged_read_count_summary.csv")

QC_plot <- read_counts %>% 
  arrange(kdhost_count) %>%
  pivot_longer(., 
             cols = c("bbduk_removed", 
             "prinseq_removed", 
             "kdtrim_removed", 
             "kdtrf_removed", 
             "kdsilva_count", 
             "kdreads_count", 
             "kdhost_count"), names_to = "fraction", values_to = "proportion") %>%
  ggplot(aes(x = factid, y = proportion, fill = fraction)) + 
  geom_col(position = "fill", width = 1) + facet_wrap(vars(species_cohort), ncol = 1) + theme_classic()

ggsave("results/Read_Proportions_QC_Plot.svg", plot = QC_plot, limitsize = FALSE, dpi = 800, units = "cm", height = 100, width = 65)

KDR_QC_plot <- read_counts %>% 
  ggplot(aes(x = factid, y = kdreads_count, fill = species_cohort)) + 
  geom_col() + 
  scale_y_log10() +
  annotation_logticks()
ggsave("results/KDRs_QC_Plot.svg", plot = KDR_QC_plot, limitsize = FALSE, dpi = 800, units = "cm", height = 20, width = 100)


host_QC_plot <- read_counts %>% 
  ggplot(aes(x = factid, y = kdhost_count, fill = species_cohort)) + 
  geom_col() + 
  scale_y_log10() +
  annotation_logticks()
ggsave("results/Host_QC_Plot.svg", plot = host_QC_plot, limitsize = FALSE, dpi = 800, units = "cm", height = 20, width = 100)


proportions_QC_boxplot <- read_counts %>% 
  pivot_longer(cols = c(host_proportion, KDR_proportion), 
               names_to = "fraction", 
               values_to = "proportion") %>% 
  ggplot(aes(x = fraction, y = proportion, color = species_cohort)) + 
  geom_boxplot() + 
  theme_minimal() + 
  scale_y_continuous(breaks = seq(0, 1, by=0.1), limits=c(0,1))
ggsave("results/Proportions_QC_Boxplot.svg", plot = proportions_QC_boxplot, limitsize = FALSE, dpi = 800, units = "cm", height = 20, width = 30)


