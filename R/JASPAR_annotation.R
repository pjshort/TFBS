# scan all of the well-covered regions with the pwm list and save output

# dependencies and libraries
source("simulate.R")
source("tfbs_core.R") # depends on TFBSTools, JASPAR2014, BioStrings
source("visualization.R")
library(optparse)

# command line options
option_list <- list(
  make_option("--verbose", action="store_true", default=TRUE,
              help="Print extra output advising the user of progression through the code."),
  make_option("--regions", default="../data/regions_annotated.txt",
              help="Pass the genomic regions that should be annotated with predicted TF binding sites."),
  make_option("--de_novos", default="../data/de_novo_filtered.txt",
              help="Pass the observed de novos - these will be analyzed in order to shorten the list of TF binding sites used for prediction."),
  make_option("--full_jaspar", default=TRUE,
              help="If TRUE annotates regions with the full JASPAR Homo Sapiens set."),
  make_option("--out", default="../data/regions_JASPAR_annotated.txt",
              help="Set location to save the (likely large) tab-delimited text file storing JASPAR annotations.")
)

args <- parse_args(OptionParser(option_list=option_list))

# load full PWM list
if ( args$verbose ) {
  write("Loading JASPAR position weight matrices from database...", stderr())
}
opts = list("species" = "Homo sapiens", "all_versions" = args$full_jaspar, "matrixtype" = "PWM")
pwm_list = getMatrixSet(JASPAR2014, opts)

# de novos are filtered for pp_dnm > 0.1 and non-coding only
de_novo_filtered <- read.table(args$de_novos, sep="\t", header=TRUE)

# regions are annotated by n_snps and n_indels - created by rupit/R/pre_process.R
well_covered_regions <- read.table(as.character(args$regions), sep="\t", header=TRUE)

if ( args$verbose ) {
  write("Processing regions input and de novos for primary scan.", stderr())
}
# make new 'de novo containing regions' (dcr) data frame
dcr = well_covered_regions[well_covered_regions$n_snp > 0 | well_covered_regions$n_indel > 0, ]
dcr = merge(dcr, de_novo_filtered[, c("region_id", "pos", "ref", "alt", "person_stable_id")], by = "region_id")
dcr$rel_pos = dcr$pos - dcr$start

# REMOVE THIS - ONLY FOR TESTING
dcr = dcr[1:100, ]

if ( args$verbose ) {
  write("Starting primary scan to identify TFBS of interest from de novo set.", stderr())
}
# get matrix (TF x dcr) with 1 if de novo in that region disrupts TF binding
ninety_five_counts = mapply(TFBS_hits, dcr$seq, dcr$rel_pos, MoreArgs = list("pwm_list" = pwm_list, "min.score" = "95%"))

# get number of TFBS disruptions per de novo
dcr$TFBS_count = colSums(ninety_five_counts)

# get number of disruptions per TFBS
counts_by_TFBS = rowSums(ninety_five_counts)

# get only the TFBS that have some annotation
recurrent_TFBS = names(counts_by_TFBS)[counts_by_TFBS > 0]
pwm_list = pwm_list[recurrent_TFBS]

# NOTE: this scanning is the longest step to complete
if ( args$verbose ) {
  write("Scanning all regions supplied against TFBS identified in primary scan.", stderr())
}

# EDIT THIS - ONLY FOR TESTING
scanned_regions = scan_regions(well_covered_regions[1:68,], pwm_list, min.score = "95%")

# make a dataframe containing region IDs and predicted TF binding sites in these regions
region_dfs = mapply(function(s, n) data.frame(region_id = rep(n, length(s)), start=s@start, stop=s@start+s@width - 1, name = s@NAMES), scanned_regions, names(scanned_regions), SIMPLIFY = FALSE)
full_df = do.call(rbind, region_dfs)
full_df$tf_name = sapply(as.character(full_df$name), function(n) pwm_list[[n]]@name)
full_df$family = sapply(as.character(full_df$name), function(n) pwm_list[[n]]@tags$family)
full_df$alias = sapply(as.character(full_df$name), function(n) pwm_list[[n]]@tags$alias)
full_df$alias = vapply(full_df$alias, paste, collapse = ", ", character(1L))  # removes the NULL entries (turning to empty)

write.table(full_df, file = args$out, row.names = FALSE, sep = "\t", col.names = TRUE)
if ( args$verbose ) {
  write(sprintf("Finished! Regions annotated by JASPAR and saved to: %s", args$out), stderr())
}


