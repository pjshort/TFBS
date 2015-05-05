# scan DDD data against all JASPAR TFs

source("tfbs_core.R") # depends on TFBSTools, JASPAR2014, BioStrings
library(optparse)

# command line options
option_list <- list(
  make_option("--verbose", action="store_true", default=TRUE,
              help="Print extra output advising the user of progression through the code."),
  make_option("--regions", default="../data/regions_annotated.txt",
              help="Pass the genomic regions that should be annotated with predicted TF binding sites."),
  make_option("--de_novos", default="../data/de_novo_filtered.txt",
              help="Pass the observed de novos - these will be analyzed in order to shorten the list of TF binding sites used for prediction."),
  make_option("--out", default="../data/TFs_in_DDD_data.txt",
              help="Set location to save the (likely large) tab-delimited text file storing JASPAR annotations.")
)

args <- parse_args(OptionParser(option_list=option_list))

# load full PWM list
if ( args$verbose ) {
  write("Loading JASPAR position weight matrices from database...", stderr())
}

# NOTE: only need to initalize DB once - this should be done on install...?
#db = "myMatrixDb.sqlite"
#initializeJASPARDB(db) 
opts = list("species" = 9606, "all_versions" = TRUE, "matrixtype" = "PWM") # 9606 = "homo sapiens"
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

# save one column file with JASPAR internal names of TFs
write.table(unique(recurrent_TFBS), file = args$out, row.names = FALSE, col.names = "jaspar_internal")

