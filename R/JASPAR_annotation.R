# scan all of the well-covered regions with the pwm list and save output

# dependencies and libraries
source("tfbs_core.R") # depends on TFBSTools, JASPAR2014, BioStrings
library(optparse)

# command line options
option_list <- list(
  make_option("--verbose", action="store_true", default=TRUE,
              help="Print extra output advising the user of progression through the code."),
  make_option("--regions", default="../data/regions_annotated.txt",
              help="Pass the genomic regions that should be annotated with predicted TF binding sites."),
  make_option("--tf_list", default="../data/TFs_in_DDD_data.txt",
              help="Pass a list of TFs to be run against regions."),
  make_option("--full_jaspar", action="store_true", default=TRUE,
              help="If TRUE annotates regions with the full JASPAR Homo Sapiens set. If false, must provide a --tflist"),
  make_option("--out", default="../data/regions_JASPAR_annotated.txt",
              help="Set location to save the (likely large) tab-delimited text file storing JASPAR annotations."),
  make_option("--min_score", default="95%",
              help="Set min score for the JASPAR annotations - less than 95% may yield unrealistic number of predicted binding sites.")
)

args <- parse_args(OptionParser(option_list=option_list))

# load full PWM list
if ( args$verbose ) {
  write("Loading JASPAR position weight matrices from database...", stderr())
}

# NOTE: only need to initalize DB once - this should be done on install...
#db = "myMatrixDb.sqlite"
#initializeJASPARDB(db)
opts = list("species" = 9606, "all_versions" = TRUE, "matrixtype" = "PWM") # 9606 = "homo sapiens"
pwm_list = getMatrixSet(JASPAR2014, opts)

# regions are annotated by n_snps and n_indels - created by rupit/R/pre_process.R
well_covered_regions <- read.table(as.character(args$regions), sep="\t", header=TRUE)
well_covered_regions = well_covered_regions[1:10,]

# scan only on the TFBS specified in input file
if (args$full_jaspar == FALSE){  # if args$full_jaspar is TRUE use the entire JASPAR list - will be huge for large set of regions
  # load the names of TFs that we want to scan for
  TFBS_to_scan = read.table(args$tf_list, header = TRUE)
  pwm_list = pwm_list[unique(TFBS_to_scan$jaspar_internal)]
}

if ( args$verbose ) {
  write(sprintf("Scanning %i regulatory regions with %i different TF binding motifs...", nrow(well_covered_regions), length(pwm_list)), stderr())
}

scanned_regions = scan_regions(well_covered_regions, pwm_list, min.score = args$min_score)

if ( args$verbose ) {
  write(sprintf("Finished scanning. Scanned regions list has %i entries.", length(scanned_regions)), stderr())
}

# make a dataframe containing region IDs and predicted TF binding sites in these regions
region_dfs = mapply(function(s, n) data.frame(region_id = rep(n, length(s)), start=s@start, stop=s@start+s@width - 1, name = s@NAMES), scanned_regions, names(scanned_regions), SIMPLIFY = FALSE)
full_df = do.call(rbind, region_dfs)
full_df$tf_name = sapply(as.character(full_df$name), function(n) pwm_list[[n]]@name)
full_df$family = sapply(as.character(full_df$name), function(n) pwm_list[[n]]@tags$family)
full_df$alias = sapply(as.character(full_df$name), function(n) pwm_list[[n]]@tags$alias)
full_df$alias = vapply(full_df$alias, paste, collapse = ", ", character(1L))  # removes the NULL entries (turning to empty)

if ( args$verbose ) {
  write(sprintf("Full dataframe for output has %i rows and %i columns.", nrow(full_df), ncol(full_df)), stderr())
}

write.table(full_df, file = args$out, row.names = FALSE, sep = "\t", col.names = TRUE)
if ( args$verbose ) {
  write(sprintf("Finished! Regions annotated by JASPAR and saved to: %s", args$out), stderr())
}





