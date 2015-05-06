# script to run large set of simulations and annotate with JASPAR pipeline

# depends on:
source("../R/simulate.R")
source("../R/tfbs_core.R")
source("../R/visualization.R")
library(optparse)

# command line options
option_list <- list(
  make_option("--n_snps", default=NULL,
              help="Pass the genomic regions that should be annotated with predicted TF binding sites."),
  make_option("--iterations", default=3,
              help="Pass a list of TFs to be run against regions."),
  make_option("--out", default="../data/TFBS_simulation_stats.txt",
              help="Set location to save the (likely large) tab-delimited text file storing JASPAR annotations."),
  make_option("--jaspar_annotated_regions", default="../data/regions_JASPAR_annotated_FULL.txt",
              help="Pass a list of TFs to be run against regions."),
  make_option("--min_score", default="95%",
              help="Set min score for the JASPAR annotations - less than 95% may yield unrealistic number of predicted binding sites."),
  make_option("--verbose", action="store_true", default=TRUE,
              help="Print extra output advising the user of progression through the code.")
)

args <- parse_args(OptionParser(option_list=option_list))


# set parameters
iterations = args$iterations
if (!is.null(args$n_snps)){
  snp_total = args$n_snps
} # else, set later using sum of n_snp in well_covered_regions


## load datasets ##

opts = list("species" = 9606, "all_versions" = TRUE, "matrixtype" = "PWM") # 9606 = "homo sapiens"
pwm_list = getMatrixSet(JASPAR2014, opts)

# regions are annotated by n_snps and n_indels - created by rupit/R/pre_process.R
well_covered_regions <- read.table("../data/regions_annotated.txt", sep="\t", header=TRUE)
if (is.null(args$n_snps)){
  snp_total = sum(well_covered_regions$n_snp)
}

# regions annotated by JASPAR transcription factor binding sites
JASPAR_annotation <- read.table(args$jaspar_annotated_regions, sep="\t", header=TRUE)
jaspar_internal = as.character(unique(JASPAR_annotation$name))
tf_names  = as.character(sapply(jaspar_internal, function(x) pwm_list[x][[1]]@name))


# saved sequence relative probabilities
#seq_probabilities = relative_seq_probabilities(well_covered_regions) # run these two lines to regenerate
#save(seq_probabilities, file = "../data/sequence_probabilities.out")
attach("../data/sequence_probabilities.out")


## simulate data anbd annotate ##
if ( args$verbose ) {
  write("Simulating de novos in regulatory regions...", stderr())
}
sim_output = simulate_de_novos(well_covered_regions, seq_probabilities, snp_total, iterations)

# annotate simulation data with JASPAR transcription factor binding sites
if ( args$verbose ) {
  write("Annotating simulated de novos with predicted TFBS...", stderr())
}
annotations = jaspar_annotate(sim_output, JASPAR_annotation) 
sim_output = annotations[[1]] # adds another layer to sim_output with number of TFBS disruptions
sim_TFBS_tables = annotations[[2]] # list of tables with TFBS hits per simulation


## calculate summary statistics across TFBS count tables ##
if ( args$verbose ) {
  write("Calculating summary statistics (mean, SD, SE)...", stderr())
}
# get mean, SD, and SE on number of hits for each TF
empty_table = table(jaspar_internal) - 1 # used to standardize names on table
standardized = do.call(rbind, lapply(sim_TFBS_tables, function(t2) table_diff(t2, empty_table)))
sim_mean = apply(standardized, MARGIN=2, mean)
sim_sd = apply(standardized, MARGIN=2, sd)
sim_se = sim_sd/sqrt(iterations)

# save all of this data as a tab-delimited text file
TFBS_binding_results = data.frame("jaspar_internal" = jaspar_internal, "tf_names" = tf_names, "sim_mean" = sim_mean, "sim_sd" = sim_sd, "sim_se" = sim_se, row.names = NULL)
write.table(TFBS_binding_results, file = args$out, row.names = FALSE, sep = "\t", col.names = TRUE)
if ( args$verbose ) {
  write(sprintf("Finished! Simulated de novos analyzed by JASPAR. Stats aggregated and saved to: %s", args$out), stderr())
}




