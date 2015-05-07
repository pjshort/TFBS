# script to run large set of simulations and annotate with JASPAR pipeline

# this script will create several files of relevant statistics and annotations
  # 1. n_TFBS X iterations dataframe saved as TFBS_simulation_annotation.txt
  # 2. summary stats with number of de novos overlapping a TFBS, number of TFBS overlapped by de novos
  # 3. n_snps X iterations dataframe with the number of TFBS overlapped by each de novo

# depends on:
source("../R/simulate.R")
source("../R/tfbs_core.R")
source("../R/visualization.R")
library(optparse)

# command line options
option_list <- list(
  make_option("--n_snps", default=NULL,
              help="Pass the genomic regions that should be annotated with predicted TF binding sites."),
  make_option("--iterations", default=100,
              help="Pass a list of TFs to be run against regions."),
  make_option("--mean_only", action="store_true", default=FALSE, help = "Output simulation TFBS binding mean, SD, SE only." ),
  make_option("--out_dir", default="../data/",
              help="Set location to save the (likely large) tab-delimited text file storing JASPAR annotations."),
  make_option("--jaspar_annotated_regions", default="../data/regions_JASPAR_annotated_FULL.txt",
              help="Pass a list of TFs to be run against regions."),
  make_option("--min_score", default="95%",
              help="Set min score for the JASPAR annotations - less than 95% may yield unrealistic number of predicted binding sites."),
  make_option("--verbose", action="store_true", default=TRUE,
              help="Print extra output advising the user of progression through the analysis.")
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


if ( args$mean_only) { # don't record result of every simulation
  
  sim_mean = apply(standardized, MARGIN=2, mean)
  sim_sd = apply(standardized, MARGIN=2, sd)
  sim_se = sim_sd/sqrt(iterations)
  
  # save all of this data as a tab-delimited text file
  TFBS_binding_results = data.frame("jaspar_internal" = colnames(standardized), "tf_names" = 
                                      as.character(sapply(colnames(standardized), function(x) pwm_list[x][[1]]@name)), 
                                    "sim_mean" = sim_mean, "sim_sd" = sim_sd, "sim_se" = sim_se, row.names = NULL)
  write.table(TFBS_binding_results, file = paste0(args$out_dir, "/TFBS_simulation_annotation.txt"), row.names = FALSE, sep = "\t", col.names = TRUE, quote = FALSE)
  if ( args$verbose ) {
    write(sprintf("Finished! Simulated de novos analyzed by JASPAR. Stats aggregated and saved to: %s", paste0(args$out_dir, "/TFBS_simulation_annotation.txt"), stderr()))
  }
  
} else { # record result of every simulation
  standardized = t(standardized)
  jaspar_internal = rownames(standardized)
  tf_names = as.character(sapply(jaspar_internal, function(x) pwm_list[x][[1]]@name))
  TFBS_binding_results = cbind(jaspar_internal, tf_names, standardized)
  write.table(TFBS_binding_results, file = paste0(args$out_dir, "/TFBS_simulation_annotation.txt"), row.names = FALSE, sep = ",", col.names = TRUE, quote = FALSE)
  if ( args$verbose ) {
    write(sprintf("Finished! Simulated de novos analyzed by JASPAR. Full results saved to: %s in COMMA SEPARATED FORMAT.", paste0(args$out_dir, "/TFBS_simulation_annotation.txt"), stderr()))
  }
}

sim_hits_per_de_novo = apply(sim_output[3,,], 2, as.numeric)  # n_snps X iterations with number of TFBS hit per de novo. user can >0 this if needed

write.table(sim_hits_per_de_novo, file = paste0(args$out_dir, "/TFBS_simulation_counts.txt"))

hits_total = colSums(sim_hits_per_de_novo)
hits_logical = colSums(sim_hits_per_de_novo > 0)

write.table(rbind(hits_total, hits_logical), file = paste0(args$out_dir, "/TFBS_global_counts.txt"), row.names = c("total_hits", "n_de_novos_with_TFBS_hit"), sep = ",", col.names = FALSE, quote = FALSE)





