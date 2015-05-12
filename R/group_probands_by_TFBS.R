# group the probands together which have de novos falling in the same predicted TF binding sites for analysis of HPO term similarity

source("tfbs_core.R")

# TO DO: add command line specifications for regions, de novos, and output path

# command line options
option_list <- list(
  make_option("--regions", default = "../data/regions_annotated.txt", help = "Set of regions de novos are drawn from"),
  make_option("--de_novos", default = "../data/de_novo_filtered.txt", help = "De novos matched with probands that should be scanned"),
  make_option("--tf_list", default=FALSE,
              help="Pass a list of TFs to be run against regions. Default is to use full set of JASPAR regions."),
  make_option("--diagnosed", default=FALSE,
              help="List of probands (IDs) that the JSON should be stratified by, not necessarily diagnosis. List of diagnosed probands, for instance, would produce three JSONs: 
              one with only diagnosed probands, one with only undiagnosed, one with full set. First row should be labelled \"person_id\""),
  make_option("--out_dir", default="../data/",
              help="Set location to save (possibly multiple) JSON files."),
  make_option("--verbose", action="store_true", default=TRUE,
              help="Print extra output advising the user of progression through the analysis.")
)

args <- parse_args(OptionParser(option_list=option_list))

args$diagnosed = "../data/ddd_likely_diagnosed.txt"

# de novos are filtered for pp_dnm > 0.1 and non-coding only
de_novo_filtered <- read.table(args$de_novos, sep="\t", header=TRUE)

# regions are annotated by n_snps and n_indels - created by rupit/R/pre_process.R
well_covered_regions <- read.table(args$regions, sep="\t", header=TRUE)

# load all JASPAR pwm for homo sapiens
opts = list("species" = "Homo sapiens", "all_versions" = TRUE, "matrixtype" = "PWM")
pwm_list = getMatrixSet(JASPAR2014, opts)

if (args$tf_list != FALSE){  # user wants to use a smaller set of TFs
  # load the names of TFs that we want to scan for (must be JASPAR internal names)
  TFBS_to_scan = read.table(args$tf_list, header = TRUE)
  pwm_list = pwm_list[unique(TFBS_to_scan$jaspar_internal)]
}

# make new 'de novo containing regions' (dcr) data frame
dcr = well_covered_regions[well_covered_regions$n_snp > 0 | well_covered_regions$n_indel > 0, ]
dcr = merge(dcr, de_novo_filtered[, c("region_id", "pos", "ref", "alt", "person_stable_id")], by = "region_id")
dcr$rel_pos = dcr$pos - dcr$start + 1

# get matrix (de_novo x TF) with 1 if de novo in that region disrupts TF binding
binding_disruption_mat = matrix(0, length(dcr_overlaps), length(names(pwm_list)))
colnames(binding_disruption_mat) <- names(pwm_list)

# analyze each de novo and fill list with TFs it disrupts (if any)
dcr_overlaps = mapply(regions_TFBS_overlap, dcr$region_id, dcr$rel_pos, MoreArgs = list("JASPAR_annotation" = JASPAR_annotation))
dcr$tfbs_ids = sapply(dcr_overlaps, function(s) ifelse(length(s[[1]]) > 0, paste(s, collapse =","), ""))

# dcr_overlaps is list (one item for each de novo) - loop over to fill de_novo X TF matrix initialized above
for (i in seq(length(dcr_overlaps))){
  idx = which(names(pwm_list) %in% dcr_overlaps[[i]])
  binding_disruption_mat[i,idx] = 1
}

# convert the list created above to a json string and save full set and diagnosed probands only

# if command line opt is true, add column specifying if proband is diagnosed or not and split into two extra jsons on condition
if (args$diagnosed != FALSE){
  diagnosed <- read.table(args$diagnosed, sep = "\t", header=TRUE)
  dcr$diagnosed = dcr$person_stable_id %in% diagnosed$person_id
  
  undiag_proband_by_tfbs = apply(binding_disruption_mat, MARGIN = 2, FUN = function(x) as.character(dcr[!dcr$diagnosed & (x>0), "person_stable_id"]))
  save_to_json(undiag_proband_by_tfbs, "../data/probands_by_tfbs_UNDIAGNOSED.json")
  
  diag_proband_by_tfbs = apply(binding_disruption_mat, MARGIN = 2, FUN = function(x) as.character(dcr[dcr$diagnosed & (x>0), "person_stable_id"]))
  save_to_json(diag_proband_by_tfbs, "../data/probands_by_tfbs_DIAGNOSED.json") 
}

# save full (unsplit by condition)
proband_by_tfbs = apply(binding_disruption_mat, MARGIN = 2, FUN = function(x) as.character(dcr[which(x > 0), "person_stable_id"]))  
save_to_json(proband_by_tfbs, "../data/probands_by_tfbs.json")
