# group the probands together which have de novos falling in the same predicted TF binding sites for analysis of HPO term similarity

source("tfbs_core.R")

# TO DO: add command line specifications for regions, de novos, and output path

# de novos are filtered for pp_dnm > 0.1 and non-coding only
de_novo_filtered <- read.table("../data/de_novo_filtered.txt", sep="\t", header=TRUE)

# regions are annotated by n_snps and n_indels - created by rupit/R/pre_process.R
well_covered_regions <- read.table("../data/regions_annotated.txt", sep="\t", header=TRUE)

# load all JASPAR pwm for homo sapiens
opts = list("species" = "Homo sapiens", "all_versions" = args$full_jaspar, "matrixtype" = "PWM")
pwm_list = getMatrixSet(JASPAR2014, opts)

# make new 'de novo containing regions' (dcr) data frame
dcr = well_covered_regions[well_covered_regions$n_snp > 0 | well_covered_regions$n_indel > 0, ]
dcr = merge(dcr, de_novo_filtered[, c("region_id", "pos", "ref", "alt", "person_stable_id")], by = "region_id")
dcr$rel_pos = dcr$pos - dcr$start

# add column specifying if proband is diagnosed or not
diagnosed <- read.table("../data/ddd_likely_diagnosed.txt", sep = "\t", header=TRUE)
dcr$diagnosed = dcr$person_stable_id %in% diagnosed$person_id

# get matrix (TF x dcr) with 1 if de novo in that region disrupts TF binding
ninety_five_counts = mapply(TFBS_hits, dcr$seq, dcr$rel_pos, MoreArgs = list("pwm_list" = pwm_list, "min.score" = "95%"))

# convert the list created above to a json string and save full set and diagnosed probands only
if (UNDIAGNOSED == TRUE){
  ud_proband_by_tfbs = apply(ninety_five_counts, MARGIN = 1, FUN = function(x) as.character(dcr[!dcr$diagnosed & x, "person_stable_id"]))
  save_to_json(ud_proband_by_tfbs, "../data/undiagnosed_probands_by_tfbs.json")
} else {
  proband_by_tfbs = apply(ninety_five_counts, MARGIN = 1, FUN = function(x) as.character(dcr[x, "person_stable_id"]))  
  save_to_json(proband_by_tfbs, "../data/probands_by_tfbs.json")
}
