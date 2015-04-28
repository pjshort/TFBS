# two ways we can simulate:
# 1. same as rupit (but alter to assign to specific position)
# 2. bootstrapping (i.e. draw random regions with replacement) and pick location of mutation based on null sequence model

source("/Users/ps14/code/rupit/R/regmut_null_model.R")

# load list of filtered de novos - run rupit_core.R function "filter_denovos" if not yet filtered
DENOVOS_PATH = "/Users/ps14/code/TFBS_analysis/data/de_novo_filtered.txt"
REGIONS_PATH = "/Users/ps14/code/TFBS_analysis/data/regions_annotated.txt"


relative_seq_probabilities <- function(regions){
  
  # get probability of mutation at each point in sequence (based on sequence context) and normalize
  seq_probabilities = sapply(as.character(regions$seq), function(s) p_position(s, normalize = TRUE), USE.NAMES = FALSE)
  names(seq_probabilities) = regions$region_id
  return(seq_probabilities)
}

assign_position <- function(regions, seq_probabilities, snp_total){
  
  # takes vector of region_ids to consider, list matching region_id to sequence relative probability, and total snps to simulate
  # samples snp_total region ids (with replacement). for each region, specific relative position will be randomly sampled
  region_ids = sample(regions$region_id, snp_total, replace = TRUE, prob = regions$p_relative)
  rel_pos = sapply(as.character(region_ids), function(id) sample(seq(1,length(seq_probabilities[id][[1]])), 1, prob = seq_probabilities[id][[1]]))
  return(rbind(names(rel_pos), as.integer(rel_pos)))  # named integer
}

simulate_de_novos <- function(regions, seq_probabilities, snp_total, iterations){
  
  # populates a sparse matrix with counts of de novos in each region for each simulation
  
  regions$p_relative = regions$p_snp_null/sum(regions$p_snp_null)
  
  # rows will be region and column will be count of de novos in each iteration
  sim_output = replicate(iterations, assign_position(regions, seq_probabilities, snp_total))
  return(sim_output)
}

scan_region <- function(pwm_list)

