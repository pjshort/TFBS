# two ways we can simulate:
# 1. same as rupit (but alter to assign to specific position)
# 2. bootstrapping (i.e. draw random regions with replacement) and pick location of mutation based on null sequence model

source("../R/regmut_null_model.R") # respects the hierarchy as organized on the farm and locally
source("../R/tfbs_core.R")
library(abind)

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

simulation_TFBS_overlap <- function(s, JASPAR_annotation){
  
  # takes slice of sim_output called s and returns TF overlaps from JASPAR annotation provided
  region_ids = s[1,]
  positions = as.integer(s[2,])
  # returns TF overlaps for every region_id, position pair that is passed
  JA <- JASPAR_annotation[JASPAR_annotation$region_id %in% region_ids,]
  JA$region_id <- factor(JA$region_id)  
  TF_overlaps = mapply(check_overlap, region_ids, positions, MoreArgs = list("JASPAR_annotation" = JA))
  counts_by_TFBS = table(unlist(TF_overlaps))
  return(TF_overlaps)
}

jaspar_annotate <- function(sim_output, JASPAR_annotation){
  
  # take sim_output matrix (2 x n_snps x iterations) and add another layer (n_snps x iterations) which describes the
  # number of TFBS each simulated snp disrupts
  
  # sim_output is a 2 x n_snps x iterations matrix. row 1 is region_id, row 2 is relative position.
  tf_out <- apply(sim_output, MARGIN = 3, simulation_TFBS_overlap, "JASPAR_annotation" = JASPAR_annotation)
  list_of_counts = sapply(tf_out, function(x) sapply(x, length))
  row.names(list_of_counts) <- NULL  
  
  # TODO: sort out strange behavior where sim_output[3,,] has named rows - possibly because matrix must have 
  # all characters or all integers ( no mixed types)?
  # add a third row to sim_output with number of TFBS disruptions per de_novo
  sim_output = abind(sim_output, list_of_counts, along = 1)
  
  # get a list of tables of counts per TFBS
  sim_TFBS_counts = lapply(tf_out, function(s) table(unlist(s)))
  
  return(list("sim_output" = sim_output, "sim_TFBS_counts" = sim_TFBS_counts))
}




