# TFBS core - functions for analyzing TFBS using JASPAR in de novo regulatory regions

# load libraries
library(TFBSTools)
library(JASPAR2014)
library(Biostrings)
library(jsonlite)

get_single_site_bp_coverage <- function(site_set){
  
  # helper function for TFBS_bp_coverage
  
  starts = site_set@views@ranges@start
  ends = site_set@views@ranges@width + starts
  x = IRanges(starts[1:length(starts) - 1], ends[1:length(ends) - 1])
  y = IRanges(starts[length(starts)], ends[length(starts)])
  u = union(x,y)
  return(u)
}

TFBS_bp_coverage <- function(seq, pwm_list, min.score = "95%"){
  
  # returns the percentage of the sequence queried that is covered by a putative TFBS
  
  subject = DNAString(seq)
  site_seq_list = searchSeq(pwm_list, subject, seqname="seq1", min.score=min.score, strand="*")
  TFBS_count = sapply(site_seq_list, get_single_site_bp_coverage)
  coverage_list = unlist(IRangesList(TFBS_count))
  x = coverage_list[1:length(coverage_list)-1]
  y = coverage_list[length(coverage_list)]
  u = union(x, y)
  return(sum(u@width)/nchar(as.character(seq)))
}

TFBS_hits <- function(seq, rel_pos, pwm_list, min.score = "95%"){
  
  # returns 0 or 1 for each TFBS in the pwm_list
  # one de novo may fall in multiple TFBS, or multiple instances of one TFBS (in this case eg highly repeated sections, 1 is still returned)
  
  subject = DNAString(seq)
  site_seq_list = searchSeq(pwm_list, subject, seqname="seq1", min.score=min.score, strand="*")
  TFBS_hits = sapply(site_seq_list, function(x) sum((rel_pos - x@views@ranges@start >= 0) & (x@views@ranges@start +x@views@ranges@width - rel_pos >= 0)) > 0)
}


## Functions to annotate regions with predicted transcription factor binding sites

single_sequence_coverage <- function(seq, pwm_list, min.score = "95%"){
  
  # returns the percentage of the sequence queried that is covered by a putative TFBS
  subject = DNAString(seq)
  site_seq_list = searchSeq(pwm_list, subject, seqname="seq1", min.score=min.score, strand="*")
  sequence_coverage = sapply(site_seq_list, get_single_site_bp_coverage)
  coverage_list = unlist(IRangesList(sequence_coverage))
  return(coverage_list)
}

scan_regions <- function(regions, pwm_list, min.score = "95%"){
  
  # scan each region with pwm_list and return list of iranges intervals
  interval_list = sapply(regions$seq, function(s) single_sequence_coverage(s, pwm_list, min.score = min.score))
  names(interval_list) = regions$region_id
  return(interval_list)
}


## Functions to check overlap with PWM seq_list object

check_overlap <- function(region_id, pos, JASPAR_annotation){
  
  # return any TF binding sites that pos coincides with
  region_slice = JASPAR_annotation[JASPAR_annotation$region_id == region_id, ]
  pos_match = region_slice[region_slice$start <= pos & region_slice$stop >= pos, ]
  return(as.character(pos_match$name))
}

regions_TFBS_overlap <- function(region_ids, positions, JASPAR_annotation){
  
  # returns TF overlaps for every region_id, position pair that is passed. must be relative position!
  JASPAR_annotation$region_id <- factor(JASPAR_annotation$region_id, levels = levels(region_ids))
  JA <- JASPAR_annotation[!is.na(JASPAR_annotation$region_id),]
  TF_overlaps = mapply(check_overlap, region_ids, positions, MoreArgs = list("JASPAR_annotation" = JA))
  return(TF_overlaps)
}

save_to_json <- function(l, fname){
  json = toJSON(l, pretty = TRUE)
  save_name = sprintf(fname)
  sink(save_name)
  cat(json)
  sink()
}


## merging and subtracting tables with different column spaces (by coercing to same column space) ##

table_merge <- function(table1, table2, row.names = c("table1", "table2")){
  all_names = unique(c(names(table1), names(table2)))
  
  # set up two new tables with shared names
  new_table1 = table(row.names = all_names) - 1
  new_table2 = table(row.names = all_names) - 1
  
  # set nonzero entries from original tables
  new_table1[all_names %in% names(table1)] = table1
  new_table2[all_names %in% names(table2)] = table2
  m = rbind(new_table1, new_table2)
  row.names(m) = row.names
  
  return(m)
}

table_diff <- function(table1, table2){
  
  # subtracts table2 from table1 and unifies the column names
  
  m = table_merge(table1, table2)
  diff = m[1,] - m[2,]
  return(diff)
}
