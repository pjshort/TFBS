# TFBS core - functions for analyzing TFBS using JASPAR in de novo regulatory regions

# load libraries
library(TFBSTools)
library(JASPAR2014)
library(Biostrings)



get_TFBS_bp_coverage <- function(site_set){
  
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
  TFBS_count = sapply(site_seq_list, get_TFBS_bp_coverage)
  coverage_list = unlist(IRangesList(TFBS_count))
  x = coverage_list[1:length(coverage_list)-1]
  y = coverage_list[length(coverage_list)]
  u = union(x, y)
  return(sum(u@width)/nchar(as.character(seq)))
}

TFBS_hits <- function(seq, rel_pos, pwm_list, min.score = "95%"){
  
  # returns 0 or 1 for each TFBS in the pwm_list
  # one de novo may fall in multiple TFBS, or multiple instances of one TFBS (in this case, 1 still returned)
  
  subject = DNAString(seq)
  site_seq_list = searchSeq(pwm_list, subject, seqname="seq1", min.score=min.score, strand="*")
  TFBS_hits = sapply(site_seq_list, function(x) sum((rel_pos - x@views@ranges@start >= 0) & (x@views@ranges@start +x@views@ranges@width - rel_pos >= 0)) > 0)
}


# Functions to annotate regions with predicted transcription factor binding sites

single_sequence_coverage <- function(seq, pwm_list, min.score = "95%"){
  
  # returns the percentage of the sequence queried that is covered by a putative TFBS
  subject = DNAString(seq)
  site_seq_list = searchSeq(pwm_list, subject, seqname="seq1", min.score=min.score, strand="*")
  sequence_coverage = sapply(site_seq_list, get_TFBS_bp_coverage)
  coverage_list = unlist(IRangesList(sequence_coverage))
  return(coverage_list)
}

scan_regions <- function(regions, pwm_list, min.score = "95%"){
  
  # scan each region with pwm_list and return list of iranges intervals
  interval_list = sapply(regions$seq, function(s) single_sequence_coverage(s, pwm_list, min.score = min.score))
  names(interval_list) = regions$region_id
  return(interval_list)
}


# Functions to check overlap with PWM seq_list object

check_overlap <- function(intervals, rel_pos){
  
  # return any TF binding sites that rel_pos coincides with
  dist = rel_pos - intervals@start
  hits = which(dist >= 0 & dist < intervals@width)
  TFs = intervals@NAMES[hits]
  return(TFs)
}

regions_TFBS_overlap <- function(region_ids, positions, interval_list){
  
  # returns TF overlaps for every region_id, position pair that is passed
  TF_overlaps = mapply(function(r, p) check_overlap(interval_list[r][[1]], p), region_ids, positions)
  return(TF_overlaps)  
}
