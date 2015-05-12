# functions to do score comparison under changes to PWM

# original sequence is given by the seqList object
# relative position and ref, alt from dcr will tell us which one to change
# pwm will tell us the score difference

binding_change <- function(jaspar_internal, rel_pos, seq, ref, alt, pwm_list, min.score = "95%"){
  
  # takes PWM name, rel position, ref, alt and returns absolute score difference
  # at the moment, works ONLY for snps (not indels)
  
  # get the pwm for single TF (jaspar_internal)
  single_pwm = pwm_list[[jaspar_internal]]

  # scan the sequence with this pwm
  subject = DNAString(seq)
  site_seq_list = searchSeq(single_pwm, subject, seqname="seq1", min.score=min.score, strand="*")
  
  # get the position within the motif
  starts = site_seq_list@views@ranges@start
  stops = site_seq_list@views@ranges@start + site_seq_list@views@ranges@width - 1
  intersects = rel_pos >= starts & rel_pos <= stops
  motif_pos = rel_pos - starts[intersects] + 1
      
  # use profileMatrix to compute score change from ref to alt at single position
  pwm = single_pwm@profileMatrix
  change = as.numeric(pwm[alt, motif_pos] - pwm[ref, motif_pos])
  
  ref_score = site_seq_list@score[intersects]
  alt_score = ref_score + change
  
  return(cbind(ref_score, alt_score))
}