# r script to split regions file into specified number of chunks
library(optparse)

option_list <- list(
  make_option("--regions", default="../data/regions_annotated.txt",
              help="Pass the genomic regions that should be annotated with predicted TF binding sites."),
  make_option("--n_chunks", default = 10, help = "Choose the number of files to split into."),
  make_option("--base_name", default="../data/region_chunk",
              help="Set location to save the (likely large) tab-delimited text file storing JASPAR annotations.")
)

args <- parse_args(OptionParser(option_list=option_list))

well_covered_regions <- read.table(as.character(args$regions), sep="\t", header=TRUE)
header = names(well_covered_regions)

break_pts = seq(1, nrow(well_covered_regions), by = floor(nrow(well_covered_regions)/(args$n_chunks)))
break_pts[length(break_pts)] = nrow(well_covered_regions) + 1

rfor (i in seq(length(break_pts)-1)){
  fname = sprintf("%s.%i.txt", args$base_name, i)
  write.table(well_covered_regions[break_pts[i]:(break_pts[i+1] - 1),], file = fname, col.names = header, row.names = FALSE, sep="\t")
}