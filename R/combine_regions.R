# combine regions after splitting and analyzing

library(optparse)

option_list <- list(
  make_option("--n_chunks", default = 20, help = "Choose the number of files to split into."),
  make_option("--base_name", default="../data/chunks/region_JASPAR_annotated_FULL",
              help="Set location to save the (likely large) tab-delimited text file storing JASPAR annotations.")
)

args <- parse_args(OptionParser(option_list=option_list))

chunks = vector("list", args$n_chunks)
for (i in seq(args$n_chunks)){
  fname = sprintf("%s.%i.txt", args$base_name, i)
  chunks[[i]] = read.table(fname, sep = "\t", header = TRUE)
}

JASPAR_annotation = do.call(rbind, chunks)
write.table(JASPAR_annotation, file = "../data/regions_JASPAR_annotated_FULL.txt", sep = "\t", col.names = TRUE, row.names = FALSE)
