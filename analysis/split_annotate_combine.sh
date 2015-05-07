#!/bin/bash

### workflow to split the well_covered_regions file into 20 different files of sub regions and annotate each one separately 
with JASPAR binding sites. this can be done now with the FULL list adding the --full_jaspar=TRUE tag.###

cd TFBS/R

# split regions into chunks. base_name is the folder to save all of the chunks to
/software/R-3.1.2/bin/Rscript split_regions.R --regions="../data/regions_annotated.txt" --n_chunks="20" --base_name="../data/region_chunk"

# annotate the regions with JASPAR predicted binding
bsub -J "jaspar_annotation[1-20]" -R'select[mem>260] rusage[mem=260]' -M260 -o "/nfs/users/nfs_p/ps14/experiments/jaspar/jaspar.log.%I" \
/software/R-3.1.2/bin/Rscript /nfs/users/nfs_p/ps14/software/TFBS/R/JASPAR_annotation.R \
--regions=/nfs/users/nfs_p/ps14/software/TFBS/data/region_chunk.\$LSB_JOBINDEX.txt \
--out=/nfs/users/nfs_p/ps14/experiments/jaspar/region_JASPAR_annotated.\$LSB_JOBINDEX.txt \
--tf_list=/nfs/users/nfs_p/ps14/software/TFBS/data/TFs_in_DDD_data.txt --verbose

# combine all of the regions back into a single file (saving into base_name as "regions_JASPAR_annotated_FULL.txt")
/software/R-3.1.2/bin/Rscript combine_regions.R --n_chunks="20" --base_name="../data/region_chunk"


