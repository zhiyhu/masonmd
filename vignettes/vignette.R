## ----example 1-----------------------------------------------------------
library(masonmd)
# an example of NMD-escape mutation from TCGA
classify.nmd(gene_id = 13, ref = 37, mut_start = 151545640, mut_end = 151545640,
ref_nt = "G",mut_nt = "T")

## ------------------------------------------------------------------------
# an example of NMD-elicit mutation from TCGA
classify.nmd(gene_id = 2, ref = 37, mut_start = 9221429, mut_end = 9221429,
ref_nt = "G", mut_nt = "A")

