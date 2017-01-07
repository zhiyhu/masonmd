# R package masonmd: making sense of nonsense-mediated decay
==============

Zhiyuan Hu

2016-12-22

## 1 Background
The package masonmd (**MA**ke **S**ense **O**f **NMD**) can be used to predict the effect of nonsense-mediated decay on mutated genes. While some mutations can introduce premature termination codons (PTCs) into genes, nonsense mediated decay (NMD) can detect these PTCs and then eliminate the abnormal transcripts. PTCs and NMD play important roles in genetic diseases and cancers.

Herein, three rules are used to predict whether a PTC-generating mutation is NMD-elicit or NMD-escape:
1. PTC is more than 50-54bp upstream of the last-exon-exon junction.
2. targeted gene is not intronless.
3. PTC is more than 200bp downstream of the start codon.

Using these rules, we can predict whether a called mutation will elicit NMD on the mRNA from the mutated gene, i.e. the NMD-elicit mutations. For now, this package can only be used for human genomes.

## 2 Installation instructions
In R or Rstudio, use the following codes to install the masonmd package directly from Github:
```{r install}
install.packages(“devtools”)
devtools::install_github("ZYBunnyHu/masonmd")
```

## 3 Examples
This package is easy to use, with its main function given by the function `classify.nmd`. It reads in the information of called mutations/variants, and return the prediction results. 

The following two example mutations are from The Cancer Genome Atlas (TCGA):

The first example is a substitution mutation on *AADAC* (Entrez ID = 13). The mutation happened on the 151545640 locus of the chromosome turning an G to a T. The genome reference is NCBI-build 37.

```{r example 1}
library(masonmd)
# an example of NMD-escape mutation from TCGA
classify.nmd(gene_id = 13, ref = 37, mut_start = 151545640, mut_end = 151545640,
ref_nt = "G",mut_nt = "T")
```

The first entry of the return object is `mut_nmd = F`, so it is not an NMD-elicit mutation. From the `note` and `have.ptc = T`. We know that the mutation caused a PTC but did not trigger NMD, because it does not fulfil the 50bp rule (Rule 1). The returned information also tells us:
1. whether the wildtype transcript is effected by NMD (`wt_nmd`)
2. and the relative position of PTC (`PTC.stop`), length of mutated coding sequence (`mutseq_length`), relative position of last exon-exon junction (`last_exon_exon_junction`) and number of exons (`n.exon`)

This information gives the underlying prediction results.

The second example is a NMD-elicit mutation in *A2M* gene (Entrez ID = 2).

```{r}
# an example of NMD-elicit mutation from TCGA
classify.nmd(gene_id = 2, ref = 37, mut_start = 9221429, mut_end = 9221429,
ref_nt = "G", mut_nt = "A")
```

