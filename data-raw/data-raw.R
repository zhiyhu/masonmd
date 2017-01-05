library(TxDb.Hsapiens.UCSC.hg19.knownGene);library('BSgenome.Hsapiens.UCSC.hg19') # NCBI build 37
library(TxDb.Hsapiens.UCSC.hg18.knownGene);library('BSgenome.Hsapiens.UCSC.hg18') # NCBI build 36

txdb37 <- TxDb.Hsapiens.UCSC.hg19.knownGene# NCBI build 37
genome37 <- BSgenome.Hsapiens.UCSC.hg19# NCBI build 37

tx37 <- select(txdb37, keys = as.character(unique(keys(txdb37))),
               columns=c("GENEID","TXID","TXSTART","TXEND","EXONID","EXONCHROM",
                         "EXONSTRAND","EXONSTART","EXONEND","CDSEND" ,"CDSID","CDSSTART"), keytype="GENEID")

gene2tx37 <- select(txdb37, keys = as.character(unique(keys(txdb37))), columns="TXID", keytype="GENEID")
cds <- cdsBy(txdb37, by="tx" )
cds_seqs37 <- extractTranscriptSeqs(genome37, cds) # all CDS in NCBI build 37

# database of build 37, have all CDS
cds37 <- data.frame(txid = names(cds_seqs37),geneid = gene2tx37$GENEID[match(names(cds_seqs37),gene2tx37$TXID)],seq = as.character(cds_seqs37),width = width(cds_seqs37))
cds37 <- cds37[!is.na(cds37$geneid),]
cds37 <- cds37[!duplicated(paste(cds37$geneid,cds37$seq)),]

# build36

txdb36 <- TxDb.Hsapiens.UCSC.hg18.knownGene# NCBI build 37
genome36 <- BSgenome.Hsapiens.UCSC.hg18# NCBI build 37

tx36 <- select(txdb36, keys = as.character(unique(keys(txdb36))),
               columns=c("GENEID","TXID","TXSTART","TXEND","EXONID","EXONCHROM",
                         "EXONSTRAND","EXONSTART","EXONEND","CDSEND" ,"CDSID","CDSSTART"), keytype="GENEID")
tx36 <- tx36[!is.na(tx36$CDSSTART),]
gene2tx36 <- select(txdb36, keys = as.character(unique(keys(txdb36))), columns="TXID", keytype="GENEID")
cds <- cdsBy(txdb36, by="tx" )
cds_seqs36 <- extractTranscriptSeqs(genome36, cds) # all CDS in NCBI build 36

# database of build 36, have all CDS
cds36 <- data.frame(txid = names(cds_seqs36),geneid = gene2tx36$GENEID[match(names(cds_seqs36),gene2tx36$TXID)],seq = as.character(cds_seqs36),width = width(cds_seqs36))
cds36 <- cds36[!is.na(cds36$geneid),]
cds36 <- cds36[!duplicated(paste(cds36$geneid,cds36$seq)),]

devtools::use_data(tx37, tx36, cds37, cds36, genome36, genome37, compress = "xz", overwrite = T, internal = TRUE)
