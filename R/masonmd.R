# R functions for predicting whether a mutation can elicit Nonsense-Mediated Decay
# zhiyuan@well.ox.ac.uk (c) 2017


##--------------------------------------------------------------------------------------------
#' @title Classify NMD-elicit mutations
#'
#' @description
#' \code{classify.nmd} can detect the premature stop codon (PTC) and predict NMD-elicit mutations
#' based on three rules.
#'
#' @details
#' classify the mutations into NMD or non-NMD, and return the position of premature stop codon in CDS
#'
#' @param gene_id numeric; gene Entrez ID
#' @param ref numeric; reference genome/NCBI build; \code{ref = 37} if the mutation is called by
#' build37/hg19, or \code{ref = 36} if the mutation is called by build36/hg18
#' @param mut_start numeric; the absolution start position of the mutations on chromosome on + strand
#' @param mut_end numeric; the absolution end position of the mutations on chromosome on + strand
#' @param ref_nt charactor; the reference nucleotides of the mutations on + strand;
#' \code{ref_nt = "-"} for insertion
#' @param mut_nt charactor; the mutation nucleotides on + strand; \code{mut_nt = "-"} for deletion
#'
#' @export
#'
#' @return A list containing eight entries:
#' \itemize{
#' \item{\code{mut_nmd}} {whether the mutation is NMD-elicit.}
#' \item{\code{note}} {extra explanation of the classification results.}
#' \item{\code{wt_nmd}} {whether the wildtype CDS harbours PTC that is NMD-elicit.}
#' \item{\code{PTC.Stop}} {the position of PTC.}
#' \item{\code{have.ptc}} {if the mutated sequence harbours a PTC.}
#' \item{\code{mutseq_length}} {length of the mutated coding sequences}
#' \item{\code{last_exon_exon_junction}} {relative position of last exon-exon junction}
#' \item{\code{n.exon}} {number of exons in the gene}
#' }
#' @examples
#' library(masonmd)
#' # an example of NMD-escape mutation from TCGA
#' classify.nmd(gene_id = 13, ref = 37, mut_start = 151545640, mut_end = 151545640,
#' ref_nt = "G",mut_nt = "T")
#' # an example of NMD-elicit mutation from TCGA
#' classify.nmd(gene_id = 2, ref = 37, mut_start = 9221429, mut_end = 9221429,
#' ref_nt = "G", mut_nt = "A")
classify.nmd <- function(gene_id, ref = 37, mut_start, mut_end, ref_nt = "-", mut_nt = "-")
{
    if (ref == 19) {
        ref <- 37
    }
    if (ref == 18) {
        ref <- 36
    }
    # get further info on the mutated sequence based on given infomation
    mut_info <- mut.info(gene_id, ref, mut_start, mut_end, ref_nt, mut_nt)
    mut_seq <- as.character(mut_info[1]) # mutated seq
    wt_seq <- as.character(mut_info[4]) # wildtype seq
    n.exon <- as.numeric(mut_info[5]) # number of exon
    last_exon_exon_junction <- as.numeric(mut_info[2]) # position of last junction
    nmd.detail <- NA # classify results: mutated seq
    nmd.detail2 <- NA # classify results: wildtype seq
    note1 <- as.character(mut_info[3])
    note <- ""
    pos.ptc <- NA # position of ptc
    have.ptc <- F
    length <- NA
    if (!is.na(mut_seq)) {
        length <- nchar(mut_seq)
        if (nchar(mut_seq) > 0) { # it has mutated sequences
            start_pos = findStartPos(seq = mut_seq) # find the start position
            stop_pos = findStopPos(seq = mut_seq) # find the stop position in the frame
            if (length(start_pos) == 0) { # cannot find start codon
                note <- "no start codon"
                warning("Can't find start codon in the mutated coding sequence")
            } else if (length(stop_pos) == 0) { # cannot find start codon
                note <- "no stop codon"
                warning("Can't find stop codon in the mutated coding sequence")
            } else {  # find start and stop codon
                start_pos = min(start_pos) # choose the cloest as the start position
                ORF <- findORF(start_pos = start_pos,
                               k = 0,
                               stop_pos = stop_pos)

                if (nrow(ORF) == 0) {
                    nmd.detail <- F
                    note <- "cannot find ORF"
                    warning("Can't find ORF in the mutated coding sequence")
                    pos.ptc <- Inf # no stop codon
                } else { # find ORF
                    ORF <-  ORF[which.min(ORF[,"Stop"]),] # the closest stop codon
                    pos.ptc <- ORF["Stop"]
                    if (pos.ptc > nchar(mut_seq) - 3) { # normal coding sequence no PTC
                        have.ptc <- F
                        note <- "normal sc"
                        nmd.detail <- F
                    } else {
                        if (pos.ptc <= nchar(mut_seq) - 3) {
                            have.ptc <- T
                        }
                        if (n.exon < 1.1) { # have only one exon
                            nmd.detail <- F
                            note <- "single exon"
                            warning("The gene only has one exon")
                        } else if (pos.ptc < (last_exon_exon_junction - 50) & pos.ptc > 200) {
                            nmd.detail <- T
                        } else if (pos.ptc >= (last_exon_exon_junction - 50) &
                                   pos.ptc < (nchar(mut_seq) - 2)) {
                            nmd.detail <- F
                            note <- "PTC doesn't not fulfill 50 bp rule"
                        } else if (pos.ptc < 200) {
                            nmd.detail <- F
                            note <- "PTC too close to start codon"
                        } else if (is.infinite(pos.ptc)) {
                            nmd.detail <- F
                            note <- "nonstop mutations; no stop codon found"
                        } else {
                            nmd.detail <- F
                            note <- "unclear; can't be classified"
                        }
                    } # have PTC
                }
            }
        }
        # find the new start and stop codons in the sequence
        if(nchar(wt_seq) > 0) {# it has wt sequences

            start_pos <- findStartPos(seq = wt_seq) # find the start position
            stop_pos <- findStopPos(seq = wt_seq)# find the stop position in the frame

            if (length(start_pos) == 0) {# cannot find start codon
                note <- paste(note,"no start codon in wt")
            } else if (length(stop_pos) == 0) {# cannot find stop codon
                note <- paste(note, "no stop codon in wt")
            } else {# find start and stop codon

                start_pos <- min(start_pos) # choose the cloest as the start position
                if (start_pos != 1) {
                    note <- paste(note, "wt cds not start with start codon")
                }
                ORF <- findORF(start_pos = start_pos,k = 0,stop_pos = stop_pos)

                if (nrow(ORF) == 0) {
                    nmd.detail2 <- F
                    note <- paste(note, "can't find ORF in wt")
                } else { # find ORF in WT
                    ORF <- ORF[which.min(ORF[, "Start"]),]

                    if (ORF["Stop"] == nchar(wt_seq)) {
                        nmd.detail2 <- T
                        note <- paste(note, "PTC in wt")
                    } else {
                        nmd.detail2 <- F
                    }
                }
            }
        }
    }
    note <- paste(note1, note)
    return(c(mut_nmd = nmd.detail,
             note = note,
             wt_nmd = nmd.detail2,
             PTC = pos.ptc,
             have.ptc = have.ptc,
             mutseq_length = length,
             last_exon_exon_junction = last_exon_exon_junction,
             n.exon = n.exon))
}

##--------------------------------------------------------------------------------------------
#' @title Get the information of mutated sequence
#'
#' @description
#' \code{mut.info} returns the sequences and positional infomation of a mutated gene.
#'
#' @details Returns the sequences and positional infomation of a mutated gene.
#'
#' @inheritParams classify.nmd
#'
#' @import Biostrings
#' @export
#' @return charactor vector containing five entries:
#' \itemize{
#' \item{\code{mut_seq}} {charactor; coding sequence of the mutated gene.}
#' \item{\code{last.exexjun}} {numeric; position of last exon-exon junction.}
#' \item{\code{note}} {extra explanation.}
#' \item{\code{wt_seq}} {charactor; coding sequence of the wildtype gene.}
#' \item{\code{n.exon}} {numeric; number of exons in the gene.}
#' }
#' @examples
#' library(masonmd)
#' # an example of NMD-escape mutation from TCGA
#' mut.info(gene_id = 13, ref = 37, mut_start = 151545640, mut_end = 151545640,
#' ref_nt = "G",mut_nt = "T")
mut.info <- function(gene_id, ref, mut_start, mut_end, ref_nt, mut_nt)
{
    cds <- get.longest.cds(gene_id,ref)
    tx <- cds$tx
    txid <- cds$txid
    wt_seq <- cds$seq
    n.exon <- nrow(tx[!is.na(tx$CDSSTART),])
    wt_pos <- get_tx_seq(tx,ncbi_build = ref)
    wt_exjun <- wt_pos$exonJunction
    wt_pos <- wt_pos$wt_pos
    strand <- as.character(unique(tx$EXONSTRAND))
    note<- ""
    mut_seq <- NA
    mut_exjun <- NA
    if (length(cds$gene) > 0) {

        if (strand == "+") { # gene on the + strand
            ref_nt <- as.character(ref_nt)
            mut_nt <- as.character(mut_nt)
            mut_start <- as.numeric(mut_start)
            mut_end <- as.numeric(mut_end)
        }
        if (strand  == "-") {
            ref_nt <- as.character(ref_nt)
            mut_nt <- as.character(mut_nt)
            mut_start <- as.numeric(mut_start)
            mut_end <- as.numeric(mut_end)
            ref_nt <- as.character(reverseComplement(DNAString(ref_nt)))
            mut_nt <- as.character(reverseComplement(DNAString(mut_nt)))
        }
        if (!any(mut_start:mut_end %in% min(wt_pos):max(wt_pos))) {
            # mutation out of the transcript
            note <- "mutation outside"
            mut_seq <- wt_seq
            mut_exjun <- wt_exjun

        } else if (all(mut_start:mut_end %in% min(wt_pos):max(wt_pos)) &
                   (!any(mut_start:mut_end %in% wt_pos))) {
            # mutation position in intron

            note <- "mutation in intron"
            mut_seq <- wt_seq
            mut_exjun <- wt_exjun
        } else if (all(mut_start:mut_end %in% wt_pos)) {
            # mutation in CDS

            mut_start <- which(wt_pos == mut_start)
            mut_end <- which(wt_pos == mut_end)
            note<-""
        } else { # part of the mutation in exon
            note  <-  "mutation partially in cds; can't classify."
            mut_seq <- NA
            mut_exjun <- NA
        }

        if (note=="" ) { # check if mutation is on splice site
            mut_start <- c(mut_start,mut_end)
            mut_end <- max(mut_start)
            mut_start <- min(mut_start)
            # use function
            mut_exjun <- get_mut_seq_exjun(wt_exjun = wt_exjun, wt_seq = wt_seq,
                                           mut_nt = mut_nt, mut_start = mut_start,
                                           mut_end = mut_end, ref_nt = ref_nt,
                                           txid = txid, note = note)
            note <- mut_exjun$note
            mut_seq <- mut_exjun$seq
            mut_exjun <- mut_exjun$exjun
        } else {
            mut_seq <- NA
            mut_exjun <- NA
        }
    }
    if (length(note) > 1) {
        note <- paste(note, collapse = "; ")
    }
    note <- gsub("  ", " ", note, fixed = TRUE)

    if (!is.na(mut_seq)) {
        last.exexjun <- last.exexjun.pos(mut_exjun)
        if (last.exexjun$note != "") {
            note <- paste(note, last.exexjun$note, sep="; ")
        }
        last.exexjun <- last.exexjun$last_exon_exon_junction
    } else {
        last.exexjun <- NA
    }
    return(c(mut_seq = mut_seq,
             last.exexjun = last.exexjun,
             note = note,
             wt_seq = wt_seq,
             n.exon = n.exon))
}

##--------------------------------------------------------------------------------------------
#' @title Get position of last exon-exon junction
#'
#' @description
#' \code{last.exexjun.pos} gets position of last exon-exon junction.
#'
#' @details Get position of last exon-exon junction.
#'
#' @param mut_exjun a numeric vector that identifies the positions of exon-exon junctions
#' in the coding sequence. It is the output of the function \code{get_mut_seq_exjun} in masonmd package.
#'
#' @import Biostrings
#' @export
#' @return A list containing two entries:
#' \itemize{
#' \item{\code{last_exon_exon_junction}} {numeric; position of last exon-exon junction.}
#' \item{\code{note}} {note on weather the mutation affects splice sites or the gene has only one exon.}
#' }
last.exexjun.pos <- function(mut_exjun)
{
    note <- ""
    if ((sum(mut_exjun == 1)) %in% ((1:1000) * 2)) {

        last_exon_exon_junction <- which(mut_exjun == 1)[length(which(mut_exjun == 1)) - 1]
        note <- ""
    } else if (sum(mut_exjun == 1) %in% ((1:1000) * 2 + 1)) {
        last_exon_exon_junction <- which(mut_exjun == 1)[length(which(mut_exjun == 1)) - 1]
        note <- "splice site mutated"
        warning("Splice site mutated")
    } else if (sum(mut_exjun == 1) == 1) {
        last_exon_exon_junction <- which(mut_exjun==1)
        note<- "splice site mut"
        warning("Splice site mutated")
    } else if (sum(mut_exjun == 1) == 0) {
        last_exon_exon_junction <- length(mut_exjun)
        note <- "single exon"
        warning("Single-exon gene")
    }
    return(list(last_exon_exon_junction = last_exon_exon_junction,
                note = note))
}

##--------------------------------------------------------------------------------------------
#' @title Get information of longest isoform for given gene
#'
#' @description
#' \code{sum} returns the sum of all the values present in its arguments.
#'
#' @details Get the longest cds if the gene ID is given
#'
#' @param gene_id numeric; Entrez gene ID
#' @param numeric; reference genome/NCBI build; \code{ref = 37} if the mutation is called by
#' build37/hg19, or \code{ref = 36} if the mutation is called by build36/hg18
#'
#' @export
#'
#' @return A list containing four entries:
#' \itemize{
#' \item{\code{gene}} {gene Entrez ID; same with the input}
#' \item{\code{txid}} {transcript ID of the longest isoform}
#' \item{\code{seq}} {charactor; coding sequence of the transcript}
#' \item{\code{tx}} {infomation of the transcript in form of Txdb}
#' }
#' @examples
#' library(masonmd)
#' get.longest.cds(gene_id = 13, ref = 37)
get.longest.cds <- function(gene_id, ref) # combine
{
    #get gene annotations
    if (ref == 37) { # get the transcripts of the genes
        tx <- tx37[tx37$GENEID == gene_id,]
        cds <- cds37[cds37$geneid == gene_id,]
    } else if (ref == 36) {
        tx <- tx36[tx36$GENEID == gene_id,]
        cds <- cds36[cds36$geneid == gene_id,]
    }
    # get the longest cds
    cds <- cds[which.max(cds$width),]
    tx <- tx[tx$TXID == cds$txid,]

    return(list(gene = as.character(cds$geneid), # Geneid
                txid = as.character(cds$txid), # Txid
                seq = as.character(cds$seq), # CDS
                tx = tx)) # transcript info
}

##--------------------------------------------------------------------------------------------
#' @title Get relative positions of exon-exon junctions in mutated gene
#'
#' @description \code{get_mut_exjun} refers relative positions of exon-exon junctions in mutated gene.
#'
#' @details Indicate the positions of exon-exon junctions in the mutated sequence based on the mutation
#' information and the positions of exon-exon junctions in wildtype sequence.
#'
#' @param wt_exjun a string that consists of 0 or 1 to indicates the positions of exon-exon junctions in
#' the coding sequence;
#' @param mut_start position of mutation start in chromosome
#' @param mut_end position of mutation end in chromosome
#' @param mut_nt nucleotides after alteration
#' @export
#'
#' @return A list containing two entries:
#' \itemize{
#' \item{\code{mut_exjun}} {charactor that consists of 0 and 1; 1 represents the position of junctions}
#' \item{\code{note}} {extra note}
#' }
get_mut_exjun <- function(wt_exjun, mut_start, mut_end, mut_nt)
{
    note <- ""
    if (any(wt_exjun[mut_start:mut_end] == 1)) {# mutation on the exon-intron junction
        note <- "mutated exon intron junction"
        warning("Exon-intron junction is mutated")
        mut_exjun <- c(wt_exjun[1:(mut_start - 1)],
                       rep(-1, length(mut_nt)),
                       wt_exjun[(mut_end + 1):length(wt_exjun)])
    } else if (mut_start == 1) {
        mut_exjun <- c(rep(-1, length(mut_nt)), wt_exjun[(mut_end + 1):length(wt_exjun)])
        note <- "mutated start"
    } else if (mut_end == length(wt_exjun)) {
        mut_exjun <- c(wt_exjun[1:(mut_start - 1)],rep(-1, length(mut_nt)))
        note <- "mutated end"
    } else if (all(mut_start:mut_end %in% 1:length(wt_exjun)) &
               all(wt_exjun[mut_start:mut_end] != 1)) {
        # mutation is not at the start or the end
        mut_exjun <- c(wt_exjun[1:(mut_start - 1)],
                       rep(-1, length(mut_nt)),
                       wt_exjun[(mut_end + 1):length(wt_exjun)])
    } else {
        mut_exjun <- ""
        note <- "unknown mutation"
    }
    return(list(mut_exjun, note))
}

##--------------------------------------------------------------------------------------------
#' @title Get mutated coding sequence and junctions
#'
#' @description Get mutated coding sequence and junctions
#'
#' @details Get the mutated sequence and junctions
#' @inheritParams get_mut_exjun
#' @param wt_seq wildtype coding sequence
#' @param ref_nt the nucleotides on the reference genome at the mutated position; the nucleotides before
#' mutation
#' @param note extra notes
#' @param txid transcript ID
#'
#' @export
#'
#' @return A list containing three entries:
#' \itemize{
#' \item{\code{mut_exjun}} {numeric string with 0 and 1; 1 represents the position of junctions}
#' \item{\code{note}} {extra note}
#' }
get_mut_seq_exjun <- function(wt_exjun, wt_seq, mut_start, mut_end, mut_nt, ref_nt, note, txid)
{
    mut_seq = NA
    mut_exjun = NA
    if (mut_nt == "-" & ref_nt == "-") {
        warning("Not mutation input")
    } else if(mut_nt == "-") { # deletion
        if(ref_nt == substr(wt_seq, mut_start, mut_end)) {
            mut_seq <- paste(c(substr(wt_seq, 1, mut_start - 1),
                               substr(wt_seq, mut_end + 1, nchar(wt_seq))),
                             collapse= "")
        } else {
            note <- "reference unmatched"
        }

        if (mut_start >1 & mut_end < length(wt_exjun)) {
            # mutation is not at the start or the end
            mut_exjun <- c(wt_exjun[1:(mut_start - 1)],
                           wt_exjun[(mut_end + 1):length(wt_exjun)])
        } else if (mut_start == 1) {
            mut_exjun <- c(wt_exjun[(mut_end + 1):length(wt_exjun)])
        } else if (mut_end == length(wt_exjun)) {
            mut_exjun <- c(wt_exjun[1:(mut_start - 1)])
        }
    } else if (ref_nt == "-") { # insertion
        mut_seq <- paste(c(substr(wt_seq, 1, mut_start),
                           mut_nt,
                           substr(wt_seq, mut_end, nchar(wt_seq))),
                         collapse= "")
        if (all(mut_start:mut_end %in% 1:length(wt_exjun))) {
            mut_exjun <- get_mut_exjun(wt_exjun = wt_exjun,
                                       mut_nt = mut_nt,
                                       mut_start = mut_start,
                                       mut_end = mut_end)
        }
        mut_exjun <- mut_exjun[[1]]
    } else if(ref_nt != "-" & mut_nt != "-") { # deletion and insertion in same site
        if(ref_nt == substr(wt_seq, mut_start, mut_end)) {
            mut_seq <- paste(c(substr(wt_seq, 1, mut_start - 1),
                               mut_nt,
                               substr(wt_seq, mut_end+1, nchar(wt_seq))),
                             collapse = "")
            mut_exjun <- get_mut_exjun(wt_exjun = wt_exjun,
                                       mut_nt = mut_nt,
                                       mut_start = mut_start,
                                       mut_end = mut_end)
            if (mut_exjun[[2]] != "") {
                note <- paste(note, mut_exjun[[2]], sep = " ")
            }
            mut_exjun <- mut_exjun[[1]]
        } else {
            note <- "reference unmatched"
            warning("ref_nt doesn't match reference genome")
        }
    }
    return(list(seq = mut_seq, exjun = mut_exjun, note = note))
}

##--------------------------------------------------------------------------------------------
#' @title Get coding sequence of a transcript
#'
#' @description Get positions and junction for the coding sequence of a transcript
#'
#' @details Get the absolute positions on chromosomes for each nucleotide on the coding sequence for
#' a given transcript; and indicate the exon-exon junctions in the coding sequence
#'
#' @param tx_pos The positional information of transcipt from TxDB
#' @param ncbi_build NCBI_build
#'
#' @export
#'
#' @return A list containing two entries:
#' \itemize{
#' \item{\code{wt_pos}} {numeric string; the absolute positions on chromosomes for each nucleotide on the coding sequence}
#' \item{\code{exonJunction}} {numeric string with 0 and 1; 1 represents the position of junctions}
#' }
get_tx_seq <- function(tx_pos, ncbi_build)
{
    pos <- numeric() # position on chromosome
    exon_jun <- numeric() # junction positions on chromosome
    if (ncbi_build == 37) {
        genome <- genome37
    } else if (ncbi_build == 36) {
        genome <- genome36
    } else {
        warning("Unrecognized reference/NCBI_build")
    }
    tx_pos <- tx_pos[!is.na(tx_pos$CDSSTART),] # some genes have no CDS
    if (nrow(tx_pos) > 0) {
        for(j in 1:nrow(tx_pos)) {
            # seq <- paste(seq,as.character(getSeq(genome, tx_pos$EXONCHROM[j], tx_pos$CDSSTART[j], tx_pos$CDSEND[j])),sep = "")
            start <- tx_pos$CDSSTART[j] #
            end <- tx_pos$CDSEND[j]
            start <- as.numeric(start)
            end <- as.numeric(end)
            pos <- c(pos,start:end)
        }

        #seq <- as.character(seq)
        # +/- strand on which the transcript is
        strand <- unique(tx_pos$EXONSTRAND)
        if (length(strand) > 1) {
            warning(paste("Strand of", tx_pos$TXID, "is problematic"))
        }
        if (strand == "+") {
            pos <- sort(pos, decreasing = F)
            exon_jun <- rep(0, length(pos))
            exon_jun[(pos[-1] - pos[-length(pos)]) > 1] <- 1
            exon_jun[which(exon_jun == 1) + 1] <- 1
        } else if(strand == "-") {
            pos <- sort(pos,decreasing = T)
            exon_jun <- rep(0, length(pos))
            exon_jun[(pos[-length(pos)] - pos[-1]) > 1] <- 1
            exon_jun[which( exon_jun == 1) + 1] <-1
        }
        return(list(wt_pos = pos, exonJunction = exon_jun))
    }
    # no cds will return NULL
}

##--------------------------------------------------------------------------------------------
#' @title Find start position
#'
#' @description Find all start codons in the sequence
#'
#' @details Find start codon in the sequence
#'
#' @param start_codons a character vector, the stop codons used; default is TGA, TAA and TAG.
#' @param seq a character vector, sequence.
#' @importFrom Biostrings matchPattern
#' @return potential start codons
findStartPos <- function(start_codons = c("ATA", "ATG"), seq)
{
    start_pos <- c()
    for (codon in start_codons){
        matches <- matchPattern(codon, seq)
        start_pos <- c(start_pos, start(matches))
    }
    start_pos <- sort(start_pos)
    return(start_pos)
}

##--------------------------------------------------------------------------------------------
#' @title Find stop position
#'
#' @description Find all stop codons in the sequence
#'
#' @details Find all stop codons in the sequence
#' @param stop_codons a character vector, the stop codons used; default is TGA, TAA and TAG.
#' @param seq a character vector, sequence.
#' @importFrom Biostrings matchPattern
#' @return {stop_pos} {found all stop positions regardless of the frames}
findStopPos <- function(stop_codons = c("TGA", "TAA", "TAG"), seq)
{
    stop_pos <- c()
    for (codon in stop_codons){
        matches <- matchPattern(codon, seq)
        stop_pos <- c(stop_pos, start(matches))
    }
    stop_pos <- sort(stop_pos)
    return(stop_pos)
}

##--------------------------------------------------------------------------------------------
#' @title Find ORF
#'
#' @description Find open reading frame (ORF)
#'
#' @details Find the ORF in '+' strand
#' @param start_pos all possible start positions
#' @param k minimal length for ORF
#' @param stop_pos all possible stop positions
#' @importFrom Biostrings matchPattern
#' @return the start positions, stop positions and lengths of all posible ORFs
findORF <- function(start_pos, k, stop_pos)
{  # k : Minimum size of Open Reading Frame
    stop_pointers<- c(0,0,0) # frame 1,2,3
    count <- 0
    ORF <- list()
    for (current_start in start_pos) {
        frame <- (current_start%%3) + 1 # decide which frame by position of start codon
        stop_pointer <- stop_pointers[frame] # 0 at first
        if (stop_pointer <= length(stop_pos) && # there are other stop positions
            (stop_pointer == 0 || stop_pos[stop_pointer] < current_start))
             # current stop position is 0 or < start position
        {
            stop_pointer <- stop_pointer + 1 # pick up a larger stop postition
            while ( (stop_pointer <= length(stop_pos))
                    && ((stop_pos[stop_pointer] <= current_start)
                        || (((stop_pos[stop_pointer]%%3) + 1) != frame)) ) # not on the same frame with start position
            {
                stop_pointer <- stop_pointer + 1 # try next stop position
            }
            stop_pointers[frame] <- stop_pointer # get the proper stop position

            if (stop_pointer <= length(stop_pos))
            {
                if ((stop_pos[stop_pointer] + 2 - current_start + 1) > k )
                {
                    count <- count + 1
                    new_ORF <- c(n = count,
                                 Frame = frame,
                                 Start = current_start,
                                 Stop = stop_pos[stop_pointer],
                                 Length = stop_pos[stop_pointer] + 2 - current_start + 1)
                    ORF[[length(ORF) + 1]] <- new_ORF
                }
            }
        }
    }
    ORF <- t(as.data.frame(ORF, as.is=T))
    rownames(ORF) <- NULL
    return(ORF)
}
