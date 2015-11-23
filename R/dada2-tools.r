# Utility Funcitons for Working with Dada2->Phyloseq Workflows
#
#
#
#

#' Write the Sequences of a seqtable to disk
#' 
#' @param seqtable is the ouput from dada2's \code{makeSequenceTable} function
#' @param outfile where to write the fna file
#' 
#' @importFrom stringr str_pad
write_FNA <- function(seqtable, outfile){
  FIRST = TRUE
  seqs <- colnames(seqtable)
  for (i in seq_along(seqs)) {
    if (i == 2)  FIRST = FALSE
    
    outdata = paste0(">Seq_", str_pad(i, 7, pad=0),"\n", seqs[[i]], "\n")
    if (FIRST == TRUE) {
      write(outdata,file = outfile, sep="", append = FALSE)
    } else {
      write(outdata,file = outfile, sep="", append = TRUE)
    } 
    
  }
  print(paste0("All records are printed to fastafile ", outfile))
}

#' Combine RDS files and write the phyloseq object (RDSfile), the OTUtable (.txt), 
#' and the sequences (.fna) to disk
#'
#'
#' @param dadalist list of \code{dada2} objects where names will be samplenames
#' @param otutableout filename for otutable
#' @param phyloseqoutfile RDS filename to store phyloseqobject
#' @param FNAoutfile outpupt file for the sequences.
#' 
#' @importFrom dada2 makeSequenceTable
#' @importFrom stringr str_pad
#' @export
combineRDSfiles <- function(dadalist, otutableout, phyloseqoutfile, FNAoutfile) {
  # load everything into a phyloseq object
  seqtab <- makeSequenceTable(dadalist)
  #seqs <- colnames(seqtab)
  otab <- otu_table(seqtab, taxa_are_rows=FALSE)
  colnames(otab) <- paste0("Seq_", str_pad(seq(ncol(otab)), 7, pad=0))
  
  taxtab <- tax_table(matrix(colnames(otab), ncol=1))
  rownames(taxtab) <- colnames(otab)
  colnames(taxtab) <- "Sequence"
  ps <- phyloseq(otab, taxtab)
  
  # save the phyloseq object
  saveRDS(ps, file=phyloseqoutfile)
  print(paste0("Saved phyloseq object as ", phyloseqoutfile ))
  
  # save OTUfile with taxa as first column
  if (taxa_are_rows(ps)) {
    write.table(otu_table(ps), file = otutableout, quote=FALSE, sep="\t")
  } else {
    write.table(t(otu_table(ps)), file = otutableout, quote=FALSE, sep="\t")
  }
  
  print(paste0("Saved otutable as ", otutableout ))
  
  # save the FNA file of sequences
  write_FNA(seqtab, outfile = FNAoutfile)
}
