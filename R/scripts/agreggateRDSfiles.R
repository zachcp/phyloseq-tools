#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(docopt))

"Usage: 
agreggateRDSfiles.R [options]

Description:   Agregate RDS files of DADA output objects
Options:
--datapath=<datapath>              Forward FastQ of paired end
--experimentname=<experimentname>  Name of the Experiment
--outdir=<outdir>                  [default: '.']  Output Directory
" -> doc
opts <- docopt(doc)

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(dada2))
suppressPackageStartupMessages(library(phyloseq))
suppressPackageStartupMessages(library(stringr))

#' Write FNA from Seqfiles
#'
write_FNA <- function(seqtable, outfile){
  FIRST = TRUE
  seqs <- colnames(seqtable)
  for (i in seq_along(seqs)) {
    if (i == 2)  FIRST = FALSE
    
    outdata = paste0(">Seq_", stringr::str_pad(i, 7, pad=0),"\n", seqs[[i]], "\n")
    if (FIRST == TRUE) {
      write(outdata,file = outfile, sep="", append = FALSE)
    } else {
      write(outdata,file = outfile, sep="", append = TRUE)
    } 
    
  }
  print(paste0("All records are printed to fastafile ", outfile))
}

#' Combine RDS files
#'
combineRDSfiles <- function(rdslist, otutableout, phyloseqoutfile, FNAoutfile) {
  # load everything into a phylsoeq object
  seqtab <- makeSequenceTable(rdslist)
  #seqs <- colnames(seqtab)
  otab <- otu_table(seqtab, taxa_are_rows=FALSE)
  colnames(otab) <- paste0("Seq", seq(ncol(otab)))
  taxtab <- tax_table(matrix(colnames(otab), ncol=1))
  rownames(taxtab) <- colnames(otab)
  colnames(taxtab) <- "Sequence"
  ps <- phyloseq(otab, taxtab)
  
  # save the phyloseq object
  saveRDS(ps, file=phyloseqoutfile)
  print(paste0("Saved phyloseq object as ", phyloseqoutfile ))
  
  # save OTUfile
  write.table(otab, file = otutableout)
  print(paste0("Saved otutable as ", otutableout ))
  
  # save the FNA file of sequences
  write_FNA(seqtab, outfile = FNAoutfile)
}


# check opts
datapath       <- opts$datapath
experimentname <- opts$experimentname
outdir         <- opts$outdir

#datapath = "/data/DFD/data/"
RDSfiles <- list.files(datapath)
RDSnames <- as.character(lapply(RDSfiles, function(x) return(strsplit(x, "_")[[1]][1])))
RDSdata <- lapply(RDSfiles, function(x) readRDS(paste0(datapath,x)))
names(RDSdata) <- RDSnames

#Combine Data adn Write out, the OTUfile, the Phyloseq Object, and the FNA file
combineRDSfiles(rdslist= RDSdata, 
                otutableout =     paste0(outdir,"/", experimentname, "_otus.txt"),
                phyloseqoutfile = paste0(outdir,"/", experimentname, "_phyloseq.RDS"),
                FNAoutfile =      paste0(outdir,"/", experimentname, "_otus.fna"))




