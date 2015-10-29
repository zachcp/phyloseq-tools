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
suppressPackageStartupMessages(library(phyloseq.tools))

# check opts
datapath       <- opts$datapath
experimentname <- opts$experimentname
outdir         <- opts$outdir

RDSfiles <- list.files(datapath)
RDSnames <- as.character(lapply(RDSfiles, function(x) return(strsplit(x, "_")[[1]][1])))
RDSdata <- lapply(RDSfiles, function(x) readRDS(paste0(datapath,x)))
names(RDSdata) <- RDSnames

#Combine Data adn Write out, the OTUfile, the Phyloseq Object, and the FNA file
combineRDSfiles(rdslist= RDSdata, 
                otutableout =     paste0(outdir,"/", experimentname, "_otus.txt"),
                phyloseqoutfile = paste0(outdir,"/", experimentname, "_phyloseq.RDS"),
                FNAoutfile =      paste0(outdir,"/", experimentname, "_otus.fna"))




