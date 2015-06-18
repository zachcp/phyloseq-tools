#'
#' Import Blast File as Taxonomy Table
#' 
#' Sometimes an adhoc taxonomic table is called for and Blast can do the 
#' trick. This funciton allows you to import your blast results as a taxonomic
#' table. In this case the 'query' and 'target' corresponds with would be your
#' OTU name and taxon. 

#' @param physeq
#' @param blasttablefile
#' @param cutoff 
#' 
#' @import phyloseq
#' @importFrom dplyr select
#' @importFrom dplyr %>%
#' @importFrom dplyr left_join
#' @importFrom magrittr %<>%
#' @export
tax_from_blast <- function(physeq, blasttablefile, cutoff=NULL){
  #some basic data checks
  if (!c('otu_table') %in% getslots.phyloseq(physeq)){
    stop("Phyloseq object must have an otu_table to import a blast taxonomy")
  }
  #load taxonomy file
  blast <- load_blast(blasttablefile) %>% 
    select(query, target, percent_ident, length, evalue, bitscore) %>%
    group_by(query) %>%  
    top_n(1, wt=bitscore) %>% #keep only the top hits
    group_by(query) %>%
    slice(1:1)                #when there are identical top hits take the first
  
  #cutoff results not passing evalue threshold
  if (!is.null(cutoff)) {
    if (!is.numeric(cutoff)){
      stop("teh evalue cutoff must be numeric")
    }
    blast %<>% filter(evalue <= cutoff)
  }
  
  #check the length of the blast table
  blastdf_dimensions = dim(as.data.frame(blast))
  if (blastdf_dimensions[0] == 0) stop("Empty Blast Dataframe. check the cutoff values")
  
  
  #make a dataframe of otunames
  otus  <- taxa_names(physeq) %>% as.data.frame()
  names(otus) <- c('query')  
  
  #merge the otudata and the blastdata together 
  #set NAs on the column appropriate to the column
  blasttax <- left_join(otus,blast)
  blasttax$target[is.na(blasttax$target)] <- "No Target"
  blasttax$percent_ident[is.na(blasttax$percent_ident)] <- 0
  blasttax$length[is.na(blasttax$length)] <- 0
  blasttax$evalue[is.na(blasttax$evalue)] <- 1
  blasttax$bitscore[is.na(blasttax$bitscore)] <- 0
  rownames(blasttax) <- blasttax$query
  
  blasttax %<>% select(-query) %>% as.matrix()
  tax_table(physeq) <- blasttax
  return(physeq)
}
#'
#' Get the  unique values from a column and return a simple phylogenetic tree 
#' 
#' This function is intened for use with a dataframe of sluter information where
#' each row is a sequence/OTU and each column is a level of similariy (percent identity)
#' this function will return a tree based on the lowest percent identity column and
#' we will build off of it.
#' @importFrom stringr str_c
#' @importFrom ape read.tree
column_to_tree <- function(df, colname){
  #get unique names, make newick string and tree
  concatnames <- unique(df[[colname]]) %>% str_c(collapse=",")
  tree = paste("(", concatnames,");")
  read.tree(text=tree)
}
#'
#' add a tree to a base tree at tipname
#' 
#' @importFrom ape bind.tree
#' @return phylo object
#' @keywords internal
add_tree_at_tip <-function(basetree, incomingtree, tipname){
  bind.tree(basetree, incomingtree, where=which(basetree$tip.label==tipname))
}
#'
#' get values from adjacent row
#' 
#' in a table representing phylogenetic relationships
#' each column of the table correspods to one more 
#' level of clustering. this will get the values of
#' the "higher" column, based on shared rows with the 
#' "lowercolumn"
#' 
#' @return phylo object
#' @importFrom ape read.tree
#' @keywords internal
get_next_order_tree <- function(df, lowidcol, highidcol, val){
  newvals <- df[df[[lowidcol]] == val,][[highidcol]]
  vals    <- str_c(unique(newvals),collapse=",")
  tree    <- paste("(", vals,");")
  return(read.tree(text=tree))
}
#'
#' make a tree from a dataframe
#' @export
make_tree <- function(df){
  numcols <- length(names(df))
  # use last column for base tree and loop backwards
  basetree <- column_to_tree(df, numcols) 
  for (i in c(numcols:1)){
    #loop through the columns and start with the second column
    if (i < numcols) {
      for (val in unique(df[[i+1]])){
        #get the values for the tree...
        incomingtree <- get_next_order_tree(df=df, 
                            lowidcol= i+1,
                            highidcol= i,
                            val=val)
        #...and tack them onto the original tree.
        basetree <- add_tree_at_tip(basetree,incomingtree,val)
      }
    }
  }
  return(basetree)  
}