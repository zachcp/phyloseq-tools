#'
#' Import Blast File as Taxonomy Table
#' 
#' Sometimes an adhoc taxonomic table is called for and Blast can do the 
#' trick. This funciton allows you to import your blast results as a taxonomic
#' table. In this case the 'query' and 'target' corresponds with would be your
#' OTU name and taxon. 
#' 
#' @import phyloseq
#' @importFrom dplyr select
#' @importFrom dplyr %>%
#' @importFrom dplyr left_join
#' 
tax_from_blast <- function(physeq, blasttablefile){
  #some basic data checks
  if (!c('otu_table') %in% getslots.phyloseq(physeq)){
    stop("Phyloseq object must have an otu_table to import a blast taxonomy")
  }
  
  #import biotables
  if(!require("biotables")){
    devtools::install_github('zachcp/biotables')
    library("biotables")
  }
  
  #load taxonomy file
  blast <- load_blast(blasttablefile)
  blast <- blast %>% select(query, target, percent_ident, length, evalue, bitscore)
  
  #make a dataframe of otunames
  otus  <- taxa_names(physeq) %>% as.data.frame()
  names(otus) <- c('query')  
  
  #merge the otudata and the blastdata together 
  blasttax <- left_join(otus,blast) %>%
    as.matrix()
  rownames(blasttax) <- blasttax$query
    
  return(blasttax)
}



library(dplyr)
library(magrittr)
library(stringr)
library(lazyeval)
library(ape)

column_to_tree <- function(df, colname){
  #get unique names
  concatnames <- unique(df[[colname]]) %>% str_c(collapse=",")
  
  #make a simple newick
  tree = paste("(", concatnames,");")
  
  #creatae a tree
  read.tree(text=tree)
}


add_to_tree <- function(df, tree, lowidcol, highidcol, val){
  # get values from the lowidcol where the highidcol is val
  # add those to the tree at the val position
  newvals <- df[df[[lowidcol]] == val,][[highidcol]]
  #concat the vals and make a tree
  concatnames <- str_c(unique(newvals),collapse=",")
  newtree = read.tree( text= paste("(", concatnames,");"))
  tree = bind.tree(tree,newtree, where=which(tree$tip.label==val))
}

df <- data.frame( otu1=c('a','b','c','d','e','f','g','h','i','j','k'), 
                  otu2=c('a','a','c','c','e','e','g','e','i','g','i'),
                  otu3=c('a','a','a','a','e','e','e','e','i','i','i'),
                  stringsAsFactors = FALSE)


#base tree of otu3
basetree <- column_to_tree(df, "otu3")

#otu3 to otu2
t <- add_to_tree(df,basetree,"otu3","otu2","a")
t <- add_to_tree(df,t,"otu3","otu2","e")
t <- add_to_tree(df,t,"otu3","otu2","i")

#otu2 to otu1
t <- add_to_tree(df,t,"otu2","otu1","a")
t <- add_to_tree(df,t,"otu2","otu1","c")
t <- add_to_tree(df,t,"otu2","otu1","e")
t <- add_to_tree(df,t,"otu2","otu1","g")
t <- add_to_tree(df,t,"otu2","otu1","i")
plot(t)


