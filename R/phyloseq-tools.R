#'
#' Import Blast File as Taxonomy Table
#' 
#' Sometimes an adhoc taxonomic table is called for and Blast can do the 
#' trick. This function allows you to import your blast results as a taxonomic
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
      stop("the evalue cutoff must be numeric")
    }
    blast <- blast %>% filter(evalue <= cutoff)
  }
  
  #check the length of the blast table
  blastdf_dimensions = dim(as.data.frame(blast))
  if (blastdf_dimensions[1] == 0) stop("Empty Blast Dataframe. check the cutoff values")
  
  
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

#' Merge Taxa Based on UC file Clustering
#'
#' @importFrom dplyr %>%
mergeTaxaUC <- function(phy, ucfile){
  phy2 <- phy
  ucdata <- load_uc(ucfile)
  ucdata <- ucdata %>% filter(rectype == "H") %>% select(query, target)
  
  splits <- split(ucdata, ucdata$target)
  
  for (df in splits){
    target = unique(df$target)
    queries = df$query
    print("Merge Ahoy!")
    print(c(target,queries))
    phy2 <- merge_taxa(phy2, eqtaxa = c(target,queries), archetype = 1)
  }
  return(phy2)
}

#' Calculate Rarefaction Curves
#'
#' Generate Rarefaction Curves from Phyloseq.
#' Adopted from @@and3k https://github.com/joey711/phyloseq/issues/143
#' 
#' @param physeq phyloseq object
#' @param measures vector of rarefaction measures
#' @param depths vector of depths to rarefy to
#' @param parallel wether to use ldplyr's parallel features
#' @param ncpus nubmer of cpus
#' 
#' @importFrom plyr ldply
#' @importFrom plyr summarise
#' @importFrom reshape2 melt
#' @import foreach
#' @import doParallel
#' 
#' @export
calculate_rarefaction_curves <- function(physeq, measures, depths, parallel=FALSE, ncpus=1) {
  estimate_rarified_richness <- function(physeq, measures, depth) {
    if(max(sample_sums(physeq)) < depth) return()
    physeq <- prune_samples(sample_sums(physeq) >= depth, physeq)
    
    rarified_physeq <- rarefy_even_depth(physeq, depth, verbose = FALSE)
    
    alpha_diversity <- estimate_richness(rarified_physeq, measures = measures)
    
    # as.matrix forces the use of melt.array, which includes the Sample names (rownames)
    molten_alpha_diversity <- melt(as.matrix(alpha_diversity), 
                                   varnames = c('Sample', 'Measure'), 
                                   value.name = 'Alpha_diversity')
    
    molten_alpha_diversity
  }
  
  if (parallel){
    #if parallel setup the cluster
    library(doParallel)
    print("Running Calculation in Parallel...")
    cl <- makeCluster(ncpus)
    registerDoParallel(cl)
  }
  
  names(depths) <- depths # this enables automatic addition of the Depth to the output by ldply
  rarefaction_curve_data <- plyr::ldply(depths, estimate_rarified_richness, 
                                        physeq = physeq, measures = measures, 
                                        .id = 'Depth', 
                                        .paropts = list(.packages = c('phyloseq', 'vegan','reshape2')),
                                        .progress = ifelse(interactive(), 'text', 'none'),
                                        .parallel = parallel)
  
  # convert Depth from factor to numeric
  rarefaction_curve_data$Depth <- as.numeric(levels(rarefaction_curve_data$Depth))[rarefaction_curve_data$Depth]
  
  #add a standard deviation column
  rarefaction_curve_data_summary <- plyr::ddply(rarefaction_curve_data, c('Depth', 'Sample', 'Measure'), 
                                                plyr::summarise, 
                                                Alpha_diversity_mean = mean(Alpha_diversity), 
                                                Alpha_diversity_sd = sd(Alpha_diversity))
  
  #add the sample data
  rarefaction_curve_data_summary_verbose <- merge(rarefaction_curve_data_summary, 
                                                  data.frame(sample_data(physeq)), 
                                                  by.x = 'Sample', 
                                                  by.y = 'row.names')
  return(rarefaction_curve_data_summary_verbose)
}
