library("phyloseq.tools")
library("testthat")
library("phyloseq")

df <- data.frame( otu1=c('a','b','c','d','e','f','g','h','i','j','k'), 
                  otu2=c('a','a','c','c','e','e','g','e','i','j','i'),
                  otu3=c('a','a','a','a','e','e','e','e','i','i','i'),
                  stringsAsFactors = FALSE)

test_that("make_tree can create a phylogenetic tree from a dataframe", {
  expect_that(make_tree(df), is_a("phylo"))
})


test_that("curve rarefaction works using single or multiple threads", {
  data('GlobalPatterns')
  psdata <- GlobalPatterns
  expect_that(calculate_rarefaction_curves(psdata,c('Observed', 'Shannon'), c(100,100)), is_a("data.frame"))
  expect_that(calculate_rarefaction_curves(psdata,c('Observed', 'Shannon'), c(100,100), parallel = TRUE), is_a("data.frame"))
  expect_that(calculate_rarefaction_curves(psdata,c('Observed', 'Shannon'), c(100,100), parallel = TRUE, ncpus=2), is_a("data.frame"))
})


