library("phyloseq.tools")
library("testthat")

df <- data.frame( otu1=c('a','b','c','d','e','f','g','h','i','j','k'), 
                  otu2=c('a','a','c','c','e','e','g','e','i','j','i'),
                  otu3=c('a','a','a','a','e','e','e','e','i','i','i'),
                  stringsAsFactors = FALSE)

test_that("make_tree can create a phylogenetic tree from a dataframe", {
  expect_that(make_tree(df), is_a("phylo"))
})