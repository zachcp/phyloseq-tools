library("phyloseq.tools")
library("testthat")

#import basic data
usearch  <- system.file("extdata", "usearch.uc", package="phyloseq-tools")
blast    <- system.file("extdata", "blast.bl6", package="phyloseq-tools")

test_that("load_blast generates a dataframe", {
  expect_that(load_blast(blast), is_a("data.frame"))
})

test_that("load_blast generates a dataframe", {
  expect_that(load_uc(usearch), is_a("data.frame"))
})