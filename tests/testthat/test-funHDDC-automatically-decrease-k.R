context("test-automatically-decrease-number-of-classes-for-method-funHDDC")

set.seed(2804)
ds <- sampleFuncy(obsNr = 50, k = 4, timeNr = 8, reg = TRUE)

# prevents funcitons from spamming the the console with unrequestted output
silent <- function(expr) 
{
  null <- capture.output(result <- expr)
  return(result)
}

test_that("k is reduced on warning 'All models converged.'", {
  
  null <- capture.output(
    x <- expect_warning(
      silent(funcit(data = Data(ds), clusters = Cluster(ds),
                    methods = "funHDDC", seed = 2404, k = 3)), 
      regexp = "Clustering with 3 classes is not possible. 2 clusters are used."))
  
  # k is reduced until k=1
  expect_equal(object = length(unique(Cluster(x))),
               expected = 1)
  
})
