<!-- NEWS.md is generated from NEWS.Rmd. Please edit that file -->


funcy version 0.8.5
===================

Changes
-------

* R/plotFSCM.R:

1. `plotOverview`: regTime is now reflected on x-axis.

2. `plotOverview`: showLegend can now be set to FALSE

3. plot type `dist2centers`: removed legend, added center
curves


* R/plot.R: changed cex of function `plotLoc`


* R/xecute.R:

1. added rownames of props (`rownames(props) <- object@methodName`)

2. changed stop("Please select one method or methods=ALL.") to
   stop("Please select one method or methods='ALL'.")


