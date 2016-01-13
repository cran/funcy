<!-- NEWS.md is generated from NEWS.Rmd. Please edit that file -->


funcy version 0.8.4
===================

Changes
-------

* R/plotFSCM.R:

1. `plotOverview`: regTime is now reflected on x-axis.

2. `plotOverview`: showLegend can now be set to FALSE

3. plot type `dist2centers`: removed legend, added center
curves


* R/plot.R: changed cex of function `plotLoc`


* R/xecute.R: added rownames of props (`rownames(props) <- object@methodName`)


* la/clapack.h, la/usepack.cc:

1. removed `f77_ret_void F77_FUNC(dgegv)`

2. removed `f77_ret_void F77_FUNC(dggsvd)`

