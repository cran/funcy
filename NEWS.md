<!-- NEWS.md is generated from NEWS.Rmd. Please edit that file -->


funcy version 0.8.6
===================

Changes (NOCH NICHT SUBMITTED)
-------

* R/plot.R:

1. added `xlab=NULL, ylab=NULL` to plotFuncy: xlab and ylab can now be
   set. They are void if nothing was set. 

2. plotOverview: changed mar=c(3,2,3,10) to mar=c(3,2,3,1)

*R/xecute.R
method "summary" for "funcyOutList": changed cat("\\n") to cat("\n")

3. functions.R
added polynomial basis (create.monomial.basis)

4. funct.Rd
for method "iterSubspace" parameter name simplify was changed to simplif (as implemented)
