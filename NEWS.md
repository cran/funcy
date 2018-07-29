<!-- NEWS.md is generated from NEWS.Rmd. Please edit that file -->


funcy version 1.0.0
===================

-------

* R/format.R:
 
For function regFuncy, method "pace" is not available any longer. 
"Pace" was dependent (via an interface funcyOctave) on the R package "RcppOctave" which had been removed from 
CRAN.

* R/functions.R:
Removed import of trapz function from caTools since the package it is scheduled for archival.

* Unit tests added with testthat

* added CITATION of JSS paper

* added citations of JSS paper in Rd files
