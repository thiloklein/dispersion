rm(list=ls())

library(dispersion)
library(RcppOctave)

##############################
## share vector of length 3 ##
##############################

## dispunif3 (in R)
dispunif3(s=c(1,3,7), lb=0.1, ub=0.9, bss=2000)

## robustmeasure function (in Matlab)
mfile <- "/home/thilo/Documents/Packages/dispersion/inst/scripts/dispd.m"
o_source(mfile)
.CallOctave("robustmeasure", s=c(1,3,7), shape=1, l1=100, n=3, bss=2000, lb=10, ub=90) 

#####################################
## share vector of length larger 3 ##
#####################################

## dispunif3 (in R)
dispunif3(s=c(1,2,3,4,5,6,7,8,9,10), lb=0.1, ub=0.9, bss=2000)

## robustmeasure function (in Matlab)
set.seed(123)
mfile <- "/home/thilo/Documents/Packages/dispersion/inst/scripts/dispd.m"
o_source(mfile)
.CallOctave("robustmeasure", s=c(1,2,3,4,5,6,7,8,9,10), shape=1, l1=100, n=3, bss=2000, lb=10, ub=90) 

## robustmeasure function (in R)
set.seed(123)
robustmeasure(s=c(1,2,3,4,5,6,7,8,9,10), shape=1, l1=100, n=3, bss=2000, lb=10, ub=90)

###################################
## share matrix using libor data ##
###################################

data(libor)
robustmeasure(s=libor, shape=1, l1=100, n=3, bss=2000, lb=10, ub=90)













