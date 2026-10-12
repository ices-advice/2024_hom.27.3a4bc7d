## Extract results of interest, write TAF output tables

## Before:
## After:

library(icesTAF)
library(FLSAM)
library(ggplot2)

mkdir("output")

#- Load the model output
load(file="model/model.RData")
source("output_diagnostics.R")

#- Generate standard set of diagnostics
run 	<- "run_with_2023ibts"

FLSAM.diagnostics("output/",run=run,NSHOM.sam,NSHOM.tun,residuals=T,retro=NSHOM.retro,loo=NSHOM.loo[-2])

load(file="model/model2.RData")

#- Generate standard set of diagnostics
run 	<- "run_wh_2023ibts"

FLSAM.diagnostics("output/",run=run,NSHOM.sam,NSHOM.tun,residuals=T,retro=NSHOM.retro,loo=NSHOM.loo[-2])

# - table of 