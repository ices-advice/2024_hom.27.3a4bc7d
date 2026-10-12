## Run analysis, write model results

## Before:
## After:
 .libPaths(c(file.path(getwd(),"boot/library"),.libPaths()))
 .libPaths()

 library(icesTAF)
 taf.library(FLCore)
 library(FLAssess)
 #library(FLash)
 library(ggplotFL)
 library(FLBRP)
 taf.library(stockassessment)
 taf.library(FLSAM)
 library(FLasher)

 rm(list = ls())

mkdir("forecast")
load(file="model/model.RData", verbose =TRUE)

# data, intermediate, advice and final YEARS
NSHOM <- NSHOM+NSHOM.sam
dy <- dims(NSHOM)$maxyear
iy <- dy + 1
ay <- iy + 1
fy <- ay + 1
ry <- (dy-9):dy #recruitment years

# TAC advice current year
tac <- FLQuant(9730, dimnames=list(age='all', year=iy), units="tonnes")

# sample from ry for recruitment 
recdr <- rec(NSHOM)[,ac(sample(ry,3))]

# GENERATE targets from reference points
refpts=FLPar(c(0, 0.034, 0.042, 236228,170000,236228,NA),
          dimnames=list(params=c("F0", "Fmsy","Flim","Btrigger","Blim","Bpa","Fpa"),year=iy,iter=1),
          units=c(NA,NA,NA,"tonnes","tonnes","tonnes",NA))
targets <- Matrix::expand(as(refpts, 'FLQuant'), year=ay:fy)

# iy assumption 

refpts_iy=FLPar(c(0, 0.034, an(fbar(NSHOM)[, ac(2023)]), 9730),
          dimnames=list(params=c("F0", "Fmsy","Fsq", "TAC"),year=iy:ay,iter=1),
          units=c(NA,NA,NA,"tonnes"))

iy_targets <- Matrix::expand(as(refpts_iy, 'FLQuant'), year=iy:ay)

# --- SETUP future
# 3 years extending, taking average of 5 years wts/selex
fut <- stf(NSHOM, nyears=3, wts.nyears=5, fbar.nyears=1)
fbar(fut) #iy is status quo F


# SET  SRR
# NSsrr <- predictModel(model=rec~a, params=FLPar(c(recdr), units="thousands",
#   dimnames=list(params="a", year=seq(iy, length=3), iter=1)))


#  A deterministic short term forecast does not use an SRR (in the traditional sense). 
# Instead, it generally assumes that recruitment in the future is some mean (e.g. geometric mean) of the historic recruitments. However, we still need to have an SRR that contains this mean value, which is what we mimic for this example. First, we estimate the geometic mean recruitment, that we then add to an SRR object
  mean_rec <- exp(mean(log(rec(trim(NSHOM, year = ry))))) 
  NSHOM_sr <- as.FLSR(NSHOM, model="geomean")
  params(NSHOM_sr)['a',] <- mean_rec
  params(NSHOM_sr)

# --- INTERMEDIATE year assumption
  discards.wt(NSHOM) <- 0
  discards.n(NSHOM) <- 0

  discards.wt(fut) <- 0
  discards.n(fut)  <- 0

  # fut@stock.wt
  # fut@m
  # fut@stock.n

# UPDATE 2024, 2025
##does not matter what to put here, gets overwritten by ftcl
#fut <- fwd(fut, sr=NSsrr, f=iy_targets["Fsq",ac(2024)])
fut <- fwd(fut, sr=NSHOM_sr, f=iy_targets["Fsq",ac(2024)])
catch(fut)

# UPDATE 2025 
#fut <- fwd(fut, sr=NSsrr, f=targets["F0",])
#ssb(fut)
#fbar(fut)
#rec(fut)
catch(fut[,ac(2024)])
# Monitoring TAC targets
Cmon_vals=FLPar(c(1200, 1800,2400),
          dimnames=list(params=c("F_MTAC_10", "F_MTAC_15","F_MTAC_20"),year=iy,iter=1),
          units=c("tonnes","tonnes","tonnes"))

Cmon_targets <- Matrix::expand(as(Cmon_vals, 'FLQuant'), year=ay:fy)

Flast_vals=FLPar(c(an(fbar(fut)[, ac(iy)])),
          dimnames=list(params=c("Flast"),year=iy,iter=1),
          units=c(NA, "tonnes"))

Flast_targets <- Matrix::expand(as(Flast_vals, 'FLQuant'), year=ay:fy)



# - calculation of advice option
  # calc SSB 2026  = SSB 2025 F 
    iters_n <-10000
    flevels <- seq(0, 1, length.out = iters_n)
    range(flevels)

    # BUILD control with targets as iters
    iy_cntrl <- FLQuant(an(iy_targets["Fsq",])[1],
        dimnames=list(age='all', year=c(iy,ay, fy)))

    iy_cntrl <- propagate(iy_cntrl, iters_n)
    iy_cntrl[, ac(2025)] <-flevels
    iy_cntrl[, ac(2026)] <-flevels

    iy_cntrl[,,,,,5]
    iy_cntrl[,,,,,6]

    control <- as(FLQuants(
      fbar=iy_cntrl),  "fwdControl")

    # RUN
    fut_exp <- propagate(fut, iter = iters_n)
    fruns_undiv <- fwd(fut_exp, sr=NSHOM_sr, control=control)

    catch(fruns_undiv[,ac(2024)])
    ssb(fruns_undiv[,ac(2024)])
    fbar(fruns_undiv[,ac(2024)])


    stored_SSB_2025<- as.data.frame(ssb(fruns_undiv[,ac(2025)]))$data
    stored_SSB_2026<- as.data.frame(ssb(fruns_undiv[,ac(2026)]))$data

    #View(data.frame(flevels = flevels,diff = stored_SSB_2025 - stored_SSB_2026, SSB_2025 = stored_SSB_2025, SSB_2026 = stored_SSB_2026))
    AD_f_ssb25_isssb26<- flevels[which(abs(stored_SSB_2025 - stored_SSB_2026) == min(abs(stored_SSB_2025 - stored_SSB_2026)))]

  # Calc advice F
    #an(iy_targets["Fsq",])[1]
  
    flevels <- seq(0, an(targets["Fmsy",])[1], length.out = 100)
  
    # BUILD control with targets as iters
    iy_cntrl <- FLQuant(an(iy_targets["Fsq",])[1],
        dimnames=list(age='all', year=c(iy,ay)))
    #propagate(iy_targets["Fsq",], 100)
    iy_cntrl <- propagate(iy_cntrl, 100)
    iy_cntrl[, ac(2025)] <-flevels
  
    control <- as(FLQuants(
      fbar=iy_cntrl),  "fwdControl")
  
    # RUN
    fut_exp <- propagate(fut, iter = 100)
    fruns_undiv <- fwd(fut_exp, sr=NSsrr, control=control)
    stored_SSB<- as.data.frame(ssb(fruns_undiv[,ac(2025)]))$data
  
    Fad<-stored_SSB/an(targets["Btrigger",][1])*an(targets["Fmsy",])[1]
    AD_f<- Fad[which(abs(flevels- Fad) == min(abs(flevels- Fad)))]
  
# --- PROJECT catch options

# DEFINE catch options
catch_options <- list(

 # advice
  advice=FLQuants(fbar=FLQuant(AD_f, dimnames=list(age='all', year=c(ay, ay+1)))),

 # ssb25 = ssb26
  rolssb=FLQuants(fbar=FLQuant(AD_f_ssb25_isssb26, dimnames=list(age='all', year=c(ay, ay+1)))),

 # FMSY
  Fmsy=FLQuants(fbar=targets["Fmsy",]),

 # F0
  F0=FLQuants(fbar=FLQuant(0, dimnames=list(age='all', year=c(ay, ay+1)))),

 # Flim
  Flim=FLQuants(fbar=targets["Flim",]),
 
# #  # Bpa
# ##  Bpa=FLQuants(ssb_flash=targets["Bpa",]),
# #
# #  # Blim
# ##  Blim=FLQuants(ssb_flash=targets["Blim",]),
# #
# #  # MSYBtrigger
# ##  MSYBtrigger=FLQuants(ssb_flash=targets["Btrigger",]),
# #

  # F iy
   lastF= FLQuants(fbar=
    Matrix::expand(Flast_targets, year=c(ay, ay+1))),

 # Monitoring TAC
  MTAC_10=FLQuants(catch=
    Matrix::expand(Cmon_targets["F_MTAC_10",], year=c(ay, ay+1))),

  MTAC_15=FLQuants(catch=
    Matrix::expand(Cmon_targets["F_MTAC_15",], year=c(ay, ay+1))),
 
  MTAC_20=FLQuants(catch=
    Matrix::expand(Cmon_targets["F_MTAC_20",], year=c(ay, ay+1))),
 
 #rollover TAC
  rotac=FLQuants(catch=Matrix::expand(tac, year=ay)),

 #SSBfy = SSBay
  #SSB25=FLQuants(SSB=Matrix::expand(tac, year=ay)),
 
 # catch rollover (catch 2025 == catch 2024)

  rocatch=FLQuants(catch=Matrix::expand(catch(fut[,ac(2024)]), year=ay)),

  rocatch_min20=FLQuants(catch=Matrix::expand(catch(fut[,ac(2024)])*0.80, year=ay)),
  rocatch_plus25=FLQuants(catch=Matrix::expand(catch(fut[,ac(2024)])*1.25, year=ay))


)

#- Setup control for forecast
fctls <- lapply(catch_options, function(x)
  as(x, "fwdControl")
)

fctls$rotac <- fwdControl(data.frame(year = c(ay, ay+1),quant = c("catch", "f"), value =  c(an(tac), 1), relYear = c(NA,ay)))

fctls$rotac_min20 <- fwdControl(data.frame(year = c(ay, ay+1),quant = c("catch", "f"), value =  c(an(tac)*0.80, 1), relYear = c(NA,ay)))
fctls$rotac_plus25 <- fwdControl(data.frame(year = c(ay, ay+1),quant = c("catch", "f"), value =  c(an(tac)*1.25, 1), relYear = c(NA,ay)))

fctls$rocatch <- fwdControl(data.frame(year = c(ay, ay+1),quant = c("catch", "f"), value =  c(an(catch(fut[,ac(2024)])), 1), relYear = c(NA,ay)))
fctls$rocatch_min20 <- fwdControl(data.frame(year = c(ay, ay+1),quant = c("catch", "f"), value =  c(an(catch(fut[,ac(2024)])*0.80), 1), relYear = c(NA,ay)))
fctls$rocatch_plus25 <- fwdControl(data.frame(year = c(ay, ay+1),quant = c("catch", "f"), value =  c(an(catch(fut[,ac(2024)])*1.25), 1), relYear = c(NA,ay)))

# RUN!
#runs <- FLStocks(lapply(fctls, function(x) fwd(fut, sr=NSsrr, control=x)))
runs <- FLStocks(lapply(fctls, function(x) fwd(fut, sr=NSHOM_sr, control=x)))

names(runs)
plot(runs)
#stock.n(window(runs[["rotac"]], start = 2015, end = 2026))

ff<- FLStocks(  F0      = window(runs[["F0"]], start = 2015, end = 2026),
                MTAC_10 = window(runs[["MTAC_10"]], start = 2015, end = 2026),
                MTAC_15 = window(runs[["MTAC_15"]], start = 2015, end = 2026),
                MTAC_20 = window(runs[["MTAC_20"]], start = 2015, end = 2026),
                rotac   = window(runs[["rotac"]], start = 2015, end = 2026)
                )
plot(ff)

# SAVE
save(runs, fruns, fctls, tac, file="forecast/runs.rda", compress="xz")
  load(file="forecast/runs.rda", verbose =TRUE)


#Options summary table
opt.sum.tbl <- function(stcks) {

    options.sum.tbl <- sapply(as.list(1:length(stcks)),function(i) {
    opt <- names(stcks)[i]
    stk <- stcks[[opt]]
    #Build up the summary
    sum.tbl <- data.frame(Rationale=ac(names(stcks)[i]),
    F.ImY=fbar(stk)[,as.character(ImY),drop=TRUE],
    Catch.ImY=computeCatch(stk)[,as.character(ImY),drop=TRUE],
    SSB.ImY=ssb(stk)[,as.character(ImY),drop=TRUE],
    
    TSB.ImY=quantSums(stock.wt(stk)*stock.n(stk))[,as.character(ImY),drop=TRUE],
    
    F.AdY=fbar(stk)[,as.character(AdY),drop=TRUE],
    Catch.AdY=computeCatch(stk)[,as.character(AdY),drop=TRUE],
    SSB.AdY=ssb(stk)[,as.character(AdY),drop=TRUE],
    
    TSB.AdY=quantSums(stock.wt(stk)*stock.n(stk))[,as.character(AdY),drop=TRUE] ,
    
    SSB.CtY=ssb(stk)[,as.character(CtY),drop=TRUE],
    TSB.CtY=quantSums(stock.wt(stk)*stock.n(stk))[,as.character(CtY),drop=TRUE] )
    
    })
    options.sum.tbl <- t(options.sum.tbl)
    options.sum.tbl[,1]<-names(stcks)
    colnames(options.sum.tbl) <- c("Rationale",
    sprintf("Fbar (%i)",ImY),sprintf("Catch (%i)",ImY),sprintf("SSB (%i)",ImY),sprintf("TSB (%i)",ImY),
    sprintf("Fbar (%i)",AdY),sprintf("Catch (%i)",AdY),sprintf("SSB (%i)",AdY),sprintf("TSB (%i)",AdY),
    sprintf("SSB (%i)",CtY),sprintf("TSB (%i)",CtY) )
    print(options.sum.tbl)
    return(options.sum.tbl)

}

ImY <- iy
AdY <- ay 
CtY <- ay+1
getwd()


opt.sum.tbl_res <- opt.sum.tbl(stcks=runs)
write.csv(opt.sum.tbl_res,file="./report/options - summary_3.csv",row.names=FALSE)

#table 1 of advice sheet

  #F sq
    fbar(runs[[1]][,ac(2024)])
  #SSB
    ssb(runs[[1]][,ac(2024)])
  #rec
    rec(runs[[1]][,ac(2024:(ay+1))])
  #catch
   catch(runs[[1]][,ac(2024)])

  #for summary table
   
   #rec
   rec(NSHOM.sam)[,c("year", "lbnd","value","ubnd")]
   rec(fut[, ac(2024)])

   #ssb
   ssb(NSHOM.sam)[,c("year", "lbnd","value","ubnd")]
   ssb(fut[, ac(2024)])
   
   #ssb
   fbar(NSHOM.sam)[,c("year", "lbnd","value","ubnd")]
   fbar(fut[, ac(2024)])

   #I show the data, instead of the model estimate
   as.data.frame(catch(fut))["data"]






# --- McMC (stochstic excersize)
#number of iterations
n <- 1000

# 3 years, 5 years wts/selex
set.seed(123)

NSHOM.mc <- simulate(NSHOM,NSHOM.tun,NSHOM.ctrl,n)
NSHOM.mc_tst <- simulate(NSHOM,NSHOM.tun,NSHOM.ctrl,n)

#- New starting conditions to projected object
futmc <- stf(propagate(NSHOM, n),
  nyears=3, wts.nyears = 5, fbar.nyears=1) #status quo fishing mortality
  #whats different between iterations is stock.n and harvest, but biology is the same

for(i in 1:n){
  dmns <- dimnames(NSHOM.mc[[i]]@stock.n)$year                   
  futmc@stock.n[,dmns,,,,i] <- NSHOM.mc[[i]]@stock.n
  futmc@harvest[,dmns,,,,i] <- NSHOM.mc[[i]]@harvest
}

#plot(futmc)

#- Same sample assumption but with n iterations
recdrmc <- rec(NSHOM)[,ac(sample(ry,n*3,replace=T))]
dimnames(recdrmc)$age <- "a"

NSsrrmc <- predictModel(model=rec~a, params=FLPar(c(recdrmc), units="thousands",
  dimnames=list(params="a", year=seq(iy, length=3), iter=1:n)))

#load(file="forecast/int_save_runsmc.rda", verbose =TRUE)
# INTERMEDIATE year assumption
discards.wt(futmc)<-0
discards.wt(futmc)<-0

# UPDATE 2025
futmc <- fwd(futmc, sr=NSsrr, f=iy_targets["Fsq",ac(2024)]) #ask thomas waht to fill in here. 

  # put in the table 1 of advice sheet:advice sheet 
    fbar(futmc[, ac(2024)])
    ssb(futmc[, ac(2024)])
    
    simrec<- NSsrrmc
    simrec_df <- as.data.frame(simrec@params)
    mean(exp(aggregate(log(data) ~ year ,data = simrec_df, FUN = mean)[,2]))
    
    
    exp(yearMeans((log(params(simrec)))))

    mean((an(rec(fut[, ac(2014:2023)]))))
    catch(futmc[, ac(2024)])
    landings(futmc[, ac(2024)])

# RUN w/McMC!
library(doParallel)
ncores <- ifelse(n<detectCores()-1,n,detectCores()-1)
cl    <- makeCluster(ncores)
clusterEvalQ(cl,library(FLasher))
registerDoParallel(cl)

mcruns <- foreach(i = 1:1000) %dopar% FLStocks(lapply(fctls, function(x) fwd(iter(futmc,i), sr=iter(NSsrrmc,i), control=x)))

# SAVE
save(mcruns, file="forecast/runsmc.rda", compress="xz")
load(file="forecast/runsmc.rda", verbose =TRUE)

#- Resuls for runs
  # - advice
    ssbproj   <- do.call(rbind,lapply(mcruns,function(x){ssb(x$advice[,ac(iy:fy)])}))
    apply(ssbproj,2,mean)     
    
    fbarproj  <- do.call(rbind,lapply(mcruns,function(x){fbar(x$advice[,ac(iy:fy)])}))
    
    apply(ssbproj,2,mean)     
    
    recproj   <- do.call(rbind,lapply(mcruns,function(x){rec(x$advice[,ac(iy:fy)])}))
    apply(ssbproj,2,mean)     

  # - F0
    ssbproj <- do.call(rbind,lapply(mcruns,function(x){ssb(x$F0[,ac(iy:fy)])}))
    ssbproj apply(ssbproj,2,median) 
    
    fbarproj<- do.call(rbind,lapply(mcruns,function(x){ssb(x$F0[,ac(iy:fy)])}))
  
  # - F0
    ssbproj <- do.call(rbind,lapply(mcruns,function(x){ssb(x$F0[,ac(iy:fy)])}))
    fbarproj<- do.call(rbind,lapply(mcruns,function(x){ssb(x$F0[,ac(iy:fy)])}))


apply(ssbproj,2,median) 
 
ssbproj <- do.call(rbind,lapply(mcruns,function(x){ssb(x$F0[,ac(iy:fy)])}))
fbarproj<- do.call(rbind,lapply(mcruns,function(x){fbar(x$F0[,ac(iy:fy)])}))
#fbarproj<- do.call(rbind,lapply(mcruns,function(x){ssb(x$F0[,ac(iy:fy)])}))

mean(fbarproj[,1])
fbar(fut)