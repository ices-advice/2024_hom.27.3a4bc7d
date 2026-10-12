## Run analysis, write model results

## Before:
## After:

# - Load the required packages
 library(FLAssess)
 library(FLash)
 library(ggplotFL)
 library(FLBRP)
 library(icesAdvice)
 library(dplyr)
 library(patchwork)

# - load model output ----
 load("model/model.RData", verbose = TRUE)
 summary(NSHOM)

# - 
 
 ## - assumptions... 
 # - you know TAC will be taken. or indication in terms of %  --> use for current year
 # - otherwise constraint of Fishing mortality 
 # - biology is last 3 years 
 # - selectivity pattern 
 # - recruitment last 10 years 
  
 
# - Extending the stock object for the projections for weights, maturity ----
 maxyr_stk <- range(NSHOM)[["maxyear"]]
 NSHOM_stf <- stf(NSHOM,nyears=200,wts.nyears=20 ,na.rm=TRUE)
  catch.n(NSHOM_stf[, ac(2003:2023)])
 maxyr_stf <- range(NSHOM_stf)[["maxyear"]]
 ImY <- maxyr_stk +1
 AdY <- ImY +1
 
 NSHOM_stf@harvest[,ac(ImY:maxyr_stf)] <- yearMeans(NSHOM_stf@harvest[,ac(2018:2022)])

 cn<- as.data.frame(catch.n(NSHOM))[,c(1,2,7)]
 cn$data[is.na(cn$data)]<-0
 cn<-tidyr::spread(cn,key= age ,value= data, -c(2))
 
 agesFbar(tail(cn, 20) , plot=TRUE, main="Horse mackerel in 4b, 4c & 7d,  80% coverage of last 2 decades")
 
 ggplot(harvest(NSHOM_stf)[,ac(2019:2030)]) + geom_line(aes(x=age, y=data)) + facet_wrap(~year) + theme_bw() + ylab("Selection")


 # sel <- sweep(NSHOM_stf@harvest,2,apply(NSHOM_stf@harvest,2,max),'/')    # selection pattern
 # selmean<-yearMeans(sel[,ac((maxyr_stk-2):maxyr_stk)])                   # mean selection pattern of the last 3 years
 # fb<-yearMeans(fbar(NSHOM_stf)[,ac((maxyr_stk-2):maxyr_stk)])            # mean of  Fbars of the last 3 years
 # harvest(NSHOM_stf)[,ac(ImY:(AdY+1))]   <-sweep(selmean,2,fb,'*') 
 
 harvest(NSHOM_stf)[,ac(2025),,,,] == yearMeans(harvest(NSHOM_stf[,ac((maxyr_stk-2):maxyr_stk)])) #default takes last 3 years
 harvest(NSHOM_stf)
# - plot of the resulting extended timeseries ----
 ggplot(harvest(NSHOM_stf)[,ac((maxyr_stf-5):maxyr_stf)]) + geom_line(aes(x=age, y=data)) + facet_wrap(~year)

# - The stock numbers at age and the catch numbers at age are not forecast yet - this is what the fwd() function will perform later ---
 fbar(NSHOM_stf)
 stock.n(NSHOM_stf)[,ac((maxyr_stf-5):maxyr_stf)]


# - The stock-recruitment relationship (SRR) ----
 # A short term forecast does not use an SRR (in the traditional sense). Instead, it generally assumes that recruitment in the future is some mean (e.g. geometric mean) of the historic recruitments. However, we still need to have an SRR that contains this mean value, which is what we mimic for this example. First, we estimate the geometic mean recruitment, that we then add to an SRR object
 plot(rec(NSHOM))
 
 rec.years <- c((maxyr_stk-10):(maxyr_stk)) #rationale: since 2012 rec is stable # or maybe try from 2006
 mean_rec <- exp(mean(log(rec(trim(NSHOM, year = rec.years)))))
 mean_rec
 # for a longer term  forecastlets set a dynamic SSR based on a segmented regression
 NSHOM_sr <- as.FLSR(NSHOM, model="segreg")#geomean
 
 #params(NSHOM_sr)['a',] <- mean_rec
 NSHOM_sr <- fmle(NSHOM_sr)
 plot(NSHOM_sr)
 params(NSHOM_sr)
 rec(NSHOM_sr)
 
 params(NSHOM_sr)['b',] <- 1.79*10^5 # point of inflection as estimated from the fitted smoother during WKBHMB
 params(NSHOM_sr)['a',] <- 6.578877 # slope 

 # predict recruitment under SR relation
 b_tst <- as.vector(params(NSHOM_sr)['b',])
 a_tst <- as.vector(params(NSHOM_sr)['a',])
 
 res_plot_df<- data.frame(ssb = seq(from = 1, to = 300000, by = 100)) %>% mutate( rec = 
 ifelse(ssb <= b_tst , a_tst * ssb, a_tst * b_tst)) # 514505.7 is highest possible recruitment with dynamic option when stock is using the NSHOM specific mat
  
# ssb_test = 46067 #2025 SSB
# ssb_test = 179000 #breakpoint SSB
# 
# ifelse(ssb_test <= b_tst , a_tst * ssb_test, a_tst * b_tst) #1371806
# 
# #slope in segreg:: 911982.9/179000 = 5.094877
# #slope in loes::  exp(predict(lo,newdata=data.frame(ssb=log(179000)))) :: 1177619/179000 = 6.578877

##  loes slope with no depensation 
# 
##  let slope start from approx 55000
#  first_break <- 55000
#  model(NSHOM_sr) <- as.formula("rec ~ ifelse(ssb <= c, 0, ifelse(ssb <= b, a * (ssb - c), d * b))")
#  params(NSHOM_sr)['c',] <- first_break
#  params(NSHOM_sr)['b',] <- 1.79*10^5 #inflection point
#  params(NSHOM_sr)['a',] <- 1177619/(179000-first_break) # slope 
#  params(NSHOM_sr)['d',] <- 6.578877
#  
#  c_tst <- as.vector(params(NSHOM_sr)['c',])
#  d_tst <- as.vector(params(NSHOM_sr)['d',])
#  b_tst <- as.vector(params(NSHOM_sr)['b',])
#  a_tst <- as.vector(params(NSHOM_sr)['a',])
#  
#  ssb_tst<- 179001
#  ifelse(ssb_tst <= c_tst, 0, ifelse(ssb_tst <= b_tst, a_tst * (ssb_tst-c_tst) , d_tst * b_tst))
#  
# res_plot_df<- data.frame(ssb = seq(from = 1, to = 300000, by = 100)) %>%
#  mutate( rec = ifelse(ssb <= c_tst, 0, ifelse(ssb <= b_tst , a_tst * (ssb-c_tst), d_tst * b_tst)))
#          #  514505.7 is highest possible recruitment with dynamic option when stock is using the NSHOM specific mat

 plot_df <- as.data.frame(sr(NSHOM))[,c("year", "data", "qname")]
 plot_df <-tidyr::spread(plot_df, key = "qname", value =  "data")
 plot_df$rec <- c(plot_df$rec[-1], NA)
 
 SRpairs     <- data.frame(rec=c(rec(NSHOM[,ac(2003:2022)])),
                           ssb=c(ssb(NSHOM[,ac(2002:2021)])))
 
 #- Fit loess smoother to ssb-rec pairs
 lo    <- loess(rec ~ ssb, data=log(SRpairs),span=1) #increasing span doesn't change anything
 plot(x=SRpairs$ssb,y=SRpairs$rec,main="Smoothed SR-pairs",pch=19,xlab="SSB (tonnes)",ylab="Recruitment (thousands)",xlim=c(0,3e5),ylim=c(0,3e6))
 lines(y=exp(predict(lo,newdata=data.frame(ssb=log(seq(0,3e5,length.out=100))))),x=seq(0,3e5,length.out=100),col=2,lwd=2)
 lo$pars
 
 pred <- exp(predict(lo,newdata=data.frame(ssb=log(seq(0,3e5,length.out=1000)))))
 ssbs <- seq(0,3e5,length.out=1000)
 poi <- ssbs[as.numeric(names(which.min(abs(diff(diff(pred[which(ssbs>1.5e5 & ssbs < 2.0e5)]))))))] #point of inflection
 abline(v=poi,col=4,lty=2,lwd=2)
 
 res_plot_df_loes<-data.frame(ssb = seq(0,3e5,length.out=1000),
                              rec = exp(predict(lo,
                                                newdata = data.frame(ssb = log(seq(0,3e5,length.out=1000))))))
 # plot both the SR relation and the Loess 

 plot_df %>% ggplot(aes(x= ssb, y= rec))+
   geom_path()+ #xlim(c(0, (250000))) +ylim(c(0, 1500000))+
   geom_point()+
   geom_line(data = res_plot_df, aes(x= ssb, y= rec) )+
   geom_line(data = res_plot_df_loes, aes(x= ssb, y= rec), col = "red" )+
   #geom_smooth(col ="blue", method = "loess")+
   theme_bw()+ theme(legend.position = "bottom")+
   geom_vline(xintercept = 179000, lty = 2, col = "blue")+
   ylab("Recruitment (tonnes)") + xlab("SSB (tonnes)")
 
   ggplot(fbar(NSHOM), aes(x=year,y=data)) + geom_line()

# - Stock vs reference points and current stock status ----
 
 #SSB
 plot(ssb(NSHOM)) +
   geom_hline(yintercept = 236228, col = "orange", lty = 2) +
   annotate("text", x=2020, y=236228+5000, label="msyBtrigger") +#msyBtrigger
   geom_hline(yintercept = 179000, col = "red", lty = 2)+
   annotate("text", x=2020, y=179000 +5000, label="Blim") +#Blim+
 #F
 plot(fbar(NSHOM)) +
   geom_hline(yintercept = 0.034, col = "orange", lty = 2) +
   annotate("text", x=2020, y=0.034, label="Fmsy") +#Fmsy
   geom_hline(yintercept = 0.042, col = "red", lty = 2)+
   annotate("text", x=2020, y=0.042, label="Flim") #Flim
 
# - Set the control object - year, quantity and value for the moment ----

  # project 2023 and 2024 with TAC catch constraints first and recruitment based on geomean
  mean_rec <- exp(mean(log(rec(trim(NSHOM, year = rec.years))))) 
  NSHOM_sr_iy <- as.FLSR(NSHOM, model="geomean")
  params(NSHOM_sr_iy)['a',] <- mean_rec
  
  ctrl_target <- data.frame(year = c(2023, 2024), quantity = c("catch"), val = c(8969,9730))
  ctrl_TAC <- fwdControl(ctrl_target)
  NSHOM_zero <- fwd(NSHOM_stf, ctrl = ctrl_TAC, sr = NSHOM_sr_iy) #
  rec(NSHOM_zero)
  
  # now project remaining years for two scenarios: Monitoring TAC, and 0 catch 
  
  ctrl_target <- data.frame(year = 2024:2222, quantity = c("f",rep("catch", 198)), val = c(0,rep(1200, 198)))
  ctrl_zero_MTAC1200 <- fwdControl(ctrl_target)
  NSHOM_zero_MTAC1200 <- fwd(NSHOM_zero, ctrl = ctrl_zero_MTAC1200, sr = NSHOM_sr)
  
  ctrl_target <- data.frame(year = 2024:2222, quantity = c("f",rep("catch", 198)), val = c(0,rep(1800, 198)))
  ctrl_zero_MTAC1800 <- fwdControl(ctrl_target)
  NSHOM_zero_MTAC1800 <- fwd(NSHOM_zero, ctrl = ctrl_zero_MTAC1800, sr = NSHOM_sr)
  
  ctrl_target <- data.frame(year = 2024:2222, quantity = c("f",rep("catch", 198)), val = c(0,rep(2400, 198)))
  ctrl_zero_MTAC2400 <- fwdControl(ctrl_target)
  NSHOM_zero_MTAC2400 <- fwd(NSHOM_zero, ctrl = ctrl_zero_MTAC2400, sr = NSHOM_sr)
  
  ctrl_target <- data.frame(year = 2024:2222, quantity = c("f",rep("catch", 198)), val = c(0,rep(1200/0.89, 198)))
  ctrl_zero_rMTAC <- fwdControl(ctrl_target)
  NSHOM_zero_rMTAC <- fwd(NSHOM_zero, ctrl = ctrl_zero_rMTAC, sr = NSHOM_sr)

  ctrl_target <- data.frame(year = 2024:2222, quantity = c("f"), val = c(rep(0, 199)))
  ctrl_zero <- fwdControl(ctrl_target)
  NSHOM_zero <- fwd(NSHOM_zero, ctrl = ctrl_zero, sr = NSHOM_sr) #
  rec(NSHOM_zero)
  
  # here is another that takes the geomean instead of dynamic recruitment
  NSHOM_zero_geomen_rec <- fwd(NSHOM_zero, ctrl = ctrl_zero, sr = NSHOM_sr_iy) 
  
  ff<- FLStocks("Zero catch" = NSHOM_zero[,ac(2015:2075)],
                "Monitoring TAC (1200t)" = NSHOM_zero_MTAC1200[,ac(2015:2075)],
                "Monitoring TAC (1800t)" = NSHOM_zero_MTAC1800[,ac(2015:2075)],
                "Monitoring TAC (2400t)" = NSHOM_zero_MTAC2400[,ac(2015:2075)],
                "Monitoring TAC + bycatch (1350t)"= NSHOM_zero_rMTAC[,ac(2015:2075)])
  ggplotFL::plot(ff)+theme_bw() + scale_colour_discrete(name = "Scenario")
  
  an(ssb(ff$MTAC))/an(ssb(ff$zero))
  
  as.data.frame(ssb(NSHOM_zero))$year[as.data.frame(ssb(NSHOM_zero))$data > 179000]
<<<<<<< HEAD
   as.data.frame(ssb(NSHOM_zero))$year[as.data.frame(ssb(NSHOM_zero))$data > 236228 ]

  as.data.frame(ssb(NSHOM_zero_MTAC1200))$year[as.data.frame(ssb(NSHOM_zero_MTAC1200))$data > 179000]
   as.data.frame(ssb(NSHOM_zero_MTAC1200))$year[as.data.frame(ssb(NSHOM_zero_MTAC1200))$data > 236228 ]
     as.data.frame(ssb(NSHOM_zero_MTAC1800))$year[as.data.frame(ssb(NSHOM_zero_MTAC1800))$data > 179000]
   as.data.frame(ssb(NSHOM_zero_MTAC1800))$year[as.data.frame(ssb(NSHOM_zero_MTAC1800))$data > 236228 ]
     as.data.frame(ssb(NSHOM_zero_MTAC2400))$year[as.data.frame(ssb(NSHOM_zero_MTAC2400))$data > 179000]
   as.data.frame(ssb(NSHOM_zero_MTAC2400))$year[as.data.frame(ssb(NSHOM_zero_MTAC2400))$data > 236228 ]
     as.data.frame(ssb(NSHOM_zero_rMTAC))$year[as.data.frame(ssb(NSHOM_zero_rMTAC))$data > 179000]
   as.data.frame(ssb(NSHOM_zero_rMTAC))$year[as.data.frame(ssb(NSHOM_zero_rMTAC))$data > 236228 ]
=======
  ssb(NSHOM_zero)[, ac(2020:2026)]

>>>>>>> d786958f017810b151c6d6f4075a7bfe7f0bfa71

# - Do some checks ----
 #recruitment in the forecast stock indeed corresponds to the mean recruitment
 mean_rec
 rec(NSHOM_zero)[,ac((maxyr_stf-5):maxyr_stf)]
 
 #ssb
 
 (plot(ssb(trim(NSHOM_zero, year=2002:2025)))+
   geom_hline(yintercept = 236228, col = "orange", lty = 2) +
   annotate("text", x=2020, y=236228+5000, label="msyBtrigger") +#msyBtrigger
   geom_hline(yintercept = 179000, col = "red", lty = 2)+
   annotate("text", x=2020, y=179000 +5000, label="Blim")+
   ggtitle(paste0("Zero catch advice shrt, terminal year ssb = ", round(ssb(NSHOM_zero)[,ac(2025)]))))+
 (plot(ssb(NSHOM_zero))+
   geom_hline(yintercept = 236228, col = "orange", lty = 2) +
   annotate("text", x=2020, y=236228+5000, label="msyBtrigger") +#msyBtrigger
   geom_hline(yintercept = 179000, col = "red", lty = 2)+
   geom_vline(xintercept = 2025, col = "purple", lty = 2)+
   ggtitle(paste0("Zero catch advice, terminal year ssb = ", round(ssb(NSHOM_zero)[,ac(2222)])),0)+
   annotate("text", x=2020, y=179000 +5000, label="Blim"))

 
 #see if fbar in the forecast stock indeed corresponds to status quo fishing mortality.
 round(fbar_SQ,3)
 round(fbar(NSHOM_zero)[,ac((maxyr_stk):maxyr_stf)],3)
 round(catch(NSHOM_zero)[,ac((maxyr_stk):maxyr_stf)],3)

 #The stock numbers are calculated using the recruitment and future mortality assumptions.
 stock.n(NSHOM_zero)[,ac((maxyr_stf-5):maxyr_stf)]

 #Note that the future harvest slot is different from the harvest from the one we set up, but that the selection pattern is the same: F at age has been multiplied by a constant value to give us the target Fbar value.
 round(harvest(NSHOM_stf)[,ac((maxyr_stf-2):maxyr_stf)],3)
 round(harvest(NSHOM_zero)[,ac((maxyr_stf-2):maxyr_stf)],3)
 harvest(NSHOM_stf)[,ac((maxyr_stf-2):maxyr_stf)] / harvest(NSHOM_zero)[,ac((maxyr_stf-2):maxyr_stf)]

 #The catch numbers come from the predicted abundance and harvest rates using the Baranov equation. The catches are then split into landings and discards using the average ratios computed by stf().
 landings.n(NSHOM_zero)[,ac((maxyr_stf-5):maxyr_stf)]
 discards.n(NSHOM_zero)[,ac((maxyr_stf-5):maxyr_stf)]

# - plot the projection ----
 plot(NSHOM_zero)

# - Save forcast object ----
 save(NSHOM_zero, file="model/LTF_HOM.RData")
 load(file="model/LTF_HOM.RData", verbose =TRUE)



 
 
 
