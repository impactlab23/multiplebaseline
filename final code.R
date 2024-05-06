library(lmerTest)
library(geepack)
library(geesmv)
library(doParallel)


## Data generator
datagen <- function(
    N.people,  # Number of patients enrolled
    N.timepoints, # Number of data points collected
    baseline, # Average baseline score
    TE.period, # time effect
    TE.adh, # Treatment Effect
    TE.ex,
    TE.adh.ex,
    #TE.ex.adh,
    i.error,  # Individual level standard deviation from baseline 
    ran.error, # Across individual difference in baseline score 
    T.error, # Across individual difference in treatment effect
    start1, # Time points at which exercise is started
    transit, # number of period to transitioning
    design2,
    mod3
    
){
  N.obs <- N.people*N.timepoints
  id <- rep(1:N.people, each=N.timepoints)
  baseline.i <- baseline + rnorm(n = N.people, mean= 0, sd = ran.error)
  mu0 <- rep(baseline.i, each=N.timepoints)
  period <- rep(0:(N.timepoints-1), times=N.people)
  seqassign <- rep(rep(0:1, each=N.people/2), each=N.timepoints)
  trtassign_1 <- rep(sample(start1, N.people, replace=TRUE), each=N.timepoints)
  trtassign_2 <- trtassign_1 + transit
  if(design2 == T){
    trtassign_2[which(trtassign_1==start1[3])] <- N.timepoints+1
  }
  #trt_adh <- as.numeric(seqassign==0 & period >= trtassign_1 & period < trtassign_2)
  trt_adhex <- as.numeric(period >= trtassign_2)
  trt_adh <- ifelse(trt_adhex == 1, 1, as.numeric(seqassign==0 & period >= trtassign_1))
  #trt_adhex <- as.numeric(seqassign==0 & period >= trtassign_2)
  trt_ex <- ifelse(trt_adhex == 1, 1, as.numeric(seqassign==1 & period >= trtassign_1))
  #trt_adhex <- as.numeric(period >= trtassign_2)
  #trt_adhex <- as.numeric(period >= trtassign_2)
  #trt_exadh <- as.numeric(seqassign==1 & period >= trtassign_2)
  control <- as.numeric(trt_adh==0 & trt_adhex==0 & trt_ex==0)
  adh_alone <- as.numeric(trt_adh==1 & trt_adhex==0 & trt_ex==0)
  ex_alone <- as.numeric(trt_adh==0 & trt_adhex==0 & trt_ex==1)
  
  
  TE.adh.i <- rep(TE.adh + rnorm(n = N.people, mean= 0, sd = T.error), each=N.timepoints)
  TE.ex.i <-  rep(TE.ex + rnorm(n = N.people, mean= 0, sd = T.error), each=N.timepoints)
  TE.adh.ex.i <-  rep(TE.adh.ex + rnorm(n = N.people, mean= 0, sd = T.error), each=N.timepoints)
  #TE.ex.adh.i <-  rep(TE.ex.adh + rnorm(n = N.people, mean= 0, sd = T.error), each=N.timepoints)
  
  
  if(mod3 == T){
    y <- mu0 + period*TE.period + 
      TE.adh.i*adh_alone + 
      TE.adh.ex.i*trt_adhex + 
      TE.ex.i*ex_alone +
      #TE.ex.adh.i*trt_exadh + 
      rnorm(N.obs, sd=i.error)
  } else {
  
  y <- mu0 + period*TE.period + 
    TE.adh.i*trt_adh + 
    TE.adh.ex.i*trt_adhex + 
    TE.ex.i*trt_ex +
    #TE.ex.adh.i*trt_exadh + 
    rnorm(N.obs, sd=i.error)}
  
  #data <- data.frame(id, y, mu0=mu0, seq = seqassign, period ,control, trt_adh, trt_adhex, trt_ex, trt_exadh, trtassign_1, trtassign_2)
  data <- data.frame(id, y, mu0=mu0, seq = seqassign, period ,control, trt_adh, trt_adhex, trt_ex, adh_alone,
                     ex_alone, trtassign_1, trtassign_2)
  return(data)
}


## FUNCTION TO FIT DATA
gee_fit <- function(data){
  gee.mod0 <- summary(geeglm(y ~ period + trt_adh + trt_ex, 
                             data=data, id=id, family=gaussian, corstr="exchangeable"))
  gee.mod1 <- summary(geeglm(y ~ period + trt_adh + trt_ex + trt_adhex, 
                             data=data, id=id, family=gaussian, corstr="exchangeable"))
  gee.mod2 <- summary(geeglm(y ~ period + adh_alone + ex_alone + trt_adhex, 
                             data=data, id=id, family=gaussian, corstr="exchangeable"))
  md.mod0 <- GEE.var.md(y ~ period + trt_adh + trt_ex, 
                        data=data, id="id", family=gaussian, corstr="exchangeable")
  md.mod1 <- GEE.var.md(y ~ period + trt_adh + trt_ex + trt_adhex, 
                        data=data, id="id", family=gaussian, corstr="exchangeable")
  md.mod2 <- GEE.var.md(y ~  period + adh_alone + ex_alone + trt_adhex,  
                        data=data, id="id", family=gaussian, corstr="exchangeable")
  md.ses.mod0 <- sqrt(md.mod0$cov.beta[-c(1:2)])
  md.ses.mod1 <- sqrt(md.mod1$cov.beta[-c(1:2)])
  md.ses.mod2 <- sqrt(md.mod2$cov.beta[-c(1:2)])
  coefs.mod0 <- gee.mod0$coefficients[3:4,1]
  coefs.mod1 <- gee.mod1$coefficients[3:5,1]
  coefs.mod2 <- gee.mod2$coefficients[3:5,1]
  ses.mod0 <-  gee.mod0$coefficients[3:4,2]
  ses.mod1 <-  gee.mod1$coefficients[3:5,2]
  ses.mod2 <-  gee.mod2$coefficients[3:5,2]
  p.values.mod0 <- gee.mod0$coefficients[3:4,4]
  p.values.mod1 <- gee.mod1$coefficients[3:5,4]
  p.values.mod2 <- gee.mod2$coefficients[3:5,4]
  return(list(coefs.mod0 = coefs.mod0, coefs.mod1 = coefs.mod1, coefs.mod2 = coefs.mod2, 
              ses.mod0 = ses.mod0, ses.mod1 = ses.mod1, ses.mod2 = ses.mod2, 
              p.values.mod0 = p.values.mod0, p.values.mod1 = p.values.mod1,  p.values.mod2 = p.values.mod2, 
              md.ses.mod0= md.ses.mod0, md.ses.mod1= md.ses.mod1, md.ses.mod2= md.ses.mod2))
}
lmm_fit <- function(data){
  fit0 <- summary(lmerTest::lmer(y ~ period + trt_adh + trt_ex  + (1|id), data))
  fit1 <- summary(lmerTest::lmer(y ~ period + trt_adh + trt_ex  + trt_adhex + (1|id), data))
  fit2 <- summary(lmerTest::lmer(y ~ period + adh_alone + ex_alone + trt_adhex  + (1|id), data))
  coefs0 <- fit0$coefficients[3:4,1]
  coefs1 <- fit1$coefficients[3:5,1]
  coefs2 <- fit2$coefficients[3:5,1]
  ses0 <-  fit0$coefficients[3:4,2]
  ses1 <-  fit1$coefficients[3:5,2]
  ses2 <-  fit2$coefficients[3:5,2]
  p.values0 <- fit0$coefficients[3:4,5]
  p.values1 <- fit1$coefficients[3:5,5]
  p.values2 <- fit2$coefficients[3:5,5]
  return(list(coefs0 = coefs0, coefs1 = coefs1, coefs2 = coefs2,
              ses0 = ses0, ses1 = ses1, ses2 = ses2,
              p.values0 = p.values0, p.values1 = p.values1, p.values2 = p.values2))
}


## OUTPUT SUMMARY
#calculate pvalue
temp <- c()
cal_ci <- function(coef, se, df, t){
  col_up <- coef + qt(df=df,t)*se
  col_low <- coef - qt(df=df,t)*se
  for (i in 1:length(col_up)){
    temp[i] <- (col_up[i] > 0 & col_low[i] < 0)
  }
  power <- 1- mean(temp)
  return (power)
}
final_summary <- function(results.lmm, results.gee, df0, df){
  coef.mat1 <- ses.mat1 <- p.value.mat1 <- coef.mat2 <- ses.mat2 <- p.value.mat2 <-matrix(NA, nrow = length(results.lmm), ncol = 3) 
  coef.mat0 <- ses.mat0 <- p.value.mat0 <-matrix(NA, nrow = length(results.lmm), ncol = 2) 
  
  for(i in 1:length(results.lmm)){
    coef.mat0[i,] <- results.lmm[[i]]$coefs0
    coef.mat1[i,] <- results.lmm[[i]]$coefs1
    coef.mat2[i,] <- results.lmm[[i]]$coefs2
    ses.mat0[i,] <- results.lmm[[i]]$ses0
    ses.mat1[i,] <- results.lmm[[i]]$ses1
    ses.mat2[i,] <- results.lmm[[i]]$ses2
    p.value.mat0[i,] <- results.lmm[[i]]$p.values0
    p.value.mat1[i,] <- results.lmm[[i]]$p.values1
    p.value.mat2[i,] <- results.lmm[[i]]$p.values2
  }
  
  coef.mat.gee1 <- coef.mat.gee2 <-  md.mat.gee1 <- md.mat.gee2  <- ses.mat.gee1 <- ses.mat.gee2  <- p.value.mat.gee1 <- p.value.mat.gee2 <- matrix(NA, nrow = length(results.gee), ncol = 3) 
  coef.mat.gee0 <- md.mat.gee0 <- ses.mat.gee0 <- p.value.mat.gee0 <- matrix(NA, nrow = length(results.gee), ncol = 2) 
  for(i in 1:length(results.gee)){
    md.mat.gee0[i,] <- results.gee[[i]]$md.ses.mod0
    md.mat.gee1[i,] <- results.gee[[i]]$md.ses.mod1
    md.mat.gee2[i,] <- results.gee[[i]]$md.ses.mod2
    coef.mat.gee0[i,] <- results.gee[[i]]$coefs.mod0
    coef.mat.gee1[i,] <- results.gee[[i]]$coefs.mod1
    coef.mat.gee2[i,] <- results.gee[[i]]$coefs.mod2
    ses.mat.gee0[i,] <- results.gee[[i]]$ses.mod0
    ses.mat.gee1[i,] <- results.gee[[i]]$ses.mod1
    ses.mat.gee2[i,] <- results.gee[[i]]$ses.mod2
    p.value.mat.gee0[i,] <- results.gee[[i]]$p.values.mod0
    p.value.mat.gee1[i,] <- results.gee[[i]]$p.values.mod1
    p.value.mat.gee2[i,] <- results.gee[[i]]$p.values.mod2
  }
  
  est0 <- c(colMeans(coef.mat0), colMeans(coef.mat.gee0))
  est1 <- c(colMeans(coef.mat1), colMeans(coef.mat.gee1))
  est2 <- c(colMeans(coef.mat2), colMeans(coef.mat.gee2))
  power.unadjust0 <- c(colMeans(p.value.mat0 <=0.05/2),colMeans(p.value.mat.gee0 <=0.05/2))
  power.unadjust1 <- c(colMeans(p.value.mat1 <=0.05/3),colMeans(p.value.mat.gee1 <=0.05/3))
  power.unadjust2 <- c(colMeans(p.value.mat2 <=0.05/3),colMeans(p.value.mat.gee2 <=0.05/3))
  power.adjust0 <- c(cal_ci(coef.mat.gee0[,1],md.mat.gee0[,1], df0, t = 1-0.05/4 ), 
                     cal_ci(coef.mat.gee0[,2],md.mat.gee0[,2], df0, t = 1-0.05/4 ))
  power.adjust1 <- c(cal_ci(coef.mat.gee1[,1],md.mat.gee1[,1], df, t = 1-0.05/6 ), 
                     cal_ci(coef.mat.gee1[,2],md.mat.gee1[,2], df, t = 1-0.05/6 ),
                     cal_ci(coef.mat.gee1[,3],md.mat.gee1[,3], df, t = 1-0.05/6 ))
  power.adjust2 <- c(cal_ci(coef.mat.gee2[,1],md.mat.gee2[,1], df, t = 1-0.05/6 ), 
                     cal_ci(coef.mat.gee2[,2],md.mat.gee2[,2], df, t = 1-0.05/6 ),
                     cal_ci(coef.mat.gee2[,3],md.mat.gee2[,3], df, t = 1-0.05/6 ))
  return(c(est0, est1, est2, 
           power.unadjust0, power.unadjust1, power.unadjust2, 
           power.adjust0, power.adjust1, power.adjust2))
}





##=======================================================================
##EVALUATION
##=======================================================================
##DATA GENERATION FOR TYPE I ERROR
## WITHOUT INTERACTION EFFECT
set.seed(2024)
type_30 <- replicate(5000, datagen(
  N.people = 30, 
  N.timepoints = 5,
  baseline = 0, 
  TE.period = 1,
  TE.adh = 0, 
  TE.ex = 0,
  TE.adh.ex = 0,
  i.error = 1, 
  ran.error = 1,
  T.error = 0,
  start1 = c(1,2,3),
  transit = 1,
  design2 = F,
  mod3 = T
), simplify = F)

set.seed(2024)
type_60 <- replicate(5000, datagen(
  N.people = 60, 
  N.timepoints = 5,
  baseline = 0, 
  TE.period = 1,
  TE.adh = 0, 
  TE.ex = 0,
  TE.adh.ex = 0,
  i.error = 1, 
  ran.error = 1,
  T.error = 0,
  start1 = c(1,2,3),
  transit = 1,
  design2 = F,
  mod3 = T
), simplify = F)

set.seed(2024)
type_120 <- replicate(5000, datagen(
  N.people = 120, 
  N.timepoints = 5,
  baseline = 0, 
  TE.period = 1,
  TE.adh = 0, 
  TE.ex = 0,
  TE.adh.ex = 0,
  i.error = 1, 
  ran.error = 1,
  T.error = 0,
  start1 = c(1,2,3),
  transit = 1,
  design2 = F,
  mod3 = T
), simplify = F)


# POWER 
## WITHOUT INTERACTION EFFECT
set.seed(2024)
no_pwr_30 <- replicate(5000, datagen(
  N.people = 30, 
  N.timepoints = 5,
  baseline = 0, 
  TE.period = 1,
  TE.adh = 0.8, 
  TE.ex = 0.8,
  TE.adh.ex = 1.6,
  i.error = 1, 
  ran.error = 1,
  T.error = 0,
  start1 = c(1,2,3),
  transit = 1,
  design2 = F,
  mod3 = T
), simplify = F)

set.seed(2024)
no_pwr_60 <- replicate(5000, datagen(
  N.people = 60, 
  N.timepoints = 5,
  baseline = 0, 
  TE.period = 1,
  TE.adh = 0.8, 
  TE.ex = 0.8,
  TE.adh.ex = 1.6,
  i.error = 1, 
  ran.error = 1,
  T.error = 0,
  start1 = c(1,2,3),
  transit = 1,
  design2 = F,
  mod3 = T
), simplify = F)

set.seed(2024)
no_pwr_120 <- replicate(5000, datagen(
  N.people = 120, 
  N.timepoints = 5,
  baseline = 0, 
  TE.period = 1,
  TE.adh = 0.8, 
  TE.ex = 0.8,
  TE.adh.ex = 1.6,
  i.error = 1, 
  ran.error = 1,
  T.error = 0,
  start1 = c(1,2,3),
  transit = 1,
  design2 = F,
  mod3 = T
), simplify = F)




#WITH INTERACTION
set.seed(2024)
int_pwr_30 <- replicate(5000, datagen(
  N.people = 30, 
  N.timepoints = 5,
  baseline = 0, 
  TE.period = 1,
  TE.adh = 0.8, 
  TE.ex = 0.8,
  TE.adh.ex = 2.0,
  i.error = 1, 
  ran.error = 1,
  T.error = 0,
  start1 = c(1,2,3),
  transit = 1,
  design2 = F,
  mod3 = T
), simplify = F)

set.seed(2024)
int_pwr_60 <- replicate(5000, datagen(
  N.people = 60, 
  N.timepoints = 5,
  baseline = 0, 
  TE.period = 1,
  TE.adh = 0.8, 
  TE.ex = 0.8,
  TE.adh.ex = 2.0,
  i.error = 1, 
  ran.error = 1,
  T.error = 0,
  start1 = c(1,2,3),
  transit = 1,
  design2 = F,
  mod3 = T
), simplify = F)


set.seed(2024)
int_pwr_120 <- replicate(5000, datagen(
  N.people = 120, 
  N.timepoints = 5,
  baseline = 0, 
  TE.period = 1,
  TE.adh = 0.8, 
  TE.ex = 0.8,
  TE.adh.ex = 2.0,
  i.error = 1, 
  ran.error = 1,
  T.error = 0,
  start1 = c(1,2,3),
  transit = 1,
  design2 = F,
  mod3 = T
), simplify = F)

# PARALLEL COMPUTING
ncores <- 4
registerDoParallel(ncores)
type.gee_30 <- foreach(i=1:5000, .packages=c("geepack", "geesmv"), .errorhandling = 'remove') %dopar% {gee_fit(type_30[[i]])}
type.lmm_30 <- foreach(i=1:5000, .packages=c("lmerTest"), .errorhandling = 'remove') %dopar% {lmm_fit(type_30[[i]])}
type.gee_60 <- foreach(i=1:5000, .packages=c("geepack", "geesmv"), .errorhandling = 'remove') %dopar% {gee_fit(type_60[[i]])}
type.lmm_60 <- foreach(i=1:5000, .packages=c("lmerTest"), .errorhandling = 'remove') %dopar% {lmm_fit(type_60[[i]])}
type.gee_120 <- foreach(i=1:5000, .packages=c("geepack", "geesmv"), .errorhandling = 'remove') %dopar% {gee_fit(type_120[[i]])}
type.lmm_120 <- foreach(i=1:5000, .packages=c("lmerTest"), .errorhandling = 'remove') %dopar% {lmm_fit(type_120[[i]])}



pwr.gee_30 <- foreach(i=1:5000, .packages=c("geepack", "geesmv"), .errorhandling = 'remove') %dopar% {gee_fit(no_pwr_30[[i]])}
pwr.lmm_30 <- foreach(i=1:5000, .packages=c("lmerTest"), .errorhandling = 'remove') %dopar% {lmm_fit(no_pwr_30[[i]])}
pwr.gee_60 <- foreach(i=1:5000, .packages=c("geepack", "geesmv"), .errorhandling = 'remove') %dopar% {gee_fit(no_pwr_60[[i]])}
pwr.lmm_60 <- foreach(i=1:5000, .packages=c("lmerTest"), .errorhandling = 'remove') %dopar% {lmm_fit(no_pwr_60[[i]])}
pwr.gee_120 <- foreach(i=1:5000, .packages=c("geepack", "geesmv"), .errorhandling = 'remove') %dopar% {gee_fit(no_pwr_120[[i]])}
pwr.lmm_120 <- foreach(i=1:5000, .packages=c("lmerTest"), .errorhandling = 'remove') %dopar% {lmm_fit(no_pwr_120[[i]])}

int.pwr.gee_30 <- foreach(i=1:5000, .packages=c("geepack", "geesmv"), .errorhandling = 'remove') %dopar% {gee_fit(int_pwr_30[[i]])}
int.pwr.lmm_30 <- foreach(i=1:5000, .packages=c("lmerTest"), .errorhandling = 'remove') %dopar% {lmm_fit(int_pwr_30[[i]])}
int.pwr.gee_60 <- foreach(i=1:5000, .packages=c("geepack", "geesmv"), .errorhandling = 'remove') %dopar% {gee_fit(int_pwr_60[[i]])}
int.pwr.lmm_60 <- foreach(i=1:5000, .packages=c("lmerTest"), .errorhandling = 'remove') %dopar% {lmm_fit(int_pwr_60[[i]])}
int.pwr.gee_120 <- foreach(i=1:5000, .packages=c("geepack", "geesmv"), .errorhandling = 'remove') %dopar% {gee_fit(int_pwr_120[[i]])}
int.pwr.lmm_120 <- foreach(i=1:5000, .packages=c("lmerTest"), .errorhandling = 'remove') %dopar% {lmm_fit(int_pwr_120[[i]])}


fac.pwr.gee_30 <- foreach(i=1:5000, .packages=c("geepack", "geesmv"), .errorhandling = 'remove') %dopar% {gee_fit(data_factorial[[i]])}
fac.pwr.lmm_30 <- foreach(i=1:5000, .packages=c("lmerTest"), .errorhandling = 'remove') %dopar% {lmm_fit(data_factorial[[i]])}

fac.type.gee_30 <- foreach(i=1:5000, .packages=c("geepack", "geesmv"), .errorhandling = 'remove') %dopar% {gee_fit(null_factorial[[i]])}
fac.typ3.lmm_30 <- foreach(i=1:5000, .packages=c("lmerTest"), .errorhandling = 'remove') %dopar% {lmm_fit(null_factorial[[i]])}


type_results <- rbind(final_summary(type.lmm_30, type.gee_30, df0=26, df= 25),
                       final_summary(type.lmm_60, type.gee_60, df0=56, df= 55),
                       final_summary(type.lmm_120, type.gee_120, df0=116, df= 115))
type_results <- as.data.frame(type_results)
names(type_results) <- c("LMM0_M1", "LMM0_M2", "GEE0_M1", "GEE0_M2",
                          "LMM1_M1", "LMM1_M2","LMM1_IN", "GEE1_M1", "GEE1_M2","GEE1_IN",
                          "LMM2_M1", "LMM2_M2","LMM2_M3", "GEE2_M1", "GEE2_M2","GEE2_M3",
                          "PWR0_LMMM1", "PWR0_LMMM2", "GEE0_LMMM1", "GEE0_LMMM2",
                          "PWR1_LMMM1", "PWR1_LMMM2", "PWR1_LMMIN", "PWR1_GEEM1", "PWR1_GEEM2", "PWR1_GEEIN",
                          "PWR2_LMMM1", "PWR2_LMMM2", "PWR2_LMMM3", "PWR2_GEEM1", "PWR_GEEM2", "PWR_GEEM3",
                          "ADPWR0_GEEM1", "ADPWR0_GEEM2",
                          "ADPWR1_GEEM1", "ADPWR1_GEEM2", "ADPWR1_GEEIN", 
                          "ADPWR2_GEEM1", "ADPWR2_GEEM2", "ADPWR2_GEEM3")


no_pwr_results <- rbind(final_summary(pwr.lmm_30, pwr.gee_30, df0=26, df= 25),
                         final_summary(pwr.lmm_60, pwr.gee_60, df0=56, df= 55),
                         final_summary(pwr.lmm_120, pwr.gee_120, df0=116, df= 115))
no_pwr_results <- as.data.frame(no_pwr_results)
names(no_pwr_results) <- c("LMM0_M1", "LMM0_M2", "GEE0_M1", "GEE0_M2",
                         "LMM1_M1", "LMM1_M2","LMM1_IN", "GEE1_M1", "GEE1_M2","GEE1_IN",
                         "LMM2_M1", "LMM2_M2","LMM2_M3", "GEE2_M1", "GEE2_M2","GEE2_M3",
                         "PWR0_LMMM1", "PWR0_LMMM2", "GEE0_LMMM1", "GEE0_LMMM2",
                         "PWR1_LMMM1", "PWR1_LMMM2", "PWR1_LMMIN", "PWR1_GEEM1", "PWR1_GEEM2", "PWR1_GEEIN",
                         "PWR2_LMMM1", "PWR2_LMMM2", "PWR2_LMMM3", "PWR2_GEEM1", "PWR_GEEM2", "PWR_GEEM3",
                         "ADPWR0_GEEM1", "ADPWR0_GEEM2",
                         "ADPWR1_GEEM1", "ADPWR1_GEEM2", "ADPWR1_GEEIN", 
                         "ADPWR2_GEEM1", "ADPWR2_GEEM2", "ADPWR2_GEEM3")


int_pwr_results <- rbind(final_summary(int.pwr.lmm_30, int.pwr.gee_30, df0=26, df= 25),
                        final_summary(int.pwr.lmm_60, int.pwr.gee_60, df0=56, df= 55),
                        final_summary(int.pwr.lmm_120, int.pwr.gee_120, df0=116, df= 115))
int_pwr_results <- as.data.frame(int_pwr_results)
names(int_pwr_results) <- c("LMM0_M1", "LMM0_M2", "GEE0_M1", "GEE0_M2",
                         "LMM1_M1", "LMM1_M2","LMM1_IN", "GEE1_M1", "GEE1_M2","GEE1_IN",
                         "LMM2_M1", "LMM2_M2","LMM2_M3", "GEE2_M1", "GEE2_M2","GEE2_M3",
                         "PWR0_LMMM1", "PWR0_LMMM2", "GEE0_LMMM1", "GEE0_LMMM2",
                         "PWR1_LMMM1", "PWR1_LMMM2", "PWR1_LMMIN", "PWR1_GEEM1", "PWR1_GEEM2", "PWR1_GEEIN",
                         "PWR2_LMMM1", "PWR2_LMMM2", "PWR2_LMMM3", "PWR2_GEEM1", "PWR2_GEEM2", "PWR2_GEEM3",
                         "ADPWR0_GEEM1", "ADPWR0_GEEM2",
                         "ADPWR1_GEEM1", "ADPWR1_GEEM2", "ADPWR1_GEEIN", 
                         "ADPWR2_GEEM1", "ADPWR2_GEEM2", "ADPWR2_GEEM3")
write.csv(type_results, "type_result_new.csv")
write.csv(no_pwr_results, "no_pwr_results_new.csv")
write.csv(int_pwr_results, "int_pwr_results_new.csv")


fac.pwr <- rbind(final_summary(fac.pwr.lmm_30, fac.pwr.gee_30, df0=26, df= 25))
fac.type <- rbind(final_summary(fac.typ3.lmm_30, fac.type.gee_30, df0=26, df= 25))
fac.pwr <- as.data.frame(fac.pwr)
names(fac.pwr) <- c("LMM0_M1", "LMM0_M2", "GEE0_M1", "GEE0_M2",
                            "LMM1_M1", "LMM1_M2","LMM1_IN", "GEE1_M1", "GEE1_M2","GEE1_IN",
                            "LMM2_M1", "LMM2_M2","LMM2_M3", "GEE2_M1", "GEE2_M2","GEE2_M3",
                            "PWR0_LMMM1", "PWR0_LMMM2", "GEE0_LMMM1", "GEE0_LMMM2",
                            "PWR1_LMMM1", "PWR1_LMMM2", "PWR1_LMMIN", "PWR1_GEEM1", "PWR1_GEEM2", "PWR1_GEEIN",
                            "PWR2_LMMM1", "PWR2_LMMM2", "PWR2_LMMM3", "PWR2_GEEM1", "PWR2_GEEM2", "PWR2_GEEM3",
                            "ADPWR0_GEEM1", "ADPWR0_GEEM2",
                            "ADPWR1_GEEM1", "ADPWR1_GEEM2", "ADPWR1_GEEIN", 
                            "ADPWR2_GEEM1", "ADPWR2_GEEM2", "ADPWR2_GEEM3")
int_pwr_results <- as.data.frame(int_pwr_results)
names(int_pwr_results) <- c("LMM0_M1", "LMM0_M2", "GEE0_M1", "GEE0_M2",
                            "LMM1_M1", "LMM1_M2","LMM1_IN", "GEE1_M1", "GEE1_M2","GEE1_IN",
                            "LMM2_M1", "LMM2_M2","LMM2_M3", "GEE2_M1", "GEE2_M2","GEE2_M3",
                            "PWR0_LMMM1", "PWR0_LMMM2", "GEE0_LMMM1", "GEE0_LMMM2",
                            "PWR1_LMMM1", "PWR1_LMMM2", "PWR1_LMMIN", "PWR1_GEEM1", "PWR1_GEEM2", "PWR1_GEEIN",
                            "PWR2_LMMM1", "PWR2_LMMM2", "PWR2_LMMM3", "PWR2_GEEM1", "PWR2_GEEM2", "PWR2_GEEM3",
                            "ADPWR0_GEEM1", "ADPWR0_GEEM2",
                            "ADPWR1_GEEM1", "ADPWR1_GEEM2", "ADPWR1_GEEIN", 
                            "ADPWR2_GEEM1", "ADPWR2_GEEM2", "ADPWR2_GEEM3")






set.seed(2024)
example <- replicate(5000, datagen(
  N.people = 24, 
  N.timepoints = 24,
  baseline = 50, 
  TE.period = 1,
  TE.adh = 10, 
  TE.ex = 10,
  TE.adh.ex = 20,
  i.error = 7.5, 
  ran.error = 12.9,
  T.error = 7.5,
  start1 = c(4,6,8),
  transit = 8,
  design2 = F,
  mod3 = T
), simplify = F)

example_res <- foreach(i=1:5000, .packages=c("geepack", "geesmv"), .errorhandling = 'remove') %dopar% {gee_fit(example[[i]])}
example_res.lmm <- foreach(i=1:5000, .packages=c("lmerTest"), .errorhandling = 'remove') %dopar% {lmm_fit(example[[i]])}
example_results <- final_summary(example_res.lmm, example_res, df0=20, df= 19)
example_results  <- as.data.frame(example_results)
names(example_results) <- c("LMM0_M1", "LMM0_M2", "GEE0_M1", "GEE0_M2",
                    "LMM1_M1", "LMM1_M2","LMM1_IN", "GEE1_M1", "GEE1_M2","GEE1_IN",
                    "LMM2_M1", "LMM2_M2","LMM2_M3", "GEE2_M1", "GEE2_M2","GEE2_M3",
                    "PWR0_LMMM1", "PWR0_LMMM2", "GEE0_LMMM1", "GEE0_LMMM2",
                    "PWR1_LMMM1", "PWR1_LMMM2", "PWR1_LMMIN", "PWR1_GEEM1", "PWR1_GEEM2", "PWR1_GEEIN",
                    "PWR2_LMMM1", "PWR2_LMMM2", "PWR2_LMMM3", "PWR2_GEEM1", "PWR2_GEEM2", "PWR2_GEEM3",
                    "ADPWR0_GEEM1", "ADPWR0_GEEM2",
                    "ADPWR1_GEEM1", "ADPWR1_GEEM2", "ADPWR1_GEEIN", 
                    "ADPWR2_GEEM1", "ADPWR2_GEEM2", "ADPWR2_GEEM3")
