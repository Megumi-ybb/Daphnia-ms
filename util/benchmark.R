# Benchmark code
library(reshape2)
library(magrittr)
library(MASS)
library(readxl)
library(glmmTMB)
library(pomp)
library(panelPomp)
library(tidyverse)
library(loo)
library(deSolve)

Mesocosm_data = read_excel("./Mesocosmdata.xls",3)

sed = 0923
set.seed(sed)

dentNoPara = Mesocosm_data[91:170, ]
dentNoPara = subset(dentNoPara, select = c(rep, day, dent.adult,dent.inf,lum.adult,lum.adult.inf))
dentNoPara = dentNoPara[80: 1, ]
#Convert sampling date to natural date. Samples are collected every 5 days after the #first 7 days.
dentNoPara$day = (dentNoPara$day - 1) * 5 + 7
data = list()
trails = c("K","L","M","N","O","P","Q","S")
for (i in 1: length(trails)){
  data[[i]]=subset(dentNoPara, select = c("day", "dent.adult","dent.inf","lum.adult","lum.adult.inf"),
                   dentNoPara$rep == trails[i])
}


##########################Negative Binomial Regression##########################

fit_dent_adult_linear <- glmmTMB(dent.adult ~ day + (1|rep),
                                 family = nbinom2,
                                 data = dentNoPara)

fit_dent_inf_linear <- glmmTMB(dent.inf ~ day + (1|rep),
                               family = nbinom2,
                               data = dentNoPara)

fit_lum_adult_linear <- glmmTMB(lum.adult ~ day + (1|rep),
                                family = nbinom2,
                                data = dentNoPara)

fit_lum_adult_inf_linear <- glmmTMB(lum.adult.inf ~ day + (1|rep),
                                    family = nbinom2,
                                    data = dentNoPara)

total_ll_linear = logLik(fit_lum_adult_inf_linear) + logLik(fit_lum_adult_linear) + 
  logLik(fit_dent_inf_linear) + logLik(fit_dent_adult_linear) 

AIC_linear <- AIC(fit_dent_adult_linear) + 
  AIC(fit_dent_inf_linear) + 
  AIC(fit_lum_adult_linear) + 
  AIC(fit_lum_adult_inf_linear)

## loglik = -1009.25, AIC = 2050.499
AIC_linear
total_ll_linear


fit_dent_adult_quad <- glmmTMB(dent.adult ~ day + I(day^2) + (1|rep),
                               family = nbinom2,
                               data = dentNoPara)

fit_dent_inf_quad <- glmmTMB(dent.inf ~ day + I(day^2) + (1|rep),
                             family = nbinom2,
                             data = dentNoPara)

fit_lum_adult_quad <- glmmTMB(lum.adult ~ day + I(day^2) + (1|rep),
                              family = nbinom2,
                              data = dentNoPara)

fit_lum_adult_inf_quad <- glmmTMB(lum.adult.inf ~ day + I(day^2) + (1|rep),
                                  family = nbinom2,
                                  data = dentNoPara)

total_ll_quad = logLik(fit_lum_adult_inf_quad) + logLik(fit_lum_adult_quad) + 
  logLik(fit_dent_inf_quad) + logLik(fit_dent_adult_quad) 

AIC_quad <- AIC(fit_dent_adult_quad) + 
  AIC(fit_dent_inf_quad) + 
  AIC(fit_lum_adult_quad) + 
  AIC(fit_lum_adult_inf_quad)

##loglik = -943.7927 AIC = 1927.585
AIC_quad
total_ll_quad


fit_dent_adult_cubic <- glmmTMB(dent.adult ~ day + I(day^2) + I(day^3) + (1|rep),
                               family = nbinom2,
                               data = dentNoPara)

fit_dent_inf_cubic <- glmmTMB(dent.inf ~ day + I(day^2) + I(day^3) + (1|rep),
                             family = nbinom2,
                             data = dentNoPara)

fit_lum_adult_cubic <- glmmTMB(lum.adult ~ day + I(day^2) + I(day^3) + (1|rep),
                              family = nbinom2,
                              data = dentNoPara)

fit_lum_adult_inf_cubic <- glmmTMB(lum.adult.inf ~ day + I(day^2) + I(day^3) + (1|rep),
                                  family = nbinom2,
                                  data = dentNoPara)

total_ll_cubic = logLik(fit_lum_adult_inf_cubic) + logLik(fit_lum_adult_cubic) + 
  logLik(fit_dent_inf_cubic) + logLik(fit_dent_adult_cubic) 

AIC_cubic <- AIC(fit_dent_adult_cubic) + 
  AIC(fit_dent_inf_cubic) + 
  AIC(fit_lum_adult_cubic) + 
  AIC(fit_lum_adult_inf_cubic)

##loglik = -943.7927 AIC = 1927.585
AIC_cubic
total_ll_cubic

##########################Searle model##########################
modelS = function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    dSn = rn*(Sn + xi*In)*(1 - (Sn + In + alpha_ni*(Si + Ii))/Kn) - pn*fSn*Sn*P - delta*Sn
    
    dIn = pn*fSn*Sn*P - dn*In - delta*In
    
    dSi = ri*(Si + xi*Ii)*(1 - (Si + Ii + alpha_in*(Sn + In))/Ki) - pi*fSi*Si*P - delta*Si
    
    dIi = pi*fSi*Si*P - di*Ii - delta*Ii
    
    dP  = beta_n*dn*In + beta_i*di*Ii - fSn*Sn*Sn*P - fIn*In*P - fSi*Si*P - fIi*Ii*P - dp*P - delta*P
    
    list(c(dSn, dIn, dSi, dIi, dP))
  })
}
# Values taken from Table... in ....
parameters = list(
  rn = 0.206, xi = 0.75, alpha_ni = 2.63, Kn = 97.5, 
  pn = 1.45e-5, fSn = 0.0348, delta = 0.013, dn = 0.05,
  ri = 0.246, alpha_in = -0.286, Ki = 12.8, pi = 4.87e-5,
  fSi = 0.0361, di = 0.05, beta_n = 120000, beta_i = 124000,
  fIn = 0.0186, fIi = 0.0171, dp = 0.5
)


state = c(Sn = 2.333, In = 0, Si = 0.667, Ii = 0, P = 0)

events = data.frame(
  var = "P",     
  time = 4,       
  value = 25e3,     
  method = "add"  
)

times = seq(0, 52, by = 1)

out = ode(y = state, 
           times = times, 
           func = modelS, 
           parms = parameters, 
           events = list(data = events))

prediction = as.data.frame(out)

names(prediction)[2:5] = c("dent.adult.pred", "dent.inf.pred", "lum.adult.pred", "lum.adult.inf.pred")


use_poisson = FALSE

if (!use_poisson) {
  size_param = 0.02
}

total_ll_all = 0 

for (i in seq_along(trails)) {
  rep_id = trails[i]
  observed = data[[i]]
  
  predictions = prediction[prediction$time %in% c(7,12,17,22,27,32,37,42,47,52),]
  if (!"rep" %in% names(observed)) {
    observed$rep = rep_id
  }
  predictions$rep = rep_id
  merged_data = merge(observed, predictions, by.x = c("rep","day"), by.y = c("rep","time"))
  
  if (!use_poisson) {
    ll_dent_adult = sum(dnbinom(x = merged_data$dent.adult,
                                 size = size_param,
                                 mu = merged_data$dent.adult.pred,
                                 log = TRUE))
    ll_dent_inf = sum(dnbinom(x = merged_data$dent.inf,
                               size = size_param,
                               mu = merged_data$dent.inf.pred,
                               log = TRUE))
    ll_lum_adult = sum(dnbinom(x = merged_data$lum.adult,
                                size = size_param,
                                mu = merged_data$lum.adult.pred,
                                log = TRUE))
    ll_lum_adult_inf = sum(dnbinom(x = merged_data$lum.adult.inf,
                                    size = size_param,
                                    mu = merged_data$lum.adult.inf.pred,
                                    log = TRUE))
  } else {
    ll_dent_adult = sum(dpois(x = merged_data$dent.adult,
                               lambda = merged_data$dent.adult.pred,
                               log = TRUE))
    ll_dent_inf = sum(dpois(x = merged_data$dent.inf,
                             lambda = merged_data$dent.inf.pred,
                             log = TRUE))
    ll_lum_adult = sum(dpois(x = merged_data$lum.adult,
                              lambda = merged_data$lum.adult.pred,
                              log = TRUE))
    ll_lum_adult_inf = sum(dpois(x = merged_data$lum.adult.inf,
                                  lambda = merged_data$lum.adult.inf.pred,
                                  log = TRUE))
  }
  
  replicate_ll = ll_dent_adult + ll_dent_inf + ll_lum_adult + ll_lum_adult_inf
  cat("Replicate:", rep_id, "Log-likelihood:", replicate_ll, "\n")
  
  total_ll_all = total_ll_all + replicate_ll
}

# -1596.406 
#  AIC = 3230.813
AIC = -2 * total_ll_all + 2 * length(parameters)
cat("Total log-likelihood across all replicates:", total_ll_all, "\n")



benchmark_model_summary = data.frame(
  NB_linear = numeric(),
  NB_quad = numeric(),
  NB_cubic = numeric(),
  ModelS = numeric()
)
benchmark_model_summary[1,] = c(total_ll_linear,total_ll_quad,total_ll_cubic,total_ll_all)
benchmark_model_summary[2,] = c(AIC_linear,AIC_quad,AIC_cubic,AIC)
rownames(benchmark_model_summary) = c("ll","AIC")
benchmark_model_summary


load('data/Benchmark_model_result/benchmark_no_juvenile.RData')
no_juv_ll = pf.loglik.of.mif.estimate
AIC = -2 * no_juv + 2 * length(mif.estimate)
benchmark_model_summary$Model_no_juv = c(no_juv_ll,AIC)

save(benchmark_model_summary,file = './daphnia-article/data/benchmark_model_summary.rda')











