library(reshape2)
library(magrittr)
library(foreach)
library(readxl)
library(doParallel)
registerDoParallel(cores=36)
library(pomp)
library(panelPomp)
library(tidyverse)

# Mesocosm_data = read_excel("~/Desktop/Research/D_P/Mesocosmdata.xlsx",2)
Mesocosm_data = read_excel("/home/ybb/D_P/Mesocosmdata.xlsx",2)

sed = 0923
set.seed(0923)

name_str = "all_shared"
run_level <- 3

dentNoPara <- Mesocosm_data[1:100, ]
dentNoPara <- subset(dentNoPara, select = c(rep, day, lum.adult))
dentNoPara <- dentNoPara[100: 1, ]
#Convert sampling date to natural date. Samples are collected every 5 days after the #first 7 days.
dentNoPara$day = (dentNoPara$day - 1) * 5 + 7
data = list()
trails = c("A","B","C","D","E","F","G","H","I","J")
for (i in 1: length(trails)){
  data[[i]]=subset(dentNoPara, select = c("day", "lum.adult"),
                   dentNoPara$rep == trails[i])
}

dyn_rpro <- Csnippet("
                      double Si_term, F_term , Ji_term;
                      double noiSi  ,noiF,noiJi;
                      double delta = 0.013; //fraction of volume replaced day-1

                      noiSi = rnorm(0, sigSi * sqrt(dt));
                      noiJi = rnorm(0, sigJi * sqrt(dt));
                      noiF = rnorm(0, sigF * sqrt(dt));


                      //------------Si-------------
                      Si_term = 0.1 * Ji * dt - theta_Si * Si * dt -  delta * Si * dt + Si * noiSi;
                        
                      //-----------Ji-------------
                      
                      Ji_term = ri * f_Si * F * Si  * dt  -  0.1 * Ji * dt - theta_Ji * Ji * dt - delta * Ji * dt + Ji * noiJi;
                      
                      //-----------F---------------

                      F_term =  F * noiF - f_Si * F * (Si + 1 * Ji) * dt  - delta * F * dt + 0.37 * dt;


                      //T_D and T_I are the S and I Daphnia sample out
                      F += F_term;
                      Si += Si_term;
                      Ji += Ji_term;
                  
                      
                      
                      if (Si < 0.0 || Si > 1e5) {
                        Si = 0.0;
                        error_count += 1;
                      }
                      if (F < 0.0 || F > 1e20) {
                        F = 0.0;
                        error_count += 1000;
                      }
                      if (Ji < 0.0 || Ji > 1e5) {
                        Ji = 0.0;
                        error_count += 0.001;
                      }
  
                      
                      T_Si = fabs(Si);
                      

                      ")

# Initial state. Assume t0 = day 4
dyn_init = Csnippet("
                     Si = 3; 
                     F = 16.667;
                     Ji = 0;
                     
                     T_Si = 0.0;
                     
                     error_count = 0.0;
                     
                     ")



dmeas = Csnippet("
                 
                 if (error_count > 0.0) {
                   lik = -150;
                  }
                 else{
                    if(give_log){
                    lik  = dnbinom_mu(lumadult,k_Si,T_Si,give_log);
                    }
                    else{
                    lik  = dnbinom_mu(lumadult,k_Si,T_Si,give_log);
                    }
                      }
                 ")

rmeas = Csnippet("
                 //double delta = 0.013;
                 lumadult = rnbinom_mu(k_Si,T_Si);
                 ")

pt <- parameter_trans(
  log = c( "sigSi","sigF","f_Si",
           "ri","k_Si","sigJi",
           "theta_Si","theta_Ji")
)

pomplist=list()

parameters=  c( "sigSi" = 0, "sigF"= 0.1,"f_Si"= 0.1,
                "ri"= 10,"k_Si"= 0.1,"sigJi"= 0.1,
                "theta_Si"= 0.1,"theta_Ji"= 0.1)


for (i in 1:10){
  colnames(data[[i]]) <- c('day', 'lumadult')
  pomp(data = data[[i]],
       times = "day",
       t0=1,
       rprocess=euler(dyn_rpro,delta.t=1/4),
       rinit = dyn_init,
       dmeasure = dmeas,
       rmeasure = rmeas,
       partrans = pt,
       obsnames = c("lumadult"),
       accumvars = c("error_count"),
       paramnames = c( "sigSi","sigF","f_Si",
                       "ri","k_Si","sigJi",
                       "theta_Si","theta_Ji"),
       statenames = c("Si","Ji" ,"error_count", "F", "T_Si")
  ) -> pomplist[[i]]
  coef(pomplist[[i]])=parameters
}
names(pomplist)=paste("u", 1:10,sep = "")



shared_parameter = c(
  ri = 708.6439, f_Si = 0.0002316239, theta_Si = 0.5255919, theta_Ji = 0.07714973,
  sigSi = 0, sigJi = 0.3588666, sigF = 0.0458815, k_Si = 2.288263
)

panelfood = panelPomp(pomplist, shared=shared_parameter)


algorithmic.params <- list(
  Np =     c(50, 500, 1e3),
  Np_rep = c( 2,  10,  20),
  Mp =     c(50, 500, 1e3),
  Nmif =   c( 2,  320, 250)
)


dent_rw.sd= 0.05
parameter_candidates = list(shared_parameter)
names(parameter_candidates)=c("shared")
U = length(panelfood)
{
  foreach(
    i = 1:(10*getDoParWorkers()),
    .packages = c("pomp", "panelPomp"),
    .inorder = FALSE,
    .options.multicore = list(set.seed = TRUE)
  ) %dopar%
    {
      guessed.parameter.values <- parameter_candidates
      mif2(
        panelfood,
        Nmif = 200,
        # shared.start = shared_start,
        shared.start = guessed.parameter.values$shared,
        # specific.start = unit_start,
        rw.sd = rw_sd(
          sigSi=0,
          sigF=dent_rw.sd,
          theta_Si=dent_rw.sd,
          k_Si = dent_rw.sd,
          f_Si=dent_rw.sd,
          ri=dent_rw.sd,
          sigJi = dent_rw.sd,
          theta_Ji = dent_rw.sd),
        cooling.type = "geometric",
        cooling.fraction.50 = 0.7,
        Np = algorithmic.params$Mp[run_level]
      ) -> m1
      
      ll <- replicate(n = algorithmic.params$Np_rep[run_level],
                      unitlogLik(pfilter(m1,
                                         Np = algorithmic.params$Np[run_level])))
      
      list(mif = m1,
           ll = panel_logmeanexp(x = ll,
                                 MARGIN = 1,
                                 se = TRUE))
    }
} -> mf1

log_list = c()
for ( i in 1:length(mf1)){
  log_list = c(log_list,mf1[[i]]$ll[1])
}
select = order(log_list,decreasing = TRUE)[1:ceiling(length(mf1)/4)]
shared_dataframe = data.frame(t(mf1[[select[1]]]$mif@shared))
for (i in 1:length(select)){
  shared_dataframe[i,] = t(mf1[[select[i]]]$mif@shared)
}
shared_dataframe <- shared_dataframe[rep(1:nrow(shared_dataframe), each = 4), ]
dent_rw.sd = 0.04
#round 2
{
  foreach(
    i = 1:(10*getDoParWorkers()),
    .packages = c("pomp", "panelPomp"),
    .inorder = FALSE,
    .options.multicore = list(set.seed = TRUE)
  ) %dopar%
    {
      share_para_temp = as.numeric(shared_dataframe[i,])
      names(share_para_temp) = colnames(shared_dataframe)
      
      mif2(
        panelfood,
        Nmif = 200,
        shared.start = share_para_temp,
        rw.sd = rw_sd(
          sigSi=0,
          sigF=dent_rw.sd,
          theta_Si=dent_rw.sd,
          k_Si = dent_rw.sd,
          f_Si=dent_rw.sd,
          ri=dent_rw.sd,
          sigJi = dent_rw.sd,
          theta_Ji = dent_rw.sd),
        cooling.type = "geometric",
        cooling.fraction.50 = 0.7,
        Np = algorithmic.params$Mp[run_level]
      ) -> m1
      
      ll <- replicate(n = algorithmic.params$Np_rep[run_level],
                      unitlogLik(pfilter(m1,
                                         Np = algorithmic.params$Np[run_level])))
      
      list(mif = m1,
           ll = panel_logmeanexp(x = ll,
                                 MARGIN = 1,
                                 se = TRUE))
    }
} -> mf2


log_list = c()
for ( i in 1:length(mf2)){
  log_list = c(log_list,mf2[[i]]$ll[1])
}
select = order(log_list,decreasing = TRUE)[1:ceiling(length(mf2)/4)]
shared_dataframe = data.frame(t(mf2[[select[1]]]$mif@shared))
for (i in 1:length(select)){
  shared_dataframe[i,] = t(mf2[[select[i]]]$mif@shared)
}
shared_dataframe <- shared_dataframe[rep(1:nrow(shared_dataframe), each = 4), ]
dent_rw.sd = 0.03

#round 3
{
  foreach(
    i = 1:(10*getDoParWorkers()),
    .packages = c("pomp", "panelPomp"),
    .inorder = FALSE,
    .options.multicore = list(set.seed = TRUE)
  ) %dopar%
    {
      
      share_para_temp = as.numeric(shared_dataframe[i,])
      names(share_para_temp) = colnames(shared_dataframe)
      
      mif2(
        panelfood,
        Nmif = 400,
        shared.start = share_para_temp,
        rw.sd = rw_sd(
          sigSi=0,
          sigF=dent_rw.sd,
          theta_Si=dent_rw.sd,
          k_Si = dent_rw.sd,
          f_Si=dent_rw.sd,
          ri=dent_rw.sd,
          sigJi = dent_rw.sd,
          theta_Ji = dent_rw.sd),
        cooling.type = "geometric",
        cooling.fraction.50 = 0.7,
        Np = algorithmic.params$Mp[run_level]
      ) -> m1
      
      ll <- replicate(n = algorithmic.params$Np_rep[run_level],
                      unitlogLik(pfilter(m1,
                                         Np = algorithmic.params$Np[run_level])))
      
      list(mif = m1,
           ll = panel_logmeanexp(x = ll,
                                 MARGIN = 1,
                                 se = TRUE))
    }
} -> mf


lls <- matrix(unlist(sapply(mf, getElement, "ll")), nrow = 2)
best <- which.max(lls[1,])
mif.estimate <- coef(mf[[best]]$mif)
pf.loglik.of.mif.estimate <- unname(mf[[best]]$ll[1])
s.e.of.pf.loglik.of.mif.estimate <- unname(mf[[best]]$ll[2])

pf.loglik.of.mif.estimate
s.e.of.pf.loglik.of.mif.estimate

if (run_level == 2){
  save(mf,mf1,mf2,lls,best,mif.estimate,pf.loglik.of.mif.estimate,
       s.e.of.pf.loglik.of.mif.estimate,
       file = paste0(name_str,".RData"))
}

if (run_level == 3){
  save(mf, lls,best,mif.estimate,pf.loglik.of.mif.estimate,
       s.e.of.pf.loglik.of.mif.estimate,
       file = paste0(name_str,".RData"))
}