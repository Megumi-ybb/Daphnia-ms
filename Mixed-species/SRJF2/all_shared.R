library(reshape2)
library(magrittr)
library(foreach)
library(readxl)
library(doParallel)
registerDoParallel(cores=36)
library(pomp)
library(panelPomp)
library(tidyverse)

Mesocosm_data = read_excel("./Mesocosmdata.xlsx",3)

name_str = "all_shared"
run_level <- 3

dentNoPara <- Mesocosm_data[1:90, ]
dentNoPara <- subset(dentNoPara, select = c(rep, day, dent.adult,lum.adult))
dentNoPara <- dentNoPara[90: 1, ]
#Convert sampling date to natural date. Samples are collected every 5 days after the #first 7 days.
dentNoPara$day = (dentNoPara$day - 1) * 5 + 7
data = list()
trails = c("A","C","D","E","F","G","H","I","J")
for (i in 1: length(trails)){
  data[[i]]=subset(dentNoPara, select = c("day", "dent.adult","lum.adult"),
                   dentNoPara$rep == trails[i])
}


dyn_rpro <- Csnippet("
                      double Sn_term, F_term , Si_term,Jn_term,Ji_term;
                      double noiSn, noiSi  ,noiF,noiJn,noiJi;
                      double delta = 0.013; //fraction of volume replaced day-1

                      noiSn = rnorm(0, sigSn * sqrt(dt));
                      noiSi = rnorm(0, sigSi * sqrt(dt));
                      noiJn = rnorm(0, sigJn * sqrt(dt));
                      noiJi = rnorm(0, sigJi * sqrt(dt));
                      noiF = rnorm(0, sigF * sqrt(dt));


                      //------------Sn-------------
                      Sn_term = 0.1 * Jn * dt - theta_Sn * Sn * dt -  delta * Sn * dt + Sn * noiSn;
                        
                      //-----------Jn-------------
                      
                      Jn_term = rn * f_Sn * F * Sn  * dt  -  0.1 * Jn * dt - theta_Jn * Jn * dt - delta * Jn * dt + Jn * noiJn;
                      

                      //------------Si-------------
                      Si_term = 0.1 * Ji * dt - theta_Si * Si * dt - delta * Si * dt + Si * noiSi;

                      //-----------Ji-------------
                      
                      Ji_term = ri * f_Si *F * Si  * dt - 0.1 * Ji * dt - theta_Ji * Ji * dt - delta * Ji * dt + Ji * noiJi;
                      
                      //-----------F---------------

                      F_term =  F * noiF - f_Sn * F * (Sn + 1 * Jn) * dt - f_Si * F * (Si + 1 * Ji) * dt - delta * F * dt + 0.37 * dt;


                      //T_D and T_I are the S and I Daphnia sample out
                      F += F_term;
                      Sn += Sn_term;
                      Si += Si_term;
                      Ji += Ji_term;
                      Jn += Jn_term;
                  
                      
                      
                      if (Sn < 0.0 || Sn > 1e5) {
                        Sn = 0.0;
                        error_count += 1;
                      }
                      if (Si < 0.0 || Si > 1e5) {
                        Si = 0.0;
                        error_count += 1000000;
                      }
                      if (F < 0.0 || F > 1e20) {
                        F = 0.0;
                        error_count += 1000;
                      }
                      if (Jn < 0.0 || Jn > 1e5) {
                        Jn = 0.0;
                        error_count += 0.001;
                      }
                      if (Ji < 0.0 || Ji > 1e5) {
                        Ji = 0.0;
                        error_count += 0.000000001;
                      }   
                      
                      T_Sn = fabs(Sn);
                      T_Si = fabs(Si);
                      

                      ")

# Initial state. Assume t0 = day 4
dyn_init = Csnippet("
                     Sn = 2.333; //2.3333 = 35/15L
                     Si = 0.667; //0.667 = 10/15L
                     F = 16.667;
                     Jn = 0;
                     Ji = 0;
                     
                     T_Sn = 0.0;
                     T_Si = 0.0;
                     
                     error_count = 0.0;
                     
                     ")



dmeas = Csnippet("
                 
                 if (error_count > 0.0) {
                   lik = -150;
                  }
                 else{
                    if(give_log){
                    lik  = dnbinom_mu(lumadult,k_Si,T_Si,give_log)  + dnbinom_mu(dentadult,k_Sn,T_Sn,give_log);
                    }
                    else{
                    lik  = dnbinom_mu(lumadult,k_Si,T_Si,give_log) * dnbinom_mu(dentadult,k_Sn,T_Sn,give_log);
                    }
                      }
                 ")

rmeas = Csnippet("
                 //double delta = 0.013;
                 dentadult = rnbinom_mu(k_Sn,T_Sn);
                 lumadult = rnbinom_mu(k_Si,T_Si);
                 ")

pt <- parameter_trans(
  log = c( "sigSn", "sigSi","sigF","f_Sn","f_Si",
           "rn","ri","k_Sn","k_Si","sigJi","sigJn",
           "theta_Sn","theta_Si","theta_Jn","theta_Ji")
)

pomplist=list()

parameters=  c( "sigSn" = 0.1, "sigSi"= 0.1,"sigF"= 0.1,"f_Sn"= 0.1,"f_Si"= 0.1,
                "rn"= 10,"ri"= 10,"k_Sn"= 0.1,"k_Si"= 0.1,"sigJi"= 0.1,"sigJn"= 0.1,
                "theta_Sn"= 0.1,"theta_Si"= 0.1,"theta_Jn"= 0.1,"theta_Ji"= 0.1)


for (i in 1:9){
  colnames(data[[i]]) <- c('day', 'dentadult', 'lumadult')
  pomp(data = data[[i]],
       times = "day",
       t0=1,
       rprocess=euler(dyn_rpro,delta.t=1/4),
       rinit = dyn_init,
       dmeasure = dmeas,
       rmeasure = rmeas,
       partrans = pt,
       obsnames = c("dentadult","lumadult"),
       accumvars = c("error_count"),
       paramnames = c( "sigSn", "sigSi","sigF","f_Sn","f_Si",
                       "rn","ri","k_Sn","k_Si","sigJi","sigJn",
                       "theta_Sn","theta_Si","theta_Jn","theta_Ji"),
       statenames = c("Sn", "Si","Jn","Ji" ,"error_count", "F", "T_Sn","T_Si")
  ) -> pomplist[[i]]
  coef(pomplist[[i]])=parameters
}
names(pomplist)=paste("u", 1:9,sep = "")



shared_parameter = c(
  ri = 100.72, rn = 100.94, f_Si = 5.178821e-05, f_Sn = 9.394652e-05, 
  theta_Sn = 0.5026099, theta_Si = 0.6418282, theta_Ji = 0.8648736, theta_Jn = 0.9614655, 
  sigSn = 0, sigSi = 0, sigJi = 0.12, sigJn = 0.1879358, 
  sigF = 0.05304384, k_Si = 2.665681, k_Sn = 14.89269
)


panelfood = panelPomp(pomplist, shared=shared_parameter)


algorithmic.params <- list(
  Np =     c(50, 500, 1e3),
  Np_rep = c( 2,  10,  20),
  Mp =     c(50, 500, 1e3),
  Nmif =   c( 2,  320, 250)
)


dent_rw.sd= 0.05
# parameter_candidates = list(shared_parameter, specific_mat)
parameter_candidates = list(shared_parameter)
# names(parameter_candidates)=c("shared","specific")
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
          sigSn=0,
          sigSi=0,
          sigF=dent_rw.sd,
          theta_Sn=dent_rw.sd,
          theta_Si=dent_rw.sd,
          k_Sn = dent_rw.sd,
          k_Si = dent_rw.sd,
          f_Sn=dent_rw.sd,
          f_Si=dent_rw.sd,
          rn=dent_rw.sd,
          ri=dent_rw.sd,
          sigJi = dent_rw.sd,
          sigJn = dent_rw.sd,
          theta_Jn = dent_rw.sd,
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
          sigSn=0,
          sigSi=0,
          sigF=dent_rw.sd,
          theta_Sn=dent_rw.sd,
          theta_Si=dent_rw.sd,
          k_Sn = dent_rw.sd,
          k_Si = dent_rw.sd,
          f_Sn=dent_rw.sd,
          f_Si=dent_rw.sd,
          rn=dent_rw.sd,
          ri=dent_rw.sd,
          sigJi = dent_rw.sd,
          sigJn = dent_rw.sd,
          theta_Jn = dent_rw.sd,
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
dent_rw.sd = 0.04

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
        Nmif = 300,
        shared.start = share_para_temp,
        rw.sd = rw_sd(
          sigSn=0,
          sigSi=0,
          sigF=dent_rw.sd,
          theta_Sn=dent_rw.sd,
          theta_Si=dent_rw.sd,
          k_Sn = dent_rw.sd,
          k_Si = dent_rw.sd,
          f_Sn=dent_rw.sd,
          f_Si=dent_rw.sd,
          rn=dent_rw.sd,
          ri=dent_rw.sd,
          sigJi = dent_rw.sd,
          sigJn = dent_rw.sd,
          theta_Jn = dent_rw.sd,
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
  save(lls,best,mif.estimate,pf.loglik.of.mif.estimate,
       s.e.of.pf.loglik.of.mif.estimate,
       file = paste0(name_str,".RData"))
}