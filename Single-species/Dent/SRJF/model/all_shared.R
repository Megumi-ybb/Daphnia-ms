library(reshape2)
library(magrittr)
library(foreach)
library(readxl)
library(doParallel)
registerDoParallel(cores=36)
library(pomp)
library(panelPomp)
library(tidyverse)

# Mesocosm_data = read_excel("~/Desktop/Research/D_P/Mesocosmdata.xlsx")
Mesocosm_data = read_excel("/home/ybb/D_P/Mesocosmdata.xlsx")

sed = 0923
set.seed(0923)

name_str = "all_shared"
run_level <- 3

dentNoPara <- Mesocosm_data[1:100, ]
dentNoPara <- subset(dentNoPara, select = c(rep, day, dent.adult))
dentNoPara <- dentNoPara[100: 1, ]
#Convert sampling date to natural date. Samples are collected every 5 days after the #first 7 days.
dentNoPara$day = (dentNoPara$day - 1) * 5 + 7
data = list()
trails = c("A","B","C","D","E","F","G","H","I","J")
for (i in 1: length(trails)){
  data[[i]]=subset(dentNoPara, select = c("day", "dent.adult"),
                   dentNoPara$rep == trails[i])
}


dyn_rpro <- Csnippet("
                      double Sn_term, F_term , Jn_term;
                      double noiSn  ,noiF,noiJn;
                      double delta = 0.013; //fraction of volume replaced day-1

                      noiSn = rnorm(0, sigSn * sqrt(dt));
                      noiJn = rnorm(0, sigJn * sqrt(dt));
                      noiF = rnorm(0, sigF * sqrt(dt));


                      //------------Sn-------------
                      Sn_term = 0.1 * Jn * dt - theta_Sn * Sn * dt -  delta * Sn * dt + Sn * noiSn;
                        
                      //-----------Jn-------------
                      
                      Jn_term = rn * f_Sn * F * Sn  * dt  -  0.1 * Jn * dt - theta_Jn * Jn * dt - delta * Jn * dt + Jn * noiJn;
                      
                      //-----------F---------------

                      F_term =  F * noiF - f_Sn * F * (Sn + 1 * Jn) * dt  - delta * F * dt + 0.37 * dt;


                      //T_D and T_I are the S and I Daphnia sample out
                      F += F_term;
                      Sn += Sn_term;
                      Jn += Jn_term;
                  
                      
                      
                      if (Sn < 0.0 || Sn > 1e5) {
                        Sn = 0.0;
                        error_count += 1;
                      }
                      if (F < 0.0 || F > 1e20) {
                        F = 0.0;
                        error_count += 1000;
                      }
                      if (Jn < 0.0 || Jn > 1e5) {
                        Jn = 0.0;
                        error_count += 0.001;
                      }
  
                      
                      T_Sn = fabs(Sn);
                      

                      ")

# Initial state. Assume t0 = day 4
dyn_init = Csnippet("
                     Sn = 3; 
                     F = 16.667;
                     Jn = 0;
                     
                     T_Sn = 0.0;
                     
                     error_count = 0.0;
                     
                     ")



dmeas = Csnippet("
                 
                 if (error_count > 0.0) {
                   lik = -150;
                  }
                 else{
                    if(give_log){
                    lik  = dnbinom_mu(dentadult,k_Sn,T_Sn,give_log);
                    }
                    else{
                    lik  = dnbinom_mu(dentadult,k_Sn,T_Sn,give_log);
                    }
                      }
                 ")

rmeas = Csnippet("
                 //double delta = 0.013;
                 dentadult = rnbinom_mu(k_Sn,T_Sn);
                 ")

pt <- parameter_trans(
  log = c( "sigSn","sigF","f_Sn",
           "rn","k_Sn","sigJn",
           "theta_Sn","theta_Jn")
)

pomplist=list()

parameters=  c(
  f_Sn     = 4.271306e-04,
  sigF     = 1.459130e-03,
  sigSn    = 0,
  rn       = 2.121818e+02,
  k_Sn     = 8.917721e+01,
  sigJn    = 2.688333e-01,
  theta_Sn = 1.514775e-01,
  theta_Jn = 1.633677e-02
)


for (i in 1:10){
  colnames(data[[i]]) <- c('day', 'dentadult', 'lumadult')
  pomp(data = data[[i]],
       times = "day",
       t0=1,
       rprocess=euler(dyn_rpro,delta.t=1/4),
       rinit = dyn_init,
       dmeasure = dmeas,
       rmeasure = rmeas,
       partrans = pt,
       obsnames = c("dentadult"),
       accumvars = c("error_count"),
       paramnames = c( "sigSn","sigF","f_Sn",
                       "rn","k_Sn","sigJn",
                       "theta_Sn","theta_Jn"),
       statenames = c("Sn","Jn" ,"error_count", "F", "T_Sn")
  ) -> pomplist[[i]]
  coef(pomplist[[i]])=parameters
}
names(pomplist)=paste("u", 1:10,sep = "")



shared_parameter = c(
  f_Sn     = 1.306857e-04,
  sigF     = 8.658726e-07,
  sigSn    = 0.000000e+00,
  rn       = 1.535539e+03,
  k_Sn     = 1.417860e+01,
  sigJn    = 3.018772e-01,
  theta_Sn = 6.353239e-01,
  theta_Jn = 1.376217e-03
)

panelfood = panelPomp(pomplist, shared=shared_parameter)


algorithmic.params <- list(
  Np =     c(50, 600, 1e3),
  Np_rep = c( 2,  10,  20),
  Mp =     c(50, 500, 1e3),
  Nmif =   c( 2,  320, 250)
)


dent_rw.sd= 0.02
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
          sigSn=0,
          sigF=dent_rw.sd,
          theta_Sn=dent_rw.sd,
          k_Sn = dent_rw.sd,
          f_Sn=dent_rw.sd,
          rn=dent_rw.sd,
          sigJn = dent_rw.sd,
          theta_Jn = dent_rw.sd),
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
dent_rw.sd = 0.02
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
          sigF=dent_rw.sd,
          theta_Sn=dent_rw.sd,
          k_Sn = dent_rw.sd,
          f_Sn=dent_rw.sd,
          rn=dent_rw.sd,
          sigJn = dent_rw.sd,
          theta_Jn = dent_rw.sd),
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
dent_rw.sd = 0.02

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
        Nmif = 200,
        shared.start = share_para_temp,
        rw.sd = rw_sd(
          sigSn=0,
          sigF=dent_rw.sd,
          theta_Sn=dent_rw.sd,
          k_Sn = dent_rw.sd,
          f_Sn=dent_rw.sd,
          rn=dent_rw.sd,
          sigJn = dent_rw.sd,
          theta_Jn = dent_rw.sd),
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