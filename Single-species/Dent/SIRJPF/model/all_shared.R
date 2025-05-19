library(reshape2)
library(magrittr)
library(foreach)
library(readxl)
library(doParallel)
registerDoParallel(cores=36)
library(pomp)
library(panelPomp)
library(tidyverse)

# Mesocosm_data = read_excel("/Users/ybb/Desktop/Research//Daphnia/Mesocosmdata.xls")
Mesocosm_data = read_excel("/home/ybb/D_P/Mesocosmdata.xlsx")

sed = 0923
set.seed(0923)

name_str = "all_shared"
run_level <- 3

dentNoPara <- Mesocosm_data[101:180, ]
dentNoPara <- subset(dentNoPara, select = c(rep, day, dent.adult,dent.inf ))
dentNoPara <- dentNoPara[80: 1, ]
#Convert sampling date to natural date. Samples are collected every 5 days after the #first 7 days.
dentNoPara$day = (dentNoPara$day - 1) * 5 + 7
data = list()
trails = c("K","L","M","N","O","P","Q","S")
for (i in 1: length(trails)){
  data[[i]]=subset(dentNoPara, select = c("day", "dent.adult","dent.inf"),
                   dentNoPara$rep == trails[i])
}


#Notice that the unit of S is individuals/L, unit of F is 1e+5 cells/L, unit of P is 1000 spores/L, unit for I is individuals/L

dyn_rpro <- Csnippet("
                      double Sn_term, In_term, F_term, P_term ,Jn_term;
                      double noiSn, noiIn,noiF, noiP,noiJn;
                      double delta = 0.013; //fraction of volume replaced day-1

                      noiSn = rnorm(0, sigSn * sqrt(dt));
                      noiIn = rnorm(0, sigIn * sqrt(dt));
                      noiJn = rnorm(0, sigJn * sqrt(dt));
                      noiF = rnorm(0, sigF * sqrt(dt));
                      noiP = rnorm(0, sigP * sqrt(dt));


                      //------------Sn-------------
                      Sn_term = 0.1 * Jn * dt - theta_Sn * Sn * dt -  probn * f_Sn * Sn * P * dt - delta * Sn * dt + Sn * noiSn;
                        
                      //-----------Jn-------------
                      
                      Jn_term = rn * f_Sn * F * Sn  * dt  -  0.1 * Jn * dt - theta_Jn * Jn * dt - delta * Jn * dt + Jn * noiJn;
                      
                      //------------In--------------

                      In_term = probn * f_Sn * Sn * P * dt - theta_In * In *dt - delta * In * dt + In * noiIn;

                      //-----------F---------------

                      F_term =  F * noiF - f_Sn * F * (Sn + xi * In + 1 * Jn) * dt - delta * F * dt + 0.37 * dt;


                      //----------P---------------

                      P_term = 30 * theta_In * In * dt - f_Sn * (Sn + xi * In) * P * dt - theta_P * P * dt - delta * P * dt + P * noiP;


                      //T_D and T_I are the S and I Daphnia sample out
                      F += F_term;
                      Sn += Sn_term;
                      In += In_term;
                      Jn += Jn_term;
                      P += P_term;
                      
                      //Trace time
                      
                     if (t - 4.0 < 0.001 && t - 4.0 > -0.001){
                      //Initial statement
                        P += 25;
                      }
                      
                      
                      if (Sn < 0.0 || Sn > 1e5) {
                        Sn = 0.0;
                        error_count += 1;
                      }
                      if (F < 0.0 || F > 1e20) {
                        F = 0.0;
                        error_count += 1000;
                      }
                      if (In < 0.0 || In > 1e5) {
                        In = 0.0;
                        error_count += 0.001;
                      }
                      if (Jn < 0.0 || Jn > 1e5) {
                        Jn = 0.0;
                        error_count += 0.001;
                      }
                      if ((P < 0.0 || P > 1e20)&& t > 3.9) {
                        P = 0.0;
                        error_count += 0.000001;
                      }
                      
                      T_Sn = fabs(Sn);
                      T_In = fabs(In);
                      

                      ")

# Initial state. Assume t0 = day 4
dyn_init = Csnippet("
                     Sn = 3;
                     F = 16.667;
                     Jn = 0;
                     
                     T_Sn = 0.0;
                     T_In = 0.0;
                     In = 0.0;
                     
                     error_count = 0.0;
                     
                     //add 25 P at day 4
                     P = 0;
                     ")



dmeas = Csnippet("
                 
                 if (error_count > 0.0) {
                   lik = -150;
                  }
                 else{
                    if(give_log){
                    lik  = dnbinom_mu(dentadult,k_Sn,T_Sn,give_log) +  dnbinom_mu(dentinf,k_In,T_In,give_log);
                    }
                    else{
                    lik  = dnbinom_mu(dentadult,k_Sn,T_Sn,give_log) *  dnbinom_mu(dentinf,k_In,T_In,give_log);
                    }
                      }
                 ")

rmeas = Csnippet("
                 //double delta = 0.013;
                 dentadult = rnbinom_mu(k_Sn,T_Sn);
                 dentinf = rnbinom_mu(k_In,T_In);
                 ")

pt <- parameter_trans(
  #Without death part
  log = c( "sigSn", "sigIn", "sigF","sigP","f_Sn",
           "rn","k_In","k_Sn","sigJn","probn",
           "theta_Sn","theta_In","theta_P","theta_Jn","xi")
)

pomplist=list()

parameters =  c( "sigSn" = 0.1, "sigIn"= 0.1, "sigF"= 0.1,"sigP"= 0.1,"f_Sn"= 0.01,
                 "rn"= 10,"k_In"= 0.1,"k_Sn"= 0.1,"sigJn"= 0.1,"probn"= 0.1,
                 "theta_Sn"= 0.1,"theta_In"= 0.1,"theta_P"= 0.1,"theta_Jn"= 0.1,"xi"= 1)


for (i in 1:8){
  colnames(data[[i]]) <- c('day', 'dentadult', 'dentinf')
  pomp(data = data[[i]],
       times = "day",
       t0=1,
       rprocess=euler(dyn_rpro,delta.t=1/4),
       rinit = dyn_init,
       dmeasure = dmeas,
       rmeasure = rmeas,
       partrans = pt,
       obsnames = c("dentadult", "dentinf"),
       accumvars = c("error_count"),
       paramnames = c( "sigSn","sigIn","sigF","sigP","f_Sn","rn","k_In","k_Sn","sigJn","probn",
                       "theta_Sn","theta_In","theta_P","theta_Jn","xi"),
       statenames = c("Sn",  "In","Jn" ,"error_count", "F", "T_Sn","T_In","P")
  ) -> pomplist[[i]]
  coef(pomplist[[i]])=parameters
}
names(pomplist)=paste("u", 1:8,sep = "")



shared_parameter = c(
  rn = 69.34499, f_Sn = 0.0007281269, probn = 0.3193476, xi = 14.43326, 
  theta_Sn = 0.01660319, theta_In = 0.542057, theta_P = 7.917627e-05, theta_Jn = 0.002125109, 
  sigSn = 0, sigIn = 0.5887923, sigJn = 0.3089306, sigF = 0.07980189, 
  sigP = 0.4788064, k_In = 1.117379, k_Sn = 12.36433
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
        rw.sd = rw_sd(xi=dent_rw.sd,
                      sigSn=0,
                      sigIn=dent_rw.sd,
                      sigF=dent_rw.sd,
                      theta_Sn=dent_rw.sd,
                      theta_In=dent_rw.sd,
                      theta_P=dent_rw.sd,
                      k_Sn = dent_rw.sd,
                      k_In = dent_rw.sd,
                      f_Sn=dent_rw.sd,
                      rn=dent_rw.sd,
                      probn=dent_rw.sd,
                      sigP = dent_rw.sd,
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
        rw.sd = rw_sd(xi=dent_rw.sd,
                      sigSn=0,
                      sigIn=dent_rw.sd,
                      sigF=dent_rw.sd,
                      theta_Sn=dent_rw.sd,
                      theta_In=dent_rw.sd,
                      theta_P=dent_rw.sd,
                      k_Sn = dent_rw.sd,
                      k_In = dent_rw.sd,
                      f_Sn=dent_rw.sd,
                      rn=dent_rw.sd,
                      probn=dent_rw.sd,
                      sigP = dent_rw.sd,
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
        Nmif = 400,
        shared.start = share_para_temp,
        rw.sd = rw_sd(xi=dent_rw.sd,
                      sigSn=0,
                      sigIn=dent_rw.sd,
                      sigF=dent_rw.sd,
                      theta_Sn=dent_rw.sd,
                      theta_In=dent_rw.sd,
                      theta_P=dent_rw.sd,
                      k_Sn = dent_rw.sd,
                      k_In = dent_rw.sd,
                      f_Sn=dent_rw.sd,
                      rn=dent_rw.sd,
                      probn=dent_rw.sd,
                      sigP = dent_rw.sd,
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

if (run_level %in% c(2,3)){
  save(mf,lls,best,mif.estimate,pf.loglik.of.mif.estimate,
       s.e.of.pf.loglik.of.mif.estimate,
       file = paste0(name_str,".RData"))
}