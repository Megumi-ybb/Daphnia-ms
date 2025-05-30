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

name_str = "benchmark"
run_level <- 3

dentNoPara <- Mesocosm_data[91:170, ]
dentNoPara <- subset(dentNoPara, select = c(rep, day, dent.adult,dent.inf,lum.adult,lum.adult.inf))
dentNoPara <- dentNoPara[80: 1, ]
#Convert sampling date to natural date. Samples are collected every 5 days after the #first 7 days.
dentNoPara$day = (dentNoPara$day - 1) * 5 + 7
data = list()
trails = c("K","L","M","N","O","P","Q","S")
for (i in 1: length(trails)){
  data[[i]]=subset(dentNoPara, select = c("day", "dent.adult","dent.inf","lum.adult","lum.adult.inf"),
                   dentNoPara$rep == trails[i])
}

dyn_rpro <- Csnippet("
                      double Sn_term, In_term, F_term, P_term , Si_term, Ii_term;
                      double noiSn, noiIn, noiSi , noiIi ,noiF, noiP;
                      double delta = 0.013; //fraction of volume replaced day-1

                      noiSn = rnorm(0, sigSn * sqrt(dt));
                      noiIn = rnorm(0, sigIn * sqrt(dt));
                      noiSi = rnorm(0, sigSi * sqrt(dt));
                      noiIi = rnorm(0, sigIi * sqrt(dt));
                      noiF = rnorm(0, sigF * sqrt(dt));
                      noiP = rnorm(0, sigP * sqrt(dt));


                      //------------Sn-------------
                      Sn_term = rn * f_Sn * F * Sn  * dt - theta_Sn * Sn * dt -  probn * f_Sn * Sn * P * dt - delta * Sn * dt + Sn * noiSn;
                      
                      //------------In--------------

                      In_term = probn * f_Sn * Sn * P * dt - theta_In * In *dt - delta * In * dt + In * noiIn;

                      //------------Si-------------
                      Si_term = ri * f_Si *F * Si  * dt - theta_Si * Si * dt -  probi * f_Si * Si * P * dt - delta * Si * dt + Si * noiSi;
                      
                      //------------Ii--------------

                      Ii_term = probi * f_Si * Si * P * dt - theta_Ii * Ii *dt - delta * Ii * dt + Ii * noiIi;


                      //-----------F---------------

                      F_term =  F * noiF - f_Sn * F * (Sn + xi * In) * dt - f_Si * F * (Si + xi * Ii) * dt - delta * F * dt + 0.37 * dt;


                      //----------P---------------

                      P_term = 30 * theta_In * In * dt + 30 * theta_Ii * Ii * dt - f_Sn * (Sn + xi * In) * P * dt - f_Si * (Si + xi * Ii) * P * dt- theta_P * P * dt - delta * P * dt + P * noiP;


                      //T_D and T_I are the S and I Daphnia sample out
                      F += F_term;
                      Sn += Sn_term;
                      In += In_term;
                      Si += Si_term;
                      Ii += Ii_term;
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
                      if (Si < 0.0 || Si > 1e5) {
                        Si = 0.0;
                        error_count += 1000000;
                      }
                      if (F < 0.0 || F > 1e20) {
                        F = 0.0;
                        error_count += 1000;
                      }
                      if (In < 0.0 || In > 1e5) {
                        In = 0.0;
                        error_count += 0.001;
                      }
                      if (Ii < 0.0 || Ii > 1e5) {
                        Ii = 0.0;
                        error_count += 0.000000001;
                      }
                      if ((P < 0.0 || P > 1e20)&& t > 3.9) {
                        P = 0.0;
                        error_count += 0.000001;
                      }
                      
                      T_Sn = fabs(Sn);
                      T_In = fabs(In);
                      T_Si = fabs(Si);
                      T_Ii = fabs(Ii);
                      

                      ")

# Initial state. Assume t0 = day 4
dyn_init = Csnippet("
                     Sn = 2.333; //2.3333 = 35/15L
                     Si = 0.667; //0.667 = 10/15L
                     F = 16.667;
                     
                     T_Sn = 0.0;
                     T_Si = 0.0;
                     T_In = 0.0;
                     T_Ii = 0.0;
                     In = 0.0;
                     Ii = 0.0;
                     
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
                    lik  = dnbinom_mu(lumadult,k_Si,T_Si,give_log) +  dnbinom_mu(luminf,k_Ii,T_Ii,give_log) + dnbinom_mu(dentadult,k_Sn,T_Sn,give_log) +  dnbinom_mu(dentinf,k_In,T_In,give_log);
                    }
                    else{
                    lik  = dnbinom_mu(lumadult,k_Si,T_Si,give_log) *  dnbinom_mu(luminf,k_Ii,T_Ii,give_log) * dnbinom_mu(dentadult,k_Sn,T_Sn,give_log) *  dnbinom_mu(dentinf,k_In,T_In,give_log);
                    }
                      }
                 ")

rmeas = Csnippet("
                 //double delta = 0.013;
                 dentadult = rnbinom_mu(k_Sn,T_Sn);
                 dentinf = rnbinom_mu(k_In,T_In);
                 lumadult = rnbinom_mu(k_Si,T_Si);
                 luminf = rnbinom_mu(k_Ii,T_Ii);
                 ")

pt <- parameter_trans(
  #Without death part
  log = c( "sigSn", "sigIn", "sigSi", "sigIi", "sigF","sigP","f_Sn","f_Si",
           "rn","ri","k_Ii","k_In","k_Sn","k_Si","probn","probi",
           "theta_Sn","theta_In","theta_Si","theta_Ii","theta_P","xi")
)

pomplist=list()

parameters = c(
  ri = 13076.83, rn = 59.04676, f_Si = 1.838259e-05, f_Sn = 0.001105668, probi = 31.10083,
  probn = 0.2565626, xi = 28.6562, theta_Sn = 0.1479834, theta_Si = 0.0318604,
  theta_Ii = 0.3531879, theta_In = 0.5489315, theta_P = 0.02024991, sigSn = 0, sigSi = 0, sigIn = 0.0003063207,
  sigIi = 0.02208698, sigF = 0.1551729, sigP = 0.238589,
  k_Ii = 1.241092, k_In = 1.005756, k_Si = 4.715556, k_Sn = 4.282648
)


for (i in 1:8){
  colnames(data[[i]]) <- c('day', 'dentadult', 'dentinf','lumadult','luminf')
  pomp(data = data[[i]],
       times = "day",
       t0=1,
       rprocess=euler(dyn_rpro,delta.t=1/4),
       rinit = dyn_init,
       dmeasure = dmeas,
       rmeasure = rmeas,
       partrans = pt,
       obsnames = c("dentadult", "dentinf","lumadult","luminf"),
       accumvars = c("error_count"),
       paramnames = c( "sigSn", "sigIn", "sigSi", "sigIi", "sigF","sigP","f_Sn","f_Si",
                       "rn","ri","k_Ii","k_In","k_Sn","k_Si","probn","probi",
                       "theta_Sn","theta_In","theta_Si","theta_Ii","theta_P","xi"),
       statenames = c("Sn",  "In", "Si","Jn","Ji" , "Ii","error_count", "F", "T_Sn","T_In","T_Si","T_Ii","P")
  ) -> pomplist[[i]]
  coef(pomplist[[i]])=parameters
}
names(pomplist)=paste("u", 1:8,sep = "")



shared_parameter = c(
  ri = 13076.83, rn = 59.04676, f_Si = 1.838259e-05, f_Sn = 0.001105668, probi = 31.10083,
  probn = 0.2565626, xi = 28.6562, theta_Sn = 0.1479834, theta_Si = 0.0318604,
  theta_Ii = 0.3531879, theta_In = 0.5489315, theta_P = 0.02024991, sigSn = 0, sigSi = 0, sigIn = 0.0003063207,
  sigIi = 0.02208698, sigF = 0.1551729, sigP = 0.238589,
  k_Ii = 1.241092, k_In = 1.005756, k_Si = 4.715556, k_Sn = 4.282648
)



panelfood = panelPomp(pomplist, shared=shared_parameter)


algorithmic.params <- list(
  Np =     c(50, 500, 1e4),
  Np_rep = c( 2,  10,  10),
  Mp =     c(50, 500, 1e4),
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
        Nmif = 150,
        shared.start = guessed.parameter.values$shared,
        rw.sd = rw_sd(xi=dent_rw.sd,
                      sigSn=0,
                      sigIn=dent_rw.sd,
                      sigSi=0,
                      sigIi=dent_rw.sd,
                      sigF=dent_rw.sd,
                      theta_Sn=dent_rw.sd,
                      theta_In=dent_rw.sd,
                      theta_Si=dent_rw.sd,
                      theta_P=dent_rw.sd,
                      theta_Ii=dent_rw.sd,
                      k_Sn = dent_rw.sd,
                      k_In = dent_rw.sd,
                      k_Si = dent_rw.sd,
                      k_Ii = dent_rw.sd,
                      f_Sn=dent_rw.sd,
                      f_Si=dent_rw.sd,
                      rn=dent_rw.sd,
                      ri=dent_rw.sd,
                      probn=dent_rw.sd,
                      probi=dent_rw.sd,
                      sigP = dent_rw.sd),
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
        Nmif = 150,
        shared.start = share_para_temp,
        rw.sd = rw_sd(xi=dent_rw.sd,
                      sigSn=0,
                      sigIn=dent_rw.sd,
                      sigSi=0,
                      sigIi=dent_rw.sd,
                      sigF=dent_rw.sd,
                      theta_Sn=dent_rw.sd,
                      theta_In=dent_rw.sd,
                      theta_Si=dent_rw.sd,
                      theta_P=dent_rw.sd,
                      theta_Ii=dent_rw.sd,
                      k_Sn = dent_rw.sd,
                      k_In = dent_rw.sd,
                      k_Si = dent_rw.sd,
                      k_Ii = dent_rw.sd,
                      f_Sn=dent_rw.sd,
                      f_Si=dent_rw.sd,
                      rn=dent_rw.sd,
                      ri=dent_rw.sd,
                      probn=dent_rw.sd,
                      probi=dent_rw.sd,
                      sigP = dent_rw.sd),
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
                      sigSi=0,
                      sigIi=dent_rw.sd,
                      sigF=dent_rw.sd,
                      theta_Sn=dent_rw.sd,
                      theta_In=dent_rw.sd,
                      theta_Si=dent_rw.sd,
                      theta_P=dent_rw.sd,
                      theta_Ii=dent_rw.sd,
                      k_Sn = dent_rw.sd,
                      k_In = dent_rw.sd,
                      k_Si = dent_rw.sd,
                      k_Ii = dent_rw.sd,
                      f_Sn=dent_rw.sd,
                      f_Si=dent_rw.sd,
                      rn=dent_rw.sd,
                      ri=dent_rw.sd,
                      probn=dent_rw.sd,
                      probi=dent_rw.sd,
                      sigP = dent_rw.sd),
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