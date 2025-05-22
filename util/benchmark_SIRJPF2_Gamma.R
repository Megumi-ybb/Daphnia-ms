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

name_str = "all_shared_gamma"
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


#Notice that the unit of S is individuals/L, unit of F is 1e+5 cells/L, unit of P is 1000 spores/L, unit for I is individuals/L

dyn_rpro <- Csnippet("
                      double Sn_term, In_term, F_term, P_term , Si_term, Ii_term,Jn_term,Ji_term;
                      double noiSIn, noiIPn, noiSIi , noiIPi ,noiF, noiP,noiSJn,noiSJi;
                      double delta = 0.013; //fraction of volume replaced day-1

                      noiSIn = rgamma(dt/ (sigSIn * sigSIn), sigSIn * sigSIn);
                      noiSJn = rgamma(dt/ (sigSJn * sigSJn), sigSJn * sigSJn); 
                      noiIPn = rgamma(dt/ (sigIPn * sigIPn), sigIPn * sigIPn);
                      noiSIi = rgamma(dt/ (sigSIi * sigSIi), sigSIi * sigSIi);
                      noiSJi = rgamma(dt/ (sigSJi * sigSJi), sigSJi * sigSJi); 
                      noiIPi = rgamma(dt/ (sigIPi * sigIPi), sigIPi * sigIPi);
                      noiF = rnorm(0, sigF * sqrt(dt));
                      noiP = rgamma(dt/ (sigP * sigP), sigP * sigP);


                      //------------Sn-------------
                      Sn_term = 0.1 * Jn * dt  - theta_Sn * Sn * dt -  probn * f_Sn * Sn * P * noiSIn - delta * Sn * dt;
                        
                      //-----------Jn-------------
                      
                      Jn_term = rn * f_Sn * F * Sn  * noiSJn  -  0.1 * Jn * dt - theta_Jn * Jn * dt - delta * Jn * dt;
                      
                      //------------In--------------
                      In_term = probn * f_Sn * P * Sn * noiSIn - theta_In * In * noiIPn - delta * In * dt;

                      //------------Si-------------
                      Si_term = 0.1 * Ji * dt  - theta_Si * Si * dt -  probi * f_Si * Si * P * noiSIi - delta * Si * dt;

                      //-----------Ji-------------
                      
                      Ji_term = ri * f_Si * F * Si  * noiSJi - 0.1 * Ji * dt - theta_Ji * Ji * dt - delta * Ji * dt;
                      
                      //------------Ii--------------

                      Ii_term = probi * f_Si * P * Si * noiSIi - theta_Ii * Ii * noiIPi - delta * Ii * dt;


                      //-----------F---------------

                      F_term =  F * noiF - f_Sn * F * (Sn + xi * In + 1 * Jn) * dt - f_Si * F * (Si + xi * Ii + 1 * Ji) * dt - delta * F * dt + 0.37 * dt;


                      //----------P---------------

                      P_term = 30 * theta_In * In * noiIPn + 30 * theta_Ii * Ii * noiIPi - f_Sn * (Sn + xi * In) * P * dt - f_Si * (Si + xi * Ii) * P * dt - theta_P * P * noiP - delta * P * dt;


                      //T_D and T_I are the S and I Daphnia sample out
                      F += F_term;
                      Sn += Sn_term;
                      In += In_term;
                      Si += Si_term;
                      Ji += Ji_term;
                      Jn += Jn_term;
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
                      if (Jn < 0.0 || Jn > 1e5) {
                        Jn = 0.0;
                        error_count += 0.001;
                      }
                      if (Ji < 0.0 || Ji > 1e5) {
                        Ji = 0.0;
                        error_count += 0.000000001;
                      }   
                      if ((P < 0.0 || P > 1e20)&& t > 3.9) {
                        P = 0.0;
                        error_count += 0.000001;
                      }
                      
                      //T_Sn = 15 * delta * Sn * dt;
                      //T_In = 15 * delta * In * dt;
                      //T_Si = 15 * delta * Si * dt;
                      //T_Ii = 15 * delta * Ii * dt;
                      
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
                     Jn = 0;
                     Ji = 0;
                     
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
  log = c( "sigSIn", "sigIPn", "sigSIi", "sigIPi", "sigF","sigP","f_Sn","f_Si",
           "rn","ri","k_Ii","k_In","k_Sn","k_Si","sigSJi","sigSJn","probn","probi",
           "theta_Sn","theta_In","theta_Si","theta_Ii","theta_P","theta_Jn","theta_Ji","xi"),
)

pomplist=list()

parameters =         c(   8.34e-06,     0.4,    0.5,    0.8,     0.2,     0.01,  0.2,    0.15,       
                          0.07 ,   0.2,       0.01,    0.096  , 0.0035 ,0.003 , 0.2  , 
                          0.4, 0.6  ,0.8 ,8,8,8,8,1,1,0.1,0.1)
names(parameters) =  c( "sigSIn", "sigIPn", "sigSIi", "sigIPi", "sigF","sigP","f_Sn","f_Si",
                        "rn","ri","k_Ii","k_In","k_Sn","k_Si","sigSJi","sigSJn","probn","probi",
                        "theta_Sn","theta_In","theta_Si","theta_Ii","theta_P","theta_Jn","theta_Ji","xi")


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
       paramnames = c( "sigSIn", "sigIPn", "sigSIi", "sigIPi", "sigF","sigP","f_Sn","f_Si",
                       "rn","ri","k_Ii","k_In","k_Sn","k_Si","sigSJi","sigSJn","probn","probi",
                       "theta_Sn","theta_In","theta_Si","theta_Ii","theta_P","theta_Jn","theta_Ji","xi"),
       statenames = c("Sn",  "In", "Si","Jn","Ji" , "Ii","error_count", "F", "T_Sn","T_In","T_Si","T_Ii","P")
  ) -> pomplist[[i]]
  coef(pomplist[[i]])=parameters
}
names(pomplist)=paste("u", 1:8,sep = "")



shared_parameter = c(
  f_Si       = 3.583694e-05,
  sigSIn     = 1.490538e+00,
  sigIPn     = 4.751448e-01,
  sigSIi     = 1.216466e+00,
  sigIPi     = 1.009119e+00,
  sigF       = 1.100299e-01,
  sigP       = 4.684769e-02,
  sigSJi     = 3.694885e-02,
  sigSJn     = 1.451527e-04,
  f_Sn       = 2.354344e-04,
  theta_Jn   = 4.453522e+00,
  theta_Ji   = 1.499845e+00,
  probi      = 9.121047e+00,
  probn      = 1.159688e+00,
  ri         = 6.210064e+03,
  rn         = 2.297429e+03,
  theta_Ii   = 1.562697e-01,
  theta_In   = 3.163630e-01,
  theta_Si   = 2.100676e-06,
  theta_Sn   = 1.252143e-03,
  theta_P    = 1.736973e-06,
  xi         = 1.469241e+02,
  k_Ii       = 1.030235e+00,
  k_In       = 2.043363e+00,
  k_Sn       = 2.257618e+00,
  k_Si       = 4.076402e+00
)


panelfood = panelPomp(pomplist, shared=shared_parameter)

N_START_POINTS <- 360
r_unit_params <- list()


shared_lb <- shared_parameter / 10
shared_ub <- shared_parameter * 10


r_shared_params <- runif_design(
  lower = shared_lb,
  upper = shared_ub,
  nseq = N_START_POINTS
)



algorithmic.params <- list(
  Np =     c(50, 800, 1e4),
  Np_rep = c( 2,  15,  10),
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
      guessed.parameter.values <- as.numeric(r_shared_params[i,])
      names(guessed.parameter.values) = names(r_shared_params)
      mif2(
        panelfood,
        Nmif = 150,
        shared.start = guessed.parameter.values,
        rw.sd = rw_sd(xi=dent_rw.sd,
                      sigSIn=dent_rw.sd,
                      sigIPn=dent_rw.sd,
                      sigSIi=dent_rw.sd,
                      sigIPi=dent_rw.sd,
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
                      sigP = dent_rw.sd,
                      sigSJi = dent_rw.sd,
                      sigSJn = dent_rw.sd,
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
        Nmif = 150,
        shared.start = share_para_temp,
        rw.sd = rw_sd(xi=dent_rw.sd,
                      sigSIn=dent_rw.sd,
                      sigIPn=dent_rw.sd,
                      sigSIi=dent_rw.sd,
                      sigIPi=dent_rw.sd,
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
                      sigP = dent_rw.sd,
                      sigSJi = dent_rw.sd,
                      sigSJn = dent_rw.sd,
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
        Nmif = 400,
        shared.start = share_para_temp,
        rw.sd = rw_sd(xi=dent_rw.sd,
                      sigSIn=dent_rw.sd,
                      sigIPn=dent_rw.sd,
                      sigSIi=dent_rw.sd,
                      sigIPi=dent_rw.sd,
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
                      sigP = dent_rw.sd,
                      sigSJi = dent_rw.sd,
                      sigSJn = dent_rw.sd,
                      theta_Jn = dent_rw.sd,
                      theta_Ji = dent_rw.sd),
        cooling.type = "geometric",
        cooling.fraction.50 = 0.7,
        Np = algorithmic.params$Mp[run_level]
      ) -> m1
      
      ll <- replicate(n = algorithmic.params$Np_rep[run_level],
                      unitlogLik(pfilter(m1,
                                         Np = algorithmic.params$Np[run_level])))
      res = pfilter(m1,
                    Np = algorithmic.params$Np[run_level],save.states = 'unweighted', filter.mean = TRUE, filter.traj = TRUE)
      list(mif = m1,
           ll = panel_logmeanexp(x = ll,
                                 MARGIN = 1,
                                 se = TRUE),
           res = res)
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