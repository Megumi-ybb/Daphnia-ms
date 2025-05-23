library(reshape2)
library(magrittr)
library(foreach)
library(readxl)
library(doParallel)
n_cores = parallel::detectCores(logical = TRUE)
registerDoParallel(cores = n_cores)
library(pomp)
library(panelPomp)
library(tidyverse)

Mesocosm_data = read_excel("./Mesocosmdata.xlsx",2)

DEBUG = FALSE

 

name_str = "f_Si"
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
  theta_Ji  = 6.366139e-02,
  sigSi     = 0,
  sigF      = 3.580001e-03,
  f_Si      = 3.973271e-04,
  ri        = 4.909026e+02,
  k_Si      = 2.726407e+00,
  sigJi     = 3.856341e-01,
  theta_Si  = 6.647845e-01
)

panelfood = panelPomp(pomplist, shared=shared_parameter)

generate_parameter_profile = function(prof_name, nprof = 50) {
  shared_ub = shared_parameter * 100
  
  shared_lb = shared_ub / 100
  
  ub_unit = log(shared_ub[prof_name])
  lb_unit = log(shared_lb[prof_name])
  
  shared_lb = shared_lb[ !(names(shared_lb) %in% c("sigSi",prof_name)) ]
  shared_ub = shared_ub[ !(names(shared_ub) %in% c("sigSi",prof_name)) ]
  
  parameter_shared = pomp::profile_design(
    temp = seq(lb_unit, ub_unit, length.out = nprof),
    lower = log(shared_lb),
    upper = log(shared_ub),
    type = 'runif',
    nprof = nprof
  )
  
  parameter_shared = parameter_shared %>% 
    rename( !!prof_name := temp)  
  
  parameter_shared = exp(parameter_shared)
  parameter_shared$sigSi = 0
  return(parameter_shared)
}

generate_sd <- function(x = 0.05, profile_name){
  sd_list = c(
    sigF      = x,
    sigSi     = 0,
    f_Si      = x,
    ri        = x,
    k_Si      = x,
    sigJi     = x,
    theta_Si  = x,
    theta_Ji  = x
  )
  
  sd_list[profile_name] = 0
  return(sd_list)
}

parameter_shared <- generate_parameter_profile(name_str)

algorithmic.params = list(
  Np =     c(50, 320, 1e3),
  Np_rep = c( 2,  10,  20),
  Mp =     c(50, 400, 1000),
  Nmif =   c( 2,  320, 250)
)


dent_rw_sd_first = generate_sd(x = 0.05,profile_name = name_str)

U = length(panelfood)
{
  foreach(
    i = 1:nrow(parameter_shared),
    .packages = c("pomp", "panelPomp"),
    .inorder = FALSE,
    .options.multicore = list(set.seed = TRUE)
  ) %dopar%
    {
      guessed.parameter.values <- as.numeric(parameter_shared[i,])
      names(guessed.parameter.values) = colnames(parameter_shared)
      mif2(
        panelfood,
        Nmif = 300,
        shared.start = guessed.parameter.values,
        rw.sd = rw_sd(
          sigSi=dent_rw_sd_first['sigSi'],
          sigF=dent_rw_sd_first['sigF'],
          theta_Si=dent_rw_sd_first['theta_Si'],
          k_Si = dent_rw_sd_first['k_Si'],
          f_Si=dent_rw_sd_first['f_Si'],
          ri=dent_rw_sd_first['ri'],
          sigJi = dent_rw_sd_first['sigJi'],
          theta_Ji = dent_rw_sd_first['theta_Ji']),
        cooling.type = "geometric",
        cooling.fraction.50 = 0.7,
        Np = algorithmic.params$Mp[run_level]
        
      ) -> m1
      
      ll <- replicate(n = algorithmic.params$Np_rep[run_level],
                      unitlogLik(pfilter(m1,
                                         Np = algorithmic.params$Np[run_level])))
      
      if(DEBUG){      
        list(mif = m1,
             ll = panel_logmeanexp(x = ll,
                                   MARGIN = 1,
                                   se = TRUE))
      }else{
        list(mif_ceof = coef(m1),
             ll = panel_logmeanexp(x = ll,
                                   MARGIN = 1,
                                   se = TRUE))
      }
    }
} -> mf1

if(DEBUG){
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
  dent_rw_sd_second = generate_sd(x = 0.04,profile_name = name_str)
}else{
  log_list = c()
  for ( i in 1:length(mf1)){
    log_list = c(log_list,mf1[[i]]$ll[1])
  }
  select = order(log_list,decreasing = TRUE)[1:ceiling(length(mf1)/4)]
  shared_dataframe = data.frame(t(mf1[[select[1]]]$mif_ceof))
  for (i in 1:length(select)){
    shared_dataframe[i,] = t(mf1[[select[i]]]$mif_ceof)
  }
  shared_dataframe <- shared_dataframe[rep(1:nrow(shared_dataframe), each = 4), ]
  dent_rw_sd_second = generate_sd(x = 0.04,profile_name = name_str)
}
#round 2
{
  foreach(
    i = 1:nrow(parameter_shared),
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
          sigSi=dent_rw_sd_first['sigSi'],
          sigF=dent_rw_sd_first['sigF'],
          theta_Si=dent_rw_sd_first['theta_Si'],
          k_Si = dent_rw_sd_first['k_Si'],
          f_Si=dent_rw_sd_first['f_Si'],
          ri=dent_rw_sd_first['ri'],
          sigJi = dent_rw_sd_first['sigJi'],
          theta_Ji = dent_rw_sd_first['theta_Ji']),
        cooling.type = "geometric",
        cooling.fraction.50 = 0.7,
        Np = algorithmic.params$Mp[run_level]
      ) -> m1
      print(i)     
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

trace <- as.data.frame(traces(mf[[1]]$mif))
trace$iter = as.numeric(rownames(trace))

for (i in 2 : length(mf)){
  tmp_df <- as.data.frame(traces(mf[[i]]$mif))
  tmp_df$iter = as.numeric(rownames(tmp_df))
  trace <- dplyr::bind_rows(trace,tmp_df)
}

final_likes <- numeric(length(mf))

for (i in 1: length(mf)){
  final_likes[i] <- mf[[i]]$ll[1]
}

final_params <- trace %>%
  dplyr::filter(iter == max(iter, na.rm = TRUE))

final_params$loglik <- final_likes

pf.loglik.of.mif.estimate
s.e.of.pf.loglik.of.mif.estimate

if (run_level %in% c(2,3)){
  save(mf,final_params,lls,best,mif.estimate,pf.loglik.of.mif.estimate,
       s.e.of.pf.loglik.of.mif.estimate,
       file = paste0(name_str,".RData"))
}