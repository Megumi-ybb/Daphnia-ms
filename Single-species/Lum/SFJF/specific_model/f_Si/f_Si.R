library(reshape2)
library(magrittr)
library(foreach)
library(readxl)
library(doParallel)
registerDoParallel(cores=36)
library(pomp)
library(panelPomp)
library(tidyverse)


create_parameters <- function(parameter_names, parameters) {
  shared_parameter <- parameters[!(names(parameters) %in% parameter_names)]
  
  specific_values <- lapply(parameter_names, function(p) rep(parameters[[p]], 10))
  specific_mat_data <- unlist(specific_values)
  
  specific_mat <- matrix(
    data = specific_mat_data,
    nrow = length(parameter_names),
    ncol = 10,
    byrow = TRUE,
    dimnames = list(
      parameter_names,
      c("u1", "u2", "u3", "u4", "u5", "u6", "u7", "u8","u9","u10")
    )
  )
  
  return(list(shared_parameter = shared_parameter, specific_mat = specific_mat))
}


create_specific_name <- function(parameter_names) {
  specific_name <- paste0("specific_", paste(parameter_names, collapse = "_"))
  return(specific_name)
}


# Mesocosm_data = read_excel("/Users/ybb/Desktop/Research//Daphnia/Mesocosmdata.xls")
Mesocosm_data = read_excel("/home/ybb/D_P/Mesocosmdata.xlsx",2)

sed = 0923
set.seed(0923)

specific_names = c('f_Si')


name_str = create_specific_name(specific_names)
run_level <- 2

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

parameters=  c(
  ri = 4.931153e+02,
  f_Si = 3.717773e-04,
  theta_Si = 4.976851e-01,
  theta_Ji = 1.800770e-01,
  sigSi = 0,
  sigJi = 4.045584e-01,
  sigF = 6.020674e-02,
  k_Si = 2.688813e+00
)


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


result <- create_parameters(specific_names, parameters)

shared_parameter <- result$shared_parameter
specific_mat <- result$specific_mat

panelfood = panelPomp(pomplist, shared=shared_parameter,specific = specific_mat)

algorithmic.params <- list(
  Np =     c(50, 500, 1e4),
  Np_rep = c( 2,  10,  10),
  Mp =     c(50, 500, 1e4),
  Nmif =   c( 2,  300, 250)
)

dent_rw.sd= 0.05
parameter_candidates = list(shared_parameter, specific_mat)
names(parameter_candidates)=c("shared","specific")
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
        Nmif = algorithmic.params$Nmif[run_level],
        shared.start = guessed.parameter.values$shared,
        specific.start = guessed.parameter.values$specific,
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


shared_dataframe <- data.frame(t(mf1[[select[1]]]$mif@shared))
for (i in 1:length(select)) {
  shared_dataframe[i, ] <- t(mf1[[select[i]]]$mif@shared)
}

replications <- 4
shared_dataframe <- shared_dataframe[rep(1:nrow(shared_dataframe), each = replications), ]

shared_dataframe$replication <- 1:(nrow(shared_dataframe))

specific_list <- list()
for (i in 1:length(select)) {
  specific_list[[i]] <- mf1[[select[i]]]$mif@specific
}

replicated_specific_list <- list()
for (i in 1:length(specific_list)) {
  replicated_matrices <- replicate(replications, specific_list[[i]], simplify = FALSE)
  replicated_specific_list <- c(replicated_specific_list, replicated_matrices)
}

specific_dataframes <- lapply(seq_along(replicated_specific_list), function(k) {
  mat <- replicated_specific_list[[k]]
  df <- as.data.frame(mat)
  return(df)
})

dent_rw.sd = 0.04

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
      
      specific_para_temp = as.matrix(specific_dataframes[[i]])
      
      mif2(
        panelfood,
        Nmif = 200,
        shared.start = share_para_temp,
        specific.start = specific_para_temp,
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
  save(mf,lls,best,mif.estimate,pf.loglik.of.mif.estimate,
       s.e.of.pf.loglik.of.mif.estimate,
       file = paste0(name_str,".RData"))
}

if (run_level == 3){
  save(lls,best,mif.estimate,pf.loglik.of.mif.estimate,
       s.e.of.pf.loglik.of.mif.estimate,
       file = paste0(name_str,".RData"))
}
