library(reshape2)
library(magrittr)
library(foreach)
library(readxl)
library(doParallel)
registerDoParallel(cores=36)
library(pomp)
library(panelPomp)
library(tidyverse)

# Mesocosm_data = read_excel("/Users/ybb/Desktop/Research//Daphnia/Mesocosmdata.xls",3)
Mesocosm_data = read_excel("/home/ybb/D_P/Mesocosmdata.xlsx",3)


DEBUG = FALSE


sed = sample(1:1000000,1)
set.seed(sed)

name_str = "sigIn"
run_level = 3

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


#Notice that the unit of S is individuals/L, unit of F is 1e+5 cells/L, unit of P is 1000 spores/L, unit for I is individuals/L

dyn_rpro = Csnippet("
                      double Sn_term, In_term, F_term, P_term , Si_term, Ii_term,Jn_term,Ji_term;
                      double noiSn, noiIn, noiSi , noiIi ,noiF, noiP,noiJn,noiJi;
                      double delta = 0.013; //fraction of volume replaced day-1

                      noiSn = rnorm(0, sigSn * sqrt(dt));
                      noiIn = rnorm(0, sigIn * sqrt(dt));
                      noiSi = rnorm(0, sigSi * sqrt(dt));
                      noiIi = rnorm(0, sigIi * sqrt(dt));
                      noiJn = rnorm(0, sigJn * sqrt(dt));
                      noiJi = rnorm(0, sigJi * sqrt(dt));
                      noiF = rnorm(0, sigF * sqrt(dt));
                      noiP = rnorm(0, sigP * sqrt(dt));


                      //------------Sn-------------
                      Sn_term = 0.1 * Jn * dt - theta_Sn * Sn * dt -  probn * f_Sn * Sn * P * dt - delta * Sn * dt + Sn * noiSn;
                        
                      //-----------Jn-------------
                      
                      Jn_term = rn * f_Sn * F * Sn  * dt  -  0.1 * Jn * dt - theta_Jn * Jn * dt - delta * Jn * dt + Jn * noiJn;
                      
                      //------------In--------------

                      In_term = probn * f_Sn * Sn * P * dt - theta_In * In *dt - delta * In * dt + In * noiIn;

                      //------------Si-------------
                      Si_term = 0.1 * Ji * dt - theta_Si * Si * dt -  probi * f_Si * Si * P * dt - delta * Si * dt + Si * noiSi;

                      //-----------Ji-------------
                      
                      Ji_term = ri * f_Si *F * Si  * dt - 0.1 * Ji * dt - theta_Ji * Ji * dt - delta * Ji * dt + Ji * noiJi;
                      
                      //------------Ii--------------

                      Ii_term = probi * f_Si * Si * P * dt - theta_Ii * Ii *dt - delta * Ii * dt + Ii * noiIi;


                      //-----------F---------------

                      F_term =  F * noiF - f_Sn * F * (Sn + xi * In + 1 * Jn) * dt - f_Si * F * (Si + xi * Ii + 1 * Ji) * dt - delta * F * dt + 0.37 * dt;


                      //----------P---------------

                      P_term = 30 * theta_In * In * dt + 30 * theta_Ii * Ii * dt - f_Sn * (Sn + xi * In) * P * dt - f_Si * (Si + xi * Ii) * P * dt- theta_P * P * dt - delta * P * dt + P * noiP;


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

pt = parameter_trans(
  #Without death part
  log = c( "sigSn", "sigIn", "sigSi", "sigIi", "sigF","sigP","f_Sn","f_Si",
           "rn","ri","k_Ii","k_In","k_Sn","k_Si","sigJi","sigJn","probn","probi",
           "theta_Sn","theta_In","theta_Si","theta_Ii","theta_P","theta_Jn","theta_Ji","xi"),
)

pomplist=list()

parameters =         c(   8.34e-06,     0.4,    0.5,    0.8,     0.2,     0.01,  0.2,    0.15,       
                          0.07 ,   0.2,       0.01,    0.096  , 0.0035 ,0.003 , 0.2  , 
                          0.4, 0.6  ,0.8 ,8,8,8,8,1,1,0.1,0.1)
names(parameters) =  c("xi", "sigSn", "sigIn", "sigSi", "sigIi", "sigF","sigP","theta_Sn",
                       "theta_In","theta_Si","theta_P","theta_Ii","f_Sn","f_Si","rn","ri","probn",
                       "probi","k_Ii","k_In","k_Sn","k_Si","sigJi","sigJn","theta_Jn","theta_Ji")


for (i in 1:8){
  colnames(data[[i]]) = c('day', 'dentadult', 'dentinf','lumadult','luminf')
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
       paramnames = c("xi","sigSn", "sigIn", "sigSi", "sigIi", "sigF","sigP","theta_Sn","theta_In","theta_Si","theta_P",
                      "theta_Ii","f_Sn","f_Si","rn","ri","probn","probi","k_Ii","k_In","k_Sn","k_Si","sigJi","sigJn","theta_Jn","theta_Ji"),
       statenames = c("Sn",  "In", "Si","Jn","Ji" , "Ii","error_count", "F", "T_Sn","T_In","T_Si","T_Ii","P")
  ) -> pomplist[[i]]
  coef(pomplist[[i]])=parameters
}
names(pomplist)=paste("u", 1:8,sep = "")



shared_parameter = c(
  f_Si     = 1.980785e-04,
  f_Sn     = 1.033392e-03,
  theta_Jn = 9.511384e-04,
  theta_Ji = 1.198650e-06,
  probi    = 1.697204e+00,
  probn    = 2.933072e-01,
  ri       = 2.561561e+02,
  rn       = 4.958747e+01,
  theta_Ii = 3.942674e-01,
  theta_In = 4.986904e-01,
  theta_Si = 7.410951e-04,
  theta_Sn = 3.369496e-02,
  theta_P  = 1.499839e-03,
  xi       = 1.889229e+01,
  sigSn    = 0,
  sigIn    = 3.327712e-05,
  sigSi    = 0,
  sigIi    = 1.323957e-05,
  sigF     = 1.233165e-01,
  k_Ii     = 1.202892e+00,
  k_In     = 8.369364e-01,
  k_Sn     = 3.923281e+00,
  k_Si     = 6.019108e+00,
  sigP     = 2.972195e-01,
  sigJi    = 2.453938e-01,
  sigJn    = 2.235086e-01
)




panelfood = panelPomp(pomplist, shared=shared_parameter)

generate_parameter_profile = function(prof_name, nprof = 80) {
  shared_ub = shared_parameter * 10
  
  shared_lb = shared_ub / 100
  
  ub_unit = log(shared_ub[prof_name])
  lb_unit = log(shared_lb[prof_name])
  
  shared_lb = shared_lb[ !(names(shared_lb) %in% c("sigSn", "sigSi",prof_name)) ]
  shared_ub = shared_ub[ !(names(shared_ub) %in% c("sigSn", "sigSi",prof_name)) ]
  
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
  parameter_shared$sigSn = 0
  parameter_shared$sigSi = 0
  return(parameter_shared)
}



generate_sd = function(x = 0.05, profile_name){
  sd_list = c(
    ri        = x,
    rn        = x,
    f_Si      = x,
    f_Sn      = x,
    probi     = x,
    probn     = x,
    xi        = x,
    theta_Sn  = x,
    theta_Si  = x,
    theta_Ii  = x,
    theta_In  = x,
    theta_P   = x,
    theta_Ji  = x,
    theta_Jn  = x,
    sigSn     = 0,
    sigSi     = 0,
    sigIn     = x,
    sigIi     = x,
    sigJi     = x,
    sigJn     = x,
    sigF      = x,
    sigP      = x,
    k_Ii      = x,
    k_In      = x,
    k_Si      = x,
    k_Sn      = x
  )
  
  sd_list[profile_name] = 0
  return(sd_list)
}

parameter_shared = generate_parameter_profile(name_str)
# profile_para = profile_design(alpha = 0.0001:0.01,lower= shared_lb,upper = shared_ub,nprof = 100)

algorithmic.params = list(
  Np =     c(50, 320, 1e4),
  Np_rep = c( 2,  10,  10),
  Mp =     c(50, 400, 1e4),
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
      guessed.parameter.values = as.numeric(parameter_shared[i,])
      names(guessed.parameter.values) = colnames(parameter_shared)
      mif2(
        panelfood,
        Nmif = 200,
        shared.start = guessed.parameter.values,
        rw.sd = rw_sd(xi=dent_rw_sd_first['xi'],
                      sigSn=dent_rw_sd_first['sigSn'],
                      sigIn=dent_rw_sd_first['sigIn'],
                      sigSi=dent_rw_sd_first['sigSi'],
                      sigIi=dent_rw_sd_first['sigIi'],
                      sigF=dent_rw_sd_first['sigF'],
                      theta_Sn=dent_rw_sd_first['theta_Sn'],
                      theta_In=dent_rw_sd_first['theta_In'],
                      theta_Si=dent_rw_sd_first['theta_Si'],
                      theta_P=dent_rw_sd_first['theta_P'],
                      theta_Ii=dent_rw_sd_first['theta_Ii'],
                      k_Sn = dent_rw_sd_first['k_Sn'],
                      k_In = dent_rw_sd_first['k_In'],
                      k_Si = dent_rw_sd_first['k_Si'],
                      k_Ii = dent_rw_sd_first['k_Ii'],
                      f_Sn=dent_rw_sd_first['f_Sn'],
                      f_Si=dent_rw_sd_first['f_Si'],
                      rn=dent_rw_sd_first['rn'],
                      ri=dent_rw_sd_first['ri'],
                      probn=dent_rw_sd_first['probn'],
                      probi=dent_rw_sd_first['probi'],
                      sigP = dent_rw_sd_first['sigP'],
                      sigJi = dent_rw_sd_first['sigJi'],
                      sigJn = dent_rw_sd_first['sigJn'],
                      theta_Jn = dent_rw_sd_first['theta_Jn'],
                      theta_Ji = dent_rw_sd_first['theta_Ji']),
        cooling.type = "geometric",
        cooling.fraction.50 = 0.7,
        Np = algorithmic.params$Mp[run_level]
        
      ) -> m1
      
      ll = replicate(n = algorithmic.params$Np_rep[run_level],
                     unitlogLik(pfilter(m1,
                                        Np = algorithmic.params$Np[run_level])))
      # No need to save complete pomp object, j
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
  shared_dataframe = shared_dataframe[rep(1:nrow(shared_dataframe), each = 4), ]
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
  shared_dataframe = shared_dataframe[rep(1:nrow(shared_dataframe), each = 4), ]
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
        rw.sd = rw_sd(xi=dent_rw_sd_second['xi'],
                      sigSn=dent_rw_sd_second['sigSn'],
                      sigIn=dent_rw_sd_second['sigIn'],
                      sigSi=dent_rw_sd_second['sigSi'],
                      sigIi=dent_rw_sd_second['sigIi'],
                      sigF=dent_rw_sd_second['sigF'],
                      theta_Sn=dent_rw_sd_second['theta_Sn'],
                      theta_In=dent_rw_sd_second['theta_In'],
                      theta_Si=dent_rw_sd_second['theta_Si'],
                      theta_P=dent_rw_sd_second['theta_P'],
                      theta_Ii=dent_rw_sd_second['theta_Ii'],
                      k_Sn = dent_rw_sd_second['k_Sn'],
                      k_In = dent_rw_sd_second['k_In'],
                      k_Si = dent_rw_sd_second['k_Si'],
                      k_Ii = dent_rw_sd_second['k_Ii'],
                      f_Sn=dent_rw_sd_second['f_Sn'],
                      f_Si=dent_rw_sd_second['f_Si'],
                      rn=dent_rw_sd_second['rn'],
                      ri=dent_rw_sd_second['ri'],
                      probn=dent_rw_sd_second['probn'],
                      probi=dent_rw_sd_second['probi'],
                      sigP = dent_rw_sd_second['sigP'],
                      sigJi = dent_rw_sd_second['sigJi'],
                      sigJn = dent_rw_sd_second['sigJn'],
                      theta_Jn = dent_rw_sd_second['theta_Jn'],
                      theta_Ji = dent_rw_sd_second['theta_Ji']),
        cooling.type = "geometric",
        cooling.fraction.50 = 0.7,
        Np = algorithmic.params$Mp[run_level]
      ) -> m1
      
      ll = replicate(n = algorithmic.params$Np_rep[run_level],
                     unitlogLik(pfilter(m1,
                                        Np = algorithmic.params$Np[run_level])))
      
      list(mif = m1,
           ll = panel_logmeanexp(x = ll,
                                 MARGIN = 1,
                                 se = TRUE))
    }
} -> mf



lls = matrix(unlist(sapply(mf, getElement, "ll")), nrow = 2)
best = which.max(lls[1,])
mif.estimate = coef(mf[[best]]$mif)
pf.loglik.of.mif.estimate = unname(mf[[best]]$ll[1])
s.e.of.pf.loglik.of.mif.estimate = unname(mf[[best]]$ll[2])

trace = as.data.frame(traces(mf[[1]]$mif))
trace$iter = as.numeric(rownames(trace))

for (i in 2 : length(mf)){
  tmp_df = as.data.frame(traces(mf[[i]]$mif))
  tmp_df$iter = as.numeric(rownames(tmp_df))
  trace = dplyr::bind_rows(trace,tmp_df)
}

final_likes = numeric(length(mf))

for (i in 1: length(mf)){
  final_likes[i] = mf[[i]]$ll[1]
}

final_params = trace %>%
  dplyr::filter(iter == max(iter, na.rm = TRUE))

final_params$loglik = final_likes

pf.loglik.of.mif.estimate
s.e.of.pf.loglik.of.mif.estimate

if (run_level %in% c(2,3)){
  save(mf,final_params,lls,best,mif.estimate,pf.loglik.of.mif.estimate,
       s.e.of.pf.loglik.of.mif.estimate,
       file = paste0(name_str,".RData"))
}

# if (run_level == 3){
#   save(lls,best,mif.estimate,pf.loglik.of.mif.estimate,
#        s.e.of.pf.loglik.of.mif.estimate,
#        file = paste0(name_str,".RData"))
# }

