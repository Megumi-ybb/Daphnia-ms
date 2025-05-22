library(ggplot2)
library(gridExtra)
library(reshape2)
library(magrittr)
library(foreach)
library(readxl)
library(doParallel)
registerDoParallel(cores=36)
library(pomp)
library(panelPomp)
library(tidyverse)

Mesocosm_data = read_excel("./Mesocosmdata.xls",3)

sed = 0923
set.seed(0416)


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

pt <- parameter_trans(
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
       paramnames = c("xi","sigSn", "sigIn", "sigSi", "sigIi", "sigF","sigP","theta_Sn","theta_In","theta_Si","theta_P",
                      "theta_Ii","f_Sn","f_Si","rn","ri","probn","probi","k_Ii","k_In","k_Sn","k_Si","sigJi","sigJn","theta_Jn","theta_Ji"),
       statenames = c("Sn",  "In", "Si","Jn","Ji" , "Ii","error_count", "F", "T_Sn","T_In","T_Si","T_Ii","P")
  ) -> pomplist[[i]]
  coef(pomplist[[i]])=parameters
}
names(pomplist)=paste("u", 1:8,sep = "")



shared_parameter = c(
  ri        = 2.152877e+05,
  rn        = 4.082784e+01,
  f_Si      = 2.419416e-07,
  f_Sn      = 1.100347e-03,
  probi     = 1.341780e+03,
  probn     = 2.715528e-01,
  xi        = 2.223829e+01,
  theta_Sn  = 8.478459e-04,
  theta_Si  = 2.524539e-03,
  theta_Ii  = 3.854897e-01,
  theta_In  = 5.837784e-01,
  theta_P   = 9.477111e-04,
  theta_Ji  = 5.620201e-04,
  theta_Jn  = 1.868651e-05,
  sigSn     = 0.000000e+00,
  sigSi     = 0.000000e+00,
  sigIn     = 2.930410e-04,
  sigIi     = 1.727540e-07,
  sigJi     = 3.021246e-01,
  sigJn     = 2.840041e-01,
  sigF      = 1.436943e-01,
  sigP      = 2.714232e-01,
  k_Ii      = 1.387153e+00,
  k_In      = 9.020138e-01,
  k_Si      = 5.262009e+00,
  k_Sn      = 4.103463e+00
)




panelfood = panelPomp(pomplist, shared=shared_parameter)
load('data/Residual_data/xi_001.RData')
pf.loglik.of.mif.estimate
s.e.of.pf.loglik.of.mif.estimate
dentNoPara = Mesocosm_data[91:170, ]
dentNoPara <- subset(dentNoPara, select = c(rep, day, dent.adult,dent.inf,lum.adult,lum.adult.inf,
                                            lum.juv,dent.juv))
dentNoPara <- dentNoPara[80: 1, ]
#Convert sampling date to natural date. Samples are collected every 5 days after the #first 7 days.
dentNoPara$day = (dentNoPara$day - 1) * 5 + 7
data = list()
trails = c("K","L","M","N","O","P","Q","S")
for (i in 1: length(trails)){
  data[[i]]=subset(dentNoPara, select = c("day", "dent.adult","dent.inf","lum.adult","lum.adult.inf",
                                          "lum.juv","dent.juv"),
                   dentNoPara$rep == trails[i])
}



coef(panelfood) <- coef(mf[[best]]$mif)

foreach(u = names(panelfood), .combine = rbind) %do% {
  unit_model <- unit_objects(panelfood)[[u]]
  shared <- coef(mf[[best]]$mif)
  pomp::simulate(
    unit_model, 
    nsim = 1000,
    format = "data.frame",
    params = c(shared)
  ) -> sims
  
  sims$unit <- u
  sims
} -> all_sims

data_df <- data[[1]]
data_df$unit <- "u1"
for (i in 2:length(data)) {
  df_unit <- data[[i]]
  df_unit$unit <- paste0('u', i)
  
  data_df <- rbind(data_df, df_unit)
}
all_sims[all_sims$error_count == 0,] -> all_sims001

load('./Mixed-species/SIRJPF2/model/best_result.RData')
pf.loglik.of.mif.estimate
s.e.of.pf.loglik.of.mif.estimate
dentNoPara = Mesocosm_data[91:170, ]
dentNoPara <- subset(dentNoPara, select = c(rep, day, dent.adult,dent.inf,lum.adult,lum.adult.inf,
                                            lum.juv,dent.juv))
dentNoPara <- dentNoPara[80: 1, ]
#Convert sampling date to natural date. Samples are collected every 5 days after the #first 7 days.
dentNoPara$day = (dentNoPara$day - 1) * 5 + 7
data = list()
trails = c("K","L","M","N","O","P","Q","S")
for (i in 1: length(trails)){
  data[[i]]=subset(dentNoPara, select = c("day", "dent.adult","dent.inf","lum.adult","lum.adult.inf",
                                          "lum.juv","dent.juv"),
                   dentNoPara$rep == trails[i])
}



coef(panelfood) <- coef(mf[[best]]$mif)

foreach(u = names(panelfood), .combine = rbind) %do% {
  unit_model <- unit_objects(panelfood)[[u]]
  shared <- coef(mf[[best]]$mif)
  pomp::simulate(
    unit_model, 
    nsim = 1000,
    format = "data.frame",
    params = c(shared)
  ) -> sims
  
  sims$unit <- u
  sims
} -> all_sims



all_sims001_summary_lum <- all_sims001 %>%
  group_by(day) %>%
  summarise(lower = quantile(Ji, 0.025),
            upper = quantile(Ji, 0.975),
            mean  = mean(Ji)) %>%
  mutate(species = "D. lumholtzi", simulation_type = "all_sims001")

all_sims001_summary_den <- all_sims001 %>%
  group_by(day) %>%
  summarise(lower = quantile(Jn, 0.025),
            upper = quantile(Jn, 0.975),
            mean  = mean(Jn)) %>%
  mutate(species = "D. dentifera", simulation_type = "all_sims001")

all_sims_summary_lum <- all_sims %>%
  group_by(day) %>%
  summarise(lower = quantile(Ji, 0.025),
            upper = quantile(Ji, 0.975),
            mean  = mean(Ji)) %>%
  mutate(species = "D. lumholtzi", simulation_type = "all_sims")

all_sims_summary_den <- all_sims %>%
  group_by(day) %>%
  summarise(lower = quantile(Jn, 0.025),
            upper = quantile(Jn, 0.975),
            mean  = mean(Jn)) %>%
  mutate(species = "D. dentifera", simulation_type = "all_sims")

all_sims_summary_combined <- bind_rows(all_sims001_summary_lum, all_sims001_summary_den,
                                       all_sims_summary_lum, all_sims_summary_den)

all_sims_summary_combined <- all_sims_summary_combined %>%
  mutate(
    species_label = case_when(
      species == "D. lumholtzi" ~ "italic(D.~lumholtzi)",
      species == "D. dentifera" ~ "italic(D.~dentifera)"
    ),
    simulation_label = case_when(
      simulation_type == "all_sims001" ~ "xi[J,u]==0.001",
      simulation_type == "all_sims" ~ "xi[J,u]==1"
    )
  )

data_df_lum_all_sims001 <- data_df %>%
  mutate(species_label = "italic(D.~lumholtzi)",
         simulation_label = "xi[J,u]==0.001",
         y = lum.juv)

data_df_lum_all_sims <- data_df %>%
  mutate(species_label = "italic(D.~lumholtzi)",
         simulation_label = "xi[J,u]==1",
         y = lum.juv)

data_df_den_all_sims001 <- data_df %>%
  mutate(species_label = "italic(D.~dentifera)",
         simulation_label = "xi[J,u]==0.001",
         y = dent.juv)

data_df_den_all_sims <- data_df %>%
  mutate(species_label = "italic(D.~dentifera)",
         simulation_label = "xi[J,u]==1",
         y = dent.juv)

data_df_faceted <- bind_rows(
  data_df_lum_all_sims001,
  data_df_lum_all_sims,
  data_df_den_all_sims001,
  data_df_den_all_sims
)


p_combined <- ggplot() +
  geom_ribbon(
    data = subset(all_sims_summary_combined, simulation_type == "all_sims001"),
    aes(x = day, ymin = lower, ymax = upper),
    fill = "#6CA6CD", alpha = 0.7
  ) +
  
  geom_ribbon(
    data = subset(all_sims_summary_combined, simulation_type == "all_sims"),
    aes(x = day, ymin = lower, ymax = upper),
    fill = "#ADD8E6", alpha = 0.5
  ) +
  

  geom_line(data = data_df_faceted,
            aes(x = day, y = y, group = unit, color = unit)) +
  

  geom_line(
    data = all_sims_summary_combined,
    aes(x = day, y = mean),
    size = 1
  ) +
  
  facet_grid(species_label ~ simulation_label, labeller = label_parsed) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 10),
    legend.position = "none"
  ) +
  labs(x = "Day", y = "Juvenile Density") +
  scale_x_continuous(limits = c(0, 52), breaks = c(0, 25, 50))+
  scale_y_continuous(transform = 'sqrt')


p_combined

