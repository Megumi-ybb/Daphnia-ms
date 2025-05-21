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
  f_Si = 3.275316e-05, f_Sn = 1.105668e-03, theta_Jn = 7.833849e-06, theta_Ji = 2.599044e-04,
  probi = 3.419323e+01, probn = 3.660066e-01, ri = 1.307683e+04, rn = 7.360543e+01,
  theta_Ii = 3.531879e-01, theta_In = 5.489315e-01, theta_Si = 3.186040e-02, theta_Sn = 1.479834e-01,
  theta_P = 2.024991e-02, xi = 2.865620e+01,sigSn = 2.073521e-02, sigIn = 5.171328e-04, sigSi = 3.631343e-06, sigIi = 6.836503e-06, 
  sigF = 1.346211e-01,k_Ii = 1.317484e+00, k_In = 1.428054e+00, k_Sn = 3.580089e+00, 
  k_Si = 5.623301e+00, sigP = 5.278157e-02, sigJi = 2.878004e-01, sigJn = 2.641643e-01
)




panelfood = panelPomp(pomplist, shared=shared_parameter)
load('./Mixed-species/SIRJPF2/model/best_result.rda')

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



coef(panelfood) <- mif.estimate

foreach(u = names(panelfood), .combine = rbind) %do% {
  unit_model <- unit_objects(panelfood)[[u]]
  shared <-  mif.estimate
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
all_sims[all_sims$error_count == 0,] -> all_sims


all_sims_summary <- all_sims %>%
  group_by(day) %>%
  summarise(lower = quantile(Sn, 0.025),   # 2.5% quantile
            upper = quantile(Sn, 0.975),
            mean = mean(Sn)) 

p1 = ggplot() +
  geom_ribbon(data = all_sims_summary, aes(x = day, ymin = lower, ymax = upper), fill = "#ADD8E6", alpha = 0.5) +
  geom_line(data = data_df, aes(x = day, y =  dent.adult, group = unit, color = unit)) +
  geom_line(data = all_sims_summary, aes(x = day, y =  mean),size = 1) +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(x = NULL, y = expression(paste("Susceptible ", italic("D. dentifera"), " density"))) +
  theme(axis.title = element_text(size = 10), axis.text.x = element_blank())+
  scale_x_continuous(limits = c(0, 52), breaks = c(0, 25, 50))


all_sims_summary <- all_sims %>%
  group_by(day) %>%
  summarise(lower = quantile(In, 0.025),   # 2.5% quantile
            upper = quantile(In, 0.975),
            mean = mean(In)) 

p2 = ggplot() +
  geom_ribbon(data = all_sims_summary, aes(x = day, ymin = lower, ymax = upper), fill = "#ADD8E6", alpha = 0.5) +
  geom_line(data = data_df, aes(x = day, y =  dent.inf, group = unit, color = unit)) +
  geom_line(data = all_sims_summary, aes(x = day, y =  mean),size = 1) +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(x = NULL, y = expression(paste("Infected ", italic("D. dentifera"), " density"))) +
  theme(axis.title = element_text(size = 10), axis.text.x = element_blank())+
  scale_x_continuous(limits = c(0, 52), breaks = c(0, 25, 50))


all_sims_summary <- all_sims %>%
  group_by(day) %>%
  summarise(lower = quantile(Si, 0.025),   # 2.5% quantile
            upper = quantile(Si, 0.975),
            mean = mean(Si)) 

p3 = ggplot() +
  geom_ribbon(data = all_sims_summary, aes(x = day, ymin = lower, ymax = upper), fill = "#ADD8E6", alpha = 0.5) +
  geom_line(data = data_df, aes(x = day, y =  lum.adult, group = unit, color = unit)) +
  geom_line(data = all_sims_summary, aes(x = day, y =  mean),size = 1) +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(x = "day", y = expression(paste("Susceptible ", italic("D. lumholtzi"), " density"))) +
  theme(axis.title = element_text(size = 10))+
  scale_x_continuous(limits = c(0, 52), breaks = c(0, 25, 50))


all_sims_summary <- all_sims %>%
  group_by(day) %>%
  summarise(lower = quantile(Ii, 0.025),   # 2.5% quantile
            upper = quantile(Ii, 0.975),
            mean = mean(Ii)) 

p4 = ggplot() +
  geom_ribbon(data = all_sims_summary, aes(x = day, ymin = lower, ymax = upper), fill = "#ADD8E6", alpha = 0.5) +
  geom_line(data = data_df, aes(x = day, y =  lum.adult.inf, group = unit, color = unit)) +
  geom_line(data = all_sims_summary, aes(x = day, y =  mean),size = 1) +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(x = "day", y = expression(paste("Infected ", italic("D. lumholtzi"), " density"))) +
  theme(axis.title = element_text(size = 10))+
  scale_x_continuous(limits = c(0, 52), breaks = c(0, 25, 50))




#Juv
all_sims_summary <- all_sims %>%
  group_by(day) %>%
  summarise(lower = quantile(Ji, 0.025),   # 2.5% quantile
            upper =  quantile(Ji,0.975),
            mean =  mean(Ji)) 

p5 = ggplot() +
  geom_ribbon(data = all_sims_summary, aes(x = day, ymin = lower, ymax = upper), fill = "#ADD8E6", alpha = 0.5) +
  geom_line(data = data_df, aes(x = day, y =  lum.juv, group = unit, color = unit)) +
  geom_line(data = all_sims_summary, aes(x = day, y =  mean),size = 1) +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(x = "day", y = expression(paste("Juvenile ", italic("D. lumholtzi"), " density"))) +
  theme(axis.title = element_text(size = 10))+
  scale_x_continuous(limits = c(0, 52), breaks = c(0, 25, 50))


all_sims_summary <- all_sims %>%
  group_by(day) %>%
  summarise(lower = quantile(Jn, 0.025),   # 2.5% quantile
            upper = quantile(Jn, 0.975),
            mean =mean(Jn)) 

p6 = ggplot() +
  geom_ribbon(data = all_sims_summary, aes(x = day, ymin = lower, ymax = upper), fill = "#ADD8E6", alpha = 0.5) +
  geom_line(data = data_df, aes(x = day, y =  dent.juv, group = unit, color = unit)) +
  geom_line(data = all_sims_summary, aes(x = day, y =  mean),size = 1) +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(x = NULL, y = expression(paste("Juvenile ", italic("D. dentifera"), " density"))) +
  theme(axis.title = element_text(size = 10), axis.text.x = element_blank())+
  scale_x_continuous(limits = c(0, 52), breaks = c(0, 25, 50))





all_sims_summary <- all_sims %>%
  group_by(day) %>%
  summarise(lower = quantile(F, 0.1),   # 2.5% quantile
            upper =  quantile(F, 0.9),
            mean =  mean(F)) 
filtered_data <- all_sims %>%
  filter(.id == 1) %>%
  select(day, F, unit)


p7 = ggplot() +
  geom_ribbon(data = all_sims_summary, aes(x = day, ymin = lower, ymax = upper), fill = "#ADD8E6", alpha = 0.5) +
  geom_line(data = all_sims_summary, aes(x = day, y =  mean),size = 1) +
  geom_line(data = filtered_data, aes(x = day, y = F, group = unit, color = unit),linetype="dashed")+
  theme_minimal() +
  theme(legend.position = "none") +
  labs(x = NULL, y = "Food") +
  theme(axis.title = element_text(size = 10), axis.text.x = element_blank())+
  scale_x_continuous(limits = c(0, 52), breaks = c(0, 25, 50))


all_sims_summary <- all_sims %>%
  group_by(day) %>%
  summarise(lower = quantile(P, 0.025),   # 2.5% quantile
            upper = quantile(P, 0.975),
            mean =mean(P)) 
filtered_data <- all_sims %>%
  filter(.id == 1) %>%
  select(day, P, unit)
p8 = ggplot() +
  geom_ribbon(data = all_sims_summary, aes(x = day, ymin = lower, ymax = upper), fill = "#ADD8E6", alpha = 0.5) +
  geom_line(data = all_sims_summary, aes(x = day, y =  mean),size = 1) +
  geom_line(data = filtered_data, aes(x = day, y = P, group = unit, color = unit),linetype="dashed")+
  theme_minimal() +
  theme(legend.position = "none") +
  labs(x = "day", y = "Parasite") +
  theme(axis.title = element_text(size = 10))+
  scale_x_continuous(limits = c(0, 52), breaks = c(0, 25, 50))



grid.arrange( p1,p2,p6,p7,p3,p4,p5,p8,  nrow = 2, ncol = 4)



all_sims$unit <- factor(all_sims$unit,
                        levels = paste0("u", 1:8),
                        labels = paste0("Mesocosm-", 1:8))

data_df$unit <- factor(data_df$unit,
                       levels = paste0("u", 1:8),
                       labels = paste0("Mesocosm-", 1:8))


set.seed(0023)
selected_ids <- sample(unique(all_sims$.id), 20)
p1 = ggplot() +
  geom_line(data = all_sims[all_sims$.id %in% selected_ids, ], 
            aes(x = day, y = Sn, group = .id), col = 'blue') +
  geom_line(data = data_df, aes(x = day, y = dent.adult), col = 'red',linewidth = 1) +
  facet_grid(. ~ as.factor(unit), scales = "free_y", switch = "y") +
  labs(y = expression("Susceptible " * italic("D. dentifera")))+
  theme_bw() +
  scale_y_continuous(trans = "sqrt") +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.y = element_text(size = 8),
    axis.text.y = element_text(size = 9),
    strip.text = element_text(size = 9),
    strip.placement = "outside",
    plot.margin = margin(0, 0.5, 0, 0.5, "cm")
  )+
  scale_x_continuous(limits = c(0, 52), breaks = c(0, 25, 50))



p2 = ggplot() +
  geom_line(data = all_sims[all_sims$.id %in% selected_ids, ], 
            aes(x = day, y = Si, group = .id), col = 'blue') +
  geom_line(data = data_df, aes(x = day, y = lum.adult), col = 'red',linewidth = 1) +
  facet_grid(. ~ as.factor(unit), scales = "free_y", switch = "y") +
  labs(y = expression("Susceptible " * italic("D. lumholtzi")))+
  theme_bw() +
  scale_y_continuous(trans = "sqrt") +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.y = element_text(size = 8),
    axis.text.y = element_text(size = 9),
    strip.text = element_blank(),
    strip.placement = "outside",
    plot.margin = margin(0, 0.5, 0, 0.5, "cm")
  )+
  scale_x_continuous(limits = c(0, 52), breaks = c(0, 25, 50))


p3 = ggplot() +
  geom_line(data = all_sims[all_sims$.id %in% selected_ids, ], 
            aes(x = day, y = In, group = .id), col = 'blue') +
  geom_line(data = data_df, aes(x = day, y = dent.inf), col = 'red',linewidth = 1) +
  facet_grid(. ~ as.factor(unit), scales = "free_y", switch = "y") +
  labs(y = expression("Infected " * italic("D. dentifera") ))+
  theme_bw() +
  scale_y_continuous(trans = "sqrt") +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.y = element_text(size = 8),
    axis.text.y = element_text(size = 9),
    strip.text = element_blank(),
    strip.placement = "outside",
    strip.background = element_blank(),
    plot.margin = margin(0, 0.5, 0, 0.5, "cm")
  )+
  scale_x_continuous(limits = c(0, 52), breaks = c(0, 25, 50))


p4 = ggplot() +
  geom_line(data = all_sims[all_sims$.id %in% selected_ids, ], 
            aes(x = day, y = Ii, group = .id), col = 'blue') +
  geom_line(data = data_df, aes(x = day, y = lum.adult.inf), col = 'red',linewidth = 1) +
  facet_grid(. ~ as.factor(unit), scales = "free_y", switch = "y") +
  labs(y = expression("Infected " * italic("D. lumholtzi") ))+
  theme_bw() +
  scale_y_continuous(trans = "sqrt") +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.y = element_text(size = 8),
    axis.text.y = element_text(size = 9),
    strip.text = element_blank(),
    strip.placement = "outside",
    strip.background = element_blank(),
    plot.margin = margin(0, 0.5, 0, 0.5, "cm")
  )+
  scale_x_continuous(limits = c(0, 52), breaks = c(0, 25, 50))


p5 = ggplot() +
  geom_line(data = all_sims[all_sims$.id %in% selected_ids, ], 
            aes(x = day, y = Jn, group = .id), col = 'blue') +
  geom_line(data = data_df, aes(x = day, y = dent.juv), col = 'red',linewidth = 1) +
  facet_grid(. ~ as.factor(unit), scales = "free_y", switch = "y") +
  labs(y = expression("Juvenile " * italic("D. dentifera") ))+
  theme_bw() +
  scale_y_continuous(trans = "sqrt") +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.y = element_text(size = 8),
    axis.text.y = element_text(size = 9),
    strip.text = element_blank(),
    strip.placement = "outside",
    strip.background = element_blank(),
    plot.margin = margin(0, 0.5, 0, 0.5, "cm")
  )+
  scale_x_continuous(limits = c(0, 52), breaks = c(0, 25, 50))


p6 = ggplot() +
  geom_line(data = all_sims[all_sims$.id %in% selected_ids, ], 
            aes(x = day, y = Ji, group = .id), col = 'blue') +
  geom_line(data = data_df, aes(x = day, y = lum.juv), col = 'red',linewidth = 1) +
  facet_grid(. ~ as.factor(unit), scales = "free_y", switch = "y") +
  labs(y = expression("Juvenile " * italic("D. lumholtzi") ))+
  theme_bw() +
  scale_y_continuous(trans = "sqrt") +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.y = element_text(size = 8),
    axis.text.y = element_text(size = 9),
    strip.text = element_blank(),
    strip.placement = "outside",
    strip.background = element_blank(),
    plot.margin = margin(0, 0.5, 0, 0.5, "cm")
  )+
  scale_x_continuous(limits = c(0, 52), breaks = c(0, 25, 50))


p7 = ggplot() +
  geom_line(data = all_sims[all_sims$.id %in% selected_ids, ], 
            aes(x = day, y = F, group = .id), col = 'blue') +
  facet_grid(. ~ as.factor(unit), scales = "free_y", switch = "y") +
  labs(y = expression("Alga"))+
  theme_bw() +
  scale_y_continuous(trans = "sqrt") +
  theme(
    axis.title.y = element_text(size = 8),
    axis.text.y = element_text(size = 9),
    strip.text = element_blank(),
    strip.placement = "outside",
    strip.background = element_blank(),
    plot.margin = margin(0, 0.5, 0, 0.5, "cm")
  )+
  scale_x_continuous(limits = c(0, 52), breaks = c(0, 25, 50))



library(cowplot)
plot_grid(p1, p2, p3,p4,p5,p6,p7, ncol = 1, align = "v", axis = "l")
