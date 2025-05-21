library(reshape2)
library(magrittr)
library(foreach)
library(readxl)
library(gridExtra)
library(doParallel)
registerDoParallel(cores=36)
library(pomp)
library(panelPomp)
library(tidyverse)

Mesocosm_data = read_excel("./Mesocosmdata.xls",2)

sed = 0923
set.seed(0923)

dentNoPara <- Mesocosm_data[101:190, ]
dentNoPara <- subset(dentNoPara, select = c(rep, day, lum.adult,lum.adult.inf))
dentNoPara <- dentNoPara[90: 1, ]
#Convert sampling date to natural date. Samples are collected every 5 days after the #first 7 days.
dentNoPara$day = (dentNoPara$day - 1) * 5 + 7
data = list()
trails = c("L","M","N","O","P","Q","R","S","T")
for (i in 1: length(trails)){
  data[[i]]=subset(dentNoPara, select = c("day", "lum.adult","lum.adult.inf"),
                   dentNoPara$rep == trails[i])
}


#Notice that the unit of S is individuals/L, unit of F is 1e+5 cells/L, unit of P is 1000 spores/L, unit for I is individuals/L

dyn_rpro <- Csnippet("
                      double Si_term, Ii_term, F_term, P_term ,Ji_term;
                      double noiSi, noiIi,noiF, noiP,noiJi;
                      double delta = 0.013; //fraction of volume replaced day-1

                      noiSi = rnorm(0, sigSi * sqrt(dt));
                      noiIi = rnorm(0, sigIi * sqrt(dt));
                      noiJi = rnorm(0, sigJi * sqrt(dt));
                      noiF = rnorm(0, sigF * sqrt(dt));
                      noiP = rnorm(0, sigP * sqrt(dt));


                      //------------Si-------------
                      Si_term = 0.1 * Ji * dt - theta_Si * Si * dt -  probi * f_Si * Si * P * dt - delta * Si * dt + Si * noiSi;
                        
                      //-----------Ji-------------
                      
                      Ji_term = ri * f_Si * F * Si  * dt  -  0.1 * Ji * dt - theta_Ji * Ji * dt - delta * Ji * dt + Ji * noiJi;
                      
                      //------------Ii--------------

                      Ii_term = probi * f_Si * Si * P * dt - theta_Ii * Ii *dt - delta * Ii * dt + Ii * noiIi;

                      //-----------F---------------

                      F_term =  F * noiF - f_Si * F * (Si + xi * Ii + 1 * Ji) * dt - delta * F * dt + 0.37 * dt;


                      //----------P---------------

                      P_term = 30 * theta_Ii * Ii * dt - f_Si * (Si + xi * Ii) * P * dt - theta_P * P * dt - delta * P * dt + P * noiP;


                      //T_D and T_I are the S and I Daphnia sample out
                      F += F_term;
                      Si += Si_term;
                      Ii += Ii_term;
                      Ji += Ji_term;
                      P += P_term;
                      
                      //Trace time
                      
                     if (t - 4.0 < 0.001 && t - 4.0 > -0.001){
                      //Initial statement
                        P += 25;
                      }
                      
                      
                      if (Si < 0.0 || Si > 1e5) {
                        Si = 0.0;
                        error_count += 1;
                      }
                      if (F < 0.0 || F > 1e20) {
                        F = 0.0;
                        error_count += 1000;
                      }
                      if (Ii < 0.0 || Ii > 1e5) {
                        Ii = 0.0;
                        error_count += 0.001;
                      }
                      if (Ji < 0.0 || Ji > 1e5) {
                        Ji = 0.0;
                        error_count += 0.001;
                      }
                      if ((P < 0.0 || P > 1e20)&& t > 3.9) {
                        P = 0.0;
                        error_count += 0.000001;
                      }
                      
                      T_Si = fabs(Si);
                      T_Ii = fabs(Ii);
                      

                      ")

# Initial state. Assume t0 = day 4
dyn_init = Csnippet("
                     Si = 3;
                     F = 16.667;
                     Ji = 0;
                     
                     T_Si = 0.0;
                     T_Ii = 0.0;
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
                    lik  = dnbinom_mu(dentadult,k_Si,T_Si,give_log) +  dnbinom_mu(dentinf,k_Ii,T_Ii,give_log);
                    }
                    else{
                    lik  = dnbinom_mu(dentadult,k_Si,T_Si,give_log) *  dnbinom_mu(dentinf,k_Ii,T_Ii,give_log);
                    }
                      }
                 ")

rmeas = Csnippet("
                 //double delta = 0.013;
                 dentadult = rnbinom_mu(k_Si,T_Si);
                 dentinf = rnbinom_mu(k_Ii,T_Ii);
                 ")

pt <- parameter_trans(
  #Without death part
  log = c( "sigSi", "sigIi", "sigF","sigP","f_Si",
           "ri","k_Ii","k_Si","sigJi","probi",
           "theta_Si","theta_Ii","theta_P","theta_Ji","xi")
)

pomplist=list()

parameters =  c( "sigSi" = 0.1, "sigIi"= 0.1, "sigF"= 0.1,"sigP"= 0.1,"f_Si"= 0.01,
                 "ri"= 10,"k_Ii"= 0.1,"k_Si"= 0.1,"sigJi"= 0.1,"probi"= 0.1,
                 "theta_Si"= 0.1,"theta_Ii"= 0.1,"theta_P"= 0.1,"theta_Ji"= 0.1,"xi"= 1)


for (i in 1:9){
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
       paramnames = c( "sigSi","sigIi","sigF","sigP","f_Si","ri","k_Ii","k_Si","sigJi","probi",
                       "theta_Si","theta_Ii","theta_P","theta_Ji","xi"),
       statenames = c("Si",  "Ii","Ji" ,"error_count", "F", "T_Si","T_Ii","P")
  ) -> pomplist[[i]]
  coef(pomplist[[i]])=parameters
}
names(pomplist)=paste("u", 1:9,sep = "")



shared_parameter = c( "sigSi" = 0.1, "sigIi"= 0.1, "sigF"= 0.1,"sigP"= 0.1,"f_Si"= 0.01,
                      "ri"= 10,"k_Ii"= 0.1,"k_Si"= 0.1,"sigJi"= 0.1,"probi"= 0.1,
                      "theta_Si"= 0.1,"theta_Ii"= 0.1,"theta_P"= 0.1,"theta_Ji"= 0.1,"xi"= 1)

panelfood = panelPomp(pomplist, shared=shared_parameter)


load("./Single-species/Lum/SIRJPF/model/best_result.rda")

dentNoPara <- Mesocosm_data[101:190, ]
dentNoPara <- subset(dentNoPara, select = c(rep, day, lum.adult,lum.adult.inf,lum.juv))
dentNoPara <- dentNoPara[90: 1, ]
#Convert sampling date to natural date. Samples are collected every 5 days after the #first 7 days.
dentNoPara$day = (dentNoPara$day - 1) * 5 + 7
data = list()
trails = c("L","M","N","O","P","Q","R","S","T")
for (i in 1: length(trails)){
  data[[i]]=subset(dentNoPara, select = c("day", "lum.adult","lum.adult.inf","lum.juv"),
                   dentNoPara$rep == trails[i])
}


coef(panelfood) <- mif.estimate

foreach(u = names(panelfood), .combine = rbind) %do% {
  unit_model <- unit_objects(panelfood)[[u]]
  shared <- mif.estimate
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
  summarise(lower = quantile(Si, 0.025),   # 2.5% quantile
            upper = quantile(Si, 0.975),
            mean = mean(Si)) 

p1 = ggplot() +
  geom_ribbon(data = all_sims_summary, aes(x = day, ymin = lower, ymax = upper), fill = "#ADD8E6", alpha = 0.5) +
  geom_line(data = data_df, aes(x = day, y =  lum.adult, group = unit, color = unit)) +
  geom_line(data = all_sims_summary, aes(x = day, y =  mean),size = 1) +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(x = NULL, y = expression(paste("Susceptible ", italic("D. lumholtzi")))) +
  theme(axis.title = element_text(size = 10), axis.text.x = element_blank())+
  scale_x_continuous(limits = c(0, 52), breaks = c(0, 25, 50))


all_sims_summary <- all_sims %>%
  group_by(day) %>%
  summarise(lower = quantile(Ji, 0.025),   # 2.5% quantile
            upper = quantile(Ji, 0.975),
            mean =mean(Ji)) 

p2 = ggplot() +
  geom_ribbon(data = all_sims_summary, aes(x = day, ymin = lower, ymax = upper), fill = "#ADD8E6", alpha = 0.5) +
  geom_line(data = data_df, aes(x = day, y =  lum.juv, group = unit, color = unit)) +
  geom_line(data = all_sims_summary, aes(x = day, y =  mean),size = 1) +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(x = "day", y = expression(paste("Juvenile ", italic("D. lumholtzi")))) +
  theme(axis.title = element_text(size = 10))+
  scale_x_continuous(limits = c(0, 52), breaks = c(0, 25, 50))



all_sims_summary <- all_sims %>%
  group_by(day) %>%
  summarise(lower = quantile(F, 0.025),   # 2.5% quantile
            upper =  quantile(F, 0.975),
            mean =  mean(F)) 

filtered_data <- all_sims %>%
  filter(.id == 3) %>%
  select(day, F, unit)


p3 = ggplot() +
  geom_ribbon(data = all_sims_summary, aes(x = day, ymin = lower, ymax = upper), fill = "#ADD8E6", alpha = 0.5) +
  geom_line(data = all_sims_summary, aes(x = day, y =  mean),size = 1) +
  geom_line(data = filtered_data, aes(x = day, y = F, group = unit, color = unit),linetype="dashed")+
  theme_minimal() +
  theme(legend.position = "none") +
  labs(x = "day", y = "Alga") +
  theme(axis.title = element_text(size = 10))+
  scale_x_continuous(limits = c(0, 52), breaks = c(0, 25, 50))



all_sims_summary <- all_sims %>%
  group_by(day) %>%
  summarise(lower = quantile(Ii, 0.025),   # 2.5% quantile
            upper = quantile(Ii, 0.975),
            mean =mean(Ii)) 

p4 = ggplot() +
  geom_ribbon(data = all_sims_summary, aes(x = day, ymin = lower, ymax = upper), fill = "#ADD8E6", alpha = 0.5) +
  geom_line(data = data_df, aes(x = day, y =  lum.adult.inf, group = unit, color = unit)) +
  geom_line(data = all_sims_summary, aes(x = day, y =  mean),size = 1) +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(x = NULL, y = expression(paste("Infected ", italic("D. lumholtzi")))) +
  theme(axis.title = element_text(size = 10), axis.text.x = element_blank())+
  scale_x_continuous(limits = c(0, 52), breaks = c(0, 25, 50))



grid.arrange( p1, p4, p2,p3,  nrow = 2, ncol = 2)


all_sims$unit <- factor(all_sims$unit,
                        levels = paste0("u", 1:9),
                        labels = paste0("Mesocosm-", 1:9))

data_df$unit <- factor(data_df$unit,
                       levels = paste0("u", 1:9),
                       labels = paste0("Mesocosm-", 1:9))

set.seed(0416)
selected_ids <- sample(unique(all_sims$.id), 20)
p1 = ggplot() +
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
    axis.title.y = element_text(size = 9),
    axis.text.y = element_text(size = 9),
    strip.text = element_text(size = 9),
    strip.placement = "outside",
    plot.margin = margin(0, 0.5, 0, 0.5, "cm")
  )+
  scale_x_continuous(limits = c(0, 52), breaks = c(0, 25, 50))


p2 = ggplot() +
  geom_line(data = all_sims[all_sims$.id %in% selected_ids, ], 
            aes(x = day, y = Ii, group = .id), col = 'blue') +
  geom_line(data = data_df, aes(x = day, y = lum.adult.inf), col = 'red',linewidth = 1) +
  facet_grid(. ~ as.factor(unit), scales = "free_y", switch = "y") +
  labs(y = expression("Infected " * italic("D. lumholtzi")))+
  theme_bw() +
  scale_y_continuous(trans = "sqrt") +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.y = element_text(size = 9),
    axis.text.y = element_text(size = 9),
    strip.text = element_blank(),
    strip.placement = "outside",
    strip.background = element_blank(),
    plot.margin = margin(0, 0.5, 0, 0.5, "cm")
  )+
  scale_x_continuous(limits = c(0, 52), breaks = c(0, 25, 50))


p3 = ggplot() +
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
    axis.title.y = element_text(size = 9),
    axis.text.y = element_text(size = 9),
    strip.text = element_blank(),
    strip.placement = "outside",
    strip.background = element_blank(),
    plot.margin = margin(0, 0.5, 0, 0.5, "cm")
  )+
  scale_x_continuous(limits = c(0, 52), breaks = c(0, 25, 50))




p4 = ggplot() +
  geom_line(data = all_sims[all_sims$.id %in% selected_ids, ], 
            aes(x = day, y = F, group = .id), col = 'blue') +
  facet_grid(. ~ as.factor(unit), scales = "free_y", switch = "y") +
  labs(y = expression("Alga"))+
  theme_bw() +
  scale_y_continuous(trans = "sqrt") +
  theme(
    axis.title.y = element_text(size = 9),
    axis.text.y = element_text(size = 9),
    strip.text = element_blank(),
    strip.placement = "outside",
    strip.background = element_blank(),
    plot.margin = margin(0, 0.5, 0, 0.5, "cm")
  )+
  scale_x_continuous(limits = c(0, 52), breaks = c(0, 25, 50))



library(cowplot)
plot_grid(p1, p2, p3,p4, ncol = 1, align = "v", axis = "l")





























