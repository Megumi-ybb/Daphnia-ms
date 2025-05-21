library(reshape2)
library(magrittr)
library(foreach)
library(gridExtra)
library(readxl)
library(doParallel)
registerDoParallel(cores=36)
library(pomp)
library(panelPomp)
library(tidyverse)

Mesocosm_data = read_excel("Mesocosmdata.xls")

sed = 0923
set.seed(0923)

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



shared_parameter = c( "sigSn" = 0.1, "sigIn"= 0.1, "sigF"= 0.1,"sigP"= 0.1,"f_Sn"= 0.01,
                      "rn"= 10,"k_In"= 0.1,"k_Sn"= 0.1,"sigJn"= 0.1,"probn"= 0.1,
                      "theta_Sn"= 0.1,"theta_In"= 0.1,"theta_P"= 0.1,"theta_Jn"= 0.1,"xi"= 1)

panelfood = panelPomp(pomplist, shared=shared_parameter)

load("./Single-species/Dent/SIRJPF/model/best_result.rda")

dentNoPara <- Mesocosm_data[101:180, ]
dentNoPara <- subset(dentNoPara, select = c(rep, day, dent.adult,dent.inf,dent.juv ))
dentNoPara <- dentNoPara[80: 1, ]
dentNoPara$day = (dentNoPara$day - 1) * 5 + 7
data = list()
trails = c("K","L","M","N","O","P","Q","S")
for (i in 1: length(trails)){
  data[[i]]=subset(dentNoPara, select = c("day", "dent.adult","dent.inf","dent.juv"),
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
  summarise(lower = quantile(Sn, 0.025),   # 2.5% quantile
            upper = quantile(Sn, 0.975),
            mean = mean(Sn)) 

p1 = ggplot() +
  geom_ribbon(data = all_sims_summary, aes(x = day, ymin = lower, ymax = upper), fill = "#ADD8E6", alpha = 0.5) +
  geom_line(data = data_df, aes(x = day, y =  dent.adult, group = unit, color = unit)) +
  geom_line(data = all_sims_summary, aes(x = day, y =  mean),size = 1) +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(x = NULL, y = expression(paste("Susceptible ", italic("D. dentifera")))) +
  scale_x_continuous(breaks = c(0, 25, 50)) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank()
  )

all_sims_summary <- all_sims %>%
  group_by(day) %>%
  summarise(lower = quantile(Jn, 0.025),   # 2.5% quantile
            upper = quantile(Jn, 0.975),
            mean =mean(Jn)) 

p2 = ggplot() +
  geom_ribbon(data = all_sims_summary, aes(x = day, ymin = lower, ymax = upper), fill = "#ADD8E6", alpha = 0.5) +
  geom_line(data = data_df, aes(x = day, y =  dent.juv, group = unit, color = unit)) +
  geom_line(data = all_sims_summary, aes(x = day, y =  mean),size = 1) +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(x = "day", y = expression(paste("Juvenile ", italic("D. dentifera")))) +
  scale_x_continuous(breaks = c(0, 25, 50)) +
  theme(axis.title = element_text(size = 10))


all_sims_summary <- all_sims %>%
  group_by(day) %>%
  summarise(lower = quantile(F, 0.025),   # 2.5% quantile
            upper =  quantile(F, 0.975),
            mean =  mean(F)) 

filtered_data <- all_sims %>%
  filter(.id == 2) %>%
  select(day, F, unit)


p3 = ggplot() +
  geom_ribbon(data = all_sims_summary, aes(x = day, ymin = lower, ymax = upper), fill = "#ADD8E6", alpha = 0.5) +
  geom_line(data = all_sims_summary, aes(x = day, y =  mean),size = 1) +
  geom_line(data = filtered_data, aes(x = day, y = F, group = unit, color = unit),linetype="dashed")+
  theme_minimal() +
  theme(legend.position = "none") +
  labs(x = "day", y = "Alga") +
  scale_x_continuous(breaks = c(0, 25, 50)) +
  theme(axis.title = element_text(size = 10))


all_sims_summary <- all_sims %>%
  group_by(day) %>%
  summarise(lower = quantile(In, 0.025),   # 2.5% quantile
            upper = quantile(In, 0.975),
            mean =mean(In)) 

p4 = ggplot() +
  geom_ribbon(data = all_sims_summary, aes(x = day, ymin = lower, ymax = upper), fill = "#ADD8E6", alpha = 0.5) +
  geom_line(data = data_df, aes(x = day, y =  dent.inf, group = unit, color = unit)) +
  geom_line(data = all_sims_summary, aes(x = day, y =  mean),size = 1) +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(x = NULL, y = expression(paste("Infected ", italic("D. dentifera")))) +
  scale_x_continuous(breaks = c(0, 25, 50)) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank()
  )


grid.arrange( p1, p4, p2,p3,  nrow = 2, ncol = 2)




all_sims$unit <- factor(all_sims$unit,
                        levels = paste0("u", 1:8),
                        labels = paste0("Mesocosm-", 1:8))

data_df$unit <- factor(data_df$unit,
                       levels = paste0("u", 1:8),
                       labels = paste0("Mesocosm-", 1:8))

set.seed(0416)
selected_ids <- sample(unique(all_sims$.id), 20)
p1 = ggplot() +
  geom_line(data = all_sims[all_sims$.id %in% selected_ids, ], 
            aes(x = day, y = Sn, group = .id), col = 'blue') +
  geom_line(data = data_df, aes(x = day, y = dent.adult), col = 'red',linewidth = 1) +
  facet_grid(. ~ as.factor(unit), scales = "free_y", switch = "y") +
  labs(y = expression("Susceptible " * italic("D. dentifera")))+
  theme_bw() +
  scale_y_continuous(trans = "sqrt", limits = c(0,600)) +
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
            aes(x = day, y = In, group = .id), col = 'blue') +
  geom_line(data = data_df, aes(x = day, y = dent.inf), col = 'red',linewidth = 1) +
  facet_grid(. ~ as.factor(unit), scales = "free_y", switch = "y") +
  labs(y = expression("Infected " * italic("D. dentifera") ))+
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
            aes(x = day, y = Jn, group = .id), col = 'blue') +
  geom_line(data = data_df, aes(x = day, y = dent.juv), col = 'red',linewidth = 1) +
  facet_grid(. ~ as.factor(unit), scales = "free_y", switch = "y") +
  labs(y = expression("Juvenile " * italic("D. dentifera") ))+
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































