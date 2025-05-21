library(reshape2)
library(magrittr)
library(foreach)
library(readxl)
library(doParallel)
registerDoParallel(cores=36)
library(pomp)
library(panelPomp)
library(tidyverse)
library(gridExtra)

Mesocosm_data = read_excel("./Mesocosmdata.xls",2)

sed = 0923
set.seed(0923)


dentNoPara <- Mesocosm_data[1:100, ]
dentNoPara <- subset(dentNoPara, select = c(rep, day, lum.adult))
dentNoPara <- dentNoPara[100: 1, ]
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

parameters=  c( "sigSi" = 0.1, "sigF"= 0.1,"f_Si"= 0.1,
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



shared_parameter = c( "sigSi" = 0.1, "sigF"= 0.1,"f_Si"= 0.1,
                      "ri"= 10,"k_Si"= 0.1,"sigJi"= 0.1,
                      "theta_Si"= 0.1,"theta_Ji"= 0.1)

panelfood = panelPomp(pomplist, shared=shared_parameter)



load("./Single-species/Lum/SRJF/model/best_result.rda")

dentNoPara <- Mesocosm_data[1:100, ]
dentNoPara <- subset(dentNoPara, select = c(rep, day, lum.adult,lum.juv))
dentNoPara <- dentNoPara[100: 1, ]
dentNoPara$day = (dentNoPara$day - 1) * 5 + 7
data = list()
trails = c("A","B","C","D","E","F","G","H","I","J")
for (i in 1: length(trails)){
  data[[i]]=subset(dentNoPara, select = c("day", "lum.adult",'lum.juv'),
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
  labs(x = "day", y = expression(paste("Susceptible ", italic("D. lumholtzi"), ""))) +
  theme(axis.title = element_text(size = 10) )+
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
  labs(x = "day", y = expression(paste("Juvenile ", italic("D. lumholtzi"), ""))) +
  theme(axis.title = element_text(size = 10) )+
  scale_x_continuous(limits = c(0, 52), breaks = c(0, 25, 50))



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
  theme(axis.title = element_text(size = 10))+
  scale_x_continuous(limits = c(0, 52), breaks = c(0, 25, 50))



grid.arrange( p1, p2,p3,  nrow = 1, ncol = 3)




all_sims$unit <- factor(all_sims$unit,
                        levels = paste0("u", 1:10),
                        labels = paste0("Mesocosm-", 1:10))

data_df$unit <- factor(data_df$unit,
                       levels = paste0("u", 1:10),
                       labels = paste0("Mesocosm-", 1:10))


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
            aes(x = day, y = Ji, group = .id), col = 'blue') +
  geom_line(data = data_df, aes(x = day, y = lum.juv), col = 'red',linewidth = 1) +
  facet_grid(. ~ as.factor(unit), scales = "free_y", switch = "y") +
  labs(y = expression("Juvenile " * italic("D. lumholtzi")))+
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
plot_grid(p1, p2, p3, ncol = 1, align = "v", axis = "l")






























