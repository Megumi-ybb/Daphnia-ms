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

# Mesocosm_data = read_excel("~/Desktop/Research/D_P/Mesocosmdata.xlsx")
Mesocosm_data = read_excel("Mesocosmdata.xls")

sed = 0923
set.seed(0923)


dentNoPara <- Mesocosm_data[1:100, ]
dentNoPara <- subset(dentNoPara, select = c(rep, day, dent.adult))
dentNoPara <- dentNoPara[100: 1, ]
#Convert sampling date to natural date. Samples are collected every 5 days after the #first 7 days.
dentNoPara$day = (dentNoPara$day - 1) * 5 + 7
data = list()
trails = c("A","B","C","D","E","F","G","H","I","J")
for (i in 1: length(trails)){
  data[[i]]=subset(dentNoPara, select = c("day", "dent.adult"),
                   dentNoPara$rep == trails[i])
}


dyn_rpro <- Csnippet("
                      double Sn_term, F_term , Jn_term;
                      double noiSn  ,noiF,noiJn;
                      double delta = 0.013; //fraction of volume replaced day-1

                      noiSn = rnorm(0, sigSn * sqrt(dt));
                      noiJn = rnorm(0, sigJn * sqrt(dt));
                      noiF = rnorm(0, sigF * sqrt(dt));


                      //------------Sn-------------
                      Sn_term = 0.1 * Jn * dt - theta_Sn * Sn * dt -  delta * Sn * dt + Sn * noiSn;
                        
                      //-----------Jn-------------
                      
                      Jn_term = rn * f_Sn * F * Sn  * dt  -  0.1 * Jn * dt - theta_Jn * Jn * dt - delta * Jn * dt + Jn * noiJn;
                      
                      //-----------F---------------

                      F_term =  F * noiF - f_Sn * F * (Sn + 1 * Jn) * dt  - delta * F * dt + 0.37 * dt;


                      //T_D and T_I are the S and I Daphnia sample out
                      F += F_term;
                      Sn += Sn_term;
                      Jn += Jn_term;
                  
                      
                      
                      if (Sn < 0.0 || Sn > 1e5) {
                        Sn = 0.0;
                        error_count += 1;
                      }
                      if (F < 0.0 || F > 1e20) {
                        F = 0.0;
                        error_count += 1000;
                      }
                      if (Jn < 0.0 || Jn > 1e5) {
                        Jn = 0.0;
                        error_count += 0.001;
                      }
  
                      
                      T_Sn = fabs(Sn);
                      

                      ")

# Initial state. Assume t0 = day 4
dyn_init = Csnippet("
                     Sn = 3; 
                     F = 16.667;
                     Jn = 0;
                     
                     T_Sn = 0.0;
                     
                     error_count = 0.0;
                     
                     ")



dmeas = Csnippet("
                 
                 if (error_count > 0.0) {
                   lik = -150;
                  }
                 else{
                    if(give_log){
                    lik  = dnbinom_mu(dentadult,k_Sn,T_Sn,give_log);
                    }
                    else{
                    lik  = dnbinom_mu(dentadult,k_Sn,T_Sn,give_log);
                    }
                      }
                 ")

rmeas = Csnippet("
                 //double delta = 0.013;
                 dentadult = rnbinom_mu(k_Sn,T_Sn);
                 ")

pt <- parameter_trans(
  log = c( "sigSn","sigF","f_Sn",
           "rn","k_Sn","sigJn",
           "theta_Sn","theta_Jn")
)

pomplist=list()

parameters=  c( "sigSn" = 0.1, "sigF"= 0.1,"f_Sn"= 0.1,
                "rn"= 10,"k_Sn"= 0.1,"sigJn"= 0.1,
                "theta_Sn"= 0.1,"theta_Jn"= 0.1)


for (i in 1:10){
  colnames(data[[i]]) <- c('day', 'dentadult', 'lumadult')
  pomp(data = data[[i]],
       times = "day",
       t0=1,
       rprocess=euler(dyn_rpro,delta.t=1/4),
       rinit = dyn_init,
       dmeasure = dmeas,
       rmeasure = rmeas,
       partrans = pt,
       obsnames = c("dentadult"),
       accumvars = c("error_count"),
       paramnames = c( "sigSn","sigF","f_Sn",
                       "rn","k_Sn","sigJn",
                       "theta_Sn","theta_Jn"),
       statenames = c("Sn","Jn" ,"error_count", "F", "T_Sn")
  ) -> pomplist[[i]]
  coef(pomplist[[i]])=parameters
}
names(pomplist)=paste("u", 1:10,sep = "")



shared_parameter = c( "sigSn" = 0.1, "sigF"= 0.1,"f_Sn"= 0.1,
                      "rn"= 10,"k_Sn"= 0.1,"sigJn"= 0.1,
                      "theta_Sn"= 0.1,"theta_Jn"= 0.1)

panelfood = panelPomp(pomplist, shared=shared_parameter)


load("./Simple_dynamics/Dent/no_para/model/best_result.RData")

dentNoPara <- Mesocosm_data[1:100, ]
dentNoPara <- subset(dentNoPara, select = c(rep, day, dent.adult,dent.juv))
dentNoPara <- dentNoPara[100: 1, ]
dentNoPara$day = (dentNoPara$day - 1) * 5 + 7
data = list()
trails = c("A","B","C","D","E","F","G","H","I","J")
for (i in 1: length(trails)){
  data[[i]]=subset(dentNoPara, select = c("day", "dent.adult","dent.juv"),
                   dentNoPara$rep == trails[i])
}

lls <- matrix(unlist(sapply(mf, getElement, "ll")), nrow = 2)
best <- which.max(lls[1,])
mif.estimate <- coef(mf[[best]]$mif)
pf.loglik.of.mif.estimate <- unname(mf[[best]]$ll[1])
s.e.of.pf.loglik.of.mif.estimate <- unname(mf[[best]]$ll[2])

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


#Remove unrealistic scenario
# bad_ids <- unique(all_sims$.id[all_sims$Jn > 1000])
# all_sims <- subset(all_sims, !(.id %in% bad_ids))



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
  labs(x = "day", y = expression(paste("Susceptible ", italic("D. dentifera"))))  +
  scale_x_continuous(breaks = c(0, 25, 50)) +
  theme(axis.title = element_text(size = 10))

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
  labs(x = "day", y = expression(paste("Juvenile ", italic("D. dentifera"))))  +
  scale_x_continuous(breaks = c(0, 25, 50)) +
  theme(axis.title = element_text(size = 10))


all_sims_summary <- all_sims %>%
  group_by(day) %>%
  summarise(lower = quantile(F, 0.025),   # 2.5% quantile
            upper =  quantile(F, 0.975),
            mean =  mean(F)) 

filtered_data <- all_sims %>%
  filter(.id == 22) %>%
  select(day, F, unit)


p3 = ggplot() +
  geom_ribbon(data = all_sims_summary, aes(x = day, ymin = lower, ymax = upper), fill = "#ADD8E6", alpha = 0.5) +
  geom_line(data = all_sims_summary, aes(x = day, y = mean), size = 1) +
  geom_line(data = filtered_data, aes(x = day, y = F, group = unit, color = unit), linetype = "dashed") +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(x = "day", y = "Alga") +
  scale_x_continuous(breaks = c(0, 25, 50)) +
  theme(axis.title = element_text(size = 10))


grid.arrange( p1, p2,p3,  nrow = 1, ncol = 3)


bad_ids <- unique(all_sims$.id[all_sims$Sn > 600])
all_sims <- subset(all_sims, !(.id %in% bad_ids))


save_all_sims = all_sims
save_data_df = data_df



all_sims = save_all_sims
data_df = save_data_df

all_sims$unit <- factor(all_sims$unit,
                        levels = paste0("u", 1:10),
                        labels = paste0("Mesocosm-", 1:10))

data_df$unit <- factor(data_df$unit,
                        levels = paste0("u", 1:10),
                        labels = paste0("Mesocosm-", 1:10))


set.seed(0923)
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
            aes(x = day, y = Jn, group = .id), col = 'blue') +
  geom_line(data = data_df, aes(x = day, y = dent.juv), col = 'red',linewidth = 1) +
  facet_grid(. ~ as.factor(unit), scales = "free_y", switch = "y") +
  labs(y = expression("Juvenile " * italic("D. dentifera")))+
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



























