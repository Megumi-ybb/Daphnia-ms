library(ggplot2)
library(gridExtra)
library(dplyr)
library(tidyr)

likelihood_df = data.frame(id = NA,
                           unit = NA,
                           time = NA,
                           cond.loglik = NA,
                           loglik_Sn = NA,
                           loglik_In = NA,
                           loglik_Si = NA,
                           loglik_Ii = NA)
temp = mf[[best]]$res
idx = best
for (inner_idx in 1:length(temp@unit_objects)){
  time_idx = 1
  for (time_idx in 1:10){
    # temp_weights = temp@unit.objects[[inner_idx]]@saved.states$weights[[time_idx]]
    tempstates = temp@unit_objects[[inner_idx]]@saved.states[[time_idx]]
    
    # tempstates = temp@unit_objects[[inner_idx]]@saved.states[[time_idx]]
    tempstates = as.data.frame(t(tempstates))
    temppara = as.data.frame(t(temp@shared))
    temp_data = temp@unit_objects[[inner_idx]]@data
    
    loglik_Sn = dnbinom(x = temp_data[1,time_idx], size = temppara$k_Sn, p = temppara$k_Sn/ (tempstates$T_Sn + temppara$k_Sn),log = TRUE)
    loglik_In = dnbinom(x = temp_data[2,time_idx], size = temppara$k_In, p = temppara$k_In/ (tempstates$T_In + temppara$k_In),log = TRUE)
    loglik_Si = dnbinom(x = temp_data[3,time_idx], size = temppara$k_Si, p = temppara$k_Si/ (tempstates$T_Si + temppara$k_Si),log = TRUE)
    loglik_Ii = dnbinom(x = temp_data[4,time_idx], size = temppara$k_Ii, p = temppara$k_Ii/ (tempstates$T_Ii + temppara$k_Ii),log = TRUE)
    
    # loglik_Sn + loglik_In + loglik_Si + loglik_Ii - temp_weights
    
    conditional_loglik = log(mean(exp(loglik_Sn + loglik_In + loglik_Si + loglik_Ii)))
    conditional_loglik_Sn = log(mean(exp(loglik_Sn)))
    conditional_loglik_In = log(mean(exp(loglik_In)))
    conditional_loglik_Si = log(mean(exp(loglik_Si)))
    conditional_loglik_Ii = log(mean(exp(loglik_Ii)))
    
    likelihood_df = rbind(likelihood_df,c(idx,inner_idx,time_idx,conditional_loglik,
                                          conditional_loglik_Sn,conditional_loglik_In,conditional_loglik_Si,conditional_loglik_Ii))   
  }
}
likelihood_df = likelihood_df[2:nrow(likelihood_df),]
gamma_object = temp
#######################xi = 0.001##############
likelihood_df0001 = data.frame(id = NA,
                               unit = NA,
                               time = NA,
                               cond.loglik = NA,
                               loglik_Sn = NA,
                               loglik_In = NA,
                               loglik_Si = NA,
                               loglik_Ii = NA)
temp = mf[[best]]$res
idx = best
for (inner_idx in 1:length(temp@unit_objects)){
  time_idx = 1
  for (time_idx in 1:10){
    # temp_weights = temp@unit.objects[[inner_idx]]@saved.states$weights[[time_idx]]
    tempstates = temp@unit_objects[[inner_idx]]@saved.states[[time_idx]]
    
    # tempstates = temp@unit_objects[[inner_idx]]@saved.states[[time_idx]]
    tempstates = as.data.frame(t(tempstates))
    temppara = as.data.frame(t(temp@shared))
    temp_data = temp@unit_objects[[inner_idx]]@data
    
    loglik_Sn = dnbinom(x = temp_data[1,time_idx], size = temppara$k_Sn, p = temppara$k_Sn/ (tempstates$T_Sn + temppara$k_Sn),log = TRUE)
    loglik_In = dnbinom(x = temp_data[2,time_idx], size = temppara$k_In, p = temppara$k_In/ (tempstates$T_In + temppara$k_In),log = TRUE)
    loglik_Si = dnbinom(x = temp_data[3,time_idx], size = temppara$k_Si, p = temppara$k_Si/ (tempstates$T_Si + temppara$k_Si),log = TRUE)
    loglik_Ii = dnbinom(x = temp_data[4,time_idx], size = temppara$k_Ii, p = temppara$k_Ii/ (tempstates$T_Ii + temppara$k_Ii),log = TRUE)
    
    # loglik_Sn + loglik_In + loglik_Si + loglik_Ii - temp_weights
    
    conditional_loglik = log(mean(exp(loglik_Sn + loglik_In + loglik_Si + loglik_Ii)))
    conditional_loglik_Sn = log(mean(exp(loglik_Sn)))
    conditional_loglik_In = log(mean(exp(loglik_In)))
    conditional_loglik_Si = log(mean(exp(loglik_Si)))
    conditional_loglik_Ii = log(mean(exp(loglik_Ii)))
    
    likelihood_df0001 = rbind(likelihood_df0001,c(idx,inner_idx,time_idx,conditional_loglik,
                                                  conditional_loglik_Sn,conditional_loglik_In,conditional_loglik_Si,conditional_loglik_Ii))   
  }
}
likelihood_df0001 = likelihood_df0001[2:nrow(likelihood_df0001),]
gaussian_object = temp

likelihood_df = rbind(likelihood_df,likelihood_df0001)
likelihood_df$model = ifelse(likelihood_df$id == 268,"Gamma","Gaussian")
# save(likelihood_df, file = 'data/residual_analysis_results.rda')

likelihood_df$time = 7 + (likelihood_df$time - 1)*5


p2 = ggplot(likelihood_df, aes(x = time, y = loglik_Sn, color = factor(id))) +  
  geom_line() +
  facet_wrap(~ unit, labeller = labeller(unit = function(x) paste0("Mesocosm-", x)), ncol = 8) +
  theme_bw() +
  labs(
       x = NULL,
       y = expression("Susecptible " * italic("D. dentifera")),
       color = "Legend") +  
  scale_color_manual(values = c("268" = "blue", "137" = "red"),  
                     labels = c("268" = "Gamma noise", "137" = "Gaussian noise")) + 
  theme(
    axis.title.y = element_text(size = 8),
    axis.text.y = element_text(size = 9),
    axis.text.x = element_blank(),
    strip.placement = "outside",
    strip.text = element_text(size = 9),
    plot.margin = margin(0, 0.5, 0, 0.5, "cm"),
    legend.position = "none"
  )+
  scale_x_continuous(limits = c(0, 52), breaks = c(0, 25, 52))

p3 = ggplot(likelihood_df, aes(x = time, y = loglik_In, color = factor(id))) +  
  geom_line() +
  facet_wrap(~ unit, labeller = labeller(unit = function(x) paste0("Mesocosm-", x)), ncol = 8) +
  theme_bw() +
  labs(
       x = NULL,
       y = expression("Infected " * italic("D. dentifera")),
       color = "Legend") +  
  scale_color_manual(values = c("268" = "blue", "137" = "red"),  
                     labels = c("268" = "xi_J = 1", "137" = "xi_J = 1000")) + 
  theme(
    axis.title.y = element_text(size = 8),
    axis.text.y = element_text(size = 9),
    axis.text.x = element_blank(),
    strip.placement = "outside",
    strip.text = element_blank(),
    plot.margin = margin(0, 0.5, 0, 0.5, "cm"),
    legend.position = "none"
  )+
  scale_x_continuous(limits = c(0, 52), breaks = c(0, 25, 52))

p4 = ggplot(likelihood_df, aes(x = time, y = loglik_Si, color = factor(id))) +  
  geom_line() +
  facet_wrap(~ unit, labeller = labeller(unit = function(x) paste0("Mesocosm-", x)), ncol = 8) +
  theme_bw() + 
  labs(
       x = NULL,
       y = expression("Susecptible " * italic("D. lumholtzi")),
       color = "Legend") +  
  scale_color_manual(values = c("268" = "blue", "137" = "red"),  
                     labels = c("268" = "xi_J = 1", "137" = "xi_J = 1000")) + 
  theme(
    axis.title.y = element_text(size = 8),
    axis.text.y = element_text(size = 9),
    axis.text.x = element_blank(),
    strip.placement = "outside",
    strip.text = element_blank(),
    plot.margin = margin(0, 0.5, 0, 0.5, "cm"),
    legend.position = "none"
  )+
  scale_x_continuous(limits = c(0, 52), breaks = c(0, 25, 52))

p5 = ggplot(likelihood_df, aes(x = time, y = loglik_Ii, color = factor(id))) +  
  geom_line() +
  facet_wrap(~ unit, labeller = labeller(unit = function(x) paste0("Mesocosm-", x)), ncol = 8) +
  theme_bw() +
  labs(title = ,
       x = "Time",
       y = expression("Infected " * italic("D. lumholtzi") ),
       color = "Legend") +  
  scale_color_manual(values = c("268" = "blue", "137" = "red"),  
                     labels = c("268" = "xi_J = 1", "137" = "xi_J = 1000")) + 
  theme(
    axis.title.y = element_text(size = 8),
    axis.text.y = element_text(size = 9),
    axis.text.x = element_text(size = 9),
    strip.placement = "outside",
    strip.text = element_blank(),
    plot.margin = margin(0, 0.5, 0, 0.5, "cm"),
    legend.position = "none"
  )+
  scale_x_continuous(limits = c(0, 52), breaks = c(0, 25, 52))


library(cowplot)
plot_grid(p2, p3, p4,p5,ncol = 1, align = "v", axis = "l")




#Gaussian - Gamma
difference = likelihood_df[likelihood_df$id == 137,][,1:ncol(likelihood_df)-1] - likelihood_df[likelihood_df$id == 268,][,1:ncol(likelihood_df)-1]
difference$unit = likelihood_df[likelihood_df$id == 137,]$unit
difference$time = likelihood_df[likelihood_df$id == 137,]$time

time_mean <- difference %>%
  group_by(time) %>%
  summarise(mean_cond_loglik = mean(cond.loglik, na.rm = TRUE))

p1 <- ggplot(difference, aes(x = time, y = cond.loglik, color = factor(unit))) +
  geom_line()+
  geom_line(data = time_mean, aes(x = time, y = mean_cond_loglik), linetype = "dashed", color = "black", inherit.aes = FALSE)+
  labs(title = "cond.loglik", x = "Time", y = "Value") +
  theme_minimal()+
  theme(legend.position = "none") 

time_mean <- difference %>%
  group_by(time) %>%
  summarise(mean_cond_loglik = mean(loglik_Sn, na.rm = TRUE))

p2 <- ggplot(difference, aes(x = time, y = loglik_Sn, color = factor(unit))) +
  geom_line() +
  geom_line(data = time_mean, aes(x = time, y = mean_cond_loglik), linetype = "dashed", color = "black", inherit.aes = FALSE)+
  labs(title = "loglik_Sn", x = "Time", y = "Value") +
  theme_minimal()+
  theme(legend.position = "none") 


time_mean <- difference %>%
  group_by(time) %>%
  summarise(mean_cond_loglik = mean(loglik_In, na.rm = TRUE))

p3 <- ggplot(difference, aes(x = time, y = loglik_In, color = factor(unit))) +
  geom_line() +
  geom_line(data = time_mean, aes(x = time, y = mean_cond_loglik), linetype = "dashed", color = "black", inherit.aes = FALSE)+
  labs(title = "loglik_In", x = "Time", y = "Value") +
  theme_minimal() +
  theme(legend.position = "none") 

time_mean <- difference %>%
  group_by(time) %>%
  summarise(mean_cond_loglik = mean(loglik_Si, na.rm = TRUE))

p4 <- ggplot(difference, aes(x = time, y = loglik_Si, color = factor(unit))) +
  geom_line() +
  geom_line(data = time_mean, aes(x = time, y = mean_cond_loglik), linetype = "dashed", color = "black", inherit.aes = FALSE)+
  labs(title = "loglik_Si", x = "Time", y = "Value") +
  theme_minimal()+
  theme(legend.position = "none") 

time_mean <- difference %>%
  group_by(time) %>%
  summarise(mean_cond_loglik = mean(loglik_Ii, na.rm = TRUE))

p5 <- ggplot(difference, aes(x = time, y = loglik_Ii, color = factor(unit))) +
  geom_line() +
  geom_line(data = time_mean, aes(x = time, y = mean_cond_loglik), linetype = "dashed", color = "black", inherit.aes = FALSE)+
  labs(title = "loglik_Ii", x = "Time", y = "Value") +
  theme_minimal()+
  theme(legend.position = "none") 

grid.arrange(p1, p2, p3, p4,p5, ncol = 3, nrow = 2)



######################Plot of filter mean######################
library(ggplot2)
library(gridExtra)
library(reshape2)
library(magrittr)
library(readxl)
library(pomp)
library(panelPomp)
library(tidyverse)

Mesocosm_data = read_excel("./Mesocosmdata.xls",3)

sed = 0923
set.seed(0416)


dentNoPara <- Mesocosm_data[91:170, ]
dentNoPara <- subset(dentNoPara, select = c(rep, day, dent.adult,dent.inf,lum.adult,lum.adult.inf, dent.juv, lum.juv))
dentNoPara <- dentNoPara[80: 1, ]
#Convert sampling date to natural date. Samples are collected every 5 days after the #first 7 days.
dentNoPara$day = (dentNoPara$day - 1) * 5 + 7
data = list()
trails = c("K","L","M","N","O","P","Q","S")
for (i in 1: length(trails)){
  data[[i]]=subset(dentNoPara, select = c("day", "dent.adult","dent.inf","lum.adult","lum.adult.inf","dent.juv","lum.juv"),
                   dentNoPara$rep == trails[i])
}





plot_observed_vs_two_models <- function(
    data, 
    gamma_object,
    gaussian_object,
    observed_col,
    gamma_modeled_name,
    gaussian_modeled_name,
    observed_label,
    gamma_label,
    gaussian_label
) {
  df_list <- vector("list", length(data))
  
  for(i in seq_along(data)) {
    tmp_df <- data.frame(
      day             = data[[i]]$day,
      observed        = data[[i]][[observed_col]],
      gamma_modeled   = gamma_object@unit_objects[[i]]@filter.mean[gamma_modeled_name, ],
      gaussian_modeled= gaussian_object@unit_objects[[i]]@filter.mean[gaussian_modeled_name, ]
    )
    
    tmp_long <- pivot_longer(
      tmp_df,
      cols      = c("observed", "gamma_modeled", "gaussian_modeled"),
      names_to  = "type",
      values_to = "value"
    )
    tmp_long$replicate <- i
    
    df_list[[i]] <- tmp_long
  }
  
  big_df <- do.call(rbind, df_list)
  
  ggplot(big_df, aes(x = day, y = value, color = type)) +
    geom_line() +
    geom_point() +
    facet_wrap(~ replicate, ncol = 4, scales = "fixed") + 
    scale_color_manual(
      name   = "Type",
      values = c("observed"         = "black", 
                 "gamma_modeled"    = "red",
                 "gaussian_modeled" = "blue"),
      labels = c("observed"         = observed_label,
                 "gamma_modeled"    = gamma_label,
                 "gaussian_modeled" = gaussian_label)
    ) +
    labs(x = "Day", y = observed_col) +
    theme_bw() +
    ggtitle(paste(observed_col, "vs", gamma_modeled_name, "and", gaussian_modeled_name))
}


p1 <- plot_observed_vs_two_models(
  data                 = data,
  gamma_object         = gamma_object,
  gaussian_object      = gaussian_object,
  observed_col         = "dent.adult",
  gamma_modeled_name   = "Sn",
  gaussian_modeled_name= "Sn",  
  observed_label       = "Observed Dent Adult",
  gamma_label          = "Gamma-modeled Dent Adult",
  gaussian_label       = "Gaussian-modeled Dent Adult"
)


p2 <- plot_observed_vs_two_models(
  data                 = data,
  gamma_object         = gamma_object,
  gaussian_object      = gaussian_object,
  observed_col         = "lum.adult",
  gamma_modeled_name   = "Si",
  gaussian_modeled_name= "Si", 
  observed_label       = "Observed Lum Adult",
  gamma_label          = "Gamma-modeled Lum Adult",
  gaussian_label       = "Gaussian-modeled Lum Adult"
)


p3 <- plot_observed_vs_two_models(
  data                 = data,
  gamma_object         = gamma_object,
  gaussian_object      = gaussian_object,
  observed_col         = "lum.adult.inf",
  gamma_modeled_name   = "Ii",
  gaussian_modeled_name= "Ii", 
  observed_label       = "Observed Lum Inf",
  gamma_label          = "Gamma-modeled Lum Inf",
  gaussian_label       = "Gaussian-modeled Lum Inf"
)

p4 <- plot_observed_vs_two_models(
  data                 = data,
  gamma_object         = gamma_object,
  gaussian_object      = gaussian_object,
  observed_col         = "dent.inf",
  gamma_modeled_name   = "In",
  gaussian_modeled_name= "In", 
  observed_label       = "Observed Dent Inf",
  gamma_label          = "Gamma-modeled Dent Inf",
  gaussian_label       = "Gaussian-modeled Dent Inf"
)

p1
p2
p3
p4











