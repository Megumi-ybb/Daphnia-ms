library(dplyr)
library(latex2exp)
library(gridExtra)
library(dplyr)
library(pomp)
library(ggplot2)

#Please set the working to be the 'Daphnia-ms' path

load("data/Simple_dynamics/Lum/no_para/profile_graph_data.rda")
load('Single-species/Lum/SRJF/model/best_result.rda')

load_option = FALSE

if(load_option){
  load("./Single-species/Lum/SRJF/profile/f_Si/f_Si.RData")
  subset_data_f_Si <- final_params %>%
    group_by(f_Si) %>%
    filter(loglik == max(loglik))
}

plot(x = log(subset_data_f_Si$f_Si), y = subset_data_f_Si$loglik)
plot(x = subset_data_f_Si$f_Si, y = subset_data_f_Si$loglik)
subset_data_f_Si$log_f_Si <- log(subset_data_f_Si$f_Si)

mcap(subset_data_f_Si$loglik, subset_data_f_Si$log_f_Si,  level = 0.95, span = 0.95, Ngrid = 1000) -> mcap_object_f_Si
mcap_object_f_Si$mle -> f_Si_mle
f_Si_p <- ggplot() +
  geom_point(data = subset_data_f_Si, aes(x = log_f_Si, y = loglik)) +
  geom_line(data = mcap_object_f_Si$fit, aes(x = parameter, y = smoothed), col = 'red') +
  # geom_vline(xintercept = mcap_object_f_Si$ci[1], linetype = 'dashed') +
  geom_vline(xintercept = mcap_object_f_Si$ci[2], linetype = 'dashed') +
  geom_vline(xintercept = mcap_object_f_Si$mle, col = 'blue') +
  geom_vline(xintercept = log(mif.estimate[['f_Si']]), col = 'red') +
  geom_hline(yintercept = pf.loglik.of.mif.estimate, col = 'red',linetype = "longdash") +
  geom_point(aes(x = log(mif.estimate[['f_Si']]), 
                 y = pf.loglik.of.mif.estimate), 
             color = "red", size = 2) + 
  labs(x =  TeX("$\\log(f^i_{S})$"), y = "log likelihood") +
  theme(axis.text = element_text(size = 11),
        axis.title = element_text(size = 11)) +
  ylim(-375, -370)+
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
f_Si_p




if(load_option){
  load("./Single-species/Lum/SRJF/profile/k_Si/k_Si.RData")
  subset_data_k_Si <- final_params %>%
    group_by(k_Si) %>%
    filter(loglik == max(loglik))
}

plot(x = log(subset_data_k_Si$k_Si), y = subset_data_k_Si$loglik)
plot(x = subset_data_k_Si$k_Si, y = subset_data_k_Si$loglik)
subset_data_k_Si$log_k_Si <- log(subset_data_k_Si$k_Si)

mcap(subset_data_k_Si$loglik, subset_data_k_Si$log_k_Si,  level = 0.95, span = 0.95, Ngrid = 1000) -> mcap_object_k_Si
mcap_object_k_Si$mle -> k_Si_mle
k_Si_p <- ggplot() +
  geom_point(data = subset_data_k_Si, aes(x = log_k_Si, y = loglik)) +
  geom_line(data = mcap_object_k_Si$fit, aes(x = parameter, y = smoothed), col = 'red') +
  geom_vline(xintercept = mcap_object_k_Si$ci[1], linetype = 'dashed') +
  geom_vline(xintercept = mcap_object_k_Si$ci[2], linetype = 'dashed') +
  geom_vline(xintercept = mcap_object_k_Si$mle, col = 'blue') +
  geom_vline(xintercept = log(mif.estimate[['k_Si']]), col = 'red') +
  geom_hline(yintercept = pf.loglik.of.mif.estimate, col = 'red',linetype = "longdash") +
  geom_point(aes(x = log(mif.estimate[['k_Si']]), 
                 y = pf.loglik.of.mif.estimate), 
             color = "red", size = 2) + 
  labs(x =  TeX("$\\log(\\tau^i_{S})$"), y = "log likelihood") +
  theme(axis.text = element_text(size = 11),
        axis.title = element_text(size = 11)) +
  ylim(-375, -370)+
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
k_Si_p


if(load_option){
  load("./Single-species/Lum/SRJF/profile/ri/ri.RData")
  subset_data_ri <- final_params %>%
    group_by(ri) %>%
    filter(loglik == max(loglik))
}

plot(x = log(subset_data_ri$ri), y = subset_data_ri$loglik)
plot(x = subset_data_ri$ri, y = subset_data_ri$loglik)
subset_data_ri$log_ri <- log(subset_data_ri$ri)

mcap(subset_data_ri$loglik, subset_data_ri$log_ri,  level = 0.95, span = 0.6, Ngrid = 1000) -> mcap_object_ri
mcap_object_ri$mle -> ri_mle
ri_p <- ggplot() +
  geom_point(data = subset_data_ri, aes(x = log_ri, y = loglik)) +
  geom_line(data = mcap_object_ri$fit, aes(x = parameter, y = smoothed), col = 'red') +
  geom_vline(xintercept = mcap_object_ri$ci[1], linetype = 'dashed') +
  # geom_vline(xintercept = mcap_object_ri$ci[2], linetype = 'dashed') +
  geom_vline(xintercept = mcap_object_ri$mle, col = 'blue') +
  geom_vline(xintercept = log(mif.estimate[['ri']]), col = 'red') +
  geom_hline(yintercept = pf.loglik.of.mif.estimate, col = 'red',linetype = "longdash") +
  geom_point(aes(x = log(mif.estimate[['ri']]), 
                 y = pf.loglik.of.mif.estimate), 
             color = "red", size = 2) + 
  labs(x =  TeX("$\\log(r^i)$"), y = "log likelihood") +
  theme(axis.text = element_text(size = 11),
        axis.title = element_text(size = 11)) +
  ylim(-375, -370)+
  theme_bw() +
  theme(axis.title.y = element_blank())
ri_p





if(load_option){
  load("./Single-species/Lum/SRJF/profile/sigF/sigF.RData")
  subset_data_sigF <- final_params %>%
    group_by(sigF) %>%
    filter(loglik == max(loglik))
}

plot(x = log(subset_data_sigF$sigF), y = subset_data_sigF$loglik)
plot(x = subset_data_sigF$sigF, y = subset_data_sigF$loglik)
subset_data_sigF$log_sigF <- log(subset_data_sigF$sigF)
subset_data_sigF = subset_data_sigF[subset_data_sigF$loglik > max(subset_data_sigF$loglik) - 15,]
mcap(subset_data_sigF$loglik, subset_data_sigF$log_sigF,  level = 0.95, span = 0.5, Ngrid = 1000) -> mcap_object_sigF
mcap_object_sigF$mle -> sigF_mle
sigF_p <- ggplot() +
  geom_point(data = subset_data_sigF, aes(x = log_sigF, y = loglik)) +
  geom_line(data = mcap_object_sigF$fit, aes(x = parameter, y = smoothed), col = 'red') +
  # geom_vline(xintercept = mcap_object_sigF$ci[1], linetype = 'dashed') +
  geom_vline(xintercept = mcap_object_sigF$ci[2], linetype = 'dashed') +
  geom_vline(xintercept = mcap_object_sigF$mle, col = 'blue') +
  geom_vline(xintercept = log(mif.estimate[['sigF']]), col = 'red') +
  geom_hline(yintercept = pf.loglik.of.mif.estimate, col = 'red',linetype = "longdash") +
  geom_point(aes(x = log(mif.estimate[['sigF']]), 
                 y = pf.loglik.of.mif.estimate), 
             color = "red", size = 2) + 
  labs(x =  TeX("$\\log(\\sigma_{F})$"), y = "log likelihood") +
  theme(axis.text = element_text(size = 11),
        axis.title = element_text(size = 11)) +
  ylim(-375, -370)+
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
sigF_p







if(load_option){
  load("./Single-species/Lum/SRJF/profile/sigJi/sigJi.RData")
  subset_data_sigJi <- final_params %>%
    group_by(sigJi) %>%
    filter(loglik == max(loglik))
}

plot(x = log(subset_data_sigJi$sigJi), y = subset_data_sigJi$loglik)
plot(x = subset_data_sigJi$sigJi, y = subset_data_sigJi$loglik)
subset_data_sigJi$log_sigJi <- log(subset_data_sigJi$sigJi)
subset_data_sigJi = subset_data_sigJi[subset_data_sigJi$loglik > -400,]
subset_data_sigJi = subset_data_sigJi[subset_data_sigJi$log_sigJi > -2.5,]
mcap(subset_data_sigJi$loglik, subset_data_sigJi$log_sigJi,  level = 0.95, span = 0.7, Ngrid = 1000) -> mcap_object_sigJi
mcap_object_sigJi$mle -> sigJi_mle
sigJi_p <- ggplot() +
  geom_point(data = subset_data_sigJi, aes(x = log_sigJi, y = loglik)) +
  geom_line(data = mcap_object_sigJi$fit, aes(x = parameter, y = smoothed), col = 'red') +
  # geom_vline(xintercept = mcap_object_sigJi$ci[1], linetype = 'dashed') +
  geom_vline(xintercept = mcap_object_sigJi$ci[2], linetype = 'dashed') +
  geom_vline(xintercept = mcap_object_sigJi$mle, col = 'blue') +
  geom_vline(xintercept = log(mif.estimate[['sigJi']]), col = 'red') +
  geom_hline(yintercept = pf.loglik.of.mif.estimate, col = 'red',linetype = "longdash") +
  geom_point(aes(x = log(mif.estimate[['sigJi']]), 
                 y = pf.loglik.of.mif.estimate), 
             color = "red", size = 2) + 
  labs(x =  TeX("$\\log(\\sigma^i_{J})$"), y = "log likelihood") +
  theme(axis.text = element_text(size = 11),
        axis.title = element_text(size = 11)) +
  ylim(-375, -370)+
  # xlim(-1.5,-0)+
  theme_bw() +
  theme(axis.title.y = element_blank()) 

sigJi_p







if(load_option){
  load("./Single-species/Lum/SRJF/profile/theta_Ji/theta_Ji.RData")
  subset_data_theta_Ji <- final_params %>%
    group_by(theta_Ji) %>%
    filter(loglik == max(loglik))
}

plot(x = log(subset_data_theta_Ji$theta_Ji), y = subset_data_theta_Ji$loglik)
plot(x = subset_data_theta_Ji$theta_Ji, y = subset_data_theta_Ji$loglik)
subset_data_theta_Ji$log_theta_Ji <- log(subset_data_theta_Ji$theta_Ji)
# 
mcap(subset_data_theta_Ji$loglik, subset_data_theta_Ji$log_theta_Ji,  level = 0.95, span = 0.95, Ngrid = 1000) -> mcap_object_theta_Ji
mcap_object_theta_Ji$mle -> theta_Ji_mle
theta_Ji_p <- ggplot() +
  geom_point(data = subset_data_theta_Ji, aes(x = log_theta_Ji, y = loglik)) +
  geom_line(data = mcap_object_theta_Ji$fit, aes(x = parameter, y = smoothed), col = 'red') +
  # geom_vline(xintercept = mcap_object_theta_Ji$ci[1], linetype = 'dashed') +
  geom_vline(xintercept = mcap_object_theta_Ji$ci[2], linetype = 'dashed') +
  geom_vline(xintercept = mcap_object_theta_Ji$mle, col = 'blue') +
  geom_vline(xintercept = log(mif.estimate[['theta_Ji']]), col = 'red') +
  geom_hline(yintercept = pf.loglik.of.mif.estimate, col = 'red',linetype = "longdash") +
  geom_point(aes(x = log(mif.estimate[['theta_Ji']]), 
                 y = pf.loglik.of.mif.estimate), 
             color = "red", size = 2) + 
  labs(x =  TeX("$\\log(\\theta^i_{J})$"), y = "log likelihood") +
  theme(axis.text = element_text(size = 11),
        axis.title = element_text(size = 11)) +
  ylim(-375, -370)+
  # xlim(-2.2,0)+
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) 
theta_Ji_p






if(load_option){
  load("./Single-species/Lum/SRJF/profile/theta_Si/theta_Si.RData")
  subset_data_theta_Si <- final_params %>%
    group_by(theta_Si) %>%
    filter(loglik == max(loglik))
}

plot(x = log(subset_data_theta_Si$theta_Si), y = subset_data_theta_Si$loglik)
plot(x = subset_data_theta_Si$theta_Si, y = subset_data_theta_Si$loglik)
subset_data_theta_Si$log_theta_Si <- log(subset_data_theta_Si$theta_Si)

mcap(subset_data_theta_Si$loglik, subset_data_theta_Si$log_theta_Si,  level = 0.95, span = 0.9, Ngrid = 1000) -> mcap_object_theta_Si
mcap_object_theta_Si$mle -> theta_Si_mle
theta_Si_p <- ggplot() +
  geom_point(data = subset_data_theta_Si, aes(x = log_theta_Si, y = loglik)) +
  geom_line(data = mcap_object_theta_Si$fit, aes(x = parameter, y = smoothed), col = 'red') +
  geom_vline(xintercept = mcap_object_theta_Si$ci[1], linetype = 'dashed') +
  geom_vline(xintercept = mcap_object_theta_Si$ci[2], linetype = 'dashed') +
  geom_vline(xintercept = mcap_object_theta_Si$mle, col = 'blue') +
  geom_vline(xintercept = log(mif.estimate[['theta_Si']]), col = 'red') +
  geom_hline(yintercept = pf.loglik.of.mif.estimate, col = 'red',linetype = "longdash") +
  geom_point(aes(x = log(mif.estimate[['theta_Si']]), 
                 y = pf.loglik.of.mif.estimate), 
             color = "red", size = 2) + 
  labs(x =  TeX("$\\log(\\theta^i_{S})$"), y = "log likelihood") +
  theme(axis.text = element_text(size = 11),
        axis.title = element_text(size = 11)) +
  ylim(-375, -370)+
  # xlim(-2.2,0)+
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) 
theta_Si_p



combine_subset = rbind(subset_data_ri[,1:9],
                       subset_data_f_Si[,1:9],
                       subset_data_theta_Si[,1:9],
                       subset_data_theta_Ji[,1:9],
                       subset_data_sigJi[,1:9],
                       subset_data_sigF[,1:9],
                       subset_data_k_Si[,1:9])

combine_subset$ri_f_Si = combine_subset$ri * combine_subset$f_Si
combine_subset$log_ri_f_Si <- log(combine_subset$ri_f_Si)
plot(x = combine_subset$log_ri_f_Si, y = combine_subset$loglik)
combine_subset_ri_f_Si <- combine_subset %>%
  ungroup() %>%
  mutate(
    bin = cut(
      log_ri_f_Si,
      breaks = seq(
        min(log_ri_f_Si, na.rm = TRUE),
        max(log_ri_f_Si, na.rm = TRUE),
        length.out = 101
      ),
      include.lowest = TRUE,
      right = FALSE
    )
  ) %>%
  group_by(bin) %>%
  slice_max(loglik, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(-bin)

mcap(combine_subset_ri_f_Si$loglik, combine_subset_ri_f_Si$log_ri_f_Si,  level = 0.95, span = 0.95, Ngrid = 1000) -> mcap_object_ri_f_Si
mcap_object_ri_f_Si$mle -> ri_f_Si_mle
ri_f_Si_p <- ggplot() +
  geom_point(data = combine_subset_ri_f_Si, aes(x = log_ri_f_Si, y = loglik)) +
  geom_line(data = mcap_object_ri_f_Si$fit, aes(x = parameter, y = smoothed), col = 'red') +
  geom_vline(xintercept = mcap_object_ri_f_Si$ci[1], linetype = 'dashed') +
  geom_vline(xintercept = mcap_object_ri_f_Si$ci[2], linetype = 'dashed') +
  geom_vline(xintercept = mcap_object_ri_f_Si$mle, col = 'blue') +
  geom_vline(xintercept = log(mif.estimate[['ri']] * mif.estimate[['f_Si']]), col = 'red') +
  geom_hline(yintercept = pf.loglik.of.mif.estimate, col = 'red',linetype = "longdash") +
  geom_point(aes(x = log(mif.estimate[['ri']] * mif.estimate[['f_Si']]), 
                 y = pf.loglik.of.mif.estimate), 
             color = "red", size = 2) + 
  labs(x =  TeX("$\\log(\\r^i \\cdot f^i_{S})$"), y = "log likelihood") +
  theme(axis.text = element_text(size = 11),
        axis.title = element_text(size = 11)) +
  ylim(-375, -370)+
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

ri_f_Si_p



grid.arrange( ri_p,  f_Si_p,
              theta_Si_p, theta_Ji_p,
             sigJi_p,sigF_p,k_Si_p,ri_f_Si_p,
              nrow = 2, ncol = 4)

g <- arrangeGrob( ri_p,  f_Si_p,
              theta_Si_p, theta_Ji_p,
              sigJi_p,sigF_p,k_Si_p,ri_f_Si_p,
              nrow = 2, ncol = 4)

ggsave(
  filename = "./daphnia-article/si/profile/Simple_dynamics/Lum/no_para/Profile_plot.png",
  plot     = g,         
  width    = 16,        
  height   = 8,         
  dpi      = 300,       
  units    = "in"
)

save(subset_data_ri,
     subset_data_f_Si,
     subset_data_theta_Si,
     subset_data_theta_Ji,
     subset_data_sigJi,
     subset_data_sigF,
     subset_data_k_Si,
     combine_subset_ri_f_Si,file = "data/Simple_dynamics/Lum/no_para/profile_graph_data.rda")


save(mcap_object_ri,
     mcap_object_f_Si,
     mcap_object_theta_Si,
     mcap_object_theta_Ji,
     mcap_object_sigJi,
     mcap_object_sigF,
     mcap_object_k_Si,
     mcap_object_ri_f_Si,file = "data/Simple_dynamics/Lum/no_para/profile_mcap_object.rda")


ci_table = cbind(mcap_object_ri$ci,
                 mcap_object_f_Si$ci,
                 mcap_object_theta_Si$ci,
                 mcap_object_theta_Ji$ci,
                 mcap_object_sigJi$ci,
                 mcap_object_sigF$ci,
                 mcap_object_k_Si$ci,
                 mcap_object_ri_f_Si$ci)

ci_table = rbind(ci_table,c(mcap_object_ri$mle,
                            mcap_object_f_Si$mle,
                            mcap_object_theta_Si$mle,
                            mcap_object_theta_Ji$mle,
                            mcap_object_sigJi$mle,
                            mcap_object_sigF$mle,
                            mcap_object_k_Si$mle,
                            mcap_object_ri_f_Si$mle))


rownames(ci_table) = c("2.5%","97.5%","MLE")
colnames(ci_table) = c('ri','f_Si','theta_Si',
                       'theta_Ji','sigJi','sigF','k_Si','ri_f_Si')
ci_table = as.data.frame(ci_table)
save(ci_table, file = 'data/Simple_dynamics/Lum/no_para/profile_ci_table.rda')
















































































