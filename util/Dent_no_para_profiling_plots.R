library(dplyr)
library(latex2exp)
library(gridExtra)
library(dplyr)
library(pomp)
library(ggplot2)

#Please set the working to be the 'Daphnia-ms' path

load("data/Simple_dynamics/Dent/no_para/profile_graph_data.rda")
load("Single-species/Dent/SRJF/model/best_result.rda")

load_option = FALSE



if(load_option){
  load("./Single-species/Dent/SRJF/profile/f_Sn/f_Sn.RData")
  subset_data_f_Sn <- final_params %>%
    group_by(f_Sn) %>%
    filter(loglik == max(loglik))
}

plot(x = log(subset_data_f_Sn$f_Sn), y = subset_data_f_Sn$loglik)
plot(x = subset_data_f_Sn$f_Sn, y = subset_data_f_Sn$loglik)
subset_data_f_Sn$log_f_Sn <- log(subset_data_f_Sn$f_Sn)

mcap(subset_data_f_Sn$loglik, subset_data_f_Sn$log_f_Sn,  level = 0.95, span = 0.95, Ngrid = 1000) -> mcap_object_f_Sn
mcap_object_f_Sn$mle -> f_Sn_mle
f_Sn_p <- ggplot() +
  geom_point(data = subset_data_f_Sn, aes(x = log_f_Sn, y = loglik)) +
  geom_line(data = mcap_object_f_Sn$fit, aes(x = parameter, y = smoothed), col = 'red') +
  geom_vline(xintercept = mcap_object_f_Sn$ci[1], linetype = 'dashed') +
  geom_vline(xintercept = mcap_object_f_Sn$ci[2], linetype = 'dashed') +
  geom_vline(xintercept = mcap_object_f_Sn$mle, col = 'blue') +
  geom_vline(xintercept = log(mif.estimate[['f_Sn']]), col = 'red') +
  geom_hline(yintercept = pf.loglik.of.mif.estimate, col = 'red',linetype = "longdash") +
  geom_point(aes(x = log(mif.estimate[['f_Sn']]), 
                 y = pf.loglik.of.mif.estimate), 
             color = "red", size = 2) + 
  labs(x =  TeX("$\\log(f^n_{S})$"), y = "log likelihood") +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 10)) +
  ylim(-505, -495)+
  # xlim(-9.5,-5.5)+
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())  + 
  annotate("text", x = mcap_object_f_Sn$mle, y = -900, label = sprintf("f_Sn_mle: %s", formatC(mcap_object_f_Sn$mle, format = 'e', digits = 3)), hjust = 1.05, vjust = -0.5, size = 3)

f_Sn_p




if(load_option){
  load("./Single-species/Dent/SRJF/profile/k_Sn/k_Sn.RData")
  subset_data_k_Sn <- final_params %>%
    group_by(k_Sn) %>%
    filter(loglik == max(loglik))
}

plot(x = log(subset_data_k_Sn$k_Sn), y = subset_data_k_Sn$loglik)
plot(x = subset_data_k_Sn$k_Sn, y = subset_data_k_Sn$loglik)
subset_data_k_Sn$log_k_Sn <- log(subset_data_k_Sn$k_Sn)

mcap(subset_data_k_Sn$loglik, subset_data_k_Sn$log_k_Sn,  level = 0.95, span = 0.95, Ngrid = 1000) -> mcap_object_k_Sn
mcap_object_k_Sn$mle -> k_Sn_mle
k_Sn_p <- ggplot() +
  geom_point(data = subset_data_k_Sn, aes(x = log_k_Sn, y = loglik)) +
  geom_line(data = mcap_object_k_Sn$fit, aes(x = parameter, y = smoothed), col = 'red') +
  geom_vline(xintercept = mcap_object_k_Sn$ci[1], linetype = 'dashed') +
  geom_vline(xintercept = mcap_object_k_Sn$ci[2], linetype = 'dashed') +
  geom_vline(xintercept = mcap_object_k_Sn$mle, col = 'blue') +
  geom_vline(xintercept = log(mif.estimate[['k_Sn']]), col = 'red') +
  geom_hline(yintercept = pf.loglik.of.mif.estimate, col = 'red',linetype = "longdash") +
  geom_point(aes(x = log(mif.estimate[['k_Sn']]), 
                 y = pf.loglik.of.mif.estimate), 
             color = "red", size = 2) + 
  labs(x =  TeX("$\\log(\\tau^n_{S})$"), y = "log likelihood") +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 10)) +
  ylim(-505, -495)+
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())  + 
  annotate("text", x = mcap_object_k_Sn$mle, y = -900, label = sprintf("k_Sn_mle: %s", formatC(mcap_object_k_Sn$mle, format = 'e', digits = 3)), hjust = 1.05, vjust = -0.5, size = 3)

k_Sn_p






if(load_option){
  load("./Single-species/Dent/SRJF/profile/rn/rn.RData")
  subset_data_rn <- final_params %>%
    group_by(rn) %>%
    filter(loglik == max(loglik))
}

plot(x = log(subset_data_rn$rn), y = subset_data_rn$loglik)
plot(x = subset_data_rn$rn, y = subset_data_rn$loglik)
subset_data_rn$log_rn <- log(subset_data_rn$rn)

mcap(subset_data_rn$loglik, subset_data_rn$log_rn,  level = 0.95, span = 0.95, Ngrid = 1000) -> mcap_object_rn
mcap_object_rn$mle -> rn_mle
rn_p <- ggplot() +
  geom_point(data = subset_data_rn, aes(x = log_rn, y = loglik)) +
  geom_line(data = mcap_object_rn$fit, aes(x = parameter, y = smoothed), col = 'red') +
  geom_vline(xintercept = mcap_object_rn$ci[1], linetype = 'dashed') +
  geom_vline(xintercept = mcap_object_rn$ci[2], linetype = 'dashed') +
  geom_vline(xintercept = mcap_object_rn$mle, col = 'blue') +
  geom_vline(xintercept = log(mif.estimate[['rn']]), col = 'red') +
  geom_hline(yintercept = pf.loglik.of.mif.estimate, col = 'red',linetype = "longdash") +
  geom_point(aes(x = log(mif.estimate[['rn']]), 
                 y = pf.loglik.of.mif.estimate), 
             color = "red", size = 2) + 
  labs(x =  TeX("$\\log(r^n)$"), y = "log likelihood") +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 10)) +
  ylim(-505, -495)+
  theme_bw() +
  theme(axis.title.y = element_blank())  + 
  annotate("text", x = mcap_object_rn$mle, y = -900, label = sprintf("rn_mle: %s", formatC(mcap_object_rn$mle, format = 'e', digits = 3)), hjust = 1.05, vjust = -0.5, size = 3)

rn_p




if(load_option){
  load("./Single-species/Dent/SRJF/profile/sigF/sigF.RData")
  subset_data_sigF <- final_params %>%
    group_by(sigF) %>%
    filter(loglik == max(loglik))
}

plot(x = log(subset_data_sigF$sigF), y = subset_data_sigF$loglik)
plot(x = subset_data_sigF$sigF, y = subset_data_sigF$loglik)
subset_data_sigF$log_sigF <- log(subset_data_sigF$sigF)

mcap(subset_data_sigF$loglik, subset_data_sigF$log_sigF,  level = 0.95, span = 0.95, Ngrid = 1000) -> mcap_object_sigF
mcap_object_sigF$mle -> sigF_mle
sigF_p <- ggplot() +
  geom_point(data = subset_data_sigF, aes(x = log_sigF, y = loglik)) +
  geom_line(data = mcap_object_sigF$fit, aes(x = parameter, y = smoothed), col = 'red') +
  geom_vline(xintercept = mcap_object_sigF$ci[1], linetype = 'dashed') +
  geom_vline(xintercept = mcap_object_sigF$ci[2], linetype = 'dashed') +
  geom_vline(xintercept = mcap_object_sigF$mle, col = 'blue') +
  geom_vline(xintercept = log(mif.estimate[['sigF']]), col = 'red') +
  geom_hline(yintercept = pf.loglik.of.mif.estimate, col = 'red',linetype = "longdash") +
  geom_point(aes(x = log(mif.estimate[['sigF']]), 
                 y = pf.loglik.of.mif.estimate), 
             color = "red", size = 2) + 
  labs(x =  TeX("$\\log(\\sigma_{F})$"), y = "log likelihood") +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 10)) +
  ylim(-505, -495)+
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())  + 
  annotate("text", x = mcap_object_sigF$mle, y = -900, label = sprintf("sigF_mle: %s", formatC(mcap_object_sigF$mle, format = 'e', digits = 3)), hjust = 1.05, vjust = -0.5, size = 3)

sigF_p







if(load_option){
  load("./Single-species/Dent/SRJF/profile/sigJn/sigJn.RData")
  subset_data_sigJn <- final_params %>%
    group_by(sigJn) %>%
    filter(loglik == max(loglik))
}

plot(x = log(subset_data_sigJn$sigJn), y = subset_data_sigJn$loglik)
plot(x = subset_data_sigJn$sigJn, y = subset_data_sigJn$loglik)
subset_data_sigJn$log_sigJn <- log(subset_data_sigJn$sigJn)

mcap(subset_data_sigJn$loglik, subset_data_sigJn$log_sigJn,  level = 0.95, span = 0.95, Ngrid = 1000) -> mcap_object_sigJn
mcap_object_sigJn$mle -> sigJn_mle
sigJn_p <- ggplot() +
  geom_point(data = subset_data_sigJn, aes(x = log_sigJn, y = loglik)) +
  geom_line(data = mcap_object_sigJn$fit, aes(x = parameter, y = smoothed), col = 'red') +
  geom_vline(xintercept = mcap_object_sigJn$ci[1], linetype = 'dashed') +
  geom_vline(xintercept = mcap_object_sigJn$ci[2], linetype = 'dashed') +
  geom_vline(xintercept = mcap_object_sigJn$mle, col = 'blue') +
  geom_vline(xintercept = log(mif.estimate[['sigJn']]), col = 'red') +
  geom_hline(yintercept = pf.loglik.of.mif.estimate, col = 'red',linetype = "longdash") +
  geom_point(aes(x = log(mif.estimate[['sigJn']]), 
                 y = pf.loglik.of.mif.estimate), 
             color = "red", size = 2) + 
  labs(x =  TeX("$\\log(\\sigma^n_{J})$"), y = "log likelihood") +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 10)) +
  ylim(-505, -495)+
  # xlim(-3,0)+
  theme_bw() +
  theme(axis.title.y = element_blank())  + 
  annotate("text", x = mcap_object_sigJn$mle, y = -900, label = sprintf("sigJn_mle: %s", formatC(mcap_object_sigJn$mle, format = 'e', digits = 3)), hjust = 1.05, vjust = -0.5, size = 3)

sigJn_p










if(load_option){
  load("./Single-species/Dent/SRJF/profile/theta_Jn/theta_Jn.RData")
  subset_data_theta_Jn <- final_params %>%
    group_by(theta_Jn) %>%
    filter(loglik == max(loglik))
}

plot(x = log(subset_data_theta_Jn$theta_Jn), y = subset_data_theta_Jn$loglik)
plot(x = subset_data_theta_Jn$theta_Jn, y = subset_data_theta_Jn$loglik)
subset_data_theta_Jn$log_theta_Jn <- log(subset_data_theta_Jn$theta_Jn)

mcap(subset_data_theta_Jn$loglik, subset_data_theta_Jn$log_theta_Jn,  level = 0.95, span = 0.95, Ngrid = 1000) -> mcap_object_theta_Jn
mcap_object_theta_Jn$mle -> theta_Jn_mle
theta_Jn_p <- ggplot() +
  geom_point(data = subset_data_theta_Jn, aes(x = log_theta_Jn, y = loglik)) +
  geom_line(data = mcap_object_theta_Jn$fit, aes(x = parameter, y = smoothed), col = 'red') +
  # geom_vline(xintercept = mcap_object_theta_Jn$ci[1], linetype = 'dashed') +
  geom_vline(xintercept = mcap_object_theta_Jn$ci[2], linetype = 'dashed') +
  geom_vline(xintercept = mcap_object_theta_Jn$mle, col = 'blue') +
  geom_vline(xintercept = log(mif.estimate[['theta_Jn']]), col = 'red') +
  geom_hline(yintercept = pf.loglik.of.mif.estimate, col = 'red',linetype = "longdash") +
  geom_point(aes(x = log(mif.estimate[['theta_Jn']]), 
                 y = pf.loglik.of.mif.estimate), 
             color = "red", size = 2) + 
  labs(x =  TeX("$\\log(\\theta^n_{J})$"), y = "log likelihood") +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 10)) +
  ylim(-505, -495)+
  # xlim(-2.2,0)+
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())  
theta_Jn_p






if(load_option){
  load("./Single-species/Dent/SRJF/profile/theta_Sn/theta_Sn.RData")
  subset_data_theta_Sn <- final_params %>%
    group_by(theta_Sn) %>%
    filter(loglik == max(loglik))
}

plot(x = log(subset_data_theta_Sn$theta_Sn), y = subset_data_theta_Sn$loglik)
plot(x = subset_data_theta_Sn$theta_Sn, y = subset_data_theta_Sn$loglik)
subset_data_theta_Sn$log_theta_Sn <- log(subset_data_theta_Sn$theta_Sn)

mcap(subset_data_theta_Sn$loglik, subset_data_theta_Sn$log_theta_Sn,  level = 0.95, span = 0.95, Ngrid = 1000) -> mcap_object_theta_Sn
mcap_object_theta_Sn$mle -> theta_Sn_mle
theta_Sn_p <- ggplot() +
  geom_point(data = subset_data_theta_Sn, aes(x = log_theta_Sn, y = loglik)) +
  geom_line(data = mcap_object_theta_Sn$fit, aes(x = parameter, y = smoothed), col = 'red') +
  geom_vline(xintercept = mcap_object_theta_Sn$ci[1], linetype = 'dashed') +
  geom_vline(xintercept = mcap_object_theta_Sn$ci[2], linetype = 'dashed') +
  geom_vline(xintercept = mcap_object_theta_Sn$mle, col = 'blue') +
  geom_vline(xintercept = log(mif.estimate[['theta_Sn']]), col = 'red') +
  geom_hline(yintercept = pf.loglik.of.mif.estimate, col = 'red',linetype = "longdash") +
  geom_point(aes(x = log(mif.estimate[['theta_Sn']]), 
                 y = pf.loglik.of.mif.estimate), 
             color = "red", size = 2) + 
  labs(x =  TeX("$\\log(\\theta^n_{S})$"), y = "log likelihood") +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 10)) +
  ylim(-505, -495)+
  # xlim(-2.2,0)+
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())  + 
  annotate("text", x = mcap_object_theta_Sn$mle, y = -900, label = sprintf("theta_Sn_mle: %s", formatC(mcap_object_theta_Sn$mle, format = 'e', digits = 3)), hjust = 1.05, vjust = -0.5, size = 3)

theta_Sn_p



grid.arrange( rn_p,f_Sn_p,
              theta_Sn_p,theta_Jn_p,sigJn_p,sigF_p,k_Sn_p,
              nrow = 2, ncol = 4)



save(subset_data_rn,
     subset_data_f_Sn,
     subset_data_theta_Sn,
     subset_data_theta_Jn,
     subset_data_sigJn,
     subset_data_sigF,
     subset_data_k_Sn, file = "data/Simple_dynamics/Dent/no_para/profile_graph_data.rda")


save(mcap_object_rn,
     mcap_object_f_Sn,
     mcap_object_theta_Sn,
     mcap_object_theta_Jn, 
     mcap_object_sigJn,
     mcap_object_sigF,
     mcap_object_k_Sn, file = "data/Simple_dynamics/Dent/no_para/profile_mcap_object.rda")


ci_table = cbind(mcap_object_rn$ci,
                 mcap_object_f_Sn$ci,
                 mcap_object_theta_Sn$ci,
                 mcap_object_theta_Jn$ci,
                 mcap_object_sigJn$ci,
                 mcap_object_sigF$ci,
                 mcap_object_k_Sn$ci)

ci_table = rbind(ci_table,c(mcap_object_rn$mle,
                            mcap_object_f_Sn$mle,
                            mcap_object_theta_Sn$mle,
                            mcap_object_theta_Jn$mle,
                            mcap_object_sigJn$mle,
                            mcap_object_sigF$mle,
                            mcap_object_k_Sn$mle))


rownames(ci_table) = c("2.5%","97.5%","MLE")
colnames(ci_table) = c('rn','f_Sn','theta_Sn','theta_Jn','sigJn','sigF','k_Sn')
ci_table = as.data.frame(ci_table)
save(ci_table, file = 'data/Simple_dynamics/Dent/no_para/profile_ci_table.rda')
















































































