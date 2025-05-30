library(dplyr)
library(latex2exp)
library(gridExtra)
library(dplyr)
library(pomp)
library(ggplot2)


#Please set the working to be the 'Daphnia-ms' path


load("data/Target_dynamics/no_para/profile_graph_data.rda")
load("./Mixed-species/SRJF2/best_result.rda")


load_option = FALSE

if(load_option){
  load("./Mixed-species/SRJF2/profile/f_Si/f_Si.RData")
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
  theme(axis.text = element_text(size =11),
        axis.title = element_text(size = 11)) +
  ylim(-710, -699)+
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
f_Si_p









if(load_option){
  load("./Mixed-species/SRJF2/profile/f_Sn/f_Sn.RData")
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
  theme(axis.text = element_text(size =11),
        axis.title = element_text(size = 11)) +
  ylim(-710, -699)+
  # xlim(-9.5,-5.5)+
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) 
f_Sn_p












if(load_option){
  load("./Mixed-species/SRJF2/profile/k_Si/k_Si.RData")
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
  theme(axis.text = element_text(size =11),
        axis.title = element_text(size = 11)) +
  ylim(-710, -699)+
  xlim(-0.5, 2.5)+
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())  + 
  annotate("text", x = mcap_object_k_Si$mle, y = -900, label = sprintf("k_Si_mle: %s", formatC(mcap_object_k_Si$mle, format = 'e', digits = 3)), hjust = 1.05, vjust = -0.5, size = 3)

k_Si_p





if(load_option){
  load("./Mixed-species/SRJF2/profile/k_Sn/k_Sn.RData")
  subset_data_k_Sn <- final_params %>%
    group_by(k_Sn) %>%
    filter(loglik == max(loglik))
}

plot(x = log(subset_data_k_Sn$k_Sn), y = subset_data_k_Sn$loglik)
plot(x = subset_data_k_Sn$k_Sn, y = subset_data_k_Sn$loglik)
subset_data_k_Sn$log_k_Sn <- log(subset_data_k_Sn$k_Sn)
subset_data_k_Sn = subset_data_k_Sn[subset_data_k_Sn$log_k_Sn <=5,]
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
  theme(axis.text = element_text(size =11),
        axis.title = element_text(size = 11)) +
  ylim(-710, -699)+
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())  + 
  annotate("text", x = mcap_object_k_Sn$mle, y = -900, label = sprintf("k_Sn_mle: %s", formatC(mcap_object_k_Sn$mle, format = 'e', digits = 3)), hjust = 1.05, vjust = -0.5, size = 3)

k_Sn_p







if(load_option){
  load("./Mixed-species/SRJF2/profile/ri/ri.RData")
  subset_data_ri <- final_params %>%
    group_by(ri) %>%
    filter(loglik == max(loglik))
}
plot(x = log(subset_data_ri$ri), y = subset_data_ri$loglik)
subset_data_ri$log_ri <- log(subset_data_ri$ri)

mcap(subset_data_ri$loglik, subset_data_ri$log_ri,  level = 0.95, span = 0.95, Ngrid = 1000) -> mcap_object_ri
mcap_object_ri$mle -> ri_mle
ri_p <- ggplot() +
  geom_point(data = subset_data_ri, aes(x = log_ri, y = loglik)) +
  geom_line(data = mcap_object_ri$fit, aes(x = parameter, y = smoothed), col = 'red') +
  geom_vline(xintercept = mcap_object_ri$ci[1], linetype = 'dashed') +
  geom_vline(xintercept = mcap_object_ri$ci[2], linetype = 'dashed') +
  geom_vline(xintercept = mcap_object_ri$mle, col = 'blue') +
  geom_vline(xintercept = log(mif.estimate[['ri']]), col = 'red') +
  geom_hline(yintercept = pf.loglik.of.mif.estimate, col = 'red',linetype = "longdash") +
  geom_point(aes(x = log(mif.estimate[['ri']]), 
                 y = pf.loglik.of.mif.estimate), 
             color = "red", size = 2) + 
  labs(x =  TeX("$r^i$"), y = "log likelihood") +
  theme(axis.text = element_text(size =11),
        axis.title = element_text(size = 11)) +
  ylim(-710, -699)+
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) 
ri_p








if(load_option){
  load("./Mixed-species/SRJF2/profile/rn/rn.RData")
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
  labs(x =  TeX("$r^n$"), y = "log likelihood") +
  theme(axis.text = element_text(size =11),
        axis.title = element_text(size = 11)) +
  ylim(-710, -699)+
  theme_bw() +
  theme(axis.title.y = element_blank())  + 
  annotate("text", x = mcap_object_rn$mle, y = -900, label = sprintf("rn_mle: %s", formatC(mcap_object_rn$mle, format = 'e', digits = 3)), hjust = 1.05, vjust = -0.5, size = 3)

rn_p






if(load_option){
  load("./Mixed-species/SRJF2/profile/sigF/sigF.RData")
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
  theme(axis.text = element_text(size =11),
        axis.title = element_text(size = 11)) +
  ylim(-710, -699)+
  theme_bw() +
  theme(axis.title.y = element_blank())  + 
  annotate("text", x = mcap_object_sigF$mle, y = -900, label = sprintf("sigF_mle: %s", formatC(mcap_object_sigF$mle, format = 'e', digits = 3)), hjust = 1.05, vjust = -0.5, size = 3)

sigF_p










if(load_option){
  load("./Mixed-species/SRJF2/profile/sigJi/sigJi.RData")
  subset_data_sigJi <- final_params %>%
    group_by(sigJi) %>%
    filter(loglik == max(loglik))
}
plot(x = log(subset_data_sigJi$sigJi), y = subset_data_sigJi$loglik)
plot(x = subset_data_sigJi$sigJi, y = subset_data_sigJi$loglik)
subset_data_sigJi$log_sigJi <- log(subset_data_sigJi$sigJi)
subset_data_sigJi = subset_data_sigJi[subset_data_sigJi$loglik > max(subset_data_sigJi$loglik) - 15,]
mcap(subset_data_sigJi$loglik, subset_data_sigJi$log_sigJi,  level = 0.95, span = 0.6, Ngrid = 1000) -> mcap_object_sigJi
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
  theme(axis.text = element_text(size =11),
        axis.title = element_text(size = 11)) +
  ylim(-710, -699)+
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
sigJi_p









if(load_option){
  load("./Mixed-species/SRJF2/profile/sigJn/sigJn.RData")
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
  theme(axis.text = element_text(size =11),
        axis.title = element_text(size = 11)) +
  ylim(-710, -699)+
  # xlim(-3,0)+
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())  + 
  annotate("text", x = mcap_object_sigJn$mle, y = -900, label = sprintf("sigJn_mle: %s", formatC(mcap_object_sigJn$mle, format = 'e', digits = 3)), hjust = 1.05, vjust = -0.5, size = 3)

sigJn_p







if(load_option){
  load("./Mixed-species/SRJF2/profile/theta_Ji/theta_Ji.RData")
  subset_data_theta_Ji <- final_params %>%
    group_by(theta_Ji) %>%
    filter(loglik == max(loglik))
}
plot(x = log(subset_data_theta_Ji$theta_Ji), y = subset_data_theta_Ji$loglik)
plot(x = subset_data_theta_Ji$theta_Ji, y = subset_data_theta_Ji$loglik)
subset_data_theta_Ji$log_theta_Ji <- log(subset_data_theta_Ji$theta_Ji)

mcap(subset_data_theta_Ji$loglik, subset_data_theta_Ji$log_theta_Ji,  level = 0.95, span = 0.95, Ngrid = 1000) -> mcap_object_theta_Ji
mcap_object_theta_Ji$mle -> theta_Ji_mle
theta_Ji_p <- ggplot() +
  geom_point(data = subset_data_theta_Ji, aes(x = log_theta_Ji, y = loglik)) +
  geom_line(data = mcap_object_theta_Ji$fit, aes(x = parameter, y = smoothed), col = 'red') +
  geom_vline(xintercept = mcap_object_theta_Ji$ci[1], linetype = 'dashed') +
  geom_vline(xintercept = mcap_object_theta_Ji$ci[2], linetype = 'dashed') +
  geom_vline(xintercept = mcap_object_theta_Ji$mle, col = 'blue') +
  geom_vline(xintercept = log(mif.estimate[['theta_Ji']]), col = 'red') +
  geom_hline(yintercept = pf.loglik.of.mif.estimate, col = 'red',linetype = "longdash") +
  geom_point(aes(x = log(mif.estimate[['theta_Ji']]), 
                 y = pf.loglik.of.mif.estimate), 
             color = "red", size = 2) + 
  labs(x =  TeX("$\\log(\\theta^i_{J})$"), y = "log likelihood") +
  theme(axis.text = element_text(size =11),
        axis.title = element_text(size = 11)) +
  ylim(-710, -699)+
  # xlim(-2.2,0)+
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) 
theta_Ji_p








if(load_option){
  load("./Mixed-species/SRJF2/profile/theta_Jn/theta_Jn.RData")
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
  geom_vline(xintercept = mcap_object_theta_Jn$ci[1], linetype = 'dashed') +
  geom_vline(xintercept = mcap_object_theta_Jn$ci[2], linetype = 'dashed') +
  geom_vline(xintercept = mcap_object_theta_Jn$mle, col = 'blue') +
  geom_vline(xintercept = log(mif.estimate[['theta_Jn']]), col = 'red') +
  geom_hline(yintercept = pf.loglik.of.mif.estimate, col = 'red',linetype = "longdash") +
  geom_point(aes(x = log(mif.estimate[['theta_Jn']]), 
                 y = pf.loglik.of.mif.estimate), 
             color = "red", size = 2) + 
  labs(x =  TeX("$\\log(\\theta^n_{J})$"), y = "log likelihood") +
  theme(axis.text = element_text(size =11),
        axis.title = element_text(size = 11)) +
  ylim(-710, -699)+
  # xlim(-2.2,0)+
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())  + 
  annotate("text", x = mcap_object_theta_Jn$mle, y = -900, label = sprintf("theta_Jn_mle: %s", formatC(mcap_object_theta_Jn$mle, format = 'e', digits = 3)), hjust = 1.05, vjust = -0.5, size = 3)

theta_Jn_p







if(load_option){
  load("./Mixed-species/SRJF2/profile/theta_Si/theta_Si.RData")
  subset_data_theta_Si <- final_params %>%
    group_by(theta_Si) %>%
    filter(loglik == max(loglik))
}

plot(x = log(subset_data_theta_Si$theta_Si), y = subset_data_theta_Si$loglik)
plot(x = subset_data_theta_Si$theta_Si, y = subset_data_theta_Si$loglik)
subset_data_theta_Si$log_theta_Si <- log(subset_data_theta_Si$theta_Si)

mcap(subset_data_theta_Si$loglik, subset_data_theta_Si$log_theta_Si,  level = 0.95, span = 0.95, Ngrid = 1000) -> mcap_object_theta_Si
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
  theme(axis.text = element_text(size =11),
        axis.title = element_text(size = 11)) +
  ylim(-710, -699)+
  # xlim(-2.2,0)+
  theme_bw() +
  theme(axis.title.y = element_blank())  + 
  annotate("text", x = mcap_object_theta_Si$mle, y = -900, label = sprintf("theta_Si_mle: %s", formatC(mcap_object_theta_Si$mle, format = 'e', digits = 3)), hjust = 1.05, vjust = -0.5, size = 3)

theta_Si_p






if(load_option){
  load("./Mixed-species/SRJF2/profile/theta_Sn/theta_Sn.RData")
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
  theme(axis.text = element_text(size =11),
        axis.title = element_text(size = 11)) +
  ylim(-710, -699)+
  # xlim(-2.2,0)+
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())  + 
  annotate("text", x = mcap_object_theta_Sn$mle, y = -900, label = sprintf("theta_Sn_mle: %s", formatC(mcap_object_theta_Sn$mle, format = 'e', digits = 3)), hjust = 1.05, vjust = -0.5, size = 3)

theta_Sn_p


combine_subset = rbind(subset_data_ri[,1:16],
                       subset_data_rn[,1:16],
                       subset_data_f_Si[,1:16],
                       subset_data_f_Sn[,1:16],
                       subset_data_theta_Sn[,1:16],
                       subset_data_theta_Si[,1:16],
                       subset_data_theta_Ji[,1:16],
                       subset_data_theta_Jn[,1:16],
                       subset_data_sigJi[,1:16],
                       subset_data_sigJn[,1:16],
                       subset_data_sigF[,1:16],
                       subset_data_k_Si[,1:16],
                       subset_data_k_Sn[,1:16])


combine_subset$rn_f_Sn = combine_subset$rn * combine_subset$f_Sn
combine_subset$log_rn_f_Sn <- log(combine_subset$rn_f_Sn)
plot(x = combine_subset$log_rn_f_Sn, y = combine_subset$loglik)

combine_subset_rn_f_Sn <- combine_subset %>%
  ungroup() %>%
  mutate(
    bin = cut(
      log_rn_f_Sn,
      breaks = seq(
        min(log_rn_f_Sn, na.rm = TRUE),
        max(log_rn_f_Sn, na.rm = TRUE),
        length.out = 51
      ),
      include.lowest = TRUE,
      right = FALSE
    )
  ) %>%
  group_by(bin) %>%
  slice_max(loglik, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(-bin)

mcap(combine_subset_rn_f_Sn$loglik, combine_subset_rn_f_Sn$log_rn_f_Sn,  level = 0.95, span = 0.95, Ngrid = 1000) -> mcap_object_rn_f_Sn
mcap_object_rn_f_Sn$mle -> rn_f_Sn_mle
rn_f_Sn_p <- ggplot() +
  geom_point(data = combine_subset_rn_f_Sn, aes(x = log_rn_f_Sn, y = loglik)) +
  geom_line(data = mcap_object_rn_f_Sn$fit, aes(x = parameter, y = smoothed), col = 'red') +
  geom_vline(xintercept = mcap_object_rn_f_Sn$ci[1], linetype = 'dashed') +
  geom_vline(xintercept = mcap_object_rn_f_Sn$ci[2], linetype = 'dashed') +
  geom_vline(xintercept = mcap_object_rn_f_Sn$mle, col = 'blue') +
  geom_vline(xintercept = log(mif.estimate[['rn']] * mif.estimate[['f_Sn']]), col = 'red') +
  geom_hline(yintercept = pf.loglik.of.mif.estimate, col = 'red',linetype = "longdash") +
  geom_point(aes(x = log(mif.estimate[['rn']] * mif.estimate[['f_Sn']]), 
                 y = pf.loglik.of.mif.estimate), 
             color = "red", size = 2) + 
  labs(x =  TeX("$\\log(\\r^n \\cdot f^n_{S})$"), y = "log likelihood") +
  theme(axis.text = element_text(size = 11),
        axis.title = element_text(size = 11)) +
  ylim(-710, -699)+
  # xlim(-2.2,0)+
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

rn_f_Sn_p




combine_subset$ri_f_Si = combine_subset$ri * combine_subset$f_Si
combine_subset$log_ri_f_Si <- log(combine_subset$ri_f_Si)
plot(x = combine_subset$log_ri_f_Si, y = combine_subset$loglik)
# combine_subset_ri_f_Si = combine_subset[combine_subset$log_ri_f_Si < -1,]
# combine_subset_ri_f_Si = combine_subset_ri_f_Si[combine_subset_ri_f_Si$log_ri_f_Si > -2.5,]
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

mcap(combine_subset_ri_f_Si$loglik, combine_subset_ri_f_Si$log_ri_f_Si,  level = 0.95, span = 0.7, Ngrid = 1000) -> mcap_object_ri_f_Si
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
  ylim(-710, -699)+
  # xlim(-2.2,0)+
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

ri_f_Si_p





grid.arrange( rn_p, ri_p, f_Sn_p, f_Si_p,theta_Sn_p,
               theta_Si_p, theta_Jn_p,theta_Ji_p,sigJn_p,sigJi_p,
              sigF_p,k_Sn_p,k_Si_p,rn_f_Sn_p,ri_f_Si_p,
              nrow = 3, ncol = 5)

g <- arrangeGrob( rn_p, ri_p, f_Sn_p, f_Si_p,theta_Sn_p,
                  theta_Si_p, theta_Jn_p,theta_Ji_p,sigJn_p,sigJi_p,
                  sigF_p,k_Sn_p,k_Si_p,rn_f_Sn_p,ri_f_Si_p,
                  nrow = 3, ncol = 5)

ggsave(
  filename = "./daphnia-article/si/profile/Target_dynamics/no_para/Profile_plot.png",
  plot     = g,         
  width    = 16,        
  height   = 8,         
  dpi      = 300,       
  units    = "in"
)

save(subset_data_ri,
     subset_data_rn,
     subset_data_f_Si,
     subset_data_f_Sn,
     subset_data_theta_Sn,
     subset_data_theta_Si,
     subset_data_theta_Ji,
     subset_data_theta_Jn,
     subset_data_sigJi,
     subset_data_sigJn,
     subset_data_sigF,
     subset_data_k_Si,
     subset_data_k_Sn,
     combine_subset_ri_f_Si,
     combine_subset_rn_f_Sn, file = "data/Target_dynamics/no_para/profile_graph_data.rda")


save(mcap_object_ri,
     mcap_object_rn,
     mcap_object_f_Si,
     mcap_object_f_Sn,
     mcap_object_theta_Sn,
     mcap_object_theta_Si,
     mcap_object_theta_Ji,
     mcap_object_theta_Jn,
     mcap_object_sigJi,
     mcap_object_sigJn,
     mcap_object_sigF,
     mcap_object_k_Si,
     mcap_object_k_Sn,
     mcap_object_ri_f_Si,
     mcap_object_rn_f_Sn, file = "data/Target_dynamics/no_para/profile_mcap_object.rda")


ci_table = cbind(mcap_object_ri$ci,
                 mcap_object_rn$ci,
                 mcap_object_f_Si$ci,
                 mcap_object_f_Sn$ci,
                 mcap_object_theta_Sn$ci,
                 mcap_object_theta_Si$ci,
                 mcap_object_theta_Ji$ci,
                 mcap_object_theta_Jn$ci,
                 mcap_object_sigJi$ci,
                 mcap_object_sigJn$ci,
                 mcap_object_sigF$ci,
                 mcap_object_k_Si$ci,
                 mcap_object_k_Sn$ci,
                 mcap_object_ri_f_Si$ci,
                 mcap_object_rn_f_Sn$ci)

ci_table = rbind(ci_table,c(mcap_object_ri$mle,
                            mcap_object_rn$mle,
                            mcap_object_f_Si$mle,
                            mcap_object_f_Sn$mle,
                            mcap_object_theta_Sn$mle,
                            mcap_object_theta_Si$mle,
                            mcap_object_theta_Ji$mle,
                            mcap_object_theta_Jn$mle,
                            mcap_object_sigJi$mle,
                            mcap_object_sigJn$mle,
                            mcap_object_sigF$mle,
                            mcap_object_k_Si$mle,
                            mcap_object_k_Sn$mle,
                            mcap_object_ri_f_Si$mle,
                            mcap_object_rn_f_Sn$mle))


rownames(ci_table) = c("2.5%","97.5%","MLE")
colnames(ci_table) = c('ri','rn','f_Si','f_Sn','theta_Sn','theta_Si',
                       'theta_Ji','theta_Jn','sigJi','sigJn','sigF','k_Si','k_Sn','ri_f_Si','rn_f_Sn')
ci_table = as.data.frame(ci_table)
save(ci_table, file = 'data/Target_dynamics/no_para/profile_ci_table.rda')
















































































