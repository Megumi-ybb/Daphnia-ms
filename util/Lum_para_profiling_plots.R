library(dplyr)
library(latex2exp)
library(gridExtra)
library(dplyr)
library(pomp)
library(ggplot2)

#Please set the working to be the 'Daphnia-ms' path

# 
load("data/Simple_dynamics/Lum/para/profile_graph_data.rda")
load("./Single-species/Lum/SIRJPF/model/best_result.rda")

load_option = FALSE

if(load_option){
  load("./Single-species/Lum/SIRJPF/profile/f_Si/f_Si.RData")
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
  geom_vline(xintercept = mcap_object_f_Si$ci[1], linetype = 'dashed') +
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
  ylim(-435, -425)+
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())  + 
  annotate("text", x = mcap_object_f_Si$mle, y = -900, label = sprintf("f_Si_mle: %s", formatC(mcap_object_f_Si$mle, format = 'e', digits = 3)), hjust = 1.05, vjust = -0.5, size = 3)

f_Si_p





if(load_option){
  load("./Single-species/Lum/SIRJPF/profile/k_Ii/k_Ii.RData")
  subset_data_k_Ii <- final_params %>%
    group_by(k_Ii) %>%
    filter(loglik == max(loglik))
}

plot(x = log(subset_data_k_Ii$k_Ii), y = subset_data_k_Ii$loglik)
plot(x = subset_data_k_Ii$k_Ii, y = subset_data_k_Ii$loglik)
subset_data_k_Ii$log_k_Ii <- log(subset_data_k_Ii$k_Ii)

mcap(subset_data_k_Ii$loglik, subset_data_k_Ii$log_k_Ii,  level = 0.95, span = 0.95, Ngrid = 1000) -> mcap_object_k_Ii
mcap_object_k_Ii$mle -> k_Ii_mle
k_Ii_p <- ggplot() +
  geom_point(data = subset_data_k_Ii, aes(x = log_k_Ii, y = loglik)) +
  geom_line(data = mcap_object_k_Ii$fit, aes(x = parameter, y = smoothed), col = 'red') +
  geom_vline(xintercept = mcap_object_k_Ii$ci[1], linetype = 'dashed') +
  geom_vline(xintercept = mcap_object_k_Ii$ci[2], linetype = 'dashed') +
  geom_vline(xintercept = mcap_object_k_Ii$mle, col = 'blue') +
  geom_vline(xintercept = log(mif.estimate[['k_Ii']]), col = 'red') +
  geom_hline(yintercept = pf.loglik.of.mif.estimate, col = 'red',linetype = "longdash") +
  geom_point(aes(x = log(mif.estimate[['k_Ii']]), 
                 y = pf.loglik.of.mif.estimate), 
             color = "red", size = 2) + 
  labs(x =  TeX("$\\log(\\tau^i_{I})$"), y = "log likelihood") +
  theme(axis.text = element_text(size = 11),
        axis.title = element_text(size = 11)) +
  ylim(-435, -425)+
  # xlim(-1,1.5)+
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())  + 
  annotate("text", x = mcap_object_k_Ii$mle, y = -900, label = sprintf("k_Ii_mle: %s", formatC(mcap_object_k_Ii$mle, format = 'e', digits = 3)), hjust = 1.05, vjust = -0.5, size = 3)

k_Ii_p






if(load_option){
  load("./Single-species/Lum/SIRJPF/profile/k_Si/k_Si.RData")
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
  # geom_vline(xintercept = mcap_object_k_Si$ci[2], linetype = 'dashed') +
  geom_vline(xintercept = mcap_object_k_Si$mle, col = 'blue') +
  geom_vline(xintercept = log(mif.estimate[['k_Si']]), col = 'red') +
  geom_hline(yintercept = pf.loglik.of.mif.estimate, col = 'red',linetype = "longdash") +
  geom_point(aes(x = log(mif.estimate[['k_Si']]), 
                 y = pf.loglik.of.mif.estimate), 
             color = "red", size = 2) + 
  labs(x =  TeX("$\\log(\\tau^i_{S})$"), y = "log likelihood") +
  theme(axis.text = element_text(size = 11),
        axis.title = element_text(size = 11)) +
  ylim(-435, -425)+
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
k_Si_p


if(load_option){
  load("./Single-species/Lum/SIRJPF/profile/probi/probi.RData")
  subset_data_probi <- final_params %>%
    group_by(probi) %>%
    filter(loglik == max(loglik))
}

plot(x = log(subset_data_probi$probi), y = subset_data_probi$loglik)
plot(x = subset_data_probi$probi, y = subset_data_probi$loglik)
subset_data_probi$log_probi <- log(subset_data_probi$probi)

mcap(subset_data_probi$loglik, subset_data_probi$log_probi,  level = 0.95, span = 0.6, Ngrid = 1000) -> mcap_object_probi
mcap_object_probi$mle -> probi_mle
probi_p <- ggplot() +
  geom_point(data = subset_data_probi, aes(x = log_probi, y = loglik)) +
  geom_line(data = mcap_object_probi$fit, aes(x = parameter, y = smoothed), col = 'red') +
  geom_vline(xintercept = mcap_object_probi$ci[1], linetype = 'dashed') +
  geom_vline(xintercept = mcap_object_probi$ci[2], linetype = 'dashed') +
  geom_vline(xintercept = mcap_object_probi$mle, col = 'blue') +
  geom_vline(xintercept = log(mif.estimate[['probi']]), col = 'red') +
  geom_hline(yintercept = pf.loglik.of.mif.estimate, col = 'red',linetype = "longdash") +
  geom_point(aes(x = log(mif.estimate[['probi']]), 
                 y = pf.loglik.of.mif.estimate), 
             color = "red", size = 2) + 
  labs(x =  TeX("$\\log(p^i)$"), y = "log likelihood") +
  theme(axis.text = element_text(size = 11),
        axis.title = element_text(size = 11)) +
  ylim(-435, -425)+
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())  + 
  annotate("text", x = mcap_object_probi$mle, y = -900, label = sprintf("probi_mle: %s", formatC(mcap_object_probi$mle, format = 'e', digits = 3)), hjust = 1.05, vjust = -0.5, size = 3)

probi_p



if(load_option){
  load("./Single-species/Lum/SIRJPF/profile/ri/ri.RData")
  subset_data_ri <- final_params %>%
    group_by(ri) %>%
    filter(loglik == max(loglik))
}

plot(x = log(subset_data_ri$ri), y = subset_data_ri$loglik)
plot(x = subset_data_ri$ri, y = subset_data_ri$loglik)
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
  labs(x =  TeX("$\\log(r^i)$"), y = "log likelihood") +
  theme(axis.text = element_text(size = 11),
        axis.title = element_text(size = 11)) +
  ylim(-435, -425)+
  theme_bw() +
  theme(axis.title.y = element_blank())  + 
  annotate("text", x = mcap_object_ri$mle, y = -900, label = sprintf("ri_mle: %s", formatC(mcap_object_ri$mle, format = 'e', digits = 3)), hjust = 1.05, vjust = -0.5, size = 3)

ri_p





if(load_option){
  load("./Single-species/Lum/SIRJPF/profile/sigF/sigF.RData")
  subset_data_sigF <- final_params %>%
    group_by(sigF) %>%
    filter(loglik == max(loglik))
}

plot(x = log(subset_data_sigF$sigF), y = subset_data_sigF$loglik)
plot(x = subset_data_sigF$sigF, y = subset_data_sigF$loglik)
subset_data_sigF$log_sigF <- log(subset_data_sigF$sigF)

mcap(subset_data_sigF$loglik, subset_data_sigF$log_sigF,  level = 0.95, span = 0.6, Ngrid = 1000) -> mcap_object_sigF
mcap_object_sigF$mle -> sigF_mle
sigF_p <- ggplot() +
  geom_point(data = subset_data_sigF, aes(x = log_sigF, y = loglik)) +
  geom_line(data = mcap_object_sigF$fit, aes(x = parameter, y = smoothed), col = 'red') +
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
  ylim(-435, -425)+
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
sigF_p






if(load_option){
  load("./Single-species/Lum/SIRJPF/profile/sigIi/sigIi.RData")
  subset_data_sigIi <- final_params %>%
    group_by(sigIi) %>%
    filter(loglik == max(loglik))
}

plot(x = log(subset_data_sigIi$sigIi), y = subset_data_sigIi$loglik)
plot(x = subset_data_sigIi$sigIi, y = subset_data_sigIi$loglik)
subset_data_sigIi$log_sigIi <- log(subset_data_sigIi$sigIi)
subset_data_sigIi = subset_data_sigIi[subset_data_sigIi$loglik > max(subset_data_sigIi$loglik) - 15,]
mcap(subset_data_sigIi$loglik, subset_data_sigIi$log_sigIi,  level = 0.95, span = 0.5, Ngrid = 1000) -> mcap_object_sigIi
mcap_object_sigIi$mle -> sigIi_mle
sigIi_p <- ggplot() +
  geom_point(data = subset_data_sigIi, aes(x = log_sigIi, y = loglik)) +
  geom_line(data = mcap_object_sigIi$fit, aes(x = parameter, y = smoothed), col = 'red') +
  # geom_vline(xintercept = mcap_object_sigIi$ci[1], linetype = 'dashed') +
  geom_vline(xintercept = mcap_object_sigIi$ci[2], linetype = 'dashed') +
  geom_vline(xintercept = mcap_object_sigIi$mle, col = 'blue') +
  geom_vline(xintercept = log(mif.estimate[['sigIi']]), col = 'red') +
  geom_hline(yintercept = pf.loglik.of.mif.estimate, col = 'red',linetype = "longdash") +
  geom_point(aes(x = log(mif.estimate[['sigIi']]), 
                 y = pf.loglik.of.mif.estimate), 
             color = "red", size = 2) + 
  labs(x =  TeX("$\\log(\\sigma^i_{I})$"), y = "log likelihood") +
  theme(axis.text = element_text(size = 11),
        axis.title = element_text(size = 11)) +
  ylim(-435, -425)+
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())  + 
  annotate("text", x = mcap_object_sigIi$mle, y = -900, label = sprintf("sigIi_mle: %s", formatC(mcap_object_sigIi$mle, format = 'e', digits = 3)), hjust = 1.05, vjust = -0.5, size = 3)

sigIi_p







if(load_option){
  load("./Single-species/Lum/SIRJPF/profile/sigJi/sigJi.RData")
  subset_data_sigJi <- final_params %>%
    group_by(sigJi) %>%
    filter(loglik == max(loglik))
}

plot(x = log(subset_data_sigJi$sigJi), y = subset_data_sigJi$loglik)
plot(x = subset_data_sigJi$sigJi, y = subset_data_sigJi$loglik)
subset_data_sigJi$log_sigJi <- log(subset_data_sigJi$sigJi)

mcap(subset_data_sigJi$loglik, subset_data_sigJi$log_sigJi,  level = 0.95, span = 0.95, Ngrid = 1000) -> mcap_object_sigJi
mcap_object_sigJi$mle -> sigJi_mle
sigJi_p <- ggplot() +
  geom_point(data = subset_data_sigJi, aes(x = log_sigJi, y = loglik)) +
  geom_line(data = mcap_object_sigJi$fit, aes(x = parameter, y = smoothed), col = 'red') +
  geom_vline(xintercept = mcap_object_sigJi$ci[1], linetype = 'dashed') +
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
  ylim(-435, -425)+
  xlim(-2,-0)+
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())  + 
  annotate("text", x = mcap_object_sigJi$mle, y = -900, label = sprintf("sigJi_mle: %s", formatC(mcap_object_sigJi$mle, format = 'e', digits = 3)), hjust = 1.05, vjust = -0.5, size = 3)

sigJi_p









if(load_option){
  load("./Single-species/Lum/SIRJPF/profile/sigP/sigP.RData")
  subset_data_sigP <- final_params %>%
    group_by(sigP) %>%
    filter(loglik == max(loglik))
}

plot(x = log(subset_data_sigP$sigP), y = subset_data_sigP$loglik)
plot(x = subset_data_sigP$sigP, y = subset_data_sigP$loglik)
subset_data_sigP$log_sigP <- log(subset_data_sigP$sigP)

mcap(subset_data_sigP$loglik, subset_data_sigP$log_sigP,  level = 0.95, span = 0.95, Ngrid = 1000) -> mcap_object_sigP
mcap_object_sigP$mle -> sigP_mle
sigP_p <- ggplot() +
  geom_point(data = subset_data_sigP, aes(x = log_sigP, y = loglik)) +
  geom_line(data = mcap_object_sigP$fit, aes(x = parameter, y = smoothed), col = 'red') +
  geom_vline(xintercept = mcap_object_sigP$ci[1], linetype = 'dashed') +
  geom_vline(xintercept = mcap_object_sigP$ci[2], linetype = 'dashed') +
  geom_vline(xintercept = mcap_object_sigP$mle, col = 'blue') +
  geom_vline(xintercept = log(mif.estimate[['sigP']]), col = 'red') +
  geom_hline(yintercept = pf.loglik.of.mif.estimate, col = 'red',linetype = "longdash") +
  geom_point(aes(x = log(mif.estimate[['sigP']]), 
                 y = pf.loglik.of.mif.estimate), 
             color = "red", size = 2) + 
  labs(x =  TeX("$\\log(\\sigma_{P})$"), y = "log likelihood") +
  theme(axis.text = element_text(size = 11),
        axis.title = element_text(size = 11)) +
  ylim(-435, -425)+
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())  + 
  annotate("text", x = mcap_object_sigP$mle, y = -900, label = sprintf("sigP_mle: %s", formatC(mcap_object_sigP$mle, format = 'e', digits = 3)), hjust = 1.05, vjust = -0.5, size = 3)

sigP_p



if(load_option){
  load("./Single-species/Lum/SIRJPF/profile/theta_Ii/theta_Ii.RData")
  subset_data_theta_Ii <- final_params %>%
    group_by(theta_Ii) %>%
    filter(loglik == max(loglik))
}

plot(x = log(subset_data_theta_Ii$theta_Ii), y = subset_data_theta_Ii$loglik)
plot(x = subset_data_theta_Ii$theta_Ii, y = subset_data_theta_Ii$loglik)
subset_data_theta_Ii$log_theta_Ii <- log(subset_data_theta_Ii$theta_Ii)

mcap(subset_data_theta_Ii$loglik, subset_data_theta_Ii$log_theta_Ii,  level = 0.95, span = 0.95, Ngrid = 1000) -> mcap_object_theta_Ii
mcap_object_theta_Ii$mle -> theta_Ii_mle
theta_Ii_p <- ggplot() +
  geom_point(data = subset_data_theta_Ii, aes(x = log_theta_Ii, y = loglik)) +
  geom_line(data = mcap_object_theta_Ii$fit, aes(x = parameter, y = smoothed), col = 'red') +
  geom_vline(xintercept = mcap_object_theta_Ii$ci[1], linetype = 'dashed') +
  geom_vline(xintercept = mcap_object_theta_Ii$ci[2], linetype = 'dashed') +
  geom_vline(xintercept = mcap_object_theta_Ii$mle, col = 'blue') +
  geom_vline(xintercept = log(mif.estimate[['theta_Ii']]), col = 'red') +
  geom_hline(yintercept = pf.loglik.of.mif.estimate, col = 'red',linetype = "longdash") +
  geom_point(aes(x = log(mif.estimate[['theta_Ii']]), 
                 y = pf.loglik.of.mif.estimate), 
             color = "red", size = 2) + 
  labs(x =  TeX("$\\log(\\theta^i_{I})$"), y = "log likelihood") +
  theme(axis.text = element_text(size = 11),
        axis.title = element_text(size = 11)) +
  ylim(-435, -425)+
  # xlim(-2.2,0)+
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
theta_Ii_p










if(load_option){
  load("./Single-species/Lum/SIRJPF/profile/theta_Ji/theta_Ji.RData")
  subset_data_theta_Ji <- final_params %>%
    group_by(theta_Ji) %>%
    filter(loglik == max(loglik))
}

plot(x = log(subset_data_theta_Ji$theta_Ji), y = subset_data_theta_Ji$loglik)
plot(x = subset_data_theta_Ji$theta_Ji, y = subset_data_theta_Ji$loglik)
subset_data_theta_Ji$log_theta_Ji <- log(subset_data_theta_Ji$theta_Ji)
subset_data_theta_Ji = subset_data_theta_Ji[subset_data_theta_Ji$loglik >= max(subset_data_theta_Ji$loglik) - 10,]
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
  ylim(-435, -425)+
  # xlim(-2.2,0)+
  theme_bw() +
  theme(axis.title.y = element_blank())
theta_Ji_p







if(load_option){
  load("./Single-species/Lum/SIRJPF/profile/theta_P/theta_P.RData")
  subset_data_theta_P <- final_params %>%
    group_by(theta_P) %>%
    filter(loglik == max(loglik))
}

plot(x = log(subset_data_theta_P$theta_P), y = subset_data_theta_P$loglik)
plot(x = subset_data_theta_P$theta_P, y = subset_data_theta_P$loglik)
subset_data_theta_P$log_theta_P <- log(subset_data_theta_P$theta_P)
# subset_data_theta_P = subset_data_theta_P[subset_data_theta_P$log_theta_P > -7,]
mcap(subset_data_theta_P$loglik, subset_data_theta_P$log_theta_P,  level = 0.95, span = 0.6, Ngrid = 1000) -> mcap_object_theta_P
mcap_object_theta_P$mle -> theta_P_mle
theta_P_p <- ggplot() +
  geom_point(data = subset_data_theta_P, aes(x = log_theta_P, y = loglik)) +
  geom_line(data = mcap_object_theta_P$fit, aes(x = parameter, y = smoothed), col = 'red') +
  geom_vline(xintercept = mcap_object_theta_P$ci[2], linetype = 'dashed') +
  geom_vline(xintercept = mcap_object_theta_P$mle, col = 'blue') +
  geom_vline(xintercept = log(mif.estimate[['theta_P']]), col = 'red') +
  geom_hline(yintercept = pf.loglik.of.mif.estimate, col = 'red',linetype = "longdash") +
  geom_point(aes(x = log(mif.estimate[['theta_P']]), 
                 y = pf.loglik.of.mif.estimate), 
             color = "red", size = 2) + 
  labs(x =  TeX("$\\log(\\theta_{P})$"), y = "log likelihood") +
  theme(axis.text = element_text(size = 11),
        axis.title = element_text(size = 11)) +
  ylim(-435, -425)+
  # xlim(-2.2,0)+
  theme_bw() +
  theme(axis.title.y = element_blank())  + 
  annotate("text", x = mcap_object_theta_P$mle, y = -900, label = sprintf("theta_P_mle: %s", formatC(mcap_object_theta_P$mle, format = 'e', digits = 3)), hjust = 1.05, vjust = -0.5, size = 3)

theta_P_p








if(load_option){
  load("./Single-species/Lum/SIRJPF/profile/theta_Si/theta_Si.RData")
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
  theme(axis.text = element_text(size = 11),
        axis.title = element_text(size = 11)) +
  ylim(-435, -425)+
  # xlim(-2.2,0)+
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())  + 
  annotate("text", x = mcap_object_theta_Si$mle, y = -900, label = sprintf("theta_Si_mle: %s", formatC(mcap_object_theta_Si$mle, format = 'e', digits = 3)), hjust = 1.05, vjust = -0.5, size = 3)

theta_Si_p






if(load_option){
  load("./Single-species/Lum/SIRJPF/profile/xi/xi.RData")
  subset_data_xi <- final_params %>%
    group_by(xi) %>%
    filter(loglik == max(loglik))
}

plot(x = log(subset_data_xi$xi), y = subset_data_xi$loglik)
plot(x = subset_data_xi$xi, y = subset_data_xi$loglik)
subset_data_xi$log_xi <- log(subset_data_xi$xi)

mcap(subset_data_xi$loglik, subset_data_xi$log_xi,  level = 0.95, span = 0.95, Ngrid = 1000) -> mcap_object_xi
mcap_object_xi$mle -> xi_mle
xi_p <- ggplot() +
  geom_point(data = subset_data_xi, aes(x = log_xi, y = loglik)) +
  geom_line(data = mcap_object_xi$fit, aes(x = parameter, y = smoothed), col = 'red') +
  geom_vline(xintercept = mcap_object_xi$ci[1], linetype = 'dashed') +
  geom_vline(xintercept = mcap_object_xi$ci[2], linetype = 'dashed') +
  geom_vline(xintercept = mcap_object_xi$mle, col = 'blue') +
  geom_vline(xintercept = log(mif.estimate[['xi']]), col = 'red') +
  geom_hline(yintercept = pf.loglik.of.mif.estimate, col = 'red',linetype = "longdash") +
  geom_point(aes(x = log(mif.estimate[['xi']]), 
                 y = pf.loglik.of.mif.estimate), 
             color = "red", size = 2) + 
  labs(x =  TeX("$\\log(\\xi)$"), y = "log likelihood") +
  theme(axis.text = element_text(size = 11),
        axis.title = element_text(size = 11)) +
  ylim(-435, -425)+
  # xlim(-2.2,0)+
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())  + 
  annotate("text", x = mcap_object_xi$mle, y = -900, label = sprintf("xi_mle: %s", formatC(mcap_object_xi$mle, format = 'e', digits = 3)), hjust = 1.05, vjust = -0.5, size = 3)

xi_p



combine_subset = rbind(subset_data_ri[,1:16],
                       subset_data_f_Si[,1:16],
                       subset_data_probi[,1:16],
                       subset_data_xi[,1:16],
                       subset_data_theta_Si[,1:16],
                       subset_data_theta_Ii[,1:16],
                       subset_data_theta_P[,1:16],
                       subset_data_theta_Ji[,1:16],
                       subset_data_sigIi[,1:16],
                       subset_data_sigJi[,1:16],
                       subset_data_sigF[,1:16],
                       subset_data_sigP[,1:16],
                       subset_data_k_Ii[,1:16],
                       subset_data_k_Si[,1:16])


combine_subset$ri_f_Si = combine_subset$ri * combine_subset$f_Si
combine_subset$log_ri_f_Si <- log(combine_subset$ri_f_Si)
plot(x = combine_subset$log_ri_f_Si, y = combine_subset$loglik)
combine_subset_ri_f_Si = combine_subset[combine_subset$log_ri_f_Si < -1,]
combine_subset_ri_f_Si = combine_subset_ri_f_Si[combine_subset_ri_f_Si$log_ri_f_Si > -2.5,]
combine_subset_ri_f_Si <- combine_subset_ri_f_Si %>%
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
  ylim(-435, -425)+
  # xlim(-2.2,0)+
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

ri_f_Si_p




combine_subset$pi_f_Si <- combine_subset$probi * combine_subset$f_Si
combine_subset$log_pi_f_Si <- log(combine_subset$pi_f_Si)

plot(x = combine_subset$log_pi_f_Si, y = combine_subset$loglik)
combine_subset_pi_f_Si <- combine_subset %>%
  # filter(log_pi_f_Si < -1, log_pi_f_Si > -2.5) %>%
  ungroup() %>%
  mutate(
    bin = cut(
      log_pi_f_Si,
      breaks = seq(
        min(log_pi_f_Si, na.rm = TRUE),
        max(log_pi_f_Si, na.rm = TRUE),
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

mcap(
  combine_subset_pi_f_Si$loglik,
  combine_subset_pi_f_Si$log_pi_f_Si,
  level = 0.95,
  span  = 0.7,
  Ngrid = 1000
) -> mcap_object_pi_f_Si

mcap_object_pi_f_Si$mle -> pi_f_Si_mle

pi_f_Si_p <- ggplot() +
  geom_point(
    data = combine_subset_pi_f_Si,
    aes(x = log_pi_f_Si, y = loglik)
  ) +
  geom_line(
    data = mcap_object_pi_f_Si$fit,
    aes(x = parameter, y = smoothed),
    col = "red"
  ) +
  geom_vline(
    xintercept = mcap_object_pi_f_Si$ci[1],
    linetype   = "dashed"
  ) +
  geom_vline(
    xintercept = mcap_object_pi_f_Si$ci[2],
    linetype   = "dashed"
  ) +
  geom_vline(
    xintercept = mcap_object_pi_f_Si$mle,
    col        = "blue"
  ) +
  geom_vline(
    xintercept = log(mif.estimate[['probi']] * mif.estimate[['f_Si']]),
    col        = "red"
  ) +
  geom_hline(
    yintercept = pf.loglik.of.mif.estimate,
    col        = "red",
    linetype   = "longdash"
  ) +
  geom_point(
    aes(
      x = log(mif.estimate[['probi']] * mif.estimate[['f_Si']]),
      y = pf.loglik.of.mif.estimate
    ),
    color = "red",
    size  = 2
  ) +
  labs(
    x = TeX("$\\log(p^i \\cdot f^i_{S})$"),
    y = "log likelihood"
  ) +
  theme(
    axis.text  = element_text(size = 11),
    axis.title = element_text(size = 11)
  ) +
  ylim(-435, -425) +
  theme_bw() +
  theme(
    axis.title.y  = element_blank()
  )

pi_f_Si_p


grid.arrange( ri_p,  f_Si_p,probi_p,
              theta_Si_p,  theta_Ii_p,  theta_Ji_p,
              sigIi_p,sigJi_p,sigF_p,sigP_p,
              theta_P_p,xi_p,k_Si_p,k_Ii_p,ri_f_Si_p,pi_f_Si_p,
              nrow = 4, ncol = 5)


g <- arrangeGrob(  ri_p,  f_Si_p,probi_p,
              theta_Si_p,  theta_Ii_p,  theta_Ji_p,
              sigIi_p,sigJi_p,sigF_p,sigP_p,
              theta_P_p,xi_p,k_Si_p,k_Ii_p,ri_f_Si_p,pi_f_Si_p,
              nrow = 4, ncol = 5)

ggsave(
  filename = "./daphnia-article/si/profile/Simple_dynamics/Lum/para/Profile_plot.png",
  plot     = g,         
  width    = 16,        
  height   = 8,         
  dpi      = 300,       
  units    = "in"
)


 save(subset_data_ri,
     subset_data_f_Si,
     subset_data_probi,
     subset_data_xi,
     subset_data_theta_Si,
     subset_data_theta_Ii,
     subset_data_theta_P,
     subset_data_theta_Ji,
     subset_data_sigIi,
     subset_data_sigJi,
     subset_data_sigF,
     subset_data_sigP,
     subset_data_k_Ii,
     subset_data_k_Si,
     combine_subset_pi_f_Si,
     combine_subset_ri_f_Si, file = "data/Simple_dynamics/Lum/para/profile_graph_data.rda")


save(mcap_object_ri,
     mcap_object_f_Si,
     mcap_object_probi,
     mcap_object_xi,
     mcap_object_theta_Si,
     mcap_object_theta_Ii,
     mcap_object_theta_P,
     mcap_object_theta_Ji,
     mcap_object_sigIi,
     mcap_object_sigJi,
     mcap_object_sigF,
     mcap_object_sigP,
     mcap_object_k_Ii,
     mcap_object_k_Si,
     mcap_object_pi_f_Si,
     mcap_object_ri_f_Si,file = "data/Simple_dynamics/Lum/para/profile_mcap_object.rda")


ci_table = cbind(mcap_object_ri$ci,
                 mcap_object_f_Si$ci,
                 mcap_object_probi$ci,
                 mcap_object_xi$ci,
                 mcap_object_theta_Si$ci,
                 mcap_object_theta_Ii$ci,
                 mcap_object_theta_P$ci,
                 mcap_object_theta_Ji$ci,
                 mcap_object_sigIi$ci,
                 mcap_object_sigJi$ci,
                 mcap_object_sigF$ci,
                 mcap_object_sigP$ci,
                 mcap_object_k_Ii$ci,
                 mcap_object_k_Si$ci,
                 mcap_object_pi_f_Si$ci,
                 mcap_object_ri_f_Si$ci)

ci_table = rbind(ci_table,c(mcap_object_ri$mle,
                            mcap_object_f_Si$mle,
                            mcap_object_probi$mle,
                            mcap_object_xi$mle,
                            mcap_object_theta_Si$mle,
                            mcap_object_theta_Ii$mle,
                            mcap_object_theta_P$mle,
                            mcap_object_theta_Ji$mle,
                            mcap_object_sigIi$mle,
                            mcap_object_sigJi$mle,
                            mcap_object_sigF$mle,
                            mcap_object_sigP$mle,
                            mcap_object_k_Ii$mle,
                            mcap_object_k_Si$mle,
                            mcap_object_pi_f_Si$mle,
                            mcap_object_ri_f_Si$mle))


rownames(ci_table) = c("2.5%","97.5%","MLE")
colnames(ci_table) = c('ri','f_Si','probi','xi','theta_Si','theta_Ii','theta_P',
                       'theta_Ji','sigIi','sigJi','sigF','sigP','k_Ii','k_Si','pi_f_Si','ri_f_Si')
ci_table = as.data.frame(ci_table)
save(ci_table, file = 'data/Simple_dynamics/Lum/para/profile_ci_table.rda')
















































































