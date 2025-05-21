library(dplyr)
library(latex2exp)
library(gridExtra)
library(dplyr)
library(pomp)
library(ggplot2)


#Please set the working to be the 'Daphnia-ms' path

load("data/Simple_dynamics/Dent/para/profile_graph_data.rda")
load("./Single-species/Dent/SIRJPF/model/best_result.rda")

load_option = FALSE


if(load_option){
  load("./Single-species/Dent/SIRJPF/profile/f_Sn/f_Sn.RData")
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
  ylim(-580, -565)+
  # xlim(-9.5,-5.5)+
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())  + 
  annotate("text", x = mcap_object_f_Sn$mle, y = -900, label = sprintf("f_Sn_mle: %s", formatC(mcap_object_f_Sn$mle, format = 'e', digits = 3)), hjust = 1.05, vjust = -0.5, size = 3)

f_Sn_p







if(load_option){
  load("./Single-species/Dent/SIRJPF/profile/k_In/k_In.RData")
  subset_data_k_In <- final_params %>%
    group_by(k_In) %>%
    filter(loglik == max(loglik))
}

plot(x = log(subset_data_k_In$k_In), y = subset_data_k_In$loglik)
plot(x = subset_data_k_In$k_In, y = subset_data_k_In$loglik)
subset_data_k_In$log_k_In <- log(subset_data_k_In$k_In)

mcap(subset_data_k_In$loglik, subset_data_k_In$log_k_In,  level = 0.95, span = 0.95, Ngrid = 1000) -> mcap_object_k_In
mcap_object_k_In$mle -> k_In_mle
k_In_p <- ggplot() +
  geom_point(data = subset_data_k_In, aes(x = log_k_In, y = loglik)) +
  geom_line(data = mcap_object_k_In$fit, aes(x = parameter, y = smoothed), col = 'red') +
  geom_vline(xintercept = mcap_object_k_In$ci[1], linetype = 'dashed') +
  geom_vline(xintercept = mcap_object_k_In$ci[2], linetype = 'dashed') +
  geom_vline(xintercept = mcap_object_k_In$mle, col = 'blue') +
  geom_vline(xintercept = log(mif.estimate[['k_In']]), col = 'red') +
  geom_hline(yintercept = pf.loglik.of.mif.estimate, col = 'red',linetype = "longdash") +
  geom_point(aes(x = log(mif.estimate[['k_In']]), 
                 y = pf.loglik.of.mif.estimate), 
             color = "red", size = 2) + 
  labs(x =  TeX("$\\log(\\tau^n_{I})$"), y = "log likelihood") +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 10)) +
  ylim(-580, -565)+
  xlim(-1,1.5)+
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())  + 
  annotate("text", x = mcap_object_k_In$mle, y = -900, label = sprintf("k_In_mle: %s", formatC(mcap_object_k_In$mle, format = 'e', digits = 3)), hjust = 1.05, vjust = -0.5, size = 3)

k_In_p






if(load_option){
  load("./Single-species/Dent/SIRJPF/profile/k_Sn/k_Sn.RData")
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
  ylim(-580, -565)+
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())  + 
  annotate("text", x = mcap_object_k_Sn$mle, y = -900, label = sprintf("k_Sn_mle: %s", formatC(mcap_object_k_Sn$mle, format = 'e', digits = 3)), hjust = 1.05, vjust = -0.5, size = 3)

k_Sn_p




if(load_option){
  load("./Single-species/Dent/SIRJPF/profile/probn/probn.RData")
  subset_data_probn <- final_params %>%
    group_by(probn) %>%
    filter(loglik == max(loglik))
}

plot(x = log(subset_data_probn$probn), y = subset_data_probn$loglik)
plot(x = subset_data_probn$probn, y = subset_data_probn$loglik)
subset_data_probn$log_probn <- log(subset_data_probn$probn)

mcap(subset_data_probn$loglik, subset_data_probn$log_probn,  level = 0.95, span = 0.95, Ngrid = 1000) -> mcap_object_probn
mcap_object_probn$mle -> probn_mle
probn_p <- ggplot() +
  geom_point(data = subset_data_probn, aes(x = log_probn, y = loglik)) +
  geom_line(data = mcap_object_probn$fit, aes(x = parameter, y = smoothed), col = 'red') +
  geom_vline(xintercept = mcap_object_probn$ci[1], linetype = 'dashed') +
  geom_vline(xintercept = mcap_object_probn$ci[2], linetype = 'dashed') +
  geom_vline(xintercept = mcap_object_probn$mle, col = 'blue') +
  geom_vline(xintercept = log(mif.estimate[['probn']]), col = 'red') +
  geom_hline(yintercept = pf.loglik.of.mif.estimate, col = 'red',linetype = "longdash") +
  geom_point(aes(x = log(mif.estimate[['probn']]), 
                 y = pf.loglik.of.mif.estimate), 
             color = "red", size = 2) + 
  labs(x =  TeX("$p^n$"), y = "log likelihood") +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 10)) +
  ylim(-580, -565)+
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())  + 
  annotate("text", x = mcap_object_probn$mle, y = -900, label = sprintf("probn_mle: %s", formatC(mcap_object_probn$mle, format = 'e', digits = 3)), hjust = 1.05, vjust = -0.5, size = 3)

probn_p








if(load_option){
  load("./Single-species/Dent/SIRJPF/profile/rn/rn.RData")
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
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 10)) +
  ylim(-580, -565)+
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())  + 
  annotate("text", x = mcap_object_rn$mle, y = -900, label = sprintf("rn_mle: %s", formatC(mcap_object_rn$mle, format = 'e', digits = 3)), hjust = 1.05, vjust = -0.5, size = 3)

rn_p




if(load_option){
  load("./Single-species/Dent/SIRJPF/profile/sigF/sigF.RData")
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
  ylim(-580, -565)+
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())  + 
  annotate("text", x = mcap_object_sigF$mle, y = -900, label = sprintf("sigF_mle: %s", formatC(mcap_object_sigF$mle, format = 'e', digits = 3)), hjust = 1.05, vjust = -0.5, size = 3)

sigF_p






if(load_option){
  load("./Single-species/Dent/SIRJPF/profile/sigIn/sigIn.RData")
  subset_data_sigIn <- final_params %>%
    group_by(sigIn) %>%
    filter(loglik == max(loglik))
}

plot(x = log(subset_data_sigIn$sigIn), y = subset_data_sigIn$loglik)
plot(x = subset_data_sigIn$sigIn, y = subset_data_sigIn$loglik)
subset_data_sigIn$log_sigIn <- log(subset_data_sigIn$sigIn)

mcap(subset_data_sigIn$loglik, subset_data_sigIn$log_sigIn,  level = 0.95, span = 0.95, Ngrid = 1000) -> mcap_object_sigIn
mcap_object_sigIn$mle -> sigIn_mle
sigIn_p <- ggplot() +
  geom_point(data = subset_data_sigIn, aes(x = log_sigIn, y = loglik)) +
  geom_line(data = mcap_object_sigIn$fit, aes(x = parameter, y = smoothed), col = 'red') +
  geom_vline(xintercept = mcap_object_sigIn$ci[1], linetype = 'dashed') +
  geom_vline(xintercept = mcap_object_sigIn$ci[2], linetype = 'dashed') +
  geom_vline(xintercept = mcap_object_sigIn$mle, col = 'blue') +
  geom_vline(xintercept = log(mif.estimate[['sigIn']]), col = 'red') +
  geom_hline(yintercept = pf.loglik.of.mif.estimate, col = 'red',linetype = "longdash") +
  geom_point(aes(x = log(mif.estimate[['sigIn']]), 
                 y = pf.loglik.of.mif.estimate), 
             color = "red", size = 2) + 
  labs(x =  TeX("$\\log(\\sigma^n_{I})$"), y = "log likelihood") +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 10)) +
  ylim(-580, -565)+
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

sigIn_p





if(load_option){
  load("./Single-species/Dent/SIRJPF/profile/sigJn/sigJn.RData")
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
  ylim(-580, -565)+
  xlim(-3,0)+
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())  + 
  annotate("text", x = mcap_object_sigJn$mle, y = -900, label = sprintf("sigJn_mle: %s", formatC(mcap_object_sigJn$mle, format = 'e', digits = 3)), hjust = 1.05, vjust = -0.5, size = 3)

sigJn_p




if(load_option){
  load("./Single-species/Dent/SIRJPF/profile/sigP/sigP.RData")
  subset_data_sigP <- final_params %>%
    group_by(sigP) %>%
    filter(loglik == max(loglik))
}

plot(x = log(subset_data_sigP$sigP), y = subset_data_sigP$loglik)
plot(x = subset_data_sigP$sigP, y = subset_data_sigP$loglik)
subset_data_sigP$log_sigP <- log(subset_data_sigP$sigP)
subset_data_sigP = subset_data_sigP[subset_data_sigP$loglik > -600,]
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
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 10)) +
  ylim(-580, -565)+
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
sigP_p





if(load_option){
  load("./Single-species/Dent/SIRJPF/profile/theta_In/theta_In.RData")
  subset_data_theta_In <- final_params %>%
    group_by(theta_In) %>%
    filter(loglik == max(loglik))
}

plot(x = log(subset_data_theta_In$theta_In), y = subset_data_theta_In$loglik)
plot(x = subset_data_theta_In$theta_In, y = subset_data_theta_In$loglik)
subset_data_theta_In$log_theta_In <- log(subset_data_theta_In$theta_In)
# 
subset_data_theta_In = subset_data_theta_In[subset_data_theta_In$log_theta_In <=1,]
mcap(subset_data_theta_In$loglik, subset_data_theta_In$log_theta_In,  level = 0.95, span = 0.95, Ngrid = 1000) -> mcap_object_theta_In
mcap_object_theta_In$mle -> theta_In_mle
theta_In_p <- ggplot() +
  geom_point(data = subset_data_theta_In, aes(x = log_theta_In, y = loglik)) +
  geom_line(data = mcap_object_theta_In$fit, aes(x = parameter, y = smoothed), col = 'red') +
  geom_vline(xintercept = mcap_object_theta_In$ci[1], linetype = 'dashed') +
  geom_vline(xintercept = mcap_object_theta_In$ci[2], linetype = 'dashed') +
  geom_vline(xintercept = mcap_object_theta_In$mle, col = 'blue') +
  geom_vline(xintercept = log(mif.estimate[['theta_In']]), col = 'red') +
  geom_hline(yintercept = pf.loglik.of.mif.estimate, col = 'red',linetype = "longdash") +
  geom_point(aes(x = log(mif.estimate[['theta_In']]), 
                 y = pf.loglik.of.mif.estimate), 
             color = "red", size = 2) + 
  labs(x =  TeX("$\\log(\\theta^n_{I})$"), y = "log likelihood") +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 10)) +
  ylim(-580, -565)+
  xlim(-2,0.5)+
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())  + 
  annotate("text", x = mcap_object_theta_In$mle, y = -900, label = sprintf("theta_In_mle: %s", formatC(mcap_object_theta_In$mle, format = 'e', digits = 3)), hjust = 1.05, vjust = -0.5, size = 3)

theta_In_p






if(load_option){
  load("./Single-species/Dent/SIRJPF/profile/theta_Jn/theta_Jn.RData")
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
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 10)) +
  ylim(-580, -565)+
  # xlim(-2.2,0)+
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())  + 
  annotate("text", x = mcap_object_theta_Jn$mle, y = -900, label = sprintf("theta_Jn_mle: %s", formatC(mcap_object_theta_Jn$mle, format = 'e', digits = 3)), hjust = 1.05, vjust = -0.5, size = 3)

theta_Jn_p






if(load_option){
  load("./Single-species/Dent/SIRJPF/profile/theta_P/theta_P.RData")
  subset_data_theta_P <- final_params %>%
    group_by(theta_P) %>%
    filter(loglik == max(loglik))
}

plot(x = log(subset_data_theta_P$theta_P), y = subset_data_theta_P$loglik)
plot(x = subset_data_theta_P$theta_P, y = subset_data_theta_P$loglik)
subset_data_theta_P$log_theta_P <- log(subset_data_theta_P$theta_P)

mcap(subset_data_theta_P$loglik, subset_data_theta_P$log_theta_P,  level = 0.95, span = 0.95, Ngrid = 1000) -> mcap_object_theta_P
mcap_object_theta_P$mle -> theta_P_mle
theta_P_p <- ggplot() +
  geom_point(data = subset_data_theta_P, aes(x = log_theta_P, y = loglik)) +
  geom_line(data = mcap_object_theta_P$fit, aes(x = parameter, y = smoothed), col = 'red') +
  geom_vline(xintercept = mcap_object_theta_P$ci[1], linetype = 'dashed') +
  geom_vline(xintercept = mcap_object_theta_P$ci[2], linetype = 'dashed') +
  geom_vline(xintercept = mcap_object_theta_P$mle, col = 'blue') +
  geom_vline(xintercept = log(mif.estimate[['theta_P']]), col = 'red') +
  geom_hline(yintercept = pf.loglik.of.mif.estimate, col = 'red',linetype = "longdash") +
  geom_point(aes(x = log(mif.estimate[['theta_P']]), 
                 y = pf.loglik.of.mif.estimate), 
             color = "red", size = 2) + 
  labs(x =  TeX("$\\log(\\theta_{P})$"), y = "log likelihood") +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 10)) +
  ylim(-580, -565)+
  # xlim(-2.2,0)+
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())  + 
  annotate("text", x = mcap_object_theta_P$mle, y = -900, label = sprintf("theta_P_mle: %s", formatC(mcap_object_theta_P$mle, format = 'e', digits = 3)), hjust = 1.05, vjust = -0.5, size = 3)

theta_P_p









if(load_option){
  load("./Single-species/Dent/SIRJPF/profile/theta_Sn/theta_Sn.RData")
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
  ylim(-580, -565)+
  # xlim(-2.2,0)+
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())  + 
  annotate("text", x = mcap_object_theta_Sn$mle, y = -900, label = sprintf("theta_Sn_mle: %s", formatC(mcap_object_theta_Sn$mle, format = 'e', digits = 3)), hjust = 1.05, vjust = -0.5, size = 3)

theta_Sn_p



if(load_option){
  load("./Single-species/Dent/SIRJPF/profile/xi/xi.RData")
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
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 10)) +
  ylim(-580, -565)+
  # xlim(-2.2,0)+
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())  + 
  annotate("text", x = mcap_object_xi$mle, y = -900, label = sprintf("xi_mle: %s", formatC(mcap_object_xi$mle, format = 'e', digits = 3)), hjust = 1.05, vjust = -0.5, size = 3)

xi_p





grid.arrange( rn_p,  f_Sn_p, probn_p, theta_Sn_p, theta_In_p,
              theta_Jn_p, sigIn_p,sigJn_p,sigF_p,sigP_p,
              theta_P_p,xi_p,k_Sn_p,k_In_p,
              nrow = 3, ncol = 5)



save( subset_data_rn,
     subset_data_f_Sn,
     subset_data_probn,
     subset_data_xi,
     subset_data_theta_Sn,
     subset_data_theta_In,
     subset_data_theta_P,
     subset_data_theta_Jn,
     subset_data_sigIn,
     subset_data_sigJn,
     subset_data_sigF,
     subset_data_sigP,
     subset_data_k_In,
     subset_data_k_Sn, file = "data/Simple_dynamics/Dent/para/profile_graph_data.rda")


save(mcap_object_rn,
     mcap_object_f_Sn,
     mcap_object_probn,
     mcap_object_xi,
     mcap_object_theta_Sn,
     mcap_object_theta_In,
     mcap_object_theta_P,
     mcap_object_theta_Jn,
     mcap_object_sigIn,
     mcap_object_sigJn,
     mcap_object_sigF,
     mcap_object_sigP,
     mcap_object_k_In,
     mcap_object_k_Sn, file = "data/Simple_dynamics/Dent/para/profile_mcap_object.rda")


ci_table = cbind(mcap_object_rn$ci,
                 mcap_object_f_Sn$ci,
                 mcap_object_probn$ci,
                 mcap_object_xi$ci,
                 mcap_object_theta_Sn$ci,
                 mcap_object_theta_In$ci,
                 mcap_object_theta_P$ci,
                 mcap_object_theta_Jn$ci,
                 mcap_object_sigIn$ci,
                 mcap_object_sigJn$ci,
                 mcap_object_sigF$ci,
                 mcap_object_sigP$ci,
                 mcap_object_k_In$ci,
                 mcap_object_k_Sn$ci)

ci_table = rbind(ci_table,c(mcap_object_rn$mle,
                            mcap_object_f_Sn$mle,
                            mcap_object_probn$mle,
                            mcap_object_xi$mle,
                            mcap_object_theta_Sn$mle,
                            mcap_object_theta_In$mle,
                            mcap_object_theta_P$mle,
                            mcap_object_theta_Jn$mle,
                            mcap_object_sigIn$mle,
                            mcap_object_sigJn$mle,
                            mcap_object_sigF$mle,
                            mcap_object_sigP$mle,
                            mcap_object_k_In$mle,
                            mcap_object_k_Sn$mle))


rownames(ci_table) = c("2.5%","97.5%","MLE")
colnames(ci_table) = c('rn','f_Sn','probn','xi','theta_Sn','theta_In','theta_P',
                       'theta_Jn','sigIn','sigJn','sigF','sigP','k_In','k_Sn')
ci_table = as.data.frame(ci_table)
save(ci_table, file = 'data/Simple_dynamics/Dent/para/profile_ci_table.rda')
















































































