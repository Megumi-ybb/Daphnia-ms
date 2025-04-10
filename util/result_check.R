true_mf = mf
true_mf1 = mf1
true_mf2 = mf2

mf = true_mf


mifLikes <- as.data.frame(traces(mf[[best]]$mif))
mifLikes$iter = as.numeric(rownames(mifLikes))
mifLikes$mif = 1
final_likes <- numeric(length(mf))
final_likes[1] <- mf[[1]]$ll[1]
for (i in 2:length(mf)) {
  tmp_df <- as.data.frame(traces(mf[[i]]$mif))
  tmp_df$iter = as.numeric(rownames(tmp_df))
  tmp_df$mif = i
  mifLikes <- dplyr::bind_rows(mifLikes, tmp_df)
  final_likes[i] <- mf[[i]]$ll[1]
}

remove_inf <- mifLikes[mifLikes$loglik != (-Inf),]
remove_inf <- remove_inf[!is.na(remove_inf$loglik),]
mifLikes = remove_inf
mifLikes%>%
  ggplot() +
  geom_point(aes(x = iter, y = loglik, group = as.factor(mif), color = as.factor(mif))) +
  theme(legend.position = 'none') +
  # geom_smooth(aes(x = iter, y = sigI))
  geom_smooth(aes(x = iter, y = loglik))


New_mf = mifLikes[mifLikes$iter >= 250,]
# New_mf = New_mf[New_mf$loglik >= -375,]

New_mf%>%
  ggplot() +
  geom_point(aes(x = iter, y = loglik, group = as.factor(mif), color = as.factor(mif))) +
  theme(legend.position = 'none') +
  geom_smooth(aes(x = iter, y = loglik))


trace <- as.data.frame(traces(mf[[1]]$mif))
trace$iter = as.numeric(rownames(trace))

for (i in 2 : length(mf)){
  tmp_df <- as.data.frame(traces(mf[[i]]$mif))
  tmp_df$iter = as.numeric(rownames(tmp_df))
  trace <- dplyr::bind_rows(trace,tmp_df)
  # print(i)
}

final_likes <- numeric(length(mf))

for (i in 1: length(mf)){
  final_likes[i] <- mf[[i]]$ll[1]
}

final_params <- trace %>%
  dplyr::filter(iter == max(iter, na.rm = TRUE))

final_params$loglik <- final_likes

GGally::ggpairs(final_params, columns = colnames(final_params)[1:27])

quantiles <- quantile(final_params$betai, probs = c(0.1, 0.9))

# Subset the data to get rows where betai is at the 2.5th or 97.5th quantiles
subset_data <- final_params[final_params$betai >= quantiles[1] & final_params$betai <= quantiles[2], ]
subset_data = subset_data[subset_data$loglik >= max(subset_data$loglik - 25),]
GGally::ggpairs(subset_data, columns = colnames(subset_data)[c(41, 1)])


lls <- matrix(unlist(sapply(mf, getElement, "ll")), nrow = 2)
best <- which.max(lls[1,])
mif.estimate <- coef(mf[[best]]$mif)
pf.loglik.of.mif.estimate <- unname(mf[[best]]$ll[1])
s.e.of.pf.loglik.of.mif.estimate <- unname(mf[[best]]$ll[2])

pf.loglik.of.mif.estimate
s.e.of.pf.loglik.of.mif.estimate

coef(panelfood) <- coef(mf[[best]]$mif)
all_params <- pparams(mf[[best]]$mif)

foreach(u = names(panelfood), .combine = rbind) %do% {
  unit_model <- unit_objects(panelfood)[[u]]
  shared <- coef(panelfood)
  # unit_specific <- all_params$specific[, u]
  # names(unit_specific) = "sigS"
  pomp::simulate(
    unit_model, 
    nsim = 20,
    format = "data.frame",
    # params = c(shared, unit_specific)
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

# all_sims %>%
#   # dplyr::filter(unit == "u4") %>%
#   ggplot(aes(x = day, y = dentinf, group = .id)) + 
#   geom_line() + 
#   facet_wrap(~unit, scales = "free_y")


ggplot() +
  geom_line(data = all_sims, aes(x = day, y = Ii, group = .id), col = 'blue') +
  geom_line(data = data_df, aes(x = day, y = luminf), col = 'red') +
  facet_wrap(~unit, scales = "free_y")
# all_sims %>%
#   filter(error_count == 0) -> large_removed 


all_sims[all_sims$error_count == 0,] -> large_removed

ggplot() +
  geom_line(data = large_removed, aes(x = day, y = F, group = .id), col = 'blue') +
  # geom_line(data = data_df, aes(x = day, y = dentadult), col = 'red') +
  facet_wrap(~unit, scales = "free_y")

mcap(final_params$loglik, final_params$k_Ii,  level = 0.95, span = 0.75, Ngrid = 10000)-> alpha_mcap
alpha_mcap$mle
ggplot() +
  # xlim(0, 100)+
  geom_point(data = final_params, aes(x = alpha, y = loglik)) +
  geom_line(data = alpha_mcap$fit, aes(x = parameter, y = smoothed), col = 'red') +
  geom_vline(xintercept = alpha_mcap$ci[1], linetype = 'dashed') +
  geom_vline(xintercept = alpha_mcap$ci[2], linetype = 'dashed') +
  geom_vline(xintercept = alpha_mcap$mle, col = 'blue') +
  labs(x = "Linear Trend in Transmission", y = 'Log Likelihood') +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 10))

library(dplyr)

midpoints <- function(x) {
  levels_numeric <- strsplit(sub("\\((.+),(.+)\\]", "\\1 \\2", as.character(x)), " ")
  levels_numeric <- do.call("rbind", lapply(levels_numeric, as.numeric))
  mean(levels_numeric[1])
  # rowMeans(levels_numeric)
}

final_params_tmp <- final_params
final_params_tmp$r = cut(final_params_tmp$r, breaks = 25)

final_params_tmp <- final_params_tmp %>%
  group_by(r) %>% 
  summarize(loglik = max(loglik),
            r = midpoints(r))
  
pomp::mcap(final_params$loglik, final_params$alpha) -> alpha_mcap
alpha_mcap$alpha
ggplot() +
  geom_point(data = final_params_tmp, aes(x = r, y = loglik)) +
  geom_line(data = alpha_mcap$fit, aes(x = parameter, y = smoothed), col = 'red') +
  geom_vline(xintercept = alpha_mcap$ci[1], linetype = 'dashed') +
  geom_vline(xintercept = alpha_mcap$ci[2], linetype = 'dashed') +
  geom_vline(xintercept = alpha_mcap$mle, col = 'blue') +
  labs(x = "Linear Trend in Transmission", y = 'Log Likelihood') +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 10))

