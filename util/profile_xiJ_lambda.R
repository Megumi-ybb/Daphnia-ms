library(dplyr)
library(latex2exp)
library(gridExtra)
library(dplyr)
library(pomp)
library(ggplot2)



subset_data_lambda_Jn <- final_params %>%
  group_by(lambda_Jn) %>%
  filter(loglik == max(loglik))

sorted_data <- subset_data_lambda_Jn %>% arrange(lambda_Jn)
rows_per_chunk <- 28

selected_rows <- data.frame()
for (i in seq(1, nrow(sorted_data), by = rows_per_chunk)) {
  chunk <- sorted_data[i:min(i + rows_per_chunk - 1, nrow(sorted_data)), ]
  if (nrow(chunk) > 0) {
    best_row <- chunk[which.max(chunk$loglik), ]
    selected_rows <- rbind(selected_rows, best_row)
  }
}
plot(x = log(selected_rows$lambda_Jn), y = selected_rows$loglik)
plot(x = selected_rows$lambda_Jn, y = selected_rows$loglik)
# selected_rows = selected_rows[selected_rows$lambda_Jn < 1e-4,]
selected_rows$log_lambda_Jn <- log(selected_rows$lambda_Jn)
log_lambda_Jn_range <- range(selected_rows$log_lambda_Jn, na.rm = TRUE)
log_lambda_Jn_breaks <- seq(from = log_lambda_Jn_range[1], to = log_lambda_Jn_range[2], length.out = 51)
max_loglik_points <- data.frame()

for (i in 1:50) {
  chunk <- selected_rows %>%
    filter(log_lambda_Jn >= log_lambda_Jn_breaks[i] & log_lambda_Jn < log_lambda_Jn_breaks[i + 1])
  if (nrow(chunk) > 0) {
    best_row <- chunk[which.max(chunk$loglik), ]
    max_loglik_points <- rbind(max_loglik_points, best_row)
  }
}
mcap(max_loglik_points$loglik, max_loglik_points$lambda_Jn,  level = 0.95, span = 0.75, Ngrid = 1000) -> mcap_object_lambda_Jn
mcap_object_lambda_Jn$mle -> lambda_Jn_mle
lambda_Jn_p <- ggplot() +
  geom_point(data = max_loglik_points, aes(x = log_lambda_Jn, y = loglik)) +
  geom_line(data = mcap_object_lambda_Jn$fit, aes(x = log(parameter), y = smoothed), col = 'red') +
  geom_vline(xintercept = log(mcap_object_lambda_Jn$ci[1]), linetype = 'dashed') +
  geom_vline(xintercept = log(mcap_object_lambda_Jn$ci[2]), linetype = 'dashed') +
  geom_vline(xintercept = log(mcap_object_lambda_Jn$mle), col = 'blue') +
  labs(x = "log(lambda_Jn)", y = "log likelihood") +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 10)) +
  ylim(-578, -565)+
  annotate("text", x = log(mcap_object_lambda_Jn$mle), y = -900, label = sprintf("lambda_Jn_mle: %s", formatC(mcap_object_lambda_Jn$mle, format = 'e', digits = 3)), hjust = 1.05, vjust = -0.5, size = 3)

lambda_Jn_p


subset_data_lambda_Ji <- final_params %>%
  group_by(lambda_Ji) %>%
  filter(loglik == max(loglik))

sorted_data <- subset_data_lambda_Ji %>% arrange(lambda_Ji)
rows_per_chunk <- 28

selected_rows <- data.frame()
for (i in seq(1, nrow(sorted_data), by = rows_per_chunk)) {
  chunk <- sorted_data[i:min(i + rows_per_chunk - 1, nrow(sorted_data)), ]
  if (nrow(chunk) > 0) {
    best_row <- chunk[which.max(chunk$loglik), ]
    selected_rows <- rbind(selected_rows, best_row)
  }
}
plot(x = log(selected_rows$lambda_Ji), y = selected_rows$loglik)
plot(x = selected_rows$lambda_Ji, y = selected_rows$loglik)
# selected_rows = selected_rows[selected_rows$lambda_Ji < 1e-4,]
selected_rows$log_lambda_Ji <- log(selected_rows$lambda_Ji)
log_lambda_Ji_range <- range(selected_rows$log_lambda_Ji, na.rm = TRUE)
log_lambda_Ji_breaks <- seq(from = log_lambda_Ji_range[1], to = log_lambda_Ji_range[2], length.out = 51)
max_loglik_points <- data.frame()

for (i in 1:50) {
  chunk <- selected_rows %>%
    filter(log_lambda_Ji >= log_lambda_Ji_breaks[i] & log_lambda_Ji < log_lambda_Ji_breaks[i + 1])
  if (nrow(chunk) > 0) {
    best_row <- chunk[which.max(chunk$loglik), ]
    max_loglik_points <- rbind(max_loglik_points, best_row)
  }
}
mcap(max_loglik_points$loglik, max_loglik_points$lambda_Ji,  level = 0.95, span = 0.5, Ngrid = 1000) -> mcap_object_lambda_Ji
mcap_object_lambda_Ji$mle -> lambda_Ji_mle
lambda_Ji_p <- ggplot() +
  geom_point(data = max_loglik_points, aes(x = log_lambda_Ji, y = loglik)) +
  geom_line(data = mcap_object_lambda_Ji$fit, aes(x = log(parameter), y = smoothed), col = 'red') +
  geom_vline(xintercept = log(mcap_object_lambda_Ji$ci[1]), linetype = 'dashed') +
  geom_vline(xintercept = log(mcap_object_lambda_Ji$ci[2]), linetype = 'dashed') +
  geom_vline(xintercept = log(mcap_object_lambda_Ji$mle), col = 'blue') +
  labs(x = "log(lambda_Ji)", y = "log likelihood") +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 10)) +
  ylim(-375, -370)+
  annotate("text", x = log(mcap_object_lambda_Ji$mle), y = -900, label = sprintf("lambda_Ji_mle: %s", formatC(mcap_object_lambda_Ji$mle, format = 'e', digits = 3)), hjust = 1.05, vjust = -0.5, size = 3)

lambda_Ji_p


subset_data_xi_J <- final_params %>%
  group_by(xi_J) %>%
  filter(loglik == max(loglik))

sorted_data <- subset_data_xi_J %>% arrange(xi_J)
rows_per_chunk <- 28

selected_rows <- data.frame()
for (i in seq(1, nrow(sorted_data), by = rows_per_chunk)) {
  chunk <- sorted_data[i:min(i + rows_per_chunk - 1, nrow(sorted_data)), ]
  if (nrow(chunk) > 0) {
    best_row <- chunk[which.max(chunk$loglik), ]
    selected_rows <- rbind(selected_rows, best_row)
  }
}
plot(x = log(selected_rows$xi_J), y = selected_rows$loglik)
plot(x = selected_rows$xi_J, y = selected_rows$loglik)
# selected_rows = selected_rows[selected_rows$xi_J < 1e-4,]
selected_rows$log_xi_J <- log(selected_rows$xi_J)
log_xi_J_range <- range(selected_rows$log_xi_J, na.rm = TRUE)
log_xi_J_breaks <- seq(from = log_xi_J_range[1], to = log_xi_J_range[2], length.out = 51)
max_loglik_points <- data.frame()

for (i in 1:50) {
  chunk <- selected_rows %>%
    filter(log_xi_J >= log_xi_J_breaks[i] & log_xi_J < log_xi_J_breaks[i + 1])
  if (nrow(chunk) > 0) {
    best_row <- chunk[which.max(chunk$loglik), ]
    max_loglik_points <- rbind(max_loglik_points, best_row)
  }
}
mcap(max_loglik_points$loglik, max_loglik_points$xi_J,  level = 0.95, span = 0.5, Ngrid = 1000) -> mcap_object_xi_J
mcap_object_xi_J$mle -> xi_J_mle
xi_J_p <- ggplot() +
  geom_point(data = max_loglik_points, aes(x = log_xi_J, y = loglik)) +
  geom_line(data = mcap_object_xi_J$fit, aes(x = log(parameter), y = smoothed), col = 'red') +
  geom_vline(xintercept = log(mcap_object_xi_J$ci[1]), linetype = 'dashed') +
  geom_vline(xintercept = log(mcap_object_xi_J$ci[2]), linetype = 'dashed') +
  geom_vline(xintercept = log(mcap_object_xi_J$mle), col = 'blue') +
  labs(x = "log(xi_J)", y = "log likelihood") +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 10)) +
  ylim(-375, -370)+
  annotate("text", x = log(mcap_object_xi_J$mle), y = -900, label = sprintf("xi_J_mle: %s", formatC(mcap_object_xi_J$mle, format = 'e', digits = 3)), hjust = 1.05, vjust = -0.5, size = 3)

xi_J_p


grid.arrange( xi_J_p,lambda_Ji_p,
              nrow = 1, ncol = 2)


