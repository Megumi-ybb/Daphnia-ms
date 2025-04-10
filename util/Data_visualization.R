library(patchwork)
library(magrittr)
library(readxl)
library(pomp)
library(panelPomp)
library(foreach)
library(doParallel)
registerDoParallel(cores=detectCores())
library(ggplot2)
library(gridExtra)
library(scales)
library(dplyr)

Mesocosm_data = read_excel("Mesocosmdata.xls",3)
# Mesocosm_data = read_excel("/home/ybb/D_P/Mesocosmdata.xlsx",3)


name_str = "all_shared"
run_level <- 2

dentNoPara <- Mesocosm_data[91:170, ]
dentNoPara <- subset(dentNoPara, select = c(rep, day, dent.adult,dent.inf,lum.adult,lum.adult.inf,
                                            lum.ephip,lum.ephip.inf,dent.ephip,dent.adult.male,lum.adult.male,lum.male.inf,dent.total,lum.total))
dentNoPara <- dentNoPara[80: 1, ]
#Convert sampling date to natural date. Samples are collected every 5 days after the #first 7 days.
dentNoPara$day = (dentNoPara$day - 1) * 5 + 7
data = list()
trails = c("K","L","M","N","O","P","Q","S")
for (i in 1: length(trails)){
  data[[i]]=subset(dentNoPara, select = c("day", "dent.adult","dent.inf","lum.adult","lum.adult.inf","lum.ephip","lum.ephip.inf","dent.ephip",
                                          "dent.adult.male","lum.adult.male","lum.male.inf","dent.total","lum.total"),
                   dentNoPara$rep == trails[i])
}
bucket_labeller <- as_labeller(function(bucket) gsub("Bucket ", "Bucket-", bucket))
dentNoPara$Bucket = paste0(rep("Bucket ",80),rep(1:8, each = 10))

ggplot(dentNoPara, aes(x = day, y = sqrt(lum.ephip))) +
  geom_line(aes(color = "Lum Egg"), linewidth = 1) +
  geom_line(aes(y = sqrt(lum.ephip.inf), color = "Infected Lum Egg"), linewidth = 1) +
  geom_line(aes(y = sqrt(dent.ephip), color = "Dent Egg"), linewidth = 1) +
  facet_wrap(~Bucket, nrow = 1, labeller = bucket_labeller) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        strip.text = element_text(size = 10),
        strip.background = element_blank(),
        strip.placement = "outside",
        plot.margin = c(.5, .5, .5, .5)) +
  ylab('Sqrt Density per liter') +
  scale_x_continuous(limits = c(0, 52), breaks = c(0, 25, 50)) +
  scale_color_manual(values = c("Lum Egg" = "blue", 
                                "Infected Lum Egg" = "red",
                                "Dent Egg" = "black")) +
  theme_minimal() +
  guides(color = guide_legend(title = "Species"))


ggplot(dentNoPara, aes(x = day, y = sqrt(dent.total))) +
  geom_line(aes(color = "Lum total"), linewidth = 1) +
  geom_line(aes(y = sqrt(lum.total), color = "Dent total"), linewidth = 1) +
  facet_wrap(~Bucket, nrow = 1, labeller = bucket_labeller) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        strip.text = element_text(size = 10),
        strip.background = element_blank(),
        strip.placement = "outside",
        plot.margin = c(.5, .5, .5, .5)) +
  ylab('Sqrt Density per liter') +
  scale_x_continuous(limits = c(0, 52), breaks = c(0, 25, 50)) +
  scale_color_manual(values = c("Lum total" = "blue", 
                                "Dent total" = "black")) +
  theme_minimal() +
  guides(color = guide_legend(title = "Species"))


ggplot(dentNoPara, aes(x = day, y = sqrt(dent.adult))) +
  geom_line(aes(color = "Susceptible female dent"), linewidth = 1) +
  geom_line(aes(y = sqrt(dent.adult.male), color = "Susceptible male dent"), linewidth = 1) +
  geom_line(aes(y = sqrt(lum.adult.male), color = "Susceptible male lum"), linewidth = 1) +
  geom_line(aes(y = sqrt(lum.adult), color = "Susceptible female lum"), linewidth = 1) +
  geom_line(aes(y = sqrt(lum.adult.inf), color = "Infected female lum"), linewidth = 1) +
  geom_line(aes(y = sqrt(lum.male.inf), color = "Infected male lum"), linewidth = 1) +
  facet_wrap(~Bucket, nrow = 1, labeller = bucket_labeller) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        strip.text = element_text(size = 10),
        strip.background = element_blank(),
        strip.placement = "outside",
        plot.margin = c(.5, .5, .5, .5)) +
  ylab('Sqrt Density per liter') +
  scale_x_continuous(limits = c(0, 52), breaks = c(0, 25, 50)) +
  scale_color_manual(values = c("Susceptible female dent" = "blue", 
                                "Susceptible male dent" = "black",
                                "Susceptible male lum" = "green",
                                "Susceptible female lum" = "orange",
                                "Infected female lum" = "red",
                                "Infected male lum" = "purple")) +
  theme_minimal() +
  guides(color = guide_legend(title = "Species"))
