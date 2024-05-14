#####################################


# Behavioural and EEG analyses of the frequency tagging 
# experiment using male and female faces images
# at a threshold and supraliminal contrasts 

# Audrey Mazancieux 2024


#####################################


## Packages ----------------------------------------------
library(tidyverse)
library(readxl)
library(magrittr)
library(reshape2)
library(broom)
library(cowplot)

# plot theme
plot_theme = theme(
  axis.title.x = element_text(size = 14),
  plot.title = element_text(size = 16),
  axis.text.x = element_text(size = 14),
  axis.title.y = element_text(size = 14))

sub_1_contrast = list(1, 4, 6, 7, 9, 10, 11, 12, 13, 16)
sub_2_contrasts = list(3, 14, 15, 17, 18, 19, 20, 24, 25, 26, 27, 29)


## Import and clean behavioral data -------------------------------------------------------------

# get data
files <-
  list.files("./Behaviour/Data", recursive = TRUE) %>% 
  as_data_frame()

files <- files %>% 
  filter(str_detect(value, regex(".csv$"))) %>% 
  filter(str_detect(value, regex("preproc")))

# load  data
resp_data <- data.frame()

for (i in files$value){
  
  pp_data <- read_csv(file.path("./Behaviour/Data", i)) %>% 
    mutate(Pp = i %>% 
             str_extract(regex("\\d+")))
  resp_data %<>% rbind(pp_data) 
}

resp_data %<>% mutate(Pp = as.numeric(Pp)) 

resp_data2 <- resp_data
for (sub in sub_1_contrast){
  resp_data2 <- resp_data2 %>% 
    filter(Pp != sub) 
}

## First-order task ------------------------------------------------------

# performance count
first_order <- resp_data2 %>% 
  mutate(count = 1) %>% 
  dcast(Pp + contrast ~ accuracy, value.var = "count", sum) %>% 
  mutate(performance = Correct/(Correct+Incorrect))

write.csv(first_order, "./Behaviour/Results/first_order.csv")

png(file="./Behaviour/Results/first_order.png", width=6, height=6, units="in", res=300)

first_order %>% 
  group_by(contrast) %>%
  summarise(VD = mean(performance),
            sd = sd(performance),
            se = sd/sqrt(length(unique(first_order$Pp))),
            CI = se * qt(.975, n() - 1)) %>%
  ggplot(aes(x = contrast, y = VD, fill=contrast)) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin = VD - CI, ymax = VD + CI), width = 0, size = 1)+
  ggtitle("Task performance per contrast") +
  scale_fill_brewer(palette="Dark2") +
  theme_bw() +
  plot_theme +
  xlab("Contrast") +
  ylab("Performance")

dev.off()


## Subjective visibility -------------------------------------------------

# PAS mean
count_pas <- resp_data2 %>% 
  mutate(count = 1) %>%
  dcast(Pp + contrast ~ pas_score+accuracy, value.var = "count", sum)


# plot mean PAS per accuracy and contrast

png(file="./Behaviour/Results/pas_acc_contrast.png", width=6, height=6, units="in", res=300)

resp_data2 %>%
  group_by(Pp, accuracy, contrast) %>%
  summarise(pas_score = mean(as.numeric(pas_score))) %>% 
  group_by(contrast, accuracy) %>%
  summarise(VD = mean(pas_score),
            sd = sd(pas_score),
            se = sd/sqrt(length(unique(resp_data2$Pp))),
            CI = se * qt(.975, n() - 1)) %>%
  ggplot(aes(x = contrast, y = VD, fill=accuracy)) +
  geom_bar(stat="identity", position=position_dodge(1)) +
  geom_errorbar(aes(ymin = VD - CI, ymax = VD + CI), position=position_dodge(1), width = 0, size = 1)+
  ggtitle("mean PAS according to accuracy") +
  scale_fill_brewer(palette="Dark2") +
  theme_bw() +
  plot_theme +
  xlab("Contrast") +
  ylab("mean PAS")

dev.off()


## Confidence -----------------------------------------------------------

# raw confidence per accuracy
png(file="./Behaviour/Results/conf_acc_contrast.png", width=6, height=6, units="in", res=300)

resp_data2 %>%
  group_by(contrast, accuracy) %>%
  summarise(VD = mean(conf_score),
            sd = sd(conf_score),
            se = sd/sqrt(length(unique(resp_data2$Pp))),
            CI = se * qt(.975, n() - 1)) %>%
  ggplot(aes(x = contrast, y = VD, fill=accuracy)) +
  geom_bar(stat="identity", position=position_dodge(0.93)) +
  geom_errorbar(aes(ymin = VD - CI, ymax = VD + CI), position=position_dodge(0.93), width = 0, size = 1)+
  ggtitle("mean confidence according to accuracy") +
  scale_fill_brewer(palette="Dark2") +
  theme_bw() +
  plot_theme +
  xlab("Contrast") +
  ylab("mean confidence")

dev.off()


# raw confidence per PAS and per contrast

png(file="./Behaviour/Results/conf_pas_contrast.png", width=6, height=6, units="in", res=300)

resp_data2 %>%
  filter(pas_score != 4) %>% 
  group_by(Pp, pas_score, contrast) %>%
  summarise(conf_score = mean(as.numeric(conf_score))) %>% 
  group_by(pas_score, contrast) %>%
  summarise(VD = mean(conf_score),
            sd = sd(conf_score),
            se = sd/sqrt(length(unique(resp_data$Pp))),
            CI = se * qt(.975, n() - 1)) %>%
  ggplot(aes(x = pas_score, y = VD, fill = contrast)) +
  geom_bar(stat="identity", position=position_dodge(0.93)) +
  geom_errorbar(aes(ymin = VD - CI, ymax = VD + CI), position=position_dodge(0.93), width = 0, size = 1)+
  ggtitle("mean confidence according to contrast and PAS") +
  scale_fill_brewer(palette="Dark2") +
  theme_bw() +
  plot_theme +
  xlab("PAS") +
  ylab("mean confidence")

dev.off()
  

## Metacognitive efficiency 

# source("Function_metad_indiv.R")
# 
# # prepare data
# nR_S1 <- resp_data_short %>% 
#   filter(correct_response == 102) %>% 
#   mutate(Count = 1,
#          Resp = case_when(
#            decision_resp == 102 ~ "R1",
#            decision_resp == 104 ~ "R2")) %>% 
#   dcast(Pp ~ conf_score + Resp, value.var = "Count", sum) %>% 
#   select(`4_R1`,
#          `3_R1`,
#          `2_R1`,
#          `1_R1`,
#          `1_R2`)
# #         `2_R2`,
# #         `3_R2`,
# #         `4_R2`)
# 
# nR_S2 <- resp_data_short %>% 
#   filter(correct_response == 104) %>% 
#   mutate(Count = 1,
#          Resp = case_when(
#            decision_resp == 102 ~ "R1",
#            decision_resp == 104 ~ "R2")) %>% 
#   dcast(Pp ~ conf_score + Resp, value.var = "Count", sum) %>% 
# #  select(`4_R1`,
#   select(`3_R1`,
# #         `2_R1`,
#          `1_R1`,
#          `1_R2`,
#          `2_R2`,
#          `3_R2`,
#          `4_R2`)
# 
# nsubj <- nrow(nR_S1)
# 
# # calculate d' and meta-d' per participant
# c <- data.frame()
# d <- data.frame()
# stats <- data.frame()
# 
# for (n in 1:(nsubj)) {
#   
#   S1 <- c(t(nR_S1[n,]))
#   S2 <- c(t(nR_S2[n,]))
#   
#   output <- metad_indiv(nR_S1 = S1, nR_S2 = S2)
#   
#   d %<>% rbind(d1) 
#   c %<>% rbind(c1) 
#   
#   Value <- summary(output)
#   
#   stat1 <- data.frame(mean = Value[["statistics"]][, "Mean"])
#   stat %<>% rbind(t(stat1)) 
#   
# }
# 
# write.csv(stat, "./Behaviour/Results/output_fit_individual_metad.csv")



## EEG RESS components and behavior ------------------------------------------------------

# get data
files <-
  list.files("./EEG_analyses/Data", recursive = TRUE) %>% 
  as_data_frame()

files <- files %>% 
  filter(str_detect(value, regex(".csv$"))) %>% 
  filter(str_detect(value, regex("ress")))

# load  data
ress_data <- data.frame()

for (i in files$value){
  pp_data <- read_csv(file.path("./EEG_analyses/Data", i))
  ress_data %<>% rbind(pp_data) 
}


# plot RESS values for participants with 1 contrast

ress_data1 <- ress_data
for (sub in sub_2_contrasts){
  ress_data1 <- ress_data1 %>% 
    filter(subject != sub) 
}

n <- length(unique(ress_data1$subject))

ress_data1 %>%
  group_by(subject, pas_score, contrast) %>%
  summarise(RESS_1_2 = mean(as.numeric(RESS_1_2))) %>% 
  group_by(contrast, pas_score) %>%
  summarise(VD = mean(RESS_1_2),
            sd = sd(RESS_1_2),
            se = sd/sqrt(n),
            CI = se * qt(.975, n() - 1)) %>%
  ggplot(aes(x = pas_score, y = VD, fill=contrast)) +
  geom_bar(stat="identity", position=position_dodge(0.93)) +
  geom_errorbar(aes(ymin = VD - CI, ymax = VD + CI), position=position_dodge(0.93), width = 0, size = 1)+
  ggtitle("RESS for 1.2 Hz according to contrast and PAS") +
  scale_fill_brewer(palette="Dark2") +
  theme_bw() +
  plot_theme +
  xlab("Contrast") +
  ylab("mean RESS for 1.2 Hz")

ress_data1 %>%
  group_by(subject, pas_score, contrast) %>%
  summarise(RESS_2_4 = mean(as.numeric(RESS_2_4))) %>% 
  group_by(contrast, pas_score) %>%
  summarise(VD = mean(RESS_2_4),
            sd = sd(RESS_2_4),
            se = sd/sqrt(n),
            CI = se * qt(.975, n() - 1)) %>%
  ggplot(aes(x = pas_score, y = VD, fill=contrast)) +
  geom_bar(stat="identity", color='black', position=position_dodge(1)) +
  geom_errorbar(aes(ymin = VD - CI, ymax = VD + CI), position=position_dodge(1), width = 0, size = 1)+
  ggtitle("RESS for 2.4 Hz according to contrast and PAS") +
  scale_fill_brewer(palette="Dark2") +
  theme_bw() +
  plot_theme +
  xlab("Contrast") +
  ylab("mean RESS for 2.4 Hz")

ress_data1 %>%
  group_by(subject, pas_score, contrast) %>%
  summarise(RESS_3_6 = mean(as.numeric(RESS_3_6))) %>% 
  group_by(contrast, pas_score) %>%
  summarise(VD = mean(RESS_3_6),
            sd = sd(RESS_3_6),
            se = sd/sqrt(n),
            CI = se * qt(.975, n() - 1)) %>%
  ggplot(aes(x = pas_score, y = VD, fill=contrast)) +
  geom_bar(stat="identity", color='black', position=position_dodge(1)) +
  geom_errorbar(aes(ymin = VD - CI, ymax = VD + CI), position=position_dodge(1), width = 0, size = 1)+
  ggtitle("RESS for 3.6 Hz according to contrast and PAS") +
  scale_fill_brewer(palette="Dark2") +
  theme_bw() +
  plot_theme +
  xlab("Contrast") +
  ylab("mean RESS for 3.6 Hz")


# plot RESS values 

ress_data2 <- ress_data
for (sub in sub_1_contrast){
  ress_data2 <- ress_data2 %>% 
    filter(subject != sub) 
}

n <- length(unique(ress_data2$subject))

png(file="./EEG_analyses/Results/ress_1_2_pas.png", width=6, height=6, units="in", res=300)

ress_data2 %>%
  group_by(subject, pas_score, contrast) %>%
  summarise(RESS_1_2 = mean(as.numeric(RESS_1_2))) %>% 
  group_by(contrast, pas_score) %>%
  summarise(VD = mean(RESS_1_2),
            sd = sd(RESS_1_2),
            se = sd/sqrt(n),
            CI = se * qt(.975, n() - 1)) %>%
  ggplot(aes(x = pas_score, y = VD, fill=contrast)) +
  geom_bar(stat="identity", position=position_dodge(0.93)) +
  geom_errorbar(aes(ymin = VD - CI, ymax = VD + CI), position=position_dodge(0.93), width = 0, size = 1)+
  ggtitle("RESS for 1.2 Hz according to contrast and PAS") +
  scale_fill_brewer(palette="Dark2") +
  theme_bw() +
  plot_theme +
  xlab("Contrast") +
  ylab("mean RESS for 1.2 Hz")

dev.off()

png(file="./EEG_analyses/Results/ress_2_4_pas.png", width=6, height=6, units="in", res=300)

ress_data2 %>%
  group_by(subject, pas_score, contrast) %>%
  summarise(RESS_2_4 = mean(as.numeric(RESS_2_4))) %>% 
  group_by(contrast, pas_score) %>%
  summarise(VD = mean(RESS_2_4),
            sd = sd(RESS_2_4),
            se = sd/sqrt(n),
            CI = se * qt(.975, n() - 1)) %>%
  ggplot(aes(x = pas_score, y = VD, fill=contrast)) +
  geom_bar(stat="identity", position=position_dodge(0.93)) +
  geom_errorbar(aes(ymin = VD - CI, ymax = VD + CI), position=position_dodge(0.93), width = 0, size = 1)+
  ggtitle("RESS for 2.4 Hz according to contrast and PAS") +
  scale_fill_brewer(palette="Dark2") +
  theme_bw() +
  plot_theme +
  xlab("Contrast") +
  ylab("mean RESS for 2.4 Hz")

dev.off()

png(file="./EEG_analyses/Results/ress_3_6_pas.png", width=6, height=6, units="in", res=300)

ress_data2 %>%
  group_by(subject, pas_score, contrast) %>%
  summarise(RESS_3_6 = mean(as.numeric(RESS_3_6))) %>% 
  group_by(contrast, pas_score) %>%
  summarise(VD = mean(RESS_3_6),
            sd = sd(RESS_3_6),
            se = sd/sqrt(n),
            CI = se * qt(.975, n() - 1)) %>%
  ggplot(aes(x = pas_score, y = VD, fill=contrast)) +
  geom_bar(stat="identity", position=position_dodge(0.93)) +
  geom_errorbar(aes(ymin = VD - CI, ymax = VD + CI), position=position_dodge(0.93), width = 0, size = 1)+
  ggtitle("RESS for 3.6 Hz according to contrast and PAS") +
  scale_fill_brewer(palette="Dark2") +
  theme_bw() +
  plot_theme +
  xlab("Contrast") +
  ylab("mean RESS for 3.6 Hz")

dev.off()


# compute linear regression for each participants between RESS and PAS

coeff_all <- data.frame()
for (contr in unique(ress_data2$contrast)){
  
  for (sub in unique(ress_data2$subject)){
    
    data_sub <- ress_data2 %>% 
      filter(subject == sub & contrast == contr)
    
    cor <- cor.test(data_sub$RESS_2_4, data_sub$pas_score)
    r <- cor[['estimate']]
    
    model <- lm(RESS_2_4 ~ pas_score, data=data_sub)
    beta <- data.frame(beta = model[['coefficients']][2])
    
    sub_coef <- data.frame(sub = sub,
                           contrast = contr,
                           r = r,
                           beta = beta)
    
    coeff_all %<>% rbind(sub_coef)
  }
}


# plot r coefficients
mean(coeff_all$r[coeff_all$contrast == '1%'], na.rm = TRUE)
sd(coeff_all$r[coeff_all$contrast == '1%'], na.rm = TRUE)

mean(coeff_all$r[coeff_all$contrast == '1.5%'], na.rm = TRUE)
sd(coeff_all$r[coeff_all$contrast == '1.5%'], na.rm = TRUE)

coeff_all %>%
  group_by(contrast) %>%
  summarise(VD = mean(r, na.rm = TRUE),
            sd = sd(r, na.rm = TRUE),
            se = sd/sqrt(n),
            CI = se * qt(.975, n() - 1)) %>%
  ggplot(aes(x = contrast, y = VD, fill=contrast)) +
  geom_point(data = coeff_all,
             aes(x = contrast, y = r, color = contrast),
             position = position_jitterdodge(0.5),
             size = 1, alpha = 0.8,
             show.legend = FALSE) +
  geom_boxplot(data = coeff_all,
               aes(x = contrast, y = r, fill = contrast),
               outlier.shape = NA,
               alpha = 1, width = .2,
               position = position_dodge(0.1),
               show.legend = FALSE) +
  geom_point(size = 1, color = 'black', position = position_dodge(0.5)) +
  geom_errorbar(aes(ymin = VD - CI, ymax = VD + CI), position=position_dodge(0.5), width = 0, size = 1)+
  ggtitle("Correlation between PAS and RESS values at 2.4 Hz for each contrast") +
  scale_fill_brewer(palette="Dark2") +
  theme_bw() +
  plot_theme +
  xlab("Contrast") +
  ylab("Correlation coefficient")


# plot beta coefficients
coeff_all %>%
  group_by(contrast) %>%
  summarise(VD = mean(beta, na.rm = TRUE),
            sd = sd(beta, na.rm = TRUE),
            se = sd/sqrt(n),
            CI = se * qt(.975, n() - 1)) %>%
  ggplot(aes(x = contrast, y = VD, fill=contrast)) +
  geom_bar(stat="identity", color='black', position=position_dodge(1)) +
  geom_errorbar(aes(ymin = VD - CI, ymax = VD + CI), position=position_dodge(1), width = 0, size = 1)+
  ggtitle("Correlation between PAS and RESS values at 2.4 Hz for each contrast") +
  scale_fill_brewer(palette="Dark2") +
  theme_bw() +
  plot_theme +
  xlab("Contrast") +
  ylab("Bate value")


## Plot RESS with confidence

ress_data2 %>%
  group_by(subject, conf_score, contrast) %>%
  summarise(RESS_2_4 = mean(as.numeric(RESS_2_4))) %>% 
  group_by(contrast, conf_score) %>%
  summarise(VD = mean(RESS_2_4),
            sd = sd(RESS_2_4),
            se = sd/sqrt(n),
            CI = se * qt(.975, n() - 1)) %>%
  ggplot(aes(x = conf_score, y = VD, fill=contrast)) +
  geom_bar(stat="identity", color='black', position=position_dodge(1)) +
  geom_errorbar(aes(ymin = VD - CI, ymax = VD + CI), position=position_dodge(1), width = 0, size = 1)+
  ggtitle("RESS for 2.4 Hz according to contrast and confidence") +
  scale_fill_brewer(palette="Dark2") +
  theme_bw() +
  plot_theme +
  xlab("Contrast") +
  ylab("mean RESS for 2.4 Hz")

coeff_conf_all <- data.frame()
for (contr in unique(ress_data2$contrast)){
  
  for (sub in unique(ress_data2$subject)){
    
    data_sub <- ress_data2 %>% 
      filter(subject == sub & contrast == contr)
    
    cor <- cor.test(data_sub$RESS_2_4, data_sub$conf_score)
    r <- cor[['estimate']]
    
    model <- lm(RESS_2_4 ~ conf_score, data=data_sub)
    beta <- data.frame(beta = model[['coefficients']][2])
    
    sub_coef <- data.frame(sub = sub,
                           contrast = contr,
                           r = r,
                           beta = beta)
    
    coeff_conf_all %<>% rbind(sub_coef)
  }
}

# plot r coefficients
coeff_conf_all %>%
  group_by(contrast) %>%
  summarise(VD = mean(r, na.rm = TRUE),
            sd = sd(r, na.rm = TRUE),
            se = sd/sqrt(n),
            CI = se * qt(.975, n() - 1)) %>%
  ggplot(aes(x = contrast, y = VD, fill=contrast)) +
  geom_bar(stat="identity", color='black', position=position_dodge(1)) +
  geom_errorbar(aes(ymin = VD - CI, ymax = VD + CI), position=position_dodge(1), width = 0, size = 1)+
  ggtitle("Correlation between PAS and RESS values at 2.4 Hz for each contrast") +
  scale_fill_brewer(palette="Dark2") +
  theme_bw() +
  plot_theme +
  xlab("Contrast") +
  ylab("Correlation coefficient")


coeff_conf_all %>%
  group_by(contrast) %>%
  summarise(VD = mean(beta, na.rm = TRUE),
            sd = sd(beta, na.rm = TRUE),
            se = sd/sqrt(n),
            CI = se * qt(.975, n() - 1)) %>%
  ggplot(aes(x = contrast, y = VD, fill=contrast)) +
  geom_bar(stat="identity", color='black', position=position_dodge(1)) +
  geom_errorbar(aes(ymin = VD - CI, ymax = VD + CI), position=position_dodge(1), width = 0, size = 1)+
  ggtitle("Correlation between PAS and RESS values at 2.4 Hz for each contrast") +
  scale_fill_brewer(palette="Dark2") +
  theme_bw() +
  plot_theme +
  xlab("Contrast") +
  ylab("Bate value")


## Plot RESS with accuracy

ress_data2 %>%
  group_by(subject, accuracy, contrast) %>%
  summarise(RESS_2_4 = mean(as.numeric(RESS_2_4))) %>% 
  group_by(contrast, accuracy) %>%
  summarise(VD = mean(RESS_2_4),
            sd = sd(RESS_2_4),
            se = sd/sqrt(n),
            CI = se * qt(.975, n() - 1)) %>%
  ggplot(aes(x = accuracy, y = VD, fill=contrast)) +
  geom_bar(stat="identity", color='black', position=position_dodge(1)) +
  geom_errorbar(aes(ymin = VD - CI, ymax = VD + CI), position=position_dodge(1), width = 0, size = 1)+
  ggtitle("RESS for 2.4 Hz according to contrast and accuracy") +
  scale_fill_brewer(palette="Dark2") +
  theme_bw() +
  plot_theme +
  xlab("Contrast") +
  ylab("mean RESS for 2.4 Hz")

ress_data2 %>%
  group_by(subject, accuracy, contrast) %>%
  summarise(RESS_1_2 = mean(as.numeric(RESS_1_2))) %>% 
  group_by(contrast, accuracy) %>%
  summarise(VD = mean(RESS_1_2),
            sd = sd(RESS_1_2),
            se = sd/sqrt(n),
            CI = se * qt(.975, n() - 1)) %>%
  ggplot(aes(x = accuracy, y = VD, fill=contrast)) +
  geom_bar(stat="identity", color='black', position=position_dodge(1)) +
  geom_errorbar(aes(ymin = VD - CI, ymax = VD + CI), position=position_dodge(1), width = 0, size = 1)+
  ggtitle("RESS for 1.2 Hz according to contrast and accuracy") +
  scale_fill_brewer(palette="Dark2") +
  theme_bw() +
  plot_theme +
  xlab("Contrast") +
  ylab("mean RESS for 1.2 Hz")


## EEG SNR per ROI and behavior --------------------------------------------------

# get data
files <-
  list.files("./EEG_analyses/", recursive = FALSE) %>% 
  as_data_frame()

files <- files %>% 
  filter(str_detect(value, regex(".csv$"))) %>% 
  filter(str_detect(value, regex("roi")))

# load  data
roi_snr_data <- data.frame()

for (i in files$value){
  pp_data <- read_csv(file.path("./EEG_analyses/", i)) %>% 
    mutate(roi = i %>% 
             str_extract(regex("\\d+")))
  roi_snr_data %<>% rbind(pp_data) 
}


# plot data at 1.2 Hz
png(file="./EEG_analyses/Results/snr_1_2_pas.png", width=6, height=6, units="in", res=300)

roi_snr_data %>%
  group_by(roi, PAS, Contrast) %>%
  summarise(VD = mean(`1_2Hz`),
            sd = sd(`1_2Hz`),
            se = sd/sqrt(n),
            CI = se * qt(.975, n() - 1)) %>%
  ggplot(aes(x = PAS, y = VD, fill=Contrast)) +
  geom_bar(stat="identity", position=position_dodge(0.93)) +
  geom_errorbar(aes(ymin = VD - CI, ymax = VD + CI), position=position_dodge(0.93), width = 0, size = 0.7)+
  facet_wrap(~ roi) +
  ggtitle("SNR at 1.2 Hz according to contrast, PAS, and ROI") +
  scale_fill_brewer(palette="Dark2") +
  theme_bw() +
  plot_theme +
  xlab("Contrast") +
  ylab("SNR at 1.2 Hz")

dev.off()


# plot data at 2.4 Hz
png(file="./EEG_analyses/Results/snr_2_4_pas.png", width=6, height=6, units="in", res=300)

roi_snr_data %>%
  group_by(roi, PAS, Contrast) %>%
  summarise(VD = mean(`2_4Hz`),
            sd = sd(`1_2Hz`),
            se = sd/sqrt(n),
            CI = se * qt(.975, n() - 1)) %>%
  ggplot(aes(x = PAS, y = VD, fill=Contrast)) +
  geom_bar(stat="identity", position=position_dodge(0.93)) +
  geom_errorbar(aes(ymin = VD - CI, ymax = VD + CI), position=position_dodge(0.93), width = 0, size = 0.7)+
  facet_wrap(~ roi) +
  ggtitle("SNR at 2.4 Hz according to contrast, PAS, and ROI") +
  scale_fill_brewer(palette="Dark2") +
  theme_bw() +
  plot_theme +
  xlab("Contrast") +
  ylab("SNR at 2.4 Hz")

dev.off()
