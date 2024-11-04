#####################################


# Behavioural and EEG analyses of the frequency 
# tagging experiment using male and female faces 
# images at threshold and supraliminal contrasts 

# Audrey Mazancieux 2024


#####################################


## Packages ----------------------------------------------
library(tidyverse)
library(readxl)
library(magrittr)
library(reshape2)
library(broom)
library(cowplot)
library(lmerTest)
library(lattice)

# plot theme
plot_theme = theme(
  axis.title.x = element_text(size = 14),
  plot.title = element_text(size = 14),
  axis.text.x = element_text(size = 14),
  axis.title.y = element_text(size = 14))


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

# gather condition 1.5% and 1.50%
resp_data2 <- resp_data2 %>% 
  mutate(contrast = ifelse(contrast == '1.5%' | contrast == '1.50%', '1.5%', '1%'))

# save clean dataset
write.csv(resp_data2, "./Behaviour/Results/behaviour_clean_data.csv")

# get number of trial per subject and PAS response
trials <- resp_data2 %>% 
  mutate(count = 1) %>% 
  dcast(Pp + contrast ~ pas_score, value.var = 'count', sum)


## First-order performance ------------------------------------------------------

# performance count
first_order <- resp_data2 %>% 
  mutate(count = 1) %>% 
  dcast(Pp + contrast ~ accuracy, value.var = "count", sum) %>% 
  mutate(performance = Correct/(Correct+Incorrect))

n <- length(unique(first_order$Pp))


png(file="./Behaviour/Results/first_order.png", width=6, height=6, units="in", res=300)

first_order %>% 
  group_by(contrast) %>%
  summarise(VD = mean(performance),
            sd = sd(performance),
            se = sd/sqrt(n),
            CI = se * qt(.975, n() - 1)) %>%
  ggplot(aes(x = contrast, y = VD, fill=contrast)) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin = VD - CI, ymax = VD + CI), width = 0, size = 1)+
  ggtitle("Task performance per contrast") +
  scale_fill_manual(values = c("#0F056B", "#9683EC")) +
  theme_classic() +
  plot_theme +
  xlab("Contrast") +
  ylab("Performance")

dev.off()

# analyses
mean(first_order$performance[first_order$contrast == '1%'])
mean(first_order$performance[first_order$contrast == '1.5%'])
sd(first_order$performance[first_order$contrast == '1%'])
sd(first_order$performance[first_order$contrast == '1.5%'])

first_order2 <- first_order %>% 
  dcast(Pp ~ contrast, value.var = 'performance') %>% 
  mutate(perf_diff = `1.5%` - `1%`,
         perf_1_to_0 = `1%` - 0.50,
         perf_1.5_to_0 = `1.5%` - 0.50)
mod <- lm(perf_diff ~ 1, data = first_order2)
mod1 <- lm(perf_1_to_0 ~ 1, data = first_order2)
mod2 <- lm(perf_1.5_to_0 ~ 1, data = first_order2)

qqnorm(residuals(mod))
qqline(residuals(mod))

data.frame(x = residuals(mod)) %>%
  ggplot(aes(x = x)) +
  geom_histogram()
shapiro.test(residuals(mod))

summary(mod)
summary(mod1)
summary(mod2)


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
            se = sd/sqrt(n),
            CI = se * qt(.975, n() - 1)) %>%
  ggplot(aes(x = contrast, y = VD, fill=accuracy)) +
  geom_bar(stat="identity", position=position_dodge(1)) +
  geom_errorbar(aes(ymin = VD - CI, ymax = VD + CI), position=position_dodge(1), width = 0, size = 1)+
  ggtitle("mean PAS according to accuracy") +
  scale_fill_manual(values = c("#003366", "#2C75FF")) +
  theme_classic() +
  plot_theme +
  xlab("Contrast") +
  ylab("mean PAS")

dev.off()

# analyses
mean(resp_data2$pas_score[resp_data2$contrast == '1%' & resp_data2$accuracy == 'Correct'])
mean(resp_data2$pas_score[resp_data2$contrast == '1.5%' & resp_data2$accuracy == 'Correct'])
mean(resp_data2$pas_score[resp_data2$contrast == '1%' & resp_data2$accuracy == 'Incorrect'])
mean(resp_data2$pas_score[resp_data2$contrast == '1.5%' & resp_data2$accuracy == 'Incorrect'])

sd(resp_data2$pas_score[resp_data2$contrast == '1%' & resp_data2$accuracy == 'Correct'])
sd(resp_data2$pas_score[resp_data2$contrast == '1.5%' & resp_data2$accuracy == 'Correct'])
sd(resp_data2$pas_score[resp_data2$contrast == '1%' & resp_data2$accuracy == 'Incorrect'])
sd(resp_data2$pas_score[resp_data2$contrast == '1.5%' & resp_data2$accuracy == 'Incorrect'])

pas_analyses <- resp_data2 %>% 
  dcast(Pp ~ accuracy + contrast, value.var = 'pas_score', mean) %>% 
  mutate(accuracy_main = ((`Correct_1%` + `Correct_1.5%`)/2) - ((`Incorrect_1%` + `Incorrect_1.5%`)/2),
         contrast_main = ((`Incorrect_1.5%` + `Correct_1.5%`)/2) - ((`Incorrect_1%` + `Correct_1%`)/2),
         interraction = ((`Correct_1%` + `Incorrect_1.5%`)/2) - ((`Incorrect_1%` + `Correct_1.5%`)/2),
         acc_diff_1 = `Correct_1%` - `Incorrect_1%`,
         acc_diff_1_5 = `Correct_1.5%` - `Incorrect_1.5%`)

mod_main_acc <- lm(accuracy_main ~ 1, data = pas_analyses)
mod_main_cont <- lm(contrast_main ~ 1, data = pas_analyses)
mod_int <- lm(interraction ~ 1, data = pas_analyses)
mod_1 <- lm(acc_diff_1 ~ 1, data = pas_analyses)
mod_1_5 <- lm(acc_diff_1_5 ~ 1, data = pas_analyses)

qqnorm(residuals(mod_main_acc))
qqline(residuals(mod_main_acc))

data.frame(x = residuals(mod_main_acc)) %>%
  ggplot(aes(x = x)) +
  geom_histogram()
shapiro.test(residuals(mod_main_acc))

summary(mod_main_acc)
summary(mod_main_cont)
summary(mod_int)
summary(mod_1)
summary(mod_1_5)


## Confidence -----------------------------------------------------------

# raw confidence per accuracy
png(file="./Behaviour/Results/conf_acc_contrast.png", width=6, height=6, units="in", res=300)

resp_data2 %>%
  group_by(contrast, accuracy) %>%
  summarise(VD = mean(conf_score),
            sd = sd(conf_score),
            se = sd/sqrt(n),
            CI = se * qt(.975, n() - 1)) %>%
  ggplot(aes(x = contrast, y = VD, fill=accuracy)) +
  geom_bar(stat="identity", position=position_dodge(0.93)) +
  geom_errorbar(aes(ymin = VD - CI, ymax = VD + CI), position=position_dodge(0.93), width = 0, size = 1)+
  ggtitle("mean confidence according to accuracy") +
  scale_fill_manual(values = c("#003366", "#2C75FF")) +
  theme_classic() +
  plot_theme +
  xlab("Contrast") +
  ylab("mean confidence")

dev.off()

# analyses
mean(resp_data2$conf_score[resp_data2$contrast == '1%' & resp_data2$accuracy == 'Correct'])
mean(resp_data2$conf_score[resp_data2$contrast == '1.5%' & resp_data2$accuracy == 'Correct'])
mean(resp_data2$conf_score[resp_data2$contrast == '1%' & resp_data2$accuracy == 'Incorrect'])
mean(resp_data2$conf_score[resp_data2$contrast == '1.5%' & resp_data2$accuracy == 'Incorrect'])

sd(resp_data2$conf_score[resp_data2$contrast == '1%' & resp_data2$accuracy == 'Correct'])
sd(resp_data2$conf_score[resp_data2$contrast == '1.5%' & resp_data2$accuracy == 'Correct'])
sd(resp_data2$conf_score[resp_data2$contrast == '1%' & resp_data2$accuracy == 'Incorrect'])
sd(resp_data2$conf_score[resp_data2$contrast == '1.5%' & resp_data2$accuracy == 'Incorrect'])

conf_analyses <- resp_data2 %>% 
  dcast(Pp ~ accuracy + contrast, value.var = 'conf_score', mean) %>% 
  mutate(accuracy_main = ((`Correct_1%` + `Correct_1.5%`)/2) - ((`Incorrect_1%` + `Incorrect_1.5%`)/2),
         contrast_main = ((`Incorrect_1.5%` + `Correct_1.5%`)/2) - ((`Incorrect_1%` + `Correct_1%`)/2),
         interraction = ((`Correct_1%` + `Incorrect_1.5%`)/2) - ((`Incorrect_1%` + `Correct_1.5%`)/2),
         acc_diff_1 = `Correct_1%` - `Incorrect_1%`,
         acc_diff_1_5 = `Correct_1.5%` - `Incorrect_1.5%`)

mod_main_acc <- lm(accuracy_main ~ 1, data = conf_analyses)
mod_main_cont <- lm(contrast_main ~ 1, data = conf_analyses)
mod_int <- lm(interraction ~ 1, data = conf_analyses)
mod_1 <- lm(acc_diff_1 ~ 1, data = conf_analyses)
mod_1_5 <- lm(acc_diff_1_5 ~ 1, data = conf_analyses)

qqnorm(residuals(mod_main_acc))
qqline(residuals(mod_main_acc))

data.frame(x = residuals(mod_main_acc)) %>%
  ggplot(aes(x = x)) +
  geom_histogram()
shapiro.test(residuals(mod_main_acc))

summary(mod_main_acc)
summary(mod_main_cont)
summary(mod_int)
summary(mod_1)
summary(mod_1_5)


# raw confidence per PAS and per contrast

png(file="./Behaviour/Results/conf_pas_contrast.png", width=6, height=6, units="in", res=300)

resp_data2 %>%
  group_by(Pp, pas_score, contrast) %>%
  summarise(conf_score = mean(as.numeric(conf_score))) %>% 
  group_by(pas_score, contrast) %>%
  summarise(VD = mean(conf_score),
            sd = sd(conf_score),
            se = sd/sqrt(n),
            CI = se * qt(.975, n() - 1)) %>%
  ggplot(aes(x = pas_score, y = VD, fill = contrast)) +
  geom_bar(stat="identity", position=position_dodge(0.93)) +
  geom_errorbar(aes(ymin = VD - CI, ymax = VD + CI), position=position_dodge(0.93), width = 0, size = 1)+
  ggtitle("mean confidence according to contrast and PAS") +
  scale_fill_manual(values = c("#0F056B", "#9683EC")) +
  theme_classic() +
  plot_theme +
  xlab("PAS") +
  ylab("mean confidence")

dev.off()
  

## Metacognitive efficiency 

source("Function_metad_indiv.R")

nR_S1 <- list(
  data.frame(),
  data.frame())

nR_S2 <- list(
  data.frame(),
  data.frame())

n_contrast <- length(unique(resp_data2$contrast))

# prepare data
for (contr in 1:(n_contrast)) {
  
  nR_S1_1 <- resp_data2 %>%
    mutate(contrast_num = ifelse(contrast == '1%', 1, 2)) %>% 
    filter(correct_response == 102 & contrast_num == contr) %>%
    mutate(Count = 1,
           Resp = case_when(
             decision_resp == 102 ~ "R1",
             decision_resp == 104 ~ "R2")) %>%
    dcast(Pp ~ conf_score + Resp, value.var = "Count", sum) %>%
    select(`4_R1`,
           `3_R1`,
           `2_R1`,
           `1_R1`,
           `1_R2`,
           `2_R2`,
           `3_R2`,
           `4_R2`)
  
  nR_S1[[contr]] %<>%  rbind(nR_S1_1)
  
  nR_S2_1 <- resp_data2 %>%
    mutate(contrast_num = ifelse(contrast == '1%', 1, 2)) %>% 
    filter(correct_response == 104 & contrast_num == contr) %>%
    mutate(Count = 1,
           Resp = case_when(
             decision_resp == 102 ~ "R1",
             decision_resp == 104 ~ "R2")) %>%
    dcast(Pp ~ conf_score + Resp, value.var = "Count", sum) %>%
    select(`4_R1`,
           `3_R1`,
           `2_R1`,
           `1_R1`,
           `1_R2`,
           `2_R2`,
           `3_R2`,
           `4_R2`)
  
  nR_S2[[contr]] %<>%  rbind(nR_S2_1)
  
}


nsubj <- nrow(nR_S1[[1]])

# calculate d' and meta-d' per participant
metad <- data.frame()

for (n in 1:(nsubj)) {
  
  for (contr in 1:(n_contrast)) {
  
    S1 <- c(t(nR_S1[[contr]][n,]))
    S2 <- c(t(nR_S2[[contr]][n,]))
  
    output <- metad_indiv(nR_S1 = S1, nR_S2 = S2)
    Value <- summary(output)
  
    stat1 <- data.frame(metad = Value[["statistics"]][, "Mean"]["meta_d"]) %>% 
      mutate(Pp = n,
             d = d1,
             c = c1,
             Mratio = metad / d,
             contrast = contr)
    metad %<>% rbind(stat1)
  }
}


write.csv(metad, "./Behaviour/Results/output_fit_individual_metad.csv")

mean(metad$Mratio[metad$contrast == 1])
mean(metad$Mratio[metad$contrast == 2])
sd(metad$Mratio[metad$contrast == 1])
sd(metad$Mratio[metad$contrast == 2])

metad2 <- metad %>% 
  dcast(Pp ~ contrast, value.var = 'Mratio') %>% 
  mutate(diff = `2` - `1`)

mod_diff <- lm(diff ~ 1, data = metad2)
mod_1 <- lm(`1` ~ 1, data = metad2)
mod_1_5 <- lm(`2` ~ 1, data = metad2)

summary(mod_diff)
summary(mod_1)
summary(mod_1_5)

cor.test(metad$Mratio[metad$contrast == 1], metad$d[metad$contrast == 1])

metad %<>%
  mutate(aware = case_when(
    Mratio < 0 ~ 0,
    Mratio > 0 ~ 1
  ))


## EEG SNR per ROI and PAS --------------------------------------------------

# get data
files <-
  list.files("./EEG_analyses/Results", recursive = FALSE) %>% 
  as_data_frame()

files <- files %>% 
  filter(str_detect(value, regex(".csv$"))) %>%
  filter(str_detect(value, regex("PAS")))

# load  data
roi_snr_data <- data.frame()

for (i in files$value){
  pp_data <- read_csv(file.path("./EEG_analyses/Results", i)) %>% 
    mutate(roi = i %>% 
             str_extract(regex("(?<=roi_)[a-zA-Z0-9]{3}")))
  roi_snr_data %<>% rbind(pp_data) 
}

roi_snr_data %<>%
  mutate(Contrast = case_when(
    Contrast == '1%' ~ '1%',
    Contrast == '1.5%' | Contrast == '1.50%' ~ '1.5%'),
    roi = case_when(
      roi == 'OCC' ~ "Occipital",
      roi == 'OT1' ~ "Left OT",
      roi == 'OT2' ~ "Right OT"))

n <- length(unique(roi_snr_data$Subject))

# plot data at 1.2 Hz
png(file="./EEG_analyses/Results/pas/snr_1_2_pas_left.png", width=6, height=4, units="in", res=300)
roi_snr_data %>%
  filter(roi == "Left OT") %>% 
  ggplot(aes(x = PAS, y = `1_2Hz`, color=Contrast)) + 
  geom_point(position = position_jitterdodge(0.3), size = 1, alpha = 0.7) + 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE, size = 1.5) + 
  ggtitle("Left OT") +
  scale_color_manual(values = c("#0F056B", "#9683EC")) +
  theme_classic() +
  plot_theme + 
  ylim(0, 8) + 
  xlab("PAS ratings") +
  ylab("SNR at 1.2 Hz")
dev.off()

png(file="./EEG_analyses/Results/pas/snr_1_2_pas_right.png", width=6, height=4, units="in", res=300)
roi_snr_data %>%
  filter(roi == "Right OT") %>% 
  ggplot(aes(x = PAS, y = `1_2Hz`, color=Contrast)) + 
  geom_point(position = position_jitterdodge(0.3), size = 1, alpha = 0.7) + 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE, size = 1.5) + 
  ggtitle("Right OT") +
  scale_color_manual(values = c("#0F056B", "#9683EC")) +
  theme_classic() +
  plot_theme + 
  ylim(0, 8) + 
  xlab("PAS ratings") +
  ylab("SNR at 1.2 Hz")
dev.off()

# plot data at 6 Hz
png(file="./EEG_analyses/Results/pas/snr_6_pas.png", width=6, height=4, units="in", res=300)
roi_snr_data %>%
  filter(roi == "Occipital") %>% 
  ggplot(aes(x = PAS, y = `6Hz`, color=Contrast)) + 
  geom_point(position = position_jitterdodge(0.3), size = 1, alpha = 0.7) + 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE, size = 1.5) + 
  ggtitle("Occipital") +
  scale_color_manual(values = c("#0F056B", "#9683EC")) +
  theme_classic() +
  plot_theme + 
  ylim(0, 80) + 
  xlab("PAS ratings") +
  ylab("SNR at 6 Hz")
dev.off()


## Mixed-models for face signal

roi_snr_data_model <- roi_snr_data %>%
  filter(roi != "Occipital") %>% 
  mutate(contrastC = ifelse(Contrast == '1%', -0.5, 0.5),
         roiC = ifelse(roi == 'Left OT', -0.5, 0.5),
         C_1 = ifelse(Contrast == '1%', 0, 1),
         C_1_5 = ifelse(Contrast == '1.5%', 0, 1))

# models for 1.2 Hz 
m_face <- lmer(`1_2Hz` ~ contrastC * PAS * roiC + (contrastC * PAS|Subject), data = roi_snr_data_model)
# use PCA to choose random effects 
summary(rePCA(m_face))

m_face <- lmer(`1_2Hz` ~ contrastC * PAS * roiC + (1|Subject), data = roi_snr_data_model)
m_face_1 <- lmer(`1_2Hz` ~ C_1 * PAS * roiC + (1|Subject), data = roi_snr_data_model)
m_face_1_5 <- lmer(`1_2Hz` ~ C_1_5 * PAS * roiC + (1|Subject), data = roi_snr_data_model)

qqnorm(residuals(m_face))
qqline(residuals(m_face))
qqmath(ranef(m_face))

data.frame(x = residuals(m_face)) %>% 
  ggplot(aes(x = x)) +
  geom_histogram()

summary(m_face)
summary(m_face_1)
summary(m_face_1_5)


## Mixed-models for image signal

roi_snr_data_model <- roi_snr_data %>%
  filter(roi == "Occipital") %>% 
  mutate(contrastC = ifelse(Contrast == '1%', -0.5, 0.5),
         Contrast_1 = ifelse(Contrast == '1%', 0, 1),
         Contrast_1_5 = ifelse(Contrast == '1.5%', 0, 1))

# models for 6 Hz
m_image <- lmer(`6Hz` ~ contrastC * PAS + (contrastC * PAS|Subject), data = roi_snr_data_model)
# use PCA to choose random effects 
summary(rePCA(m_image))

m_image <- lmer(`6Hz` ~ contrastC * PAS + (1|Subject), data = roi_snr_data_model)
m_image_1 <- lmer(`6Hz` ~ Contrast_1 * PAS + (1|Subject), data = roi_snr_data_model)
m_image_1_5 <- lmer(`6Hz` ~ Contrast_1_5 * PAS + (1|Subject), data = roi_snr_data_model)

qqnorm(residuals(m_image))
qqline(residuals(m_image))
qqmath(ranef(m_image))

data.frame(x = residuals(m_image)) %>% 
  ggplot(aes(x = x)) +
  geom_histogram()

summary(m_image)
summary(m_image_1)
summary(m_image_1_5)


## EEG SNR per ROI and accuracy --------------------------------------------------

# get data
files <-
  list.files("./EEG_analyses/Results", recursive = FALSE) %>% 
  as_data_frame()

files <- files %>% 
  filter(str_detect(value, regex(".csv$"))) %>% 
  filter(str_detect(value, regex("acc")))

# load  data
roi_snr_data <- data.frame()

for (i in files$value){
  pp_data <- read_csv(file.path("./EEG_analyses/Results", i)) %>% 
    mutate(roi = i %>% 
             str_extract(regex("(?<=roi_)[a-zA-Z0-9]{3}")))
  roi_snr_data %<>% rbind(pp_data) 
}

roi_snr_data %<>%
  mutate(Contrast = case_when(
    Contrast == '1%' ~ '1%',
    Contrast == '1.5%' | Contrast == '1.50%' ~ '1.5%'),
    roi = case_when(
      roi == 'OCC' ~ "Occipital",
      roi == 'OT1' ~ "Left OT",
      roi == 'OT2' ~ "Right OT")) 

# plot data at 1.2 Hz
png(file="./EEG_analyses/Results/accuracy/snr_1_2_acc.png", width=6, height=6, units="in", res=300)
roi_snr_data %>%
  filter(roi != "Occipital" ) %>% 
  group_by(roi, Accuracy, Contrast) %>%
  summarise(VD = mean(`1_2Hz`),
            sd = sd(`1_2Hz`),
            se = sd/sqrt(n),
            CI = se * qt(.975, n() - 1)) %>%
  ggplot(aes(x = Accuracy, y = VD, fill=Contrast)) +
  geom_bar(stat="identity", position=position_dodge(0.93)) +
  geom_errorbar(aes(ymin = VD - CI, ymax = VD + CI), position=position_dodge(0.93), width = 0, size = 0.7)+
  facet_wrap(~ roi) +
  ggtitle("SNR at 1.2 Hz according to contrast, accuracy, and ROI") +
  scale_fill_manual(values = c("#0F056B", "#9683EC")) +
  theme_classic() +
  plot_theme +
  xlab("Accuracy") +
  ylab("SNR at 1.2 Hz")
dev.off()

# plot for 6 Hz 
png(file="./EEG_analyses/Results/accuracy/snr_6_acc.png", width=6, height=6, units="in", res=300)
roi_snr_data %>%
  filter(roi == "Occipital" ) %>% 
  group_by(roi, Accuracy, Contrast) %>%
  summarise(VD = mean(`6Hz`),
            sd = sd(`6Hz`),
            se = sd/sqrt(n),
            CI = se * qt(.975, n() - 1)) %>%
  ggplot(aes(x = Accuracy, y = VD, fill=Contrast)) +
  geom_bar(stat="identity", position=position_dodge(0.93)) +
  geom_errorbar(aes(ymin = VD - CI, ymax = VD + CI), position=position_dodge(0.93), width = 0, size = 0.7)+
  ggtitle("SNR at 6 Hz according to contrast and accuracy in the occipital ROI") +
  scale_fill_manual(values = c("#0F056B", "#9683EC")) +
  theme_classic() +
  plot_theme +
  xlab("Accuracy") +
  ylab("SNR at 6 Hz")
dev.off()


## Mixed-models for face signal

# exclude outliers (sub-39) and create contrasts
roi_snr_data_model <- roi_snr_data %>%
  filter(roi != "Occipital") %>% 
  mutate(accC = ifelse(Accuracy == 'Incorrect', -0.5, 0.5),
         contrastC = ifelse(Contrast == '1%', -0.5, 0.5),
         roiC = ifelse(roi == 'Left OT', -0.5, 0.5),
         C_1 = ifelse(Contrast == '1%', 0, 1),
         C_1_5 = ifelse(Contrast == '1.5%', 0, 1))

# models for 1.2 Hz 
m_face <- lmer(`1_2Hz` ~ contrastC * accC * roiC + (contrastC * accC|Subject), data = roi_snr_data_model)
# use PCA to choose random effects 
summary(rePCA(m_face))

m_face <- lmer(`1_2Hz` ~ contrastC * accC * roiC + (1|Subject), data = roi_snr_data_model)
m_face_1 <- lmer(`1_2Hz` ~ C_1 * accC * roiC + (1|Subject), data = roi_snr_data_model)
m_face_1_5 <- lmer(`1_2Hz` ~ C_1_5 * accC * roiC + (1|Subject), data = roi_snr_data_model)

qqnorm(residuals(m_face))
qqline(residuals(m_face))
qqmath(ranef(m_face))

data.frame(x = residuals(m_face)) %>% 
  ggplot(aes(x = x)) +
  geom_histogram()

summary(m_face)
summary(m_face_1)
summary(m_face_1_5)


## Mixed-models for image signal

roi_snr_data_model <- roi_snr_data %>%
  filter(roi == "Occipital") %>% 
  mutate(accC = ifelse(Accuracy == 'Incorrect', -0.5, 0.5),
         contrastC = ifelse(Contrast == '1%', -0.5, 0.5),
         C_1 = ifelse(Contrast == '1%', 0, 1),
         C_1_5 = ifelse(Contrast == '1.5%', 0, 1))
         
# models for 6 Hz 
m_image <- lmer(`6Hz` ~ contrastC * accC +(contrastC|Subject), data = roi_snr_data_model)
# use PCA to choose random effects 
summary(rePCA(m_image))

m_image <- lmer(`6Hz` ~ contrastC * accC +(1|Subject), data = roi_snr_data_model)
m_image_1 <- lmer(`6Hz` ~ C_1 * accC + (1|Subject), data = roi_snr_data_model)
m_image_1_5 <- lmer(`6Hz` ~ C_1_5 * accC + (1|Subject), data = roi_snr_data_model)

qqnorm(residuals(m_image))
qqline(residuals(m_image))
qqmath(ranef(m_image))

data.frame(x = residuals(m_image)) %>% 
  ggplot(aes(x = x)) +
  geom_histogram()

summary(m_image)
summary(m_image_1)
summary(m_image_1_5)


## EEG SNR per ROI and confidence ------------------------------------------

# get data
files <-
  list.files("./EEG_analyses/Results", recursive = FALSE) %>% 
  as_data_frame()

files <- files %>% 
  filter(str_detect(value, regex(".csv$"))) %>%
  filter(str_detect(value, regex("conf")))

# load  data
roi_snr_data <- data.frame()

for (i in files$value){
  pp_data <- read_csv(file.path("./EEG_analyses//Results", i)) %>% 
    mutate(roi = i %>% 
             str_extract(regex("(?<=roi_)[a-zA-Z0-9]{3}")))
  roi_snr_data %<>% rbind(pp_data) 
}

roi_snr_data %<>%
  mutate(Contrast = case_when(
    Contrast == '1%' ~ '1%',
    Contrast == '1.5%' | Contrast == '1.50%' ~ '1.5%'),
    roi = case_when(
      roi == 'OCC' ~ "Occipital",
      roi == 'OT1' ~ "Left OT",
      roi == 'OT2' ~ "Right OT"))

# plot data at 1.2 Hz
png(file="./EEG_analyses/Results/confidence/snr_1_2_pas_left.png", width=6, height=4, units="in", res=300)
roi_snr_data %>%
  filter(roi == "Left OT") %>% 
  ggplot(aes(x = Confidence, y = `1_2Hz`, color=Contrast)) + 
  geom_point(position = position_jitterdodge(0.3), size = 1, alpha = 0.7) + 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE, size = 1.5) + 
  ggtitle("Left OT") +
  scale_color_manual(values = c("#0F056B", "#9683EC")) +
  theme_classic() +
  plot_theme + 
  ylim(0, 10) + 
  xlab("PAS ratings") +
  ylab("SNR at 1.2 Hz")
dev.off()

png(file="./EEG_analyses/Results/confidence/snr_1_2_pas_right.png", width=6, height=4, units="in", res=300)
roi_snr_data %>%
  filter(roi == "Right OT") %>% 
  ggplot(aes(x = Confidence, y = `1_2Hz`, color=Contrast)) + 
  geom_point(position = position_jitterdodge(0.3), size = 1, alpha = 0.7) + 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE, size = 1.5) + 
  ggtitle("Right OT") +
  scale_color_manual(values = c("#0F056B", "#9683EC")) +
  theme_classic() +
  plot_theme + 
  ylim(0, 10) + 
  xlab("PAS ratings") +
  ylab("SNR at 1.2 Hz")
dev.off()

# plot data at 6 Hz
png(file="./EEG_analyses/Results/confidence/snr_6_pas.png", width=6, height=4, units="in", res=300)
roi_snr_data %>%
  filter(roi == "Occipital") %>% 
  ggplot(aes(x = Confidence, y = `6Hz`, color=Contrast)) + 
  geom_point(position = position_jitterdodge(0.3), size = 1, alpha = 0.7) + 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE, size = 1.5) + 
  ggtitle("Occipital") +
  scale_color_manual(values = c("#0F056B", "#9683EC")) +
  theme_classic() +
  plot_theme + 
  # ylim(0, 80) + 
  xlab("PAS ratings") +
  ylab("SNR at 6 Hz")
dev.off()


## Mixed-models for face signal

roi_snr_data_model <- roi_snr_data %>%
  filter(roi != "Occipital") %>% 
  mutate(contrastC = ifelse(Contrast == '1%', -0.5, 0.5),
         roiC = ifelse(roi == 'Left OT', -0.5, 0.5),
         C_1 = ifelse(Contrast == '1%', 0, 1),
         C_1_5 = ifelse(Contrast == '1.5%', 0, 1))

# models for 1.2 Hz 
m_face <- lmer(`1_2Hz` ~ contrastC * Confidence * roiC + (contrastC * Confidence|Subject), data = roi_snr_data_model)
# use PCA to choose random effects 
summary(rePCA(m_face))

m_face <- lmer(`1_2Hz` ~ contrastC * Confidence * roiC + (1|Subject), data = roi_snr_data_model)
m_face_1 <- lmer(`1_2Hz` ~ C_1 * Confidence * roiC + (1|Subject), data = roi_snr_data_model)
m_face_1_5 <- lmer(`1_2Hz` ~ C_1_5 * Confidence * roiC + (1|Subject), data = roi_snr_data_model)

qqnorm(residuals(m_face))
qqline(residuals(m_face))
qqmath(ranef(m_face))

data.frame(x = residuals(m_face)) %>% 
  ggplot(aes(x = x)) +
  geom_histogram()

summary(m_face)
summary(m_face_1)
summary(m_face_1_5)


## Mixed-models for image signal

roi_snr_data_model <- roi_snr_data %>%
  filter(roi == "Occipital") %>% 
  mutate(contrastC = ifelse(Contrast == '1%', -0.5, 0.5),
         Contrast_1 = ifelse(Contrast == '1%', 0, 1),
         Contrast_1_5 = ifelse(Contrast == '1.5%', 0, 1))

# models for 6 Hz
m_image <- lmer(`6Hz` ~ contrastC * Confidence + (contrastC * Confidence|Subject), data = roi_snr_data_model)
# use PCA to choose random effects 
summary(rePCA(m_image))

m_image <- lmer(`6Hz` ~ contrastC * Confidence + (1|Subject), data = roi_snr_data_model)
m_image_1 <- lmer(`6Hz` ~ Contrast_1 * Confidence + (1|Subject), data = roi_snr_data_model)
m_image_1_5 <- lmer(`6Hz` ~ Contrast_1_5 * Confidence + (1|Subject), data = roi_snr_data_model)

qqnorm(residuals(m_image))
qqline(residuals(m_image))
qqmath(ranef(m_image))

data.frame(x = residuals(m_image)) %>% 
  ggplot(aes(x = x)) +
  geom_histogram()

summary(m_image)
summary(m_image_1)
summary(m_image_1_5)

