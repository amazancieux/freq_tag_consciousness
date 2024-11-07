#####################################


# Compute detection rate for the fixation
# cross change task of the frequency 
# tagging experiment 

# Audrey Mazancieux 2024


#####################################


## Packages and load data ----------------------------------------------

library(tidyverse)
library(magrittr)
library(reshape2)
library(lmerTest)
library(lattice)

# plot theme
plot_theme = theme(
  axis.title.x = element_text(size = 12),
  plot.title = element_text(size = 14),
  axis.text.x = element_text(size = 12),
  axis.title.y = element_text(size = 12))

# load data
fixCross <- read.csv2("./Behaviour/Results/result_fixCross.csv", sep = ',', dec = '.')
behaviour_dataset <- read.table("./Behaviour/Results/behaviour_clean_data.csv", header = TRUE, sep=",", dec=".", fill  = TRUE)


## Plot data ---------------------------------------------------------------

# Get performance according to PAS
fixCross2 <- fixCross %<>%
  mutate(correct = `X32`,
         block_num = contrast_type) %>% 
  dcast(sub ~ block_num, value.var = 'correct', mean) %>% 
  gather(block_num, perf, -sub)

fixCross2 %<>%
  mutate(sub2 = sub %>% 
           str_extract(regex("\\d+")),
         subject = as.numeric(sub2)) %>% 
  select(subject, block_num, perf) %>% 
  arrange(subject, block_num)

data_cross_behaviour <- behaviour_dataset %>% 
  mutate(subject = as.numeric(Pp)) %>% 
  select(subject, block_num, contrast, pas_score, conf_score, accuracy) %>% 
  arrange(subject)

data_cross_behaviour <- merge(data_cross_behaviour, fixCross2, by=c('subject', 'block_num'))


# plot mean detection per PAS
png(file="./Behaviour/Results/fixCross_contrast_pas.png", width=6, height=4, units="in", res=300)
data_cross_behaviour %>%
  mutate(pas = as.factor(pas_score)) %>% 
  group_by(subject, pas, contrast) %>%
  summarise(VD = mean(perf, na.rm = TRUE)) %>% 
  ggplot(aes(x = pas, y = VD, color=contrast)) + 
  geom_point(position = position_jitterdodge(0.3), size = 1, alpha = 0.5) + 
  geom_smooth(data = data_cross_behaviour,
              aes(x = pas_score, y = perf, color=contrast), 
              method=lm, se=FALSE, fullrange=TRUE, size = 1.5, formula = y ~ x) + 
  ggtitle("mean cross detection ratio per confidence rating") +
  scale_color_manual(values = c("#0F056B", "#9683EC")) +
  theme_classic() +
  plot_theme + 
  xlab("PAS ratings") +
  ylab("mean performance for cross detection")
dev.off()


# plot mean detection per accuracy
png(file="./Behaviour/Results/fixCross_contrast_acc.png", width=6, height=6, units="in", res=300)
data_cross_behaviour %>%
  group_by(accuracy, contrast) %>%
  summarise(VD = mean(perf, na.rm = TRUE),
            sd = sd(perf, na.rm = TRUE),
            se = sd/sqrt(length(unique(fixCross$sub))),
            CI = se * qt(.975, n() - 1)) %>%
  ggplot(aes(x = accuracy, y = VD, fill=contrast)) +
  geom_bar(stat="identity", color='black', position = position_dodge(0.9)) +
  geom_errorbar(aes(ymin = VD - CI, ymax = VD + CI), width = 0, size = 0.5, position = position_dodge(0.9))+
  ggtitle("mean cross detection ratio per accuracy") +
  scale_fill_manual(values = c("#0F056B", "#9683EC")) +
  theme_classic() +
  plot_theme + 
  xlab("accuracy") +
  ylab("mean performance for cross detection")
dev.off()

# plot mean detection per confidence
png(file="./Behaviour/Results/fixCross_contrast_conf.png", width=6, height=4, units="in", res=300)
data_cross_behaviour %>%
  mutate(conf = as.factor(conf_score)) %>% 
  group_by(subject, conf, contrast) %>%
  summarise(VD = mean(perf, na.rm = TRUE)) %>% 
  ggplot(aes(x = conf, y = VD, color=contrast)) + 
  geom_point(position = position_jitterdodge(0.3), size = 1, alpha = 0.5) + 
  geom_smooth(data = data_cross_behaviour,
              aes(x = conf_score, y = perf, color=contrast), 
              method=lm, se=FALSE, fullrange=TRUE, size = 1.5, formula = y ~ x) + 
  ggtitle("mean cross detection ratio per confidence rating") +
  scale_color_manual(values = c("#0F056B", "#9683EC")) +
  theme_classic() +
  plot_theme + 
  xlab("confidence ratings") +
  ylab("mean performance for cross detection")
dev.off()


## Mixed-effect models  ----------------------------------------------------

# create coding 
data_cross_behaviour %<>%
  mutate(accC = ifelse(accuracy == 'Incorrect', -0.5, 0.5),
         contrastC = ifelse(contrast == '1%', -0.5, 0.5),
         C_1 = ifelse(contrast == '1%', 0, 1),
         C_1_5 = ifelse(contrast == '1.5%', 0, 1))

# models 
m_pas <- lmer(perf ~ contrastC * pas_score + (pas_score|subject), data = data_cross_behaviour)
m_acc <- lmer(perf ~ contrastC * accC + (accC|subject), data = data_cross_behaviour)
m_conf <- lmer(perf ~ contrastC * conf_score + (contrastC|subject), data = data_cross_behaviour)
# use PCA to choose random effects 
summary(rePCA(m_pas))
summary(rePCA(m_acc))
summary(rePCA(m_conf))

m_pas <- lmer(perf ~ contrastC * pas_score + (1|subject), data = data_cross_behaviour)
m_acc <- lmer(perf ~ contrastC * accC + (1|subject), data = data_cross_behaviour)
m_conf <- lmer(perf ~ contrastC * conf_score + (1|subject), data = data_cross_behaviour)

qqnorm(residuals(m_pas))
qqline(residuals(m_pas))
qqmath(ranef(m_pas))

data.frame(x = residuals(m_pas)) %>% 
  ggplot(aes(x = x)) +
  geom_histogram()

summary(m_pas)
summary(m_acc)
summary(m_conf)

m_conf_1 <- lmer(perf ~ C_1 * conf_score + (1|subject), data = data_cross_behaviour)
m_conf_1_5 <- lmer(perf ~ C_1_5 * conf_score + (1|subject), data = data_cross_behaviour)

summary(m_conf_1)
summary(m_conf_1_5)

