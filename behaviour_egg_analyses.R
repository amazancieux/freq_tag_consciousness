#####################################


# Behavioural analyses of the frequency tagging 
# experiment using male and female faces with
# at a supraliminal contrast (session 2) 

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
  axis.title.x = element_text(size = 16),
  plot.title = element_text(size = 18),
  axis.text.x = element_text(size = 14),
  axis.title.y = element_text(size = 16))


## Import and clean data -------------------------------------------------------------

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


## First-order task ------------------------------------------------------

# performance count
first_order <- resp_data %>% 
  mutate(count = 1) %>% 
  dcast(Pp ~ accuracy, value.var = "count", sum) %>% 
  mutate(performance = Correct/(Correct+Incorrect))


## Subjective visibility -------------------------------------------------

# PAS mean
count_pas <- resp_data %>% 
  mutate(count = 1) %>%
  dcast(Pp ~ pas_score+accuracy, value.var = "count", sum)


# plot mean contrast per PAS
resp_data %>%
  group_by(Pp, accuracy) %>%
  summarise(pas_score = mean(as.numeric(pas_score))) %>% 
  group_by(accuracy) %>%
  summarise(VD = mean(pas_score),
            sd = sd(pas_score),
            se = sd/sqrt(length(unique(resp_data$Pp))),
            CI = se * qt(.975, n() - 1)) %>%
  ggplot(aes(x = accuracy, y = VD, fill=accuracy)) +
  geom_bar(stat="identity", color='black') +
  geom_errorbar(aes(ymin = VD - CI, ymax = VD + CI), width = 0, size = 1)+
  ggtitle("mean PAS according to accuracy") +
  scale_fill_brewer(palette="Dark2") +
  theme_bw() +
  plot_theme +
  xlab("Accuracy") +
  ylab("mean PAS")



## Confidence -----------------------------------------------------------

# raw confidence per accuracy
resp_data %>%
  group_by(Pp, accuracy) %>%
  summarise(conf_score = mean(as.numeric(conf_score))) %>% 
  group_by(accuracy) %>%
  summarise(VD = mean(conf_score),
            sd = sd(conf_score),
            se = sd/sqrt(length(unique(resp_data$Pp))),
            CI = se * qt(.975, n() - 1)) %>%
  ggplot(aes(x = accuracy, y = VD, fill=accuracy)) +
  geom_bar(stat="identity", color='black') +
  geom_errorbar(aes(ymin = VD - CI, ymax = VD + CI), width = 0, size = 1)+
  ggtitle("mean confidence according to accuracy") +
  scale_fill_brewer(palette="Dark2") +
  theme_bw() +
  plot_theme +
  xlab("Accuracy") +
  ylab("mean confidence")


## Metacognitive efficiency 

source("Function_metad_indiv.R")

# prepare data
nR_S1 <- resp_data_short %>% 
  filter(correct_response == 102) %>% 
  mutate(Count = 1,
         Resp = case_when(
           decision_resp == 102 ~ "R1",
           decision_resp == 104 ~ "R2")) %>% 
  dcast(Pp ~ conf_score + Resp, value.var = "Count", sum) %>% 
  select(`4_R1`,
         `3_R1`,
         `2_R1`,
         `1_R1`,
         `1_R2`)
#         `2_R2`,
#         `3_R2`,
#         `4_R2`)

nR_S2 <- resp_data_short %>% 
  filter(correct_response == 104) %>% 
  mutate(Count = 1,
         Resp = case_when(
           decision_resp == 102 ~ "R1",
           decision_resp == 104 ~ "R2")) %>% 
  dcast(Pp ~ conf_score + Resp, value.var = "Count", sum) %>% 
#  select(`4_R1`,
  select(`3_R1`,
#         `2_R1`,
         `1_R1`,
         `1_R2`,
         `2_R2`,
         `3_R2`,
         `4_R2`)

nsubj <- nrow(nR_S1)

# calculate d' and meta-d' per participant
c <- data.frame()
d <- data.frame()
stats <- data.frame()

for (n in 1:(nsubj)) {
  
  S1 <- c(t(nR_S1[n,]))
  S2 <- c(t(nR_S2[n,]))
  
  output <- metad_indiv(nR_S1 = S1, nR_S2 = S2)
  
  d %<>% rbind(d1) 
  c %<>% rbind(c1) 
  
  Value <- summary(output)
  
  stat1 <- data.frame(mean = Value[["statistics"]][, "Mean"])
  stat %<>% rbind(t(stat1)) 
  
}

write.csv(stat, "./Behaviour/Results/output_fit_individual_metad.csv")


