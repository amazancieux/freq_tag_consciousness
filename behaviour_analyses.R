#####################################


# Behavioural analyses of the frequency tagging 
# experiment using male and female faces with
# different contrasts (session 1)

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

# get data for all subject
files <-
  list.files("./Behaviour/Data", recursive = TRUE) %>% 
  as_data_frame()

files_cross <- files %>% 
  filter(str_detect(value, regex(".csv$"))) %>% 
  filter(str_detect(value, regex("fix_staircase")))

files_resp <- files %>% 
  filter(str_detect(value, regex(".csv$"))) %>% 
  filter(str_detect(value, regex("Stim_staircase")))


# load fix cross data
fix_cross <- data.frame()

for (i in files_cross$value){
  
  pp_data <- read_csv(file.path("./Behaviour/Data", i)) %>% 
    mutate(Pp = i %>% 
             str_extract(regex("\\d+")))
  fix_cross %<>% rbind(pp_data) 
}

# load response data
resp_data <- data.frame()

for (i in files_resp$value){
  
  pp_data <- read_csv(file.path("./Behaviour/Data", i)) %>% 
    mutate(Pp = i %>% 
             str_extract(regex("\\d+")))
  resp_data %<>% rbind(pp_data) 
}


# get relevant variables
resp_data_short <- resp_data %>% 
  filter(decision_resp != "NaN") %>% 
  mutate(accuracy = case_when(
    decision_resp == 102 & correct_response == 102 ~ "H",   # signal = 102 (female)
    decision_resp == 102 & correct_response == 104 ~ "FA",
    decision_resp == 104 & correct_response == 102 ~ "O",
    decision_resp == 104 & correct_response == 104 ~ "CR"),
    contrast = contrast_type %>% 
      str_extract(regex("\\d+.\\d+")),
    pas_score = case_when(
      pas_resp == 49 ~ "1",
      pas_resp == 50 ~ "2",
      pas_resp == 51 ~ "3",
      pas_resp == 52 ~ "4"),
    conf_score = case_when(
      conf_resp == 49 ~ "1",
      conf_resp == 50 ~ "2",
      conf_resp == 51 ~ "3",
      conf_resp == 52 ~ "4",
      conf_resp == 53 ~ "5",
      conf_resp == 54 ~ "6"))

n <- length(unique(resp_data_short$Pp))


## First-order task ------------------------------------------------------

# performance count
count_first_order <- resp_data_short %>% 
  mutate(count = 1) %>% 
  dcast(Pp ~ accuracy, value.var = "count", sum)


# get more frequent contrast
contrasts <- resp_data_short %>% 
  mutate(count = 1,
         contrast = contrast_type %>% 
           str_extract(regex("\\d+.\\d+"))) 

contrast_resume <- contrasts %>%
  dcast(Pp ~ contrast, value.var = "count", sum)

# get performance per contrast per participants
perf_contrast_pp <- resp_data_short %>%
  mutate(acc = ifelse(decision_resp == correct_response, 1, 0)) %>% 
  dcast(Pp ~ contrast, value.var = "acc", mean)

# plot performance across contrast 
resp_data_short %>%
  mutate(acc = ifelse(decision_resp == correct_response, 1, 0)) %>% 
  group_by(Pp, contrast) %>%
  summarise(acc = mean(acc)) %>% 
  group_by(contrast) %>%
  summarise(VD = mean(acc),
            sd = sd(acc),
            se = sd/sqrt(n),
            CI = se * qt(.975, n() - 1)) %>%
  ggplot(aes(x = contrast, y = VD, fill=contrast)) +
  geom_bar(stat="identity", color='black') +
  geom_errorbar(aes(ymin = VD - CI, ymax = VD + CI), width = 0, size = 1)+
  ggtitle("mean performance across contrast levels") +
  theme_bw() +
  plot_theme +
  xlab("Contrast level") +
  ylab("mean performance")


## Subjective visibility -------------------------------------------------

# plot mean contrast per PAS
resp_data_short %>%
  group_by(Pp, pas_score) %>%
  summarise(contrast = mean(as.numeric(contrast))) %>% 
  group_by(pas_score) %>%
  summarise(VD = mean(contrast),
            sd = sd(contrast),
            se = sd/sqrt(length(unique(resp_data_short$Pp))),
            CI = se * qt(.975, n() - 1)) %>%
  ggplot(aes(x = pas_score, y = VD, fill=pas_score)) +
  geom_bar(stat="identity", color='black') +
  geom_errorbar(aes(ymin = VD - CI, ymax = VD + CI), width = 0, size = 1)+
  ggtitle("mean contrast per PAS rating") +
  scale_fill_brewer(palette="Dark2") +
  theme_bw() +
  plot_theme +
  xlab("PAS") +
  ylab("mean contrast")



## Confidence -----------------------------------------------------------

# raw confidence per contrast
resp_data_short %>%
  group_by(Pp, conf_score) %>%
  summarise(contrast = mean(as.numeric(contrast))) %>% 
  group_by(conf_score) %>% 
  summarise(VD = mean(contrast),
            sd = sd(contrast),
            se = sd/sqrt(n),
            CI = se * qt(.975, n() - 1)) %>%
  ggplot(aes(x = conf_score, y = VD, fill=conf_score)) +
  geom_bar(stat="identity", color='black') +
  geom_errorbar(aes(ymin = VD - CI, ymax = VD + CI), width = 0, size = 1)+
  ggtitle("mean contrast per confidence rating") +
  scale_fill_brewer(palette="Dark2") +
  theme_bw() +
  plot_theme +
  xlab("Confidence") +
  ylab("mean contrast")

# raw confidence per contrast
resp_data_short %>%
  mutate(acc = ifelse(decision_resp == correct_response, "Correct", "Incorrect")) %>% 
  group_by(Pp, acc) %>%
  summarise(conf_score = mean(as.numeric(conf_score))) %>% 
  group_by(acc) %>%
  summarise(VD = mean(conf_score),
            sd = sd(conf_score),
            se = sd/sqrt(n),
            CI = se * qt(.975, n() - 1)) %>%
  ggplot(aes(x = acc, y = VD, fill=acc)) +
  geom_point(size = 3, color = "purple") +
  geom_errorbar(aes(ymin = VD - CI, ymax = VD + CI), width = 0, size = 1, color = "purple")+
  ggtitle("mean confidence for correct and incorrect responses") +
  theme_bw() +
  plot_theme +
  xlab("Discrimination accuracy") +
  ylab("mean confidence")


