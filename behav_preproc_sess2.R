#####################################


# Behavioural data preprocessing for session 2
# of the frequency tagging experiment

# Audrey Mazancieux 04/2024


#####################################


# Packages
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


# get data for this subject
file <- './Behaviour/Data/sub-1/FreqTagStim_sub-1.csv'
pp_data <- read_csv(file)

# get relevant variables
pp_data %<>% 
  filter(decision_resp != "NaN") %>% 
  mutate(sdt = case_when(
    decision_resp == 102 & correct_response == 102 ~ "H",   # signal = 102 (female)
    decision_resp == 102 & correct_response == 104 ~ "FA",
    decision_resp == 104 & correct_response == 102 ~ "O",
    decision_resp == 104 & correct_response == 104 ~ "CR"),
    accuracy = ifelse(decision_resp == correct_response, 'Correct', 'Incorrect'), 
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
      conf_resp == 54 ~ "6")) %>% 
  select(-contrast)

# task performance
table(pp_data$accuracy)

# save data
write.csv(pp_data, "./Behaviour/Data/sub-1/FreqTagStim_sub-1_preproc.csv")

