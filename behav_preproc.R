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


subjects = list(3, 14, 15, 17, 18, 19, 20, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 34, 35, 37, 38, 39, 41, 42, 43, 44, 45, 47, 48, 49, 50)

for (sub in subjects){
  
  # get data for this subject
  file <- sprintf('./Behaviour/Data/sub-%s/FreqTagStim_sub-%s.csv', sub, sub)
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
        conf_resp == 54 ~ "6"),
      contrast = case_when(
        contrast_type == 'Contrast_1_female' ~ '1%',
        contrast_type == 'Contrast_1_male' ~ '1%',
        contrast_type == 'Contrast_1.5_female' ~ '1.5%',
        contrast_type == 'Contrast_1.5_male' ~ '1.5%')) 
  
  # task performance
  table(pp_data$accuracy)
  
  # save data
  new_file <- sprintf('./Behaviour/Data/sub-%s/FreqTagStim_sub-%s_preproc.csv', sub, sub)
  write.csv(pp_data, new_file)
  
}


