#####################################


# Compute detection rate for the fixation
# cross change task of the frequency 
# tagging experiment 

# Audrey Mazancieux 2024


#####################################


## Packages ----------------------------------------------

library(tidyverse)
library(magrittr)
library(reshape2)

# plot theme
plot_theme = theme(
  axis.title.x = element_text(size = 12),
  plot.title = element_text(size = 14),
  axis.text.x = element_text(size = 12),
  axis.title.y = element_text(size = 12))


## Get percentage of fixation cross change detection ---------------------

# get data
files <-
  list.files("./Behaviour/Data", recursive = TRUE) %>% 
  as_data_frame()

files <- files %>% 
  filter(str_detect(value, regex(".csv$"))) %>% 
  filter(str_detect(value, regex("fixCross_sub"))) %>% 
  filter(str_detect(value, regex("fixCross_sub"))) 

files <- data.frame(files) %>% 
  mutate(tokeep = ifelse(str_detect(value, regex("noreport")), 0, 1)) %>% 
  filter(tokeep == 1)

fixCross <- data.frame()
for (i in files$value){
  
  # import data
  pp_data <- read.table(file.path("./Behaviour/Data", i), header = TRUE, sep=",", dec=".", fill  = TRUE) %>% 
    filter(fixcross_vector == 1) %>% 
    mutate(sub = i %>% 
             str_extract(regex("sub-\\d+")),
           contrast = block_num %>% 
             str_extract(regex("(?<=_)...")),
           contrast = ifelse(contrast == '1.5', '1.5%', '1%'))
  
  # find the number of cross changes
  cross_change <- data.frame()
  for (l in (1:nrow(pp_data))){
    
    if (l == 1){
      cross_change %<>%
        rbind(1)
    } else if (l == nrow(pp_data)){
      last <- nrow(cross_change)
      cross_change %<>%
        rbind(cross_change[last,])
    } else if (pp_data$frame[l+1] - pp_data$frame[l] > 1){
      last <- nrow(cross_change)
      cross_change %<>%
        rbind(cross_change[last,]+1)
    } else {
      last <- nrow(cross_change)
      cross_change %<>%
        rbind(cross_change[last,])
    }
  }
  pp_data %<>%
    cbind(cross_change) %>% 
    mutate(cross_num = X1)
  
  # get whether the cross change has been detected 
  detect <- pp_data %>% 
    mutate(count = 1) %>% 
    dcast(sub + contrast + contrast_type + cross_num ~ fixcross_resp, value.var = "count", sum)  
  
  # add this subject to dataframe
  fixCross %<>%
    rbind(detect)
  
}

# save
write.csv(fixCross, "./Behaviour/Results/result_fixCross.csv")

