#####################################


# Compute detection rate for the fixation
# cross change task of the frequency 
# tagging experiment 

# Audrey Mazancieux 2024


#####################################


library(tidyverse)
library(magrittr)
library(reshape2)


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

fixCross_seq <- data.frame()
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
  fixCross_seq %<>%
    rbind(detect)
  
  # calculate ratio for each contrast (when column 32 is 0 then there is no detection for this cross)
  detect_1 <- detect %>% 
    filter(contrast == '1%')
  ratio_1 = 1 - table(c(detect_1$`32`))["0"]/nrow(detect_1) # 32 corresponds to the ASCII code for space bar 
  
  detect_1.5 <- detect %>% 
    filter(contrast == '1.5%')
  ratio_1.5 = 1 - table(c(detect_1.5$`32`))["0"]/nrow(detect_1.5) 
  
  data <- data.frame(sub = detect$sub[1] %>% 
                       str_extract(regex("\\d+")), 
                     ratio_1 = ratio_1,
                     ratio_1.5 = ratio_1.5)

  # add this subject to dataframe
  fixCross %<>% 
    rbind(data)
  
}

# save
write.csv(fixCross, "./result_fixCross.csv")


## Perform analyses and plots --------------------------------------------

# plot mean detection per contrast
fixCross %>%
  gather(contrast, perf, -sub) %>% 
  mutate(contrast = ifelse(str_detect(contrast, regex("1.5")), '1.5%', '1%')) %>% 
  group_by(contrast) %>%
  summarise(VD = mean(perf),
            sd = sd(perf),
            se = sd/sqrt(length(unique(fixCross$sub))),
            CI = se * qt(.975, n() - 1)) %>%
  ggplot(aes(x = contrast, y = VD, fill=contrast)) +
  geom_bar(stat="identity", color='black') +
  geom_errorbar(aes(ymin = VD - CI, ymax = VD + CI), width = 0, size = 0.5)+
  ggtitle("mean cross detection ratio per contrast") +
  scale_fill_brewer(palette="Dark2") +
  theme_classic() +
  xlab("Contrast") +
  ylab("mean performance")

# statistics
mean(fixCross$ratio_1)
mean(fixCross$ratio_1.5)
sd(fixCross$ratio_1)
sd(fixCross$ratio_1.5)

fixCross %<>%
  mutate(perf_diff = ratio_1 - ratio_1.5) 
mod <- lm(perf_diff ~ 1, data = fixCross)

qqnorm(residuals(mod))
qqline(residuals(mod))

data.frame(x = residuals(mod)) %>%
  ggplot(aes(x = x)) +
  geom_histogram()
shapiro.test(residuals(mod))

summary(mod)


# get performance according to behaviour for contrast 1.5%
behaviour_dataset <- read.table("./Behaviour/Results/behaviour_clean_dataset.csv", header = TRUE, sep=",", dec=".", fill  = TRUE)

fixCross_seq2 <- fixCross_seq %<>%
  mutate(correct = `32`,
         block_num = contrast_type) %>% 
  dcast(sub ~ block_num, value.var = 'correct', mean) %>% 
  gather(block_num, perf, -sub)

fixCross_seq2 %<>%
  mutate(sub2 = sub %>% 
           str_extract(regex("(?<=-)..")),
         subject = as.numeric(sub2)) %>% 
  select(subject, block_num, perf) %>% 
  arrange(subject, block_num)

data_cross_behaviour <- behaviour_dataset %>% 
  select(subject, block_num, contrast, pas_score, conf_score, accuracy)

data_cross_behaviour <- merge(data_cross_behaviour, fixCross_seq2, by=c('subject', 'block_num'))


# plot mean detection per PAS
data_cross_behaviour %>%
  mutate(pas = as.factor(pas_score)) %>% 
  group_by(pas, contrast) %>%
  summarise(VD = mean(perf, na.rm = TRUE),
            sd = sd(perf, na.rm = TRUE),
            se = sd/sqrt(length(unique(fixCross$sub))),
            CI = se * qt(.975, n() - 1)) %>%
  ggplot(aes(x = pas, y = VD, fill=contrast)) +
  geom_bar(stat="identity", color='black', position = position_dodge(0.9)) +
  geom_errorbar(aes(ymin = VD - CI, ymax = VD + CI), width = 0, size = 0.5, position = position_dodge(0.9))+
  ggtitle("mean cross detection ratio per PAS score") +
  scale_fill_brewer(palette="Dark2") +
  theme_classic() +
  xlab("PAS score") +
  ylab("mean performance")

