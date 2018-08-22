# Calculate correlation between clocklike mutation signature 1 and 5 with age.
# 2018.07.15 CJY
# Run @ my laptop
# Did not have age info for 8 samples from batch2
# Will need to rerun when age info of these 8 additional samples become available. DONE 2018.07.16

library(tidyverse)
library(ggplot2)
library(readr)
library(readxl)
library(stringr)

# import excel sheet
sample_summary = read_excel("~/Documents/julab/projects/myeloma/sample_summary.xlsx")
sample_summary_perGB = sample_summary %>% mutate(ClocklikeMutationsPerGB = ClocklikeMutations/3) # Human Genome is 3GB

# linear regression
age_clockwise = lm(ClocklikeMutationsPerGB ~ Age, data=sample_summary_perGB)
summary(age_clockwise)
# Call:
#   lm(formula = ClocklikeMutationsPerGB ~ Age, data = sample_summary_perGB)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -965.47 -192.04  -83.78  183.92 1530.96 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept)   203.50     784.42   0.259    0.797
# Age            12.18      11.94   1.020    0.315
# 
# Residual standard error: 582.2 on 34 degrees of freedom
# (2 observations deleted due to missingness)
# Multiple R-squared:  0.02968,	Adjusted R-squared:  0.001138 
# F-statistic:  1.04 on 1 and 34 DF,  p-value: 0.3151


# draw scatter plot
ggplot(sample_summary_perGB, aes(x=Age, y=ClocklikeMutationsPerGB)) + geom_point() +   xlim(c(0, 90)) +
  geom_smooth(method = "lm", se = T) 
ggsave('mm37_clocklikemutations.png')



# only with signature 1
sig1_with_dataavailable = sample_summary_perGB %>% filter(Sig1 != 0 )
age_sig1 = lm(Sig1 ~ Age, data=sig1_with_dataavailable)

ggplot(sig1_with_dataavailable, aes(x=Age, y=Sig1)) + geom_point()  +
  geom_smooth(method = "lm", se = T)
