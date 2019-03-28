###########################################################################################
#########                                                                         #########
#########                     Daniel Elias Martin Herranz                         #########
#########                             15/03/2019                                  #########
#########                              EMBL-EBI                                   #########
#########                           Thornton group                                #########
#########                                                                         #########
###########################################################################################

###########################################################################################
#####                                 PhD Thesis                                       ####
###########################################################################################
##### Additional code for Chapter 3.                                                   ####
###########################################################################################
##### USAGE: manual                                                                    ####
###########################################################################################

setwd('~/Documents/Thesis/PhD-Thesis/additional_code/Chapter3/');

library(data.table);
library(ggplot2);

#### 1. Basic description of the cases dataset. ####

raw_cases <- as.data.frame(fread('cases_data_downstream.tsv'));
table(raw_cases$Batch, raw_cases$Sex);
median(raw_cases$Age_years[raw_cases$Batch=='Europe']); # Change the batch as necessary
min(raw_cases$Age_years); max(raw_cases$Age_years);

# Plot the chronological age distribution.

plot_age_distributions_cases <- ggplot(data=raw_cases, aes(x=Age_years)) + 
  geom_histogram(aes(y=..density..),colour="black", fill='grey') + stat_density(geom="line", col='blue', size=1.5) +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) +
  xlab("Chronological age (years)") + ylab("Density") +
  labs(title = paste0("Cases: N=", nrow(raw_cases)));
ggsave("plots/chronological_age_dist_cases.pdf", height=5, width=5);

# Find the peak. 

D <- density(raw_cases$Age_years);
D$x[which.max(D$y)]; # 4.5 years

#### End of the script. ####


