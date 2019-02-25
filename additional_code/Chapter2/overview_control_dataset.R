###########################################################################################
#########                                                                         #########
#########                     Daniel Elias Martin Herranz                         #########
#########                             20/02/2019                                  #########
#########                              EMBL-EBI                                   #########
#########                           Thornton group                                #########
#########                                                                         #########
###########################################################################################

###########################################################################################
#####                                 PhD Thesis                                       ####
###########################################################################################
##### Overview on the control blood dataset.                                           ####
###########################################################################################
##### USAGE: manual                                                                    ####
###########################################################################################

setwd('/Users/dem44/Documents/Thesis/PhD-Thesis/additional_code/Chapter2/');

library(data.table);
library(ggplot2);
library(flexmix);


#### 1. Basic description of the control dataset. ####

raw_controls <- as.data.frame(fread('final_control_data.tsv'));
table(raw_controls$Batch, raw_controls$Sex);
median(raw_controls$Age_years[raw_controls$Batch=='Europe']); # Change the batch as necessary
min(raw_controls$Age_years); max(raw_controls$Age_years);

# Plot the chronological age distribution.

plot_age_distributions_controls <- ggplot(data=raw_controls, aes(x=Age_years)) + 
  geom_histogram(aes(y=..density..),colour="black", fill='grey') + stat_density(geom="line", col='blue', size=1.5) +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) +
  xlab("Chronological age (years)") + ylab("Density") +
  labs(title = paste0("Control: N=", nrow(raw_controls)));
ggsave("plots/chronological_age_dist_controls.pdf", height=5, width=5);

# Find the peaks. 

D <- density(raw_controls$Age_years);
D$x[which.max(D$y)];
D$x[which.max(D$y[1:150])];

