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
library(RColorBrewer);
library(gridExtra);
library(cowplot);
library(IlluminaHumanMethylation450kanno.ilmn12.hg19);
library(ggthemes);
library(ggpubr);


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


#### 2. Changes of blood cell composition with age. ####

ct_palette <- brewer.pal(6,"Dark2")[order(colnames(raw_controls)[11:16])];
plot_cts <- list();
i <- 1;

for(c in colnames(raw_controls)[11:16]){
  
  to_plot_i <- data.frame(Age=raw_controls$Age_years, CT=as.numeric(raw_controls[,i+10])*100);
  sp <- as.numeric(cor.test(to_plot_i[,1], to_plot_i[,2], method='spearman')$estimate);
  pv <- as.numeric(cor.test(to_plot_i[,1], to_plot_i[,2], method='spearman')$p.value);
  lmi_int <- as.numeric(lm(CT~Age, data=to_plot_i)$coefficients)[1];
  lmi_slope <- as.numeric(lm(CT~Age, data=to_plot_i)$coefficients)[2];
  plot_cts[[i]] <- ggplot(data=to_plot_i, aes(x=Age, y=CT)) + geom_point(col=ct_palette[i]) +
    geom_smooth(method="lm", formula=y~x, show.legend=F, alpha = 0.3, col='black', fill='darkgrey') +
    theme_classic() +
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14,face="bold"),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          plot.margin = unit(c(1,1,1,1), "cm")) +
    xlab("Chronological age (years)") + ylab("Density") +
    labs(title = paste0("SCC: ", round(sp, digits=4), 
                        ifelse(pv < 2.2e-16, '; p-value < 2.2e-16', paste0('; p-value = ', format(pv, digits=4))),
                        paste0("\n Model: y = ", round(lmi_int,digits=4), 
                               ifelse(lmi_slope > 0, " + ", " - "),
                               abs(round(lmi_slope, digits=4)), "*x"))) +
    ylim(c(-1,ifelse(c=='Gran', 100, 50))) + xlab('Chronological age (years)') + ylab(paste0(c, ' (%)'));
  i <- i+1;
    
}

exp <- do.call("arrangeGrob", c(plot_cts, nrow=3, ncol=2));
ggsave(file="plots/cc_during_ageing.pdf", exp, height=15, width=10);


#### 3. aDMPs (full lifespan analysis). ####

## Read the calculated aDMPs.

aDMPs_with_CCC <- as.data.frame(fread('aDMPs_final_full_lifespan_with_CCC.csv'));
aDMPs_without_CCC <- as.data.frame(fread('aDMPs_final_full_lifespan_without_CCC.csv'));

## Order them by p-value and effect size and add info about Horvath sites. 

thr <- 0.01 / nrow(aDMPs_with_CCC); # p-value threshold (after Bonferroni correction)
horvath_sites <- as.data.frame(fread('AdditionalFile3.csv'));
horvath_sites <- horvath_sites[-1,];

aDMPs_with_CCC <- aDMPs_with_CCC[order(as.numeric(aDMPs_with_CCC$pval), -as.numeric(abs(aDMPs_with_CCC$t))),];
aDMPs_with_CCC$Direction <- ifelse(aDMPs_with_CCC$pval>=thr, 'No change', ifelse(
  aDMPs_with_CCC$beta>0, 'Hypermethylated', 'Hypomethylated'));
aDMPs_with_CCC$shape <- NA;
aDMPs_with_CCC$shape[aDMPs_with_CCC$ProbeID %in% horvath_sites$CpGmarker] <- 'In Horvath model';
aDMPs_without_CCC <- aDMPs_without_CCC[order(as.numeric(aDMPs_without_CCC$pval), -as.numeric(abs(aDMPs_without_CCC$t))),];
aDMPs_without_CCC$Direction <- ifelse(aDMPs_without_CCC$pval>=thr, 'No change', ifelse(
  aDMPs_without_CCC$beta>0, 'Hypermethylated', 'Hypomethylated'));
aDMPs_without_CCC$shape <- NA;
aDMPs_without_CCC$shape[aDMPs_without_CCC$ProbeID %in% horvath_sites$CpGmarker] <- 'In Horvath model';

## Export top 100 aDMPs for donwstream use.

# List of ProbeIDs.
top100_probes_with_CCC <- aDMPs_with_CCC$ProbeID[1:100];
write(top100_probes_with_CCC, file='top100_aDMPs_with_CCC.txt');
top100_probes_without_CCC <- aDMPs_without_CCC$ProbeID[1:100];
write(top100_probes_without_CCC, file='top100_aDMPs_without_CCC.txt');

# Table.
table_aDMPs_with_CCC <- aDMPs_with_CCC[1:100,];
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19);
aDMPs_ann <- data.frame(ProbeID=ann450k$Name[ann450k$Name %in% table_aDMPs_with_CCC$ProbeID],
                        Chr=ann450k$chr[ann450k$Name %in% table_aDMPs_with_CCC$ProbeID],
                        Coord=ann450k$pos[ann450k$Name %in% table_aDMPs_with_CCC$ProbeID],
                        Genes=ann450k$UCSC_RefGene_Name[ann450k$Name %in% table_aDMPs_with_CCC$ProbeID]);
table_aDMPs_with_CCC <- merge(table_aDMPs_with_CCC, aDMPs_ann, by='ProbeID');
table_aDMPs_with_CCC <- table_aDMPs_with_CCC[order(as.numeric(table_aDMPs_with_CCC$pval), -as.numeric(abs(table_aDMPs_with_CCC$t))),];                        
table_aDMPs_with_CCC <- table_aDMPs_with_CCC[,-6];
colnames(table_aDMPs_with_CCC) <- c('ProbeID', 'Intercept', 'Slope', 'T statistic', 'p-value', 'Methylation change', 
                                    'In Horvath model', 'Chromosome', 'Coordinate', 'Gene(s)');
table_aDMPs_with_CCC <- table_aDMPs_with_CCC[,c(1,8,9,2:7,10)];
table_aDMPs_with_CCC$`In Horvath model` <- ifelse(is.na(table_aDMPs_with_CCC$`In Horvath model`), 'No', 'Yes');
write.table(table_aDMPs_with_CCC, file='top100_aDMPs_with_CCC_table.csv', quote=F, sep=',', row.names=F);


## Create volcano-like plots.

vp_with_CCC <- ggplot(data=aDMPs_with_CCC, aes(x=beta, y=-log10(pval), col=Direction)) +
  geom_point(alpha=0.3) + scale_colour_manual(values=c('Red', 'Blue', 'Grey')) +
  geom_point(data=aDMPs_with_CCC[aDMPs_with_CCC$ProbeID %in% horvath_sites$CpGmarker,], aes(x=beta, y=-log10(pval), shape=shape), col='black', size=2) +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        plot.margin=margin(-3,0,0,0),
        legend.text=element_text(size=11),
        legend.title=element_text(size=12,face="bold")) +
  xlab(bquote(bold('Age slope ('*beta*'-value change/year)'))) + ylab(expression(bold(paste(-log[10](P-value))))) +
  guides(colour=guide_legend(title="Methylation change"), shape=guide_legend(title="")) + 
  geom_abline(slope=0, intercept=-log10(thr), linetype=2, col='forestgreen', size=0.8);
d <-ggplot(data=aDMPs_with_CCC[aDMPs_with_CCC$Direction!='No change',], aes(x=beta, fill=Direction)) + 
  geom_density() + scale_fill_manual(values=c('Red', 'Blue')) +
  theme_classic() +
  theme(axis.line=element_blank(),
        axis.text.x=element_blank(), axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(), axis.title.y=element_blank());
comb_with_CCC <- plot_grid(plot_grid(d+theme(legend.position = "none"), vp_with_CCC+theme(legend.position = "none"), ncol = 1, align = "h",rel_heights=c(1,5)), 
                           plot_grid(ggplot(), get_legend(vp_with_CCC), ncol=1, rel_heights=c(1,5)), rel_widths=c(3,1));
ggsave('plots/vp_aDMPs_with_CCC.png', height=7, width=7);

vp_without_CCC <- ggplot(data=aDMPs_without_CCC, aes(x=beta, y=-log10(pval), col=Direction)) +
  geom_point(alpha=0.3) + scale_colour_manual(values=c('Red', 'Blue', 'Grey')) +
  geom_point(data=aDMPs_without_CCC[aDMPs_without_CCC$ProbeID %in% horvath_sites$CpGmarker,], aes(x=beta, y=-log10(pval), shape=shape), col='black', size=2) +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        plot.margin=margin(-3,0,0,0),
        legend.text=element_text(size=11),
        legend.title=element_text(size=12,face="bold")) +
  xlab(bquote(bold('Age slope ('*beta*'-value change/year)'))) + ylab(expression(bold(paste(-log[10](P-value))))) +
  guides(colour=guide_legend(title="Methylation change"), shape=guide_legend(title="")) + 
  geom_abline(slope=0, intercept=-log10(thr), linetype=2, col='forestgreen', size=0.8);
d <-ggplot(data=aDMPs_without_CCC[aDMPs_without_CCC$Direction!='No change',], aes(x=beta, fill=Direction)) + 
  geom_density() + scale_fill_manual(values=c('Red', 'Blue')) +
  theme_classic() +
  theme(axis.line=element_blank(),
        axis.text.x=element_blank(), axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(), axis.title.y=element_blank());
comb_without_CCC <- plot_grid(plot_grid(d+theme(legend.position = "none"), vp_without_CCC+theme(legend.position = "none"), ncol = 1, align = "h",rel_heights=c(1,5)), 
                           plot_grid(ggplot(), get_legend(vp_without_CCC), ncol=1, rel_heights=c(1,5)), rel_widths=c(3,1));
ggsave('plots/vp_aDMPs_without_CCC.png', height=7, width=7);


## Create barplots with the overall picture. 

bp_with_CCC_df <- data.frame(N=c(as.numeric(table(aDMPs_with_CCC$Direction)['No change']),
                                 as.numeric(table(aDMPs_with_CCC$Direction)['Hypermethylated']),
                                 as.numeric(table(aDMPs_with_CCC$Direction)['Hypomethylated'])),
                             perc=c(paste0(round(as.numeric(table(aDMPs_with_CCC$Direction)['No change'])*100/nrow(aDMPs_with_CCC), digits=2), '%'),
                                    paste0(round(as.numeric(table(aDMPs_with_CCC$Direction)['Hypermethylated'])*100/nrow(aDMPs_with_CCC), digits=2), '%'),
                                    paste0(round(as.numeric(table(aDMPs_with_CCC$Direction)['Hypomethylated'])*100/nrow(aDMPs_with_CCC), digits=2), '%')),
                             Direction=c('No change', 'Hypermethlated', 'Hypomethylated'));
bp_with_CCC <- ggplot(data=bp_with_CCC_df, aes(x=Direction, y=N, fill=Direction)) +
  geom_bar(stat="identity") + scale_fill_manual(values=c('Red', 'Blue', 'Grey')) +
  geom_text(data=bp_with_CCC_df,aes(x=Direction, y=N,label=perc),vjust=-0.5, size=4) + 
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.text=element_text(size=11),
        legend.title=element_text(size=12,face="bold"),
        axis.text.x=element_blank()) +
  xlab("") + ylab("Number of aDMPs") +
  guides(fill=guide_legend(title="Methylation change")) +
  ylim(c(0, 303000));
ggsave('plots/barplots_aDMPs_with_CCC.pdf', height=5, width=5);


bp_without_CCC_df <- data.frame(N=c(as.numeric(table(aDMPs_without_CCC$Direction)['No change']),
                                 as.numeric(table(aDMPs_without_CCC$Direction)['Hypermethylated']),
                                 as.numeric(table(aDMPs_without_CCC$Direction)['Hypomethylated'])),
                             perc=c(paste0(round(as.numeric(table(aDMPs_without_CCC$Direction)['No change'])*100/nrow(aDMPs_without_CCC), digits=2), '%'),
                                    paste0(round(as.numeric(table(aDMPs_without_CCC$Direction)['Hypermethylated'])*100/nrow(aDMPs_without_CCC), digits=2), '%'),
                                    paste0(round(as.numeric(table(aDMPs_without_CCC$Direction)['Hypomethylated'])*100/nrow(aDMPs_without_CCC), digits=2), '%')),
                             Direction=c('No change', 'Hypermethlated', 'Hypomethylated'));
bp_without_CCC <- ggplot(data=bp_without_CCC_df, aes(x=Direction, y=N, fill=Direction)) +
  geom_bar(stat="identity") + scale_fill_manual(values=c('Red', 'Blue', 'Grey')) +
  geom_text(data=bp_without_CCC_df,aes(x=Direction, y=N,label=perc),vjust=-0.5, size=4) + 
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.text=element_text(size=11),
        legend.title=element_text(size=12,face="bold"),
        axis.text.x=element_blank()) +
  xlab("") + ylab("Number of aDMPs") +
  guides(fill=guide_legend(title="Methylation change")) +
  ylim(c(0, 303000));
ggsave('plots/barplots_aDMPs_without_CCC.pdf', height=5, width=5);

## Plot top aDMPs (with CCC). ##

top100_aDMPs_with_CCC <- as.data.frame(fread('betas_for_aDMPs_with_CCC.csv'));
rownames(top100_aDMPs_with_CCC) <- top100_aDMPs_with_CCC$ProbeID;
top100_aDMPs_with_CCC <- as.data.frame(t(top100_aDMPs_with_CCC[,-1]));
top100_aDMPs_with_CCC$GEO_sample <- rownames(top100_aDMPs_with_CCC);
rownames(top100_aDMPs_with_CCC) <- NULL;
final_aDMPs <- merge(raw_controls, top100_aDMPs_with_CCC, by='GEO_sample');
cgs <- c(aDMPs_with_CCC$ProbeID[head(which(aDMPs_with_CCC$beta>0), 2)], # Top hypermethylated
         aDMPs_with_CCC$ProbeID[head(which(aDMPs_with_CCC$beta<0), 2)]) # Top hypomethylated
cg_plots <- list();

for(c in cgs){
  
  cg_df <- data.frame(Age=as.numeric(final_aDMPs$Age_years), Beta=as.numeric(final_aDMPs[,which(colnames(final_aDMPs)==c)]));
  cg_plots[[c]] <- ggplot(data=cg_df, aes(x=Age, y=Beta)) + geom_point(col='grey') +
    theme_classic() +
    theme(axis.text=element_text(size=12, angle=90),
          axis.title=element_text(size=14,face="bold"),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          plot.margin=margin(3,0,0,0),) +
    labs(title=c) +
    xlab('Chronological age (years)') + ylab(bquote(bold(beta*'-value'))) +
    geom_smooth(method='lm',formula=y~x, show.legend = F, col='gray40') +
    ylim(c(0, 1));
  
}

display_cg_plots <- plot_grid(cg_plots[[1]], cg_plots[[2]], cg_plots[[3]], cg_plots[[4]], 
                              ncol = 2, nrow=2, align = "h", scale = 0.9);
ggsave('plots/aDMPs_beta_values_with_age.pdf', height=8, width=8);


#### 4. Methylation Shannon entropy. ####

## Understanding Shannon entropy.

calculate_entropy_ind <- function(betas){
  betas[betas==0] <- 0.000000000000001;
  betas[betas==1] <- 0.999999999999999;
  return(-(betas*log2(betas) + (1-betas)*log2(1-betas)));
}

entropy_df <- data.frame(Beta=seq(0,1,by=0.001), entropy=calculate_entropy_ind(seq(0,1,by=0.001)));
entropy_plot <- ggplot(entropy_df, aes(x=Beta, y=entropy)) + geom_line(col='grey', size=1.5) +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) +
  xlab(bquote(bold(beta*'-value (of a CpG site)'))) + ylab("Shannon entropy (of a CpG site)");
ggsave('plots/entropy_explanation.pdf', height=5, width=5);

## Genome-wide Shannon entropy in the full lifespan control. 

raw_entropy <- as.data.frame(fread('entropy_results.csv'));
final_entropy_df <- merge(raw_controls, raw_entropy, by='GEO_sample');

entropy_control_scatterplot <- ggplot(data=final_entropy_df, aes(x=Age_years, y=entropy)) +
  geom_point(col='grey') + 
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) +
  xlab("Chronological age (years)") + ylab("Genome-wide Shannon entropy") +
  geom_smooth(method='lm',formula=y~x, show.legend = F, col='gray40') +
  ylim(c(0.30, 0.60));
ggsave("plots/entropy_control_scatterplot.pdf", height=6, width=6);

# Calculate correlation.

cor(final_entropy_df$Age_years, final_entropy_df$entropy, method="spearman"); # 0.1984988
cor.test(final_entropy_df$Age_years, final_entropy_df$entropy, method="spearman")$p.value; # 3.82813e-21
mean(final_entropy_df$entropy); # 0.3943528

## Batch effect. 

entropy_batch_scatterplot <- ggplot(data=final_entropy_df, aes(x=Age_years, y=entropy, col=Batch)) +
  geom_point() + scale_colour_stata() + 
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) +
  xlab("Chronological age (years)") + ylab("Genome-wide Shannon entropy") +
  guides(col=guide_legend(title="Batch")) +
  ylim(c(0.30, 0.60));
ggsave("plots/entropy_batch_scatterplot.pdf", height=6, width=6);

# Correlation when removing 'GSE41273' batch --> same conclusion.

cor(final_entropy_df$Age_years[final_entropy_df$Batch!='GSE41273'], 
final_entropy_df$entropy[final_entropy_df$Batch!='GSE41273'], method="spearman");
cor.test(final_entropy_df$Age_years[final_entropy_df$Batch!='GSE41273'], 
final_entropy_df$entropy[final_entropy_df$Batch!='GSE41273'], method="spearman")$p.value;

# Correlation when removing 'GSE59065' batch --> same conclusion.

cor(final_entropy_df$Age_years[final_entropy_df$Batch!='GSE59065'], 
final_entropy_df$entropy[final_entropy_df$Batch!='GSE59065'], method="spearman");
cor.test(final_entropy_df$Age_years[final_entropy_df$Batch!='GSE59065'], 
final_entropy_df$entropy[final_entropy_df$Batch!='GSE59065'], method="spearman")$p.value;

# Correlation when removing 'GSE97362' batch --> same conclusion.

cor(final_entropy_df$Age_years[final_entropy_df$Batch!='GSE97362'], 
    final_entropy_df$entropy[final_entropy_df$Batch!='GSE97362'], method="spearman");
cor.test(final_entropy_df$Age_years[final_entropy_df$Batch!='GSE97362'], 
         final_entropy_df$entropy[final_entropy_df$Batch!='GSE97362'], method="spearman")$p.value;


#### 5. Horvath's epigenetic clock results. ####

## Plot for the DNAmAge pipeline. 

F_transf <- function(c, a=20){
  if(c <= a){
    F_c <- log((c+1)/(a+1));
  }
  if(c > a){
    F_c <- (c-a) / (a+1);
  }
  return(F_c);
}

aget_df <- data.frame(Age=seq(0,100,by=0.1), TransAge=(sapply(seq(0,100,by=0.1), F_transf)));
aget_plot <- ggplot(aget_df, aes(x=Age, y=TransAge)) + geom_line(col='grey', size=1.5) +
  geom_segment(x=-10, y=0, xend=20, yend=0, linetype='dashed') +
  geom_segment(x=20, y=-10, xend=20, yend=0, linetype='dashed') + 
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) +
  xlab('Chronological age (years)') + ylab("Transformed chronological age");
ggsave('plots/age_transformation_explanation.pdf', height=5, width=5);

## Fit linear models for full lifespan model.

lm_formula_ext_int <- paste0('DNAmAge_noob~Age_years+Sex+', paste(colnames(raw_controls)[grep('PC',colnames(raw_controls))], collapse='+'));
lm_formula_int <- paste0('DNAmAge_noob~Age_years+Sex+Gran+CD4T+CD8T+B+Mono+NK+', paste(colnames(raw_controls)[grep('PC',colnames(raw_controls))], collapse='+'));

lm_control_ext_int <- lm(lm_formula_ext_int, data=raw_controls);
lm_control_int <- lm(lm_formula_int, data=raw_controls);
raw_controls$delta_ext_int <- lm_control_ext_int$residuals;
raw_controls$delta_int <- lm_control_int$residuals;
raw_controls$Age_group <- ifelse(raw_controls$Age_years <= 20, 'Young age', ifelse(
  raw_controls$Age_years <= 55, 'Middle age', 'Old age'));
cor(raw_controls$Age_years, raw_controls$DNAmAge_noob, method="pearson"); # PCC: 0.9671172
median(abs(raw_controls$delta_ext_int)); # MAE without CCC: 2.821098
median(abs(raw_controls$delta_int)); # MAE with CCC: 2.711711

## Plots for full lifespan model.

horvath_control_scatterplot <- ggplot(data=raw_controls, aes(x=Age_years, y=DNAmAge_noob)) +
  geom_point(col='grey') +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) +
  xlab("Chronological age (years)") + ylab("DNAmAge (years)") +
  xlim(c(-1,101)) + ylim(c(-1,101)) + geom_abline(slope=1, intercept=0, linetype=2) +
  geom_smooth(method="lm", formula=y~x, show.legend=F, alpha = 0.3, col='salmon4', fill='darkgrey') +
  labs(title =paste0("Full lifespan control: N = ", nrow(raw_controls)));
ggsave("plots/horvath_control_scatterplot_full_lifespan.pdf", height=5, width=5);

full_lifespan_df <- data.frame(Age_group=factor(rep(raw_controls$Age_group, 2), levels=c('Young age', 'Middle age', 'Old age')),
                               delta=c(raw_controls$delta_ext_int, raw_controls$delta_int),
                               Model=rep(c("Without CCC", "With CCC"), each=nrow(raw_controls)));
full_lifespan_bias <- ggplot(data=full_lifespan_df, aes(x=Age_group, y=delta, fill=Model)) +
  geom_boxplot() + scale_fill_manual(values=c('red', 'blue')) +
  theme_classic() +
  scale_y_continuous(breaks=seq(-40,40,10)) +
  theme(axis.text=element_text(size=12),
        axis.text.x=element_text(angle = -90, hjust = 0, size=12),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(hjust = 0.5)) +
  xlab("Age group") + ylab("Epigenetic age acceleration (years)") +
  guides(fill=guide_legend(title="EAA model")) +
  geom_abline(slope=0, intercept=0, linetype=2, size=0.25) +
  labs(title =paste0("Full lifespan control: N = ", nrow(raw_controls)));
ggsave("plots/horvath_control_bias_full_lifespan.pdf", height=5, width=5);

## Fit linear models for 0-55 years model.

raw_controls <- as.data.frame(fread('final_control_data.tsv'));
controls_reduced <- raw_controls[raw_controls$Age_years <= 55,]; # N=1128
lm_control_ext_int_r <- lm(lm_formula_ext_int, data=controls_reduced);
lm_control_int_r <- lm(lm_formula_int, data=controls_reduced);
controls_reduced$delta_ext_int <- lm_control_ext_int_r$residuals;
controls_reduced$delta_int <- lm_control_int_r$residuals;
controls_reduced$Age_group <- ifelse(controls_reduced$Age_years <= 20, 'Young age', 'Middle age');
cor(controls_reduced$Age_years, controls_reduced$DNAmAge_noob, method="pearson"); # PCC: 0.96427
median(abs(controls_reduced$delta_ext_int)); # MAE without CCC: 2.323729
median(abs(controls_reduced$delta_int)); # MAE with CCC: 2.274239

## Plots for 0-55 years model.

horvath_control_scatterplot_half <- ggplot(data=controls_reduced, aes(x=Age_years, y=DNAmAge_noob)) +
  geom_point(col='grey') +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) +
  xlab("Chronological age (years)") + ylab("DNAmAge (years)") +
  xlim(c(-1,101)) + ylim(c(-1,101)) + geom_abline(slope=1, intercept=0, linetype=2) +
  geom_smooth(method="lm", formula=y~x, show.legend=F, alpha = 0.3, col='seagreen4', fill='darkgrey') +
  labs(title =paste0("0-55 years control: N = ", nrow(controls_reduced)));
ggsave("plots/horvath_control_scatterplot_0_55.pdf", height=5, width=5);

half_lifespan_df <- data.frame(Age_group=factor(rep(controls_reduced$Age_group, 2), levels=c('Young age', 'Middle age')),
                               delta=c(controls_reduced$delta_ext_int, controls_reduced$delta_int),
                               Model=rep(c("Without CCC", "With CCC"), each=nrow(controls_reduced)));
half_lifespan_bias <- ggplot(data=half_lifespan_df, aes(x=Age_group, y=delta, fill=Model)) +
  geom_boxplot() + scale_fill_manual(values=c('red', 'blue')) +
  theme_classic() +
  scale_y_continuous(breaks=seq(-40,40,10)) +
  theme(axis.text=element_text(size=12),
        axis.text.x=element_text(angle = -90, hjust = 0, size=12),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(hjust = 0.5)) +
  xlab("Age group") + ylab("Epigenetic age acceleration (years)") +
  guides(fill=guide_legend(title="EAA model")) +
  geom_abline(slope=0, intercept=0, linetype=2, size=0.25) +
  labs(title =paste0("0-55 years control: N = ", nrow(controls_reduced)));
ggsave("plots/horvath_control_bias_0_55.pdf", height=5, width=5);


#### 6. Other epigenetic clock results. ####

## Hannum clock. 

hannum_betas <- as.data.frame(fread('betas_for_probes_hannum_clock.csv'));
na_hannum <- which(is.na(hannum_betas$ProbeID) | hannum_betas$ProbeID=='NA.1' | hannum_betas$ProbeID=='NA.2');
hannum_betas <- hannum_betas[-na_hannum,];
hannum_coeffs <- as.data.frame(fread('hannum_clock_coefs.tsv'));
hannum_coeffs <- hannum_coeffs[-na_hannum,];
if(!all(hannum_betas$ProbeID==hannum_coeffs$ProbeID)){stop('The probeIDs do not match')};

hannum_betas_f <- t(hannum_betas[,-1]);
colnames(hannum_betas_f) <- hannum_betas$ProbeID;
hannum_coeffs_f <- matrix(hannum_coeffs$Coefficient, nrow=nrow(hannum_coeffs), ncol=1);
rownames(hannum_coeffs_f) <- hannum_coeffs$ProbeID;
hannum_predictions <- as.data.frame(hannum_betas_f %*% hannum_coeffs_f);
hannum_predictions$GEO_sample <- rownames(hannum_predictions);
colnames(hannum_predictions)[1] <- 'hannum_age';
raw_controls_f <- merge(raw_controls, hannum_predictions, by='GEO_sample');

## Plot Hannum for full lifespan model.

hannum_control_scatterplot <- ggplot(data=raw_controls_f, aes(x=Age_years, y=hannum_age)) +
  geom_point(col='grey') +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) +
  xlab("Chronological age (years)") + ylab("HannumAge (years)") +
  xlim(c(-15,105)) + ylim(c(-15,105)) + geom_abline(slope=1, intercept=0, linetype=2) +
  geom_smooth(method="lm", formula=y~x, show.legend=F, alpha = 0.3, col='salmon4', fill='darkgrey') +
  labs(title =paste0("Full lifespan control: N = ", nrow(raw_controls)));
ggsave("plots/hannum_control_scatterplot_full_lifespan.pdf", height=5, width=5);

hannum_horvath_control_scatterplot <- ggplot(data=raw_controls_f, aes(x=DNAmAge_noob, y=hannum_age)) +
  geom_point(col='grey') +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) +
  xlab("DNAmAge (years)") + ylab("HannumAge (years)") +
  xlim(c(-15,105)) + ylim(c(-15,105)) + geom_abline(slope=1, intercept=0, linetype=2) +
  geom_smooth(method="lm", formula=y~x, show.legend=F, alpha = 0.3, col='salmon4', fill='darkgrey') +
  labs(title =paste0("Full lifespan control: N = ", nrow(raw_controls)));
ggsave("plots/horvath_vs_hannum_control_scatterplot_full_lifespan.pdf", height=5, width=5);

## Fit linear models to Hannum model for full lifespan model.

lm_formula_ext_int_h <- paste0('hannum_age~Age_years+Sex+', paste(colnames(raw_controls)[grep('PC',colnames(raw_controls))], collapse='+'));
lm_formula_int_h <- paste0('hannum_age~Age_years+Sex+Gran+CD4T+CD8T+B+Mono+NK+', paste(colnames(raw_controls)[grep('PC',colnames(raw_controls))], collapse='+'));

lm_control_ext_int_h <- lm(lm_formula_ext_int_h, data=raw_controls_f);
lm_control_int_h <- lm(lm_formula_int_h, data=raw_controls_f);
raw_controls_f$hannum_ext_int <- lm_control_ext_int_h$residuals;
raw_controls_f$hannum_int <- lm_control_int_h$residuals;
cor(raw_controls_f$Age_years, raw_controls_f$hannum_age, method="pearson"); # PCC: 0.975631
median(abs(raw_controls_f$hannum_ext_int)); # MAE without CCC: 2.94844
median(abs(raw_controls_f$hannum_int)); # MAE with CCC: 2.842229
cor(raw_controls_f$DNAmAge_noob, raw_controls_f$hannum_age, method="pearson"); # PCC: 0.9778428

## Plot EAA for Horvath vs Hannum

hannum_horvath_EAA_with_CCC_scatterplot <- ggplot(data=raw_controls_f, aes(x=delta_int, y=hannum_int)) +
  geom_point(col='grey') +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) +
  xlab("Horvath EAA with CCC (years)") + ylab("Hannum EAA with CCC (years)") +
  geom_abline(slope=1, intercept=0, linetype=2) +
  geom_smooth(method="lm", formula=y~x, show.legend=F, alpha = 0.3, col='salmon4', fill='darkgrey') +
  labs(title=paste0("Full lifespan control: N = ", nrow(raw_controls)),
       subtitle=paste0('PCC: ', round(cor(raw_controls_f$delta_int, raw_controls_f$hannum_int), digits=4),
                       ifelse(cor.test(raw_controls_f$delta_int, raw_controls_f$hannum_int)$p.value < 2.2e-16, 
                              '; p-value < 2.2e-16', paste0('; p-value = ', format(cor.test(raw_controls_f$delta_int, raw_controls_f$hannum_int)$p.value, digits=4)))));
ggsave("plots/horvath_vs_hannum_EAA_with_CCC_full_lifespan.pdf", height=5, width=5);

hannum_horvath_EAA_without_CCC_scatterplot <- ggplot(data=raw_controls_f, aes(x=delta_ext_int, y=hannum_ext_int)) +
  geom_point(col='grey') +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) +
  xlab("Horvath EAA without CCC (years)") + ylab("Hannum EAA without CCC (years)") +
  geom_abline(slope=1, intercept=0, linetype=2) +
  geom_smooth(method="lm", formula=y~x, show.legend=F, alpha = 0.3, col='salmon4', fill='darkgrey') +
  labs(title=paste0("Full lifespan control: N = ", nrow(raw_controls)),
       subtitle=paste0('PCC: ', round(cor(raw_controls_f$delta_ext_int, raw_controls_f$hannum_ext_int), digits=4),
                       ifelse(cor.test(raw_controls_f$delta_ext_int, raw_controls_f$hannum_ext_int)$p.value < 2.2e-16, 
                              '; p-value < 2.2e-16', paste0('; p-value = ', format(cor.test(raw_controls_f$delta_ext_int, raw_controls_f$hannum_ext_int)$p.value, digits=4)))));
ggsave("plots/horvath_vs_hannum_EAA_without_CCC_full_lifespan.pdf", height=5, width=5);


## Mitotic clock.

mitotic_raw <- as.data.frame(fread('pcgtAge_results.csv'));
raw_controls_m <- merge(raw_controls, mitotic_raw, by='GEO_sample');

## Plot mitotic clock for full lifespan model.

pcgtAge_control_scatterplot <- ggplot(data=raw_controls_m, aes(x=Age_years, y=pcgtAge)) +
  geom_point(col='grey') +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) +
  xlab("Chronological age (years)") + ylab("pcgtAge") +
  xlim(c(0,105)) +
  geom_smooth(method="lm", formula=y~x, show.legend=F, alpha = 0.3, col='salmon4', fill='darkgrey') +
  labs(title =paste0("Full lifespan control: N = ", nrow(raw_controls)));
ggsave("plots/pcgtAge_control_scatterplot_full_lifespan.pdf", height=5, width=5);

pcgtAge_horvath_control_scatterplot <- ggplot(data=raw_controls_m, aes(x=DNAmAge_noob, y=pcgtAge)) +
  geom_point(col='grey') +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) +
  xlab("DNAmAge (years)") + ylab("pcgtAge") +
  xlim(c(0,105)) + 
  geom_smooth(method="lm", formula=y~x, show.legend=F, alpha = 0.3, col='salmon4', fill='darkgrey') +
  labs(title =paste0("Full lifespan control: N = ", nrow(raw_controls)));
ggsave("plots/horvath_vs_pcgtAge_control_scatterplot_full_lifespan.pdf", height=5, width=5);

## Fit linear models to mitotic clock for full lifespan model.

lm_formula_ext_int_m <- paste0('pcgtAge~Age_years+Sex+', paste(colnames(raw_controls)[grep('PC',colnames(raw_controls))], collapse='+'));
lm_formula_int_m <- paste0('pcgtAge~Age_years+Sex+Gran+CD4T+CD8T+B+Mono+NK+', paste(colnames(raw_controls)[grep('PC',colnames(raw_controls))], collapse='+'));

lm_control_ext_int_m <- lm(lm_formula_ext_int_m, data=raw_controls_m);
lm_control_int_m <- lm(lm_formula_int_m, data=raw_controls_m);
raw_controls_m$pcgtAge_ext_int <- lm_control_ext_int_m$residuals;
raw_controls_m$pcgtAge_int <- lm_control_int_m$residuals;
cor(raw_controls_m$Age_years, raw_controls_m$pcgtAge, method="pearson"); # PCC: 0.5130522
cor(raw_controls_m$DNAmAge_noob, raw_controls_m$pcgtAge, method="pearson"); # PCC: 0.5601623

## Plot EAA for Horvath vs pcgtAge

pcgtAge_horvath_EAA_with_CCC_scatterplot <- ggplot(data=raw_controls_m, aes(x=delta_int, y=pcgtAge_int)) +
  geom_point(col='grey') +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) +
  xlab("Horvath EAA with CCC (years)") + ylab("pcgtAge EAA with CCC") +
  geom_smooth(method="lm", formula=y~x, show.legend=F, alpha = 0.3, col='salmon4', fill='darkgrey') +
  labs(title=paste0("Full lifespan control: N = ", nrow(raw_controls)),
       subtitle=paste0('PCC: ', round(cor(raw_controls_m$delta_int, raw_controls_m$pcgtAge_int), digits=4),
                       ifelse(cor.test(raw_controls_m$delta_int, raw_controls_m$pcgtAge_int)$p.value < 2.2e-16, 
                              '; p-value < 2.2e-16', paste0('; p-value = ', format(cor.test(raw_controls_m$delta_int, raw_controls_m$pcgtAge_int)$p.value, digits=4)))));
ggsave("plots/horvath_vs_pcgtAge_EAA_with_CCC_full_lifespan.pdf", height=5, width=5);

pcgtAge_horvath_EAA_without_CCC_scatterplot <- ggplot(data=raw_controls_m, aes(x=delta_ext_int, y=pcgtAge_ext_int)) +
  geom_point(col='grey') +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) +
  xlab("Horvath EAA without CCC (years)") + ylab("pcgtAge EAA without CCC") +
  geom_smooth(method="lm", formula=y~x, show.legend=F, alpha = 0.3, col='salmon4', fill='darkgrey') +
  labs(title=paste0("Full lifespan control: N = ", nrow(raw_controls)),
       subtitle=paste0('PCC: ', round(cor(raw_controls_m$delta_ext_int, raw_controls_m$pcgtAge_ext_int), digits=4),
                       ifelse(cor.test(raw_controls_m$delta_ext_int, raw_controls_m$pcgtAge_ext_int)$p.value < 2.2e-16, 
                              '; p-value < 2.2e-16', paste0('; p-value = ', format(cor.test(raw_controls_m$delta_ext_int, raw_controls_m$pcgtAge_ext_int)$p.value, digits=4)))));
ggsave("plots/horvath_vs_pcgtAge_EAA_without_CCC_full_lifespan.pdf", height=5, width=5);

#### End of the script. ####