###########################################################################################
#########                                                                         #########
#########                     Daniel Elias Martin Herranz                         #########
#########                             22/02/2019                                  #########
#########                              EMBL-EBI                                   #########
#########                           Thornton group                                #########
#########                                                                         #########
###########################################################################################

###########################################################################################
#####                                 PhD Thesis                                       ####
###########################################################################################
##### Create plots for the different steps in the pre-processing DNA methylation data  #### 
##### pipeline.                                                                        ####
###########################################################################################
##### USAGE: manual.                                                                   ####
###########################################################################################

###########################################################
##################### Dependencies ########################
###########################################################

library(minfi);
library(RPMM);
library(IlluminaHumanMethylation27kmanifest);
library(ggplot2);
library(tidyr);


###########################################################
#####################  Arguments ##########################
###########################################################

print('Getting the input arguments ...');

## Input, output and annotation paths

path_to_raw_idat <- "/Users/dem44/Documents/Thesis/PhD-Thesis/additional_code/Chapter2/GSE41273/raw_idat/";
setwd('/Users/dem44/Documents/Thesis/PhD-Thesis/additional_code/Chapter2/');
output_path <- "/Users/dem44/Documents/Thesis/PhD-Thesis/additional_code/Chapter2/";
ann_path <- '/Users/dem44/Desktop/methylation_clock/polycomb_hypothesis/epigenetic_ageing_clock/utils/'

folder_name <- strsplit(path_to_raw_idat, '/')[[1]][length(strsplit(path_to_raw_idat, '/')[[1]])-1];


################################################################
################## Functions ###################################
################################################################

#### Function: reorganise paths of IDAT files so they can be used by the minfi package.

# path_to_raw: path to the input folder.

move_idat_files <- function(path_to_raw){
  
  all_idat_files_paths <- list.files(path_to_raw, recursive = TRUE, full.names = TRUE);
  
  for(f in all_idat_files_paths){
    
    file_name <- strsplit(f, '/')[[1]][length(strsplit(f, '/')[[1]])];
    slide_f <- strsplit(file_name, '_')[[1]][2];
    array_f <- strsplit(file_name, '_')[[1]][3];
    dir.create(file.path(path_to_raw, slide_f), showWarnings = FALSE);
    
    new_path <- file.path(path_to_raw,slide_f,file_name);
    file.rename(from=f, to=new_path);
    
  }
}


##### Function: BMIQ normalisation function (as implemented in the 'wateRmelon' package).

BMIQ <- function (beta.v, design.v, nL = 3, doH = TRUE, nfit = 50000, 
          th1.v = c(0.2, 0.75), th2.v = NULL, niter = 5, tol = 0.001, 
          plots = TRUE, sampleID = 1, pri = TRUE) 
{
  if (!library(RPMM, logical.return = TRUE, quietly = TRUE)) {
    stop("need RPMM package")
  }
  good <- !is.na(beta.v)
  out <- beta.v
  beta.v <- beta.v[good]
  design.v <- design.v[good]
  print <- function(x) {
    if (pri) 
      base::print(x)
  }
  type1.idx <- which(design.v == 1)
  type2.idx <- which(design.v == 2)
  beta1.v <- beta.v[type1.idx]
  beta2.v <- beta.v[type2.idx]
  if (min(beta1.v) == 0) {
    beta1.v[beta1.v == 0] <- min(setdiff(beta1.v, 0))
  }
  if (min(beta2.v) == 0) {
    beta2.v[beta2.v == 0] <- min(setdiff(beta2.v, 0))
  }
  if (max(beta1.v) == 1) {
    beta1.v[beta1.v == 1] <- max(setdiff(beta1.v, 1))
  }
  if (max(beta2.v) == 1) {
    beta2.v[beta2.v == 1] <- max(setdiff(beta2.v, 1))
  }
  w0.m <- matrix(0, nrow = length(beta1.v), ncol = nL)
  w0.m[which(beta1.v <= th1.v[1]), 1] <- 1
  w0.m[intersect(which(beta1.v > th1.v[1]), which(beta1.v <= 
                                                    th1.v[2])), 2] <- 1
  w0.m[which(beta1.v > th1.v[2]), 3] <- 1
  print("Fitting EM beta mixture to type1 probes")
  rand.idx <- sample(1:length(beta1.v), nfit, replace = FALSE)
  em1.o <- blc(matrix(beta1.v[rand.idx], ncol = 1), w = w0.m[rand.idx, 
                                                             ], maxiter = niter, tol = tol, verbose = pri)
  subsetclass1.v <- apply(em1.o$w, 1, which.max)
  subsetth1.v <- c(mean(max(beta1.v[rand.idx[subsetclass1.v == 
                                               1]]), min(beta1.v[rand.idx[subsetclass1.v == 2]])), mean(max(beta1.v[rand.idx[subsetclass1.v == 
                                                                                                                               2]]), min(beta1.v[rand.idx[subsetclass1.v == 3]])))
  class1.v <- rep(2, length(beta1.v))
  class1.v[which(beta1.v < subsetth1.v[1])] <- 1
  class1.v[which(beta1.v > subsetth1.v[2])] <- 3
  nth1.v <- subsetth1.v
  print("Done")
  if (plots) {
    print("Check")
    tmpL.v <- as.vector(rmultinom(1:nL, length(beta1.v), 
                                  prob = em1.o$eta))
    tmpB.v <- vector()
    for (l in 1:nL) {
      tmpB.v <- c(tmpB.v, rbeta(tmpL.v[l], em1.o$a[l, 1], 
                                em1.o$b[l, 1]))
    }
    pdf(paste("Type1fit-", sampleID, ".pdf", sep = ""), width = 6, 
        height = 4)
    plot(density(beta1.v))
    d.o <- density(tmpB.v)
    points(d.o$x, d.o$y, col = "green", type = "l")
    legend(x = 0.5, y = 3, legend = c("obs", "fit"), fill = c("black", 
                                                              "green"), bty = "n")
    dev.off()
  }
  d1U.o <- density(beta1.v[class1.v == 1])
  d1M.o <- density(beta1.v[class1.v == 3])
  mod1U <- d1U.o$x[which.max(d1U.o$y)]
  mod1M <- d1M.o$x[which.max(d1M.o$y)]
  d2U.o <- density(beta2.v[which(beta2.v < 0.4)])
  d2M.o <- density(beta2.v[which(beta2.v > 0.6)])
  mod2U <- d2U.o$x[which.max(d2U.o$y)]
  mod2M <- d2M.o$x[which.max(d2M.o$y)]
  th2.v <- vector()
  th2.v[1] <- nth1.v[1] + (mod2U - mod1U)
  th2.v[2] <- nth1.v[2] + (mod2M - mod1M)
  w0.m <- matrix(0, nrow = length(beta2.v), ncol = nL)
  w0.m[which(beta2.v <= th2.v[1]), 1] <- 1
  w0.m[intersect(which(beta2.v > th2.v[1]), which(beta2.v <= 
                                                    th2.v[2])), 2] <- 1
  w0.m[which(beta2.v > th2.v[2]), 3] <- 1
  print("Fitting EM beta mixture to type2 probes")
  rand.idx <- sample(1:length(beta1.v), nfit, replace = FALSE)
  em2.o <- blc(matrix(beta2.v[rand.idx], ncol = 1), w = w0.m[rand.idx, 
                                                             ], maxiter = niter, tol = tol, verbose = pri)
  print("Done")
  subsetclass2.v <- apply(em2.o$w, 1, which.max)
  subsetth2.v <- c(mean(max(beta2.v[rand.idx[subsetclass2.v == 
                                               1]]), min(beta2.v[rand.idx[subsetclass2.v == 2]])), mean(max(beta2.v[rand.idx[subsetclass2.v == 
                                                                                                                               2]]), min(beta2.v[rand.idx[subsetclass2.v == 3]])))
  class2.v <- rep(2, length(beta2.v))
  class2.v[which(beta2.v < subsetth2.v[1])] <- 1
  class2.v[which(beta2.v > subsetth2.v[2])] <- 3
  if (plots) {
    tmpL.v <- as.vector(rmultinom(1:nL, length(beta2.v), 
                                  prob = em2.o$eta))
    tmpB.v <- vector()
    for (lt in 1:nL) {
      tmpB.v <- c(tmpB.v, rbeta(tmpL.v[lt], em2.o$a[lt, 
                                                    1], em2.o$b[lt, 1]))
    }
    pdf(paste("Type2fit-", sampleID, ".pdf", sep = ""), width = 6, 
        height = 4)
    plot(density(beta2.v))
    d.o <- density(tmpB.v)
    points(d.o$x, d.o$y, col = "green", type = "l")
    legend(x = 0.5, y = 3, legend = c("obs", "fit"), fill = c("black", 
                                                              "green"), bty = "n")
    dev.off()
  }
  classAV1.v <- vector()
  classAV2.v <- vector()
  for (l in 1:nL) {
    classAV1.v[l] <- em1.o$mu[l, 1]
    classAV2.v[l] <- em2.o$mu[l, 1]
  }
  print("Start normalising type 2 probes")
  nbeta2.v <- beta2.v
  lt <- 1
  selU.idx <- which(class2.v == lt)
  selUR.idx <- selU.idx[which(beta2.v[selU.idx] > classAV2.v[lt])]
  selUL.idx <- selU.idx[which(beta2.v[selU.idx] < classAV2.v[lt])]
  p.v <- pbeta(beta2.v[selUR.idx], em2.o$a[lt, 1], em2.o$b[lt, 
                                                           1], lower.tail = FALSE)
  q.v <- qbeta(p.v, em1.o$a[lt, 1], em1.o$b[lt, 1], lower.tail = FALSE)
  nbeta2.v[selUR.idx] <- q.v
  p.v <- pbeta(beta2.v[selUL.idx], em2.o$a[lt, 1], em2.o$b[lt, 
                                                           1], lower.tail = TRUE)
  q.v <- qbeta(p.v, em1.o$a[lt, 1], em1.o$b[lt, 1], lower.tail = TRUE)
  nbeta2.v[selUL.idx] <- q.v
  lt <- 3
  selM.idx <- which(class2.v == lt)
  selMR.idx <- selM.idx[which(beta2.v[selM.idx] > classAV2.v[lt])]
  selML.idx <- selM.idx[which(beta2.v[selM.idx] < classAV2.v[lt])]
  p.v <- pbeta(beta2.v[selMR.idx], em2.o$a[lt, 1], em2.o$b[lt, 
                                                           1], lower.tail = FALSE)
  q.v <- qbeta(p.v, em1.o$a[lt, 1], em1.o$b[lt, 1], lower.tail = FALSE)
  nbeta2.v[selMR.idx] <- q.v
  if (doH) {
    lt <- 2
    selH.idx <- c(which(class2.v == lt), selML.idx)
    minH <- min(beta2.v[selH.idx])
    maxH <- max(beta2.v[selH.idx])
    deltaH <- maxH - minH
    deltaUH <- -max(beta2.v[selU.idx]) + min(beta2.v[selH.idx])
    deltaHM <- -max(beta2.v[selH.idx]) + min(beta2.v[selMR.idx])
    nmaxH <- min(nbeta2.v[selMR.idx]) - deltaHM
    nminH <- max(nbeta2.v[selU.idx]) + deltaUH
    ndeltaH <- nmaxH - nminH
    hf <- ndeltaH/deltaH
    nbeta2.v[selH.idx] <- nminH + hf * (beta2.v[selH.idx] - 
                                          minH)
  }
  pnbeta.v <- beta.v
  pnbeta.v[type1.idx] <- beta1.v
  pnbeta.v[type2.idx] <- nbeta2.v
  if (plots) {
    print("Generating final plot")
    d1.o <- density(beta1.v)
    d2.o <- density(beta2.v)
    d2n.o <- density(nbeta2.v)
    ymax <- max(d2.o$y, d1.o$y, d2n.o$y)
    pdf(paste("CheckBMIQ-", sampleID, ".pdf", sep = ""), 
        width = 6, height = 4)
    plot(density(beta2.v), type = "l", ylim = c(0, ymax), 
         xlim = c(0, 1))
    points(d1.o$x, d1.o$y, col = "red", type = "l")
    points(d2n.o$x, d2n.o$y, col = "blue", type = "l")
    legend(x = 0.5, y = ymax, legend = c("type1", "type2", 
                                         "type2-BMIQ"), bty = "n", fill = c("red", "black", 
                                                                            "blue"))
    dev.off()
  }
  print(paste("Finished for sample ", sampleID, sep = ""))
  out[good] <- pnbeta.v
  pnbeta.v <- out
  return(list(nbeta = pnbeta.v, class1 = class1.v, class2 = class2.v, 
              av1 = classAV1.v, av2 = classAV2.v, hf = hf, th1 = nth1.v, 
              th2 = th2.v))
}



############################################################
################## Running the pipeline ####################
############################################################

#### 1. Process the input dataset.

## Check if the folder has any IDAT files.

if(!file.exists(path_to_raw_idat)){stop(paste0('The project ', folder_name, ' has no IDAT files.'))};

## Move the IDAT files to a correct directory tree. ##

print('Rearranging IDAT files in a new directory tree ...');
move_idat_files(path_to_raw_idat);

## Classify the files into 27K / 450K. ##

print('Classifying the array platforms ...');
new_files_paths <- list.files(path_to_raw_idat, full.names = TRUE, recursive=TRUE);
array_class <- sapply(new_files_paths, function(f){
  ifelse(file.size(f) < 800000, "27K", ifelse(file.size(f) < 9000000, "450K", NA));
});

if(sum(is.na(array_class)) > 0){
  stop("The array platform could not be identified in some samples !!");
}

if(sum(array_class=="27K")){
  stop('There are 27K samples among the input IDAT files. This version of the cell-type composition estimation can only handle 450K data.');
}

## Create sample annotation files. ##

print('Creating MINFI annotation for samples ...');
sample_ann <- unique(data.frame(
  Array=as.character(sapply(names(which(array_class == "450K")), function(x){
    a <- strsplit(strsplit(x, '/')[[1]][length(strsplit(names(which(array_class == "450K")), '/')[[1]])], '_')[[1]][3];
    return(a);
  })),
  Slide=as.character(sapply(names(which(array_class == "450K")), function(x){
    s <- strsplit(x, '/')[[1]][length(strsplit(names(which(array_class == "450K")), '/')[[1]])-1];
    return(s);
  })),
  Basename=gsub('_Red.idat', '', gsub('_Grn.idat', '', names(which(array_class == "450K"))))));

print(paste0('The folder ', folder_name, ' contains ', length(list.files(path_to_raw_idat, recursive = TRUE, full.names = TRUE)),
             ' IDAT files.'));

## Get only the control samples. 

raw_control_metadata <- as.data.frame(fread('final_control_data.tsv'));
raw_control_metadata <- raw_control_metadata[raw_control_metadata$Batch=='GSE41273',];
raw_control_metadata_tempID <- paste0(raw_control_metadata$Slide_ID, '_', raw_control_metadata$Array_ID);
raw_control_metadata_tempID <- c(raw_control_metadata_tempID, '7497398049_R02C01', '7512560130_R02C02');
sample_ann <- sample_ann[which(paste0(sample_ann$Slide, '_', sample_ann$Array) %in% raw_control_metadata_tempID),];

## QC plot. 

input_rg_all <-read.metharray.exp(targets = sample_ann);
ms_all <- preprocessNoob(input_rg_all);
qc_res <- as.data.frame(getQC(ms_all));
qc_res$QC <- ifelse(as.numeric(apply(qc_res, 1, mean)) <10.5, 'Failed QC', 'Passed QC');

qc_plot <- ggplot(data=qc_res, aes(x=mMed, y=uMed, col=QC)) + 
  geom_point() + geom_abline(intercept = 21, slope = -1, linetype=2) + 
  theme_classic() + scale_colour_manual(values=c("red", "grey"), ) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = c(0.3, 0.85),
        legend.title = element_blank(),
        legend.box.background = element_rect(colour = "black", size=0.8)) +
  xlab(bquote(bold("median{"~log[2]~M[i]~"}"))) + ylab(bquote(bold("median{"~log[2]~U[i]~"}"))) +
  ylim(c(9,13)) + xlim(c(9,13));
ggsave("plots/QC_samples_failed.pdf", height=5, width=5);

## Beta-values distributions after background correction.

betas_all_bc <- as.data.frame(getBeta(ms_all, type="Illumina"));
betas_all_bc$ProbeID <- rownames(betas_all_bc);
betas_all_bc_gather <- betas_all_bc %>% gather(Sample, Betas, -ProbeID);
betas_all_bc_gather$QC <- ifelse(betas_all_bc_gather$Sample %in% c('sample15_7497398049_R02C01', 'sample76_7512560130_R02C02'),
                                  'Failed QC', 'Passed QC');
plot_betas_all_bc <- ggplot() + 
  geom_line(data=betas_all_bc_gather[betas_all_bc_gather$QC=='Passed QC',], aes(x=Betas, group=Sample, col=QC), stat="density") +
  geom_line(data=betas_all_bc_gather[betas_all_bc_gather$QC=='Failed QC',], aes(x=Betas, group=Sample, col=QC), stat="density") + 
  theme_classic() + scale_color_manual(labels = c("Failed QC", "Passed QC"), values = c("red", "grey")) + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.title = element_blank(),
        legend.box.background = element_rect(colour = "black", size=0.8)) +
  xlab(bquote(bold(beta*'-value'))) + ylab("Density") +
  xlim(c(-0.1,1.1)) + scale_x_continuous(breaks=c(0,0.5,1));
ggsave("plots/betas_distributions_all_bc.pdf", height=5, width=5);


## Beta-values distributions before any pre-processing.

ms_all_raw <- preprocessRaw(input_rg_all);
betas_all_raw <- as.data.frame(getBeta(ms_all_raw, type="Illumina"));
betas_all_raw$ProbeID <- rownames(betas_all_raw);
betas_all_raw_gather <- betas_all_raw %>% gather(Sample, Betas, -ProbeID);
betas_all_raw_gather$QC <- ifelse(betas_all_raw_gather$Sample %in% c('sample15_7497398049_R02C01', 'sample76_7512560130_R02C02'),
                                  'Failed QC', 'Passed QC');
plot_betas_all_raw <- ggplot() + 
  geom_line(data=betas_all_raw_gather[betas_all_raw_gather$QC=='Passed QC',], aes(x=Betas, group=Sample, col=QC), stat="density") +
  geom_line(data=betas_all_raw_gather[betas_all_raw_gather$QC=='Failed QC',], aes(x=Betas, group=Sample, col=QC), stat="density") + 
  theme_classic() + scale_color_manual(labels = c("Failed QC", "Passed QC"), values = c("red", "grey")) + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.title = element_blank(),
        legend.box.background = element_rect(colour = "black", size=0.8)) +
  xlab(bquote(bold(beta*'-value'))) + ylab("Density") +
  xlim(c(-0.1,1.1)) + scale_x_continuous(breaks=c(0,0.5,1));
ggsave("plots/betas_distributions_all_raw.pdf", height=5, width=5);


## Background correction, filtering and BMIQ. 

# Filter out the probes associated with SNPs.
# See https://www.bioconductor.org/help/course-materials/2014/BioC2014/minfi_BioC2014.pdf

input_ms <- dropLociWithSnps(mapToGenome(ms_all));

# Obtain beta-values.

input_betas <- as.data.frame(getBeta(input_ms, type="Illumina"));
print(paste0('After removing probes associated with SNPs, there are ', nrow(input_betas), ' probes left.'));

# Filter out cross-reactive probes.

path_cr <- paste0(ann_path, '/cross_reactive_probes_Chen_2013.txt');
cr_probes <- readLines(path_cr);
input_betas <- input_betas[!(rownames(input_betas) %in% cr_probes),];
print(paste0('After removing cross-reactive probes, there are ', nrow(input_betas),
             ' probes left.'));

# Filter out probes from sex chromosomes.

ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19);
sex_probes <- ann450k$Name[ann450k$chr %in% c("chrX","chrY")];
input_betas <- input_betas[!(rownames(input_betas) %in% sex_probes),];
print(paste0('After removing probes in sex chromosomes, there are ', nrow(input_betas), ' probes left.'));

# BMIQ.

print('Performing BMIQ normalisation ...');
info_450K_II <- getProbeInfo(IlluminaHumanMethylation450kmanifest, type = c("II"));
type_II <- which(rownames(input_betas) %in% info_450K_II$Name);
design_vector <- rep(1,nrow(input_betas));
design_vector[type_II] <- 2;
failed_BMIQ_samples <- c();
input_final_df <- as.data.frame(matrix(NA, ncol=ncol(input_betas), nrow=nrow(input_betas)));
rownames(input_final_df) <- rownames(input_betas);
colnames(input_final_df) <- colnames(input_betas);

for(i in 1:ncol(input_betas)){
  
  print(paste0('Processing sample ', colnames(input_final_df)[i], ' ...'));
  current_betas <- try(BMIQ(beta.v=input_betas[,i], design.v=design_vector, plots=F, pri=T));
  if(class(current_betas)=="try-error"){
    print(paste0('Sample ', colnames(input_final_df)[i], ' could not be BMIQ normalised.'));
    failed_BMIQ_samples <- c(colnames(input_final_df)[i], failed_BMIQ_samples);
  }else{
    input_final_df[,i] <- current_betas$nbeta;
  }
}

if(length(failed_BMIQ_samples)>0){ # Remove samples that failed BMIQ
  rm_indexes <- match(failed_BMIQ_samples, colnames(input_final_df))
  input_final_df <- input_final_df[,-rm_indexes];
  print(paste0('Number of samples that did not completed BMIQ: ', length(failed_BMIQ_samples)));
  print(paste0('Sample names: ', paste(failed_BMIQ_samples, collapse=',')));
}else{
  print('All the samples completed BMIQ.');
}

# Plot beta-values distributions after all the pre-processing.

input_final_df$ProbeID <- rownames(input_final_df);
#saveRDS(input_final_df, 'betas_after_all_preprocessing');
input_final_df_gather <- input_final_df[,-which(colnames(input_final_df) %in% c('sample15_7497398049_R02C01', 'sample76_7512560130_R02C02'))];
input_final_df_gather <- input_final_df_gather %>% gather(Sample, Betas, -ProbeID);
input_final_df_gather$QC <- ifelse(input_final_df_gather$Sample %in% c('sample15_7497398049_R02C01', 'sample76_7512560130_R02C02'),
                                  'Failed QC', 'Passed QC');

plot_betas_after_all_preprocessing <- ggplot() + 
  geom_line(data=input_final_df_gather[input_final_df_gather$QC=='Passed QC',], aes(x=Betas, group=Sample, col=QC), stat="density") +
  geom_line(data=input_final_df_gather[input_final_df_gather$QC=='Failed QC',], aes(x=Betas, group=Sample, col=QC), stat="density") + 
  theme_classic() + scale_color_manual(labels = c("Passed QC"), values = c("grey")) + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.title = element_blank(),
        legend.box.background = element_rect(colour = "black", size=0.8)) +
  xlab(bquote(bold(beta*'-value'))) + ylab("Density") +
  xlim(c(-0.1,1.1)) + scale_x_continuous(breaks=c(0,0.5,1));
ggsave("plots/betas_distributions_after_all_preprocessing.pdf", height=5, width=5);


## Supplementary plots: background correction.

raw_m_intensities <- as.data.frame(getMeth(ms_all_raw));
raw_m_intensities$ProbeID <- rownames(raw_m_intensities);
raw_m_intensities_gather <- raw_m_intensities %>% gather(Sample, Betas, -ProbeID);
raw_m_intensities_gather$QC <- ifelse(raw_m_intensities_gather$Sample %in% c('sample15_7497398049_R02C01', 'sample76_7512560130_R02C02'),
                                  'Failed QC', 'Passed QC');
plot_raw_m_intensities <- ggplot() + 
  geom_line(data=raw_m_intensities_gather[raw_m_intensities_gather$QC=='Passed QC',], aes(x=Betas, group=Sample, col=QC), stat="density") +
  geom_line(data=raw_m_intensities_gather[raw_m_intensities_gather$QC=='Failed QC',], aes(x=Betas, group=Sample, col=QC), stat="density") + 
  theme_classic() + scale_color_manual(labels = c("Failed QC", "Passed QC"), values = c("red", "grey")) + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = c(0.5, 0.5),
        legend.title = element_blank(),
        legend.box.background = element_rect(colour = "black", size=0.8)) +
  xlab(bquote(bold("Raw"~M[i]))) + ylab("Density") +
  xlim(c(0,30000));
ggsave("plots/raw_m_intensities_dist.pdf", height=5, width=5);

raw_u_intensities <- as.data.frame(getUnmeth(ms_all_raw));  
raw_u_intensities$ProbeID <- rownames(raw_u_intensities);
raw_u_intensities_gather <- raw_u_intensities %>% gather(Sample, Betas, -ProbeID);
raw_u_intensities_gather$QC <- ifelse(raw_u_intensities_gather$Sample %in% c('sample15_7497398049_R02C01', 'sample76_7512560130_R02C02'),
                                      'Failed QC', 'Passed QC');
plot_raw_u_intensities <- ggplot() + 
  geom_line(data=raw_u_intensities_gather[raw_u_intensities_gather$QC=='Passed QC',], aes(x=Betas, group=Sample, col=QC), stat="density") +
  geom_line(data=raw_u_intensities_gather[raw_u_intensities_gather$QC=='Failed QC',], aes(x=Betas, group=Sample, col=QC), stat="density") + 
  theme_classic() + scale_color_manual(labels = c("Failed QC", "Passed QC"), values = c("red", "grey")) + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = c(0.5, 0.5),
        legend.title = element_blank(),
        legend.box.background = element_rect(colour = "black", size=0.8)) +
  xlab(bquote(bold("Raw"~U[i]))) + ylab("Density") +
  xlim(c(0,30000));
ggsave("plots/raw_u_intensities_dist.pdf", height=5, width=5);

rm(raw_m_intensities,raw_m_intensities_gather,raw_u_intensities,raw_u_intensities_gather);

bc_m_intensities <- as.data.frame(getMeth(ms_all));
bc_m_intensities$ProbeID <- rownames(bc_m_intensities);
bc_m_intensities_gather <- bc_m_intensities %>% gather(Sample, Betas, -ProbeID);
bc_m_intensities_gather$QC <- ifelse(bc_m_intensities_gather$Sample %in% c('sample15_7497398049_R02C01', 'sample76_7512560130_R02C02'),
                                      'Failed QC', 'Passed QC');
plot_bc_m_intensities <- ggplot() + 
  geom_line(data=bc_m_intensities_gather[bc_m_intensities_gather$QC=='Passed QC',], aes(x=Betas, group=Sample, col=QC), stat="density") +
  geom_line(data=bc_m_intensities_gather[bc_m_intensities_gather$QC=='Failed QC',], aes(x=Betas, group=Sample, col=QC), stat="density") + 
  theme_classic() + scale_color_manual(labels = c("Failed QC", "Passed QC"), values = c("red", "grey")) + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = c(0.5, 0.5),
        legend.title = element_blank(),
        legend.box.background = element_rect(colour = "black", size=0.8)) +
  xlab(bquote(bold("Background-corrected "~M[i]))) + ylab("Density") +
  xlim(c(0,30000));
ggsave("plots/bc_m_intensities_dist.pdf", height=5, width=5);

bc_u_intensities <- as.data.frame(getUnmeth(ms_all));
bc_u_intensities$ProbeID <- rownames(bc_u_intensities);
bc_u_intensities_gather <- bc_u_intensities %>% gather(Sample, Betas, -ProbeID);
bc_u_intensities_gather$QC <- ifelse(bc_u_intensities_gather$Sample %in% c('sample15_7497398049_R02C01', 'sample76_7512560130_R02C02'),
                                     'Failed QC', 'Passed QC');
plot_bc_u_intensities <- ggplot() + 
  geom_line(data=bc_u_intensities_gather[bc_u_intensities_gather$QC=='Passed QC',], aes(x=Betas, group=Sample, col=QC), stat="density") +
  geom_line(data=bc_u_intensities_gather[bc_u_intensities_gather$QC=='Failed QC',], aes(x=Betas, group=Sample, col=QC), stat="density") + 
  theme_classic() + scale_color_manual(labels = c("Failed QC", "Passed QC"), values = c("red", "grey")) + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = c(0.5, 0.5),
        legend.title = element_blank(),
        legend.box.background = element_rect(colour = "black", size=0.8)) +
  xlab(bquote(bold("Background-corrected "~U[i]))) + ylab("Density") +
  xlim(c(0,30000));
ggsave("plots/bc_u_intensities_dist.pdf", height=5, width=5);

rm(bc_m_intensities,bc_m_intensities_gather,bc_u_intensities,bc_u_intensities_gather);


## Supplementary plots: BMIQ.

chosen_sample <- 'sample1_7497398005_R01C01';
info_450K_II <- getProbeInfo(IlluminaHumanMethylation450kmanifest, type = c("II"));
info_450K_I <- getProbeInfo(IlluminaHumanMethylation450kmanifest, type = c("I"));
BMIQ_df <- data.frame(ProbeID=input_final_df$ProbeID, Betas=input_final_df[,which(colnames(input_final_df)==chosen_sample)],
                      Probe_type=NA);
BMIQ_df$Probe_type[which(BMIQ_df$ProbeID %in% info_450K_I$Name)] <- 'Infinium I';
BMIQ_df$Probe_type[which(BMIQ_df$ProbeID %in% info_450K_II$Name)] <- 'Infinium II \nwith BMIQ';
betas_all_bc_sub <- betas_all_bc[match(BMIQ_df$ProbeID,betas_all_bc$ProbeID),];
ind <- which(betas_all_bc_sub$ProbeID %in% info_450K_II$Name);
typeiipre_df <- data.frame(ProbeID=betas_all_bc_sub$ProbeID[ind],
                           Betas=betas_all_bc_sub[ind,which(colnames(betas_all_bc_sub)==chosen_sample)],
                           Probe_type='Infinium II \nwithout BMIQ');
final_BMIQ_df <- rbind(BMIQ_df,typeiipre_df);

plot_bmiq <- ggplot() + 
  geom_line(data=final_BMIQ_df[final_BMIQ_df$Probe_type=='Infinium I',], aes(x=Betas, group=Probe_type, col=Probe_type), stat="density") +
  geom_line(data=final_BMIQ_df[final_BMIQ_df$Probe_type=='Infinium II \nwith BMIQ',], aes(x=Betas, group=Probe_type, col=Probe_type), stat="density") +
  geom_line(data=final_BMIQ_df[final_BMIQ_df$Probe_type=='Infinium II \nwithout BMIQ',], aes(x=Betas, group=Probe_type, col=Probe_type), stat="density") +
  theme_classic() + scale_color_manual(labels = c("Infinium I\n", "Infinium II \nwith BMIQ\n", "Infinium II \nwithout BMIQ"), values = c("cyan", "magenta", "deeppink4")) + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = c(0.5, 0.5),
        legend.title = element_blank(),
        legend.box.background = element_rect(colour = "black", size=0.8)) +
  xlab(bquote(bold(beta*'-value'))) + ylab("Density") +
  xlim(c(0,1));
ggsave("plots/BMIQ_probe_types.pdf", height=5, width=5);


## Plot with M-values after pre-processing. 

input_final_df <- readRDS('betas_after_all_preprocessing');
m_values <- log2(input_final_df[,-54]/(1-input_final_df[,-54]));
m_values$ProbeID <- rownames(m_values);
m_values_gather <- m_values[,-which(colnames(m_values) %in% c('sample15_7497398049_R02C01', 'sample76_7512560130_R02C02'))];
m_values_gather <- m_values_gather %>% gather(Sample, Betas, -ProbeID);
m_values_gather$QC <- ifelse(m_values_gather$Sample %in% c('sample15_7497398049_R02C01', 'sample76_7512560130_R02C02'),
                                   'Failed QC', 'Passed QC');

plot_m_values_after_all_preprocessing <- ggplot() + 
  geom_line(data=m_values_gather[m_values_gather$QC=='Passed QC',], aes(x=Betas, group=Sample, col=QC), stat="density") +
  geom_line(data=m_values_gather[m_values_gather$QC=='Failed QC',], aes(x=Betas, group=Sample, col=QC), stat="density") + 
  theme_classic() + scale_color_manual(labels = c("Passed QC"), values = c("grey")) + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.title = element_blank(),
        legend.box.background = element_rect(colour = "black", size=0.8)) +
  xlab('M-value') + ylab("Density");
ggsave("plots/m_values_distributions_after_all_preprocessing.pdf", height=5, width=5);

#########################################################
########## End of the script ############################
#########################################################
