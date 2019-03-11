###########################################################################################
#########                                                                         #########
#########                     Daniel Elias Martin Herranz                         #########
#########                             07/03/2019                                  #########
#########                              EMBL-EBI                                   #########
#########                           Thornton group                                #########
#########                                                                         #########
###########################################################################################

###########################################################################################
#####              Biological insights into the epigenetic ageing clock           #########
###########################################################################################
##### Extract beta-values from the controls for a given set of array probes.          #####
###########################################################################################
##### USAGE: manual                                                                   #####
###########################################################################################

###########################################################
##################### Dependencies ########################
###########################################################

library(data.table);


################################################################
################## Arguments ###################################
################################################################

args <- commandArgs(trailingOnly=TRUE);

input_path <- "/nfs/nobackup/thornton/dem44/methylation_clock/polycomb_hypothesis/epigenetic_syndromes/blood_control/beta_values/"; # Path where all the beta-values matrices are located
output_path <- "/nfs/research1/thornton/dem44/methylation_clock/raw_data_human/epigenetic_syndromes/aDMPs_full_lifespan/";
metadata_path <- "/nfs/research1/thornton/dem44/methylation_clock/raw_data_human/epigenetic_syndromes/aDMPs_full_lifespan/final_control_data.tsv";
probes_path <- as.character(args[1]); # Text file with cg ProbeIDs to be selected (one per line).
ccc_s <- strsplit(strsplit(probes_path, 'top100_aDMPs_')[[1]][2], '_')[[1]][1];

################################################################
################## Run the pipeline ############################
################################################################

#### 1. Read all the input data. ####

final_metadata <- as.data.frame(fread(metadata_path));
print(paste0('There are ', nrow(final_metadata), ' control samples in the original metadata.'));

beta_files <- list.files(input_path, full.names = T);

for(i in 1:length(beta_files)){
  
  batch <- strsplit(strsplit(beta_files[i], '/')[[1]][length(strsplit(beta_files[i], '/')[[1]])], '_')[[1]][1];
  print(paste0('Processing batch ', batch, ' ...'));
  temp_betas <- as.data.frame(fread(beta_files[i]));
  rownames(temp_betas) <- temp_betas$ProbeID;
  cgs <- rownames(temp_betas);
  temp_betas <- temp_betas[,-1];
  if(!(batch %in% c('GSE40279', 'GSE41273', 'Europe'))){colnames(temp_betas) <- sapply(strsplit(colnames(temp_betas), '_'), function(x){x[1]})};
  temp_betas <- temp_betas[,which(colnames(temp_betas) %in% final_metadata$GEO_sample)]; # Select only those samples that are in metadata (and therefore passed QC)
  
  if(i==1){
    
    final_betas <- temp_betas;
    
  }else{
    
    if(!all(cgs==rownames(final_betas))){stop('The probe IDs in this matrix is not in the same order as the previous ones.')};
    final_betas <- cbind(final_betas, temp_betas);
  }
}

print(paste0('There are ', ncol(final_betas), ' control samples in the final beta-values matrix.'));

# Edit metadata so it matches the samples we have. 

final_metadata <- final_metadata[match(colnames(final_betas), final_metadata$GEO_sample),];


#### 2. Extract selected probes. ####

print('Extracting probes ...');
probes_raw <- readLines(probes_path);
selected_betas <- final_betas[rownames(final_betas) %in% probes_raw,];
selected_betas <- selected_betas[match(probes_raw, rownames(selected_betas)),];
selected_betas <- cbind(rownames(selected_betas), selected_betas);
rownames(selected_betas) <- NULL;
colnames(selected_betas)[1] <- 'ProbeID';

write.table(selected_betas, file=paste0(output_path, '/betas_for_aDMPs_', ccc_s, '_CCC.csv'),
            sep=',', row.names=F, quote=F);

print('The script finished correctly.');

################################################################
################## End of the script ###########################
################################################################
