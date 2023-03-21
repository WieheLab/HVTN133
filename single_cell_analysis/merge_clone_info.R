library(stringr)
library(dplyr)
library(optparse)
library(gtools)
library(Biostrings)

option_list = list(
  make_option(c("-h", "--heavy_clone_assignments"), type="character", default="Heavy/CloneAssignments.txt",
              help="Heavy CloneAssignments.txt file from Cloanalyst", metavar="character"),
  make_option(c("-H", "--heavy_clones"), type="character", default="Heavy/Clones.txt",
              help="Heavy Clones.txt from Cloanalyst", metavar="character"),
  make_option(c("--heavy_fxnl"), type="character", default="../intermediate_files/heavies.fxnl.ids.txt",
              help="List of fxnl heavy ids", metavar="character"),
  make_option(c('--heavy_rs'), type = 'character', default='../unfiltered/Heavy/heavies.RecombinationSummaries.txt',
            help='Heavy RS Cloanalyst output', metavar='character'),
  make_option(c("-k", "--kappa_clone_assignments"), type="character", default="Kappa/CloneAssignments.txt",
              help="Kappa CloneAssignments.txt from Cloanalyst", metavar="character"),
  make_option(c("-K", "--kappa_clones"), type="character", default="Kappa/Clones.txt",
              help="Kappa Clones.txt file from Cloanalyst", metavar="character"),
  make_option(c("--kappa_fxnl"), type="character", default="../intermediate_files/kappas.fxnl.ids.txt",
              help="List of fxnl kappa ids", metavar="character"),
  make_option(c('--kappa_rs'), type = 'character', default='../unfiltered/Kappa/kappas.RecombinationSummaries.txt',
            help='Kappa RS Cloanalyst output', metavar='character'),
  make_option(c("-l", "--lambda_clone_assignments"), type="character", default="Lambda/CloneAssignments.txt",
              help="Lambda CloneAssignments.txt file from Cloanalyst", metavar="character"),
  make_option(c("-L", "--lambda_clones"), type="character", default="Lambda/Clones.txt",
              help="Lambda Clones.txt file from Cloanalyst", metavar="character"),
  make_option(c("--lambda_fxnl"), type="character", default="../intermediate_files/lambdas.fxnl.ids.txt",
              help="List of fxnl lambda ids", metavar="character"),
  make_option(c('--lambda_rs'), type = 'character', default='../unfiltered/Lambda/lambdas.RecombinationSummaries.txt',
            help='Lambda RS Cloanalyst output', metavar='character'),
  make_option(c("-c", "--chain_counts_csv"), type = "character", default = "../intermediate_files/chain_counts.csv",
              help="CSV file with chain counts by cell", metavar="character"),
  make_option(c("-I", "--allcontig_annotations"), type="character", default="../cellranger_output/all_contig_annotations_umi_filtered.csv",
              help="All contig annotations csv file output from Cell Ranger (with or without UMI filtering)", metavar="character"),
  make_option(c("-o", "--output_file"), type = "character", default = "10x_merged_clones.csv",
              help="Merged clone info output file", metavar="character"),
  make_option(c("-O", "--output_file_fxnl_1to1"), type = "character", default = "10x_merged_clones.fxnl_1to1.csv",
            help="Merged clone info output file (Fxnl and 1:1 only)", metavar="character")
);

opt_parser = OptionParser(option_list=option_list,add_help_option=FALSE);
opt = parse_args(opt_parser);

chain_counts <- read.csv(opt$chain_counts_csv)
heavy_fxnl <- read.csv(opt$heavy_fxnl, header = F, stringsAsFactors = F)
kappa_fxnl <- read.csv(opt$kappa_fxnl, header = F, stringsAsFactors = F)
lambda_fxnl <- read.csv(opt$lambda_fxnl, header = F, stringsAsFactors = F)
annots <- read.csv(opt$all_contig_annotations, header = T, stringsAsFactors = F)
isos <- annots[,c('contig_id', 'c_gene')]
colnames(isos) <- c('ReadID_Heavy', 'Isotype')
umis <- annots[,c('contig_id', 'umis')]

### Heavy

heavy_clone_assigns <- read.csv(opt$heavy_clone_assignments, sep = '\t', header = T, stringsAsFactors = F)
heavy_clone_assigns$ID <- gsub('_contig_[0-9]{1,2}', '', heavy_clone_assigns$ReadID)
heavy_clone_assigns <- heavy_clone_assigns[,c('ID', 'ReadID', 'CloneID')]

heavy_clones <- read.csv(opt$heavy_clones, sep = '\t', header = T, stringsAsFactors = F)
heavy_clones <- heavy_clones[,c('CloneID', 'X.Members')]

heavy <- merge(heavy_clone_assigns, heavy_clones, by = 'CloneID')
heavy$Heavy_Fxnl <- heavy$ReadID %in% heavy_fxnl$V1

heavy_rs <- read.csv(opt$heavy_rs, sep = '\t', header = T, stringsAsFactors = F)
heavy_rs <- select(heavy_rs, UID:VGene, DGene, JGene, CDR3, CDR3Length, MuFreq)

heavy <- merge(heavy, heavy_rs, by.x = 'ReadID', by.y = 'UID', all = F)

### Kappa
kappa_clone_assigns <- read.csv(opt$kappa_clone_assignments, sep = '\t', header = T, stringsAsFactors = F)
kappa_clone_assigns$ID <- gsub('_contig_[0-9]{1,2}', '', kappa_clone_assigns$ReadID)
kappa_clone_assigns <- kappa_clone_assigns[,c('ID', 'ReadID', 'CloneID')]

kappa_clones <- read.csv(opt$kappa_clones, sep = '\t', header = T, stringsAsFactors = F)
kappa_clones <- kappa_clones[,c('CloneID', 'X.Members')]

kappa <- merge(kappa_clone_assigns, kappa_clones, by = 'CloneID')
kappa$Chain <- 'Kappa'
kappa$Light_Fxnl <- kappa$ReadID %in% kappa_fxnl$V1

kappa_rs <- read.csv(opt$kappa_rs, sep = '\t', header = T, stringsAsFactors = F)
kappa_rs <- select(kappa_rs, UID:VGene, JGene, CDR3, CDR3Length, MuFreq)

kappa <- merge(kappa, kappa_rs, by.x = 'ReadID', by.y = 'UID', all = F)

### Lambda
lambda_clone_assigns <- read.csv(opt$lambda_clone_assignments, sep = '\t', header = T, stringsAsFactors = F)
lambda_clone_assigns$ID <- gsub('_contig_[0-9]{1,2}', '', lambda_clone_assigns$ReadID)
lambda_clone_assigns <- lambda_clone_assigns[,c('ID', 'ReadID', 'CloneID')]

lambda_clones <- read.csv(opt$lambda_clones, sep = '\t', header = T, stringsAsFactors = F)
lambda_clones <- lambda_clones[,c('CloneID', 'X.Members')]

lambda <- merge(lambda_clone_assigns, lambda_clones, by = 'CloneID')
lambda$Chain <- 'Lambda'
lambda$Light_Fxnl <- lambda$ReadID %in% lambda_fxnl$V1

lambda_rs <- read.csv(opt$lambda_rs, sep = '\t', header = T, stringsAsFactors = F)
lambda_rs <- select(lambda_rs, UID:VGene, JGene, CDR3, CDR3Length, MuFreq)

lambda <- merge(lambda, lambda_rs, by.x = 'ReadID', by.y = 'UID', all = F)

### Lights 
lights <- smartbind(kappa, lambda)

all_merged <- merge(heavy, lights, by = 'ID', suffixes = c('_Heavy', '_Light'), all = T)
all_merged <- merge(all_merged, chain_counts, by = 'ID', all.x = T)
all_merged <- merge(all_merged, isos, by = 'ReadID_Heavy', all.x = T)
all_merged <- merge(all_merged, umis, by.x = 'ReadID_Heavy', by.y = 'contig_id', all.x = T)
all_merged <- merge(all_merged, umis, by.x = 'ReadID_Light', by.y = 'contig_id', all.x = T, suffixes = c('_Heavy', '_Light'))
all_merged <- all_merged[order(all_merged$CloneID_Heavy, all_merged$CloneID_Light),]

no_na <- subset(all_merged, !is.na(CDR3_Heavy))
na <- subset(all_merged, is.na(CDR3_Heavy))
no_na$CDR3_Heavy_AA <- Biostrings::translate(DNAStringSet(no_na$CDR3_Heavy))
na$CDR3_Heavy_AA <- NA
all_merged <- smartbind(no_na, na)

no_na <- subset(all_merged, !is.na(CDR3_Light))
na <- subset(all_merged, is.na(CDR3_Light))
no_na$CDR3_Light_AA <- Biostrings::translate(DNAStringSet(no_na$CDR3_Light))
na$CDR3_Light_AA <- NA
all_merged <- smartbind(no_na, na)

all_merged$VGene_Heavy <- gsub('\\*[0-9]{1,}', '', all_merged$VGene_Heavy)
all_merged$VGene_Heavy <- gsub('IGHV', 'VH', all_merged$VGene_Heavy)
all_merged$DGene <- gsub('\\*[0-9]{1,}', '', all_merged$DGene)
all_merged$DGene <- gsub('IGHD', 'DH', all_merged$DGene)
all_merged$JGene_Heavy <- gsub('\\*[0-9]{1,}', '', all_merged$JGene_Heavy)
all_merged$JGene_Heavy <- gsub('IGHJ', 'JH', all_merged$JGene_Heavy)

all_merged$VGene_Light <- gsub('\\*[0-9]{1,}', '', all_merged$VGene_Light)
all_merged$VGene_Light <- gsub('IGLV', 'VL', all_merged$VGene_Light)
all_merged$VGene_Light <- gsub('IGKV', 'VK', all_merged$VGene_Light)
all_merged$JGene_Light <- gsub('\\*[0-9]{1,}', '', all_merged$JGene_Light)
all_merged$JGene_Light <- gsub('IGLJ', 'JL', all_merged$JGene_Light)
all_merged$JGene_Light <- gsub('IGKJ', 'JK', all_merged$JGene_Light)
all_merged$Fxnl_1to1 <- ifelse(all_merged$N_Heavy == 1 & all_merged$N_Light == 1 & all_merged$Heavy_Fxnl == 1 & all_merged$Light_Fxnl == 1, T, F)

all_merged <- select(all_merged, ID, Fxnl_1to1, CloneID_Heavy, X.Members_Heavy, ReadID_Heavy, VGene_Heavy, DGene, JGene_Heavy, CDR3Length_Heavy, Isotype, MuFreq_Heavy,
                     CloneID_Light, X.Members_Light, ReadID_Light, Chain, VGene_Light, JGene_Light, CDR3Length_Light, MuFreq_Light,
                     Heavy_Fxnl, Light_Fxnl, N_Heavy, N_Kappa, N_Lambda, N_Light, umis_Heavy, umis_Light, CDR3_Heavy, CDR3_Heavy_AA, CDR3_Light, CDR3_Light_AA)

all_merged_subset <- subset(all_merged, N_Heavy == 1 & N_Light == 1 & Heavy_Fxnl == 1 & Light_Fxnl == 1)

write.csv(all_merged, file = opt$output_file, row.names = F, quote = F)
write.csv(all_merged_subset, file = opt$output_file_fxnl_1to1, row.names = F, quote = F)

