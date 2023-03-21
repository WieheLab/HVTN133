library(stringr)
library(optparse)
library(gtools)

option_list = list(
  make_option(c("-h", "--heavy_RS"), type="character", default="Heavy/heavies.RecombinationSummaries.txt",
              help="Heavy Recombination Summaries txt file from Cloanalyst", metavar="character"),
  make_option(c("-k", "--kappa_RS"), type="character", default="Kappa/kappas.RecombinationSummaries.txt",
              help="Kappa Recombination Summaries txt from Cloanalyst", metavar="character"),
  make_option(c("-l", "--lambda_RS"), type="character", default="Lambda/lambdas.RecombinationSummaries.txt",
              help="Lambda CloneAssignments.txt file from Cloanalyst", metavar="character"),
  make_option(c("-o", "--output_file"), type = "character", default = "chain_counts.csv",
              help="Output file containig chain counts by barcode", metavar="character")
);

opt_parser = OptionParser(option_list=option_list,add_help_option=FALSE);
opt = parse_args(opt_parser);

heavy_rs <- read.csv(opt$heavy_RS, sep = '\t', header = T, stringsAsFactors = F)
kappa_rs <- read.csv(opt$kappa_RS, sep = '\t', header = T, stringsAsFactors = F)
lambda_rs <- read.csv(opt$lambda_RS, sep = '\t', header = T, stringsAsFactors = F)

heavy_rs$ID <- gsub('_contig_[0-9]{1,2}', '', heavy_rs$UID)
kappa_rs$ID <- gsub('_contig_[0-9]{1,2}', '', kappa_rs$UID)
lambda_rs$ID <- gsub('_contig_[0-9]{1,2}', '', lambda_rs$UID)

heavy_counts <- data.frame(table(heavy_rs$ID))
colnames(heavy_counts) <- c('ID', 'N_Heavy')
kappa_counts <- data.frame(table(kappa_rs$ID))
colnames(kappa_counts) <- c('ID', 'N_Kappa')
lambda_counts <- data.frame(table(lambda_rs$ID))
colnames(lambda_counts) <- c('ID', 'N_Lambda')

lights <- merge(kappa_counts, lambda_counts, by = 'ID', all = T)
all <- merge(heavy_counts, lights, by = 'ID', all = T)
all$N_Heavy <- ifelse(is.na(all$N_Heavy), 0, all$N_Heavy)
all$N_Kappa <- ifelse(is.na(all$N_Kappa), 0, all$N_Kappa)
all$N_Lambda <- ifelse(is.na(all$N_Lambda), 0, all$N_Lambda)
all$N_Light <- all$N_Kappa + all$N_Lambda

write.csv(all, opt$output_file, quote = F, row.names = F)
