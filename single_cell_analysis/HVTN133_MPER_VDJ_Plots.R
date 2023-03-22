.libPaths('/datacommons/dhvi/mb488/R_Libs/R_v4.2.2/')
library(dplyr)
library(tidyverse)
library(gtools)
library(stringr)
library(ggplot2)
library(ggpubr)

setwd('/datacommons/dhvi/10X_Genomics_data/HVTN133/Analysis/Manuscript_Figures_Data/')

# HVTN133
v1 <- read.csv('/datacommons/dhvi/10X_Genomics_data/HVTN133/V1/visit2/CR6/VDJ/analysis/clonal_analysis/10x_merged_clones.fxnl_1to1.csv')
v2 <- read.csv('/datacommons/dhvi/10X_Genomics_data/HVTN133/V2/visit2/CR6/VDJ/analysis/clonal_analysis/10x_merged_clones.fxnl_1to1.csv')
v3 <- read.csv('/datacommons/dhvi/10X_Genomics_data/HVTN133/V3/visit2/CR6/VDJ/analysis/clonal_analysis/10x_merged_clones.fxnl_1to1.csv')
v4 <- read.csv('/datacommons/dhvi/10X_Genomics_data/HVTN133/V4/visit2/CR6/VDJ/analysis/clonal_analysis/10x_merged_clones.fxnl_1to1.csv')
v6 <- read.csv('/datacommons/dhvi/10X_Genomics_data/HVTN133/V6/visit2/CR6/VDJ/analysis/clonal_analysis/10x_merged_clones.fxnl_1to1.csv')

v1$PTID <- '133-23'
v2$PTID <- '133-39'
v3$PTID <- '133-35'
v4$PTID <- '133-33'
v6$PTID <- '133-30'

all_list <- list(v1, v2, v3, v4, v6)

all <- smartbind(list = all_list)

all_vgene_counts <- all %>% group_by(PTID, VGene_Heavy) %>% summarize(n=n()) %>% mutate(freq=n/sum(n)) %>% data.frame()
pooled_vgene_counts <- all %>% group_by(VGene_Heavy) %>% summarize(n=n()) %>% mutate(freq=n/sum(n)) %>% data.frame()
pooled_vgene_counts$PTID <- 'Pooled'

all2 <- smartbind(all_vgene_counts, pooled_vgene_counts)

vh_table <- pivot_wider(all2, id_cols = 'VGene_Heavy', names_from = 'PTID', values_from = 'freq', values_fill = 0) %>% data.frame()
vh_table$gene <- gsub('VH', '', vh_table$VGene_Heavy)
vh_table <- cbind(vh_table, str_split_fixed(vh_table$gene, pattern = '-', n = Inf))
vh_table$char <- str_extract(vh_table$gene,'[A-Z]{1,}')
vh_table$`2` <- gsub('[A-Z]{1,}', '', vh_table$`2`)
vh_table <- vh_table %>% mutate_at(c('1', '2', '3'), as.numeric)
vh_table$char <- ifelse(is.na(vh_table$char), '', vh_table$char )
vh_table$`3` <- ifelse(is.na(vh_table$`3`), 0, vh_table$`3` )

vh_table <- vh_table[order(vh_table$`1`, vh_table$`2`, vh_table$`3`, vh_table$char, decreasing = c(F, F, F, F)),]
vh_table <- select(vh_table, VGene_Heavy:Pooled)

colnames(vh_table) <- gsub('X', '', colnames(vh_table))

# Kappa

all_vklgene_counts <- all %>% group_by(PTID, VGene_Light) %>% summarize(n=n()) %>% mutate(freq=n/sum(n)) %>% data.frame()
pooled_vklgene_counts <- all %>% group_by(VGene_Light) %>% summarize(n=n()) %>% mutate(freq=n/sum(n)) %>% data.frame()
pooled_vklgene_counts$PTID <- 'Pooled'

all2 <- smartbind(all_vklgene_counts, pooled_vklgene_counts)

vkl_table <- pivot_wider(all2, id_cols = 'VGene_Light', names_from = 'PTID', values_from = 'freq', values_fill = 0) %>% data.frame()
vkl_table$chain <- ifelse(grepl(vkl_table$VGene_Light, pattern = 'K'), 'K', 'L')
vkl_table$gene <- gsub('V[KL]', '', vkl_table$VGene_Light)
vkl_table <- cbind(vkl_table, str_split_fixed(vkl_table$gene, pattern = '-', n = Inf))
vkl_table$char <- str_extract(vkl_table$gene,'[A-Z]{1,}')
vkl_table$`1` <- gsub('[A-Z]{1,}', '', vkl_table$`1`)
vkl_table$`2` <- gsub('[A-Z]{1,}', '', vkl_table$`2`)
vkl_table <- vkl_table %>% mutate_at(c('1', '2'), as.numeric)
vkl_table$char <- ifelse(is.na(vkl_table$char), '', vkl_table$char )

vkl_table <- vkl_table[order(vkl_table$chain,vkl_table$`1`, vkl_table$`2`, vkl_table$char, decreasing = c(F, F, F, F)),]

vk_table <- subset(vkl_table, chain == 'K')
vl_table <- subset(vkl_table, chain == 'L')

vk_table <- select(vk_table, VGene_Light:Pooled)
vl_table <- select(vl_table, VGene_Light:Pooled)

colnames(vk_table) <- gsub('X', '', colnames(vk_table))
colnames(vl_table) <- gsub('X', '', colnames(vl_table))

mper_pre <- read.csv('2_Baseline_MPER_Abs.csv', header = T, nrows = 2)
mper_post <- read.csv('Total_MPERAbs_Post_IMM_WW120622.csv', header = T)

mper_pre$H_VGENE <- gsub('\\*[0-9]{1,}', '', mper_pre$H_VGENE) %>% gsub('IGHV', 'VH', .)
mper_post$H_VGENE <- gsub('\\*[0-9]{1,}', '', mper_post$H_VGENE) %>% gsub('IGHV', 'VH', .)


mper_pre_vgene_counts <- mper_pre %>% group_by(H_VGENE) %>% summarize(n=n()) %>% mutate(freq=n/sum(n)) %>% data.frame()
mper_post_vgene_counts <- mper_post %>% group_by(H_VGENE) %>% summarize(n=n()) %>% mutate(freq=n/sum(n)) %>% data.frame()

mper_pre_vgene_counts$PTID <- 'MPER_Pre'
mper_post_vgene_counts$PTID <- 'MPER_Post'

colnames(mper_pre_vgene_counts)[1] <- c('VGene_Heavy')
colnames(mper_post_vgene_counts)[1] <- c('VGene_Heavy')

all_list <- list(pooled_vgene_counts, mper_pre_vgene_counts, mper_post_vgene_counts)

all_mper_vgene_counts <- smartbind(list = all_list)

all_mper_vgene_counts$gene <- gsub('VH', '', all_mper_vgene_counts$VGene_Heavy)
all_mper_vgene_counts <- cbind(all_mper_vgene_counts, str_split_fixed(all_mper_vgene_counts$gene, pattern = '-', n = Inf))
all_mper_vgene_counts$char <- str_extract(all_mper_vgene_counts$gene,'[A-Z]{1,}')
all_mper_vgene_counts$`2` <- gsub('[A-Z]{1,}', '', all_mper_vgene_counts$`2`)
all_mper_vgene_counts <- all_mper_vgene_counts %>% mutate_at(c('1', '2', '3'), as.numeric)
all_mper_vgene_counts$char <- ifelse(is.na(all_mper_vgene_counts$char), '', all_mper_vgene_counts$char )
all_mper_vgene_counts$`3` <- ifelse(is.na(all_mper_vgene_counts$`3`), 0, all_mper_vgene_counts$`3` )

all_mper_vgene_counts <- all_mper_vgene_counts[order(all_mper_vgene_counts$`1`, all_mper_vgene_counts$`2`, all_mper_vgene_counts$`3`, all_mper_vgene_counts$char, decreasing = c(F, F, F, F)),]
all_mper_vgene_counts$VGene_Heavy <- factor(all_mper_vgene_counts$VGene_Heavy, levels = unique(all_mper_vgene_counts$VGene_Heavy))

all_mper_vgene_counts <- select(all_mper_vgene_counts, VGene_Heavy:PTID)

vhgenes <- c('IGHV1-69', 'IGHV2-5', 'IGHV3-15', 'IGHV3-30-3','IGHV3-49', 'IGHV3-73', 'IGHV5-51', 'IGHV7-4-1')
vhgenes <- gsub('IGHV', 'VH', vhgenes)

vh_table_subset <- subset(all_mper_vgene_counts, VGene_Heavy %in% vhgenes)
vh_table_subset$PTID2 <- case_when(vh_table_subset$PTID == 'Pooled' ~ 'Pooled Baseline (N = 209,523)', 
                                   vh_table_subset$PTID == 'MPER_Pre' ~ 'MPER Abs Baseline (N = 2)', 
                                   vh_table_subset$PTID == 'MPER_Post' ~ 'MPER Abs Post-Immunization (N = 87)')

vh_table_subset$PTID2 <- factor(vh_table_subset$PTID2, levels = c('Pooled Baseline (N = 209,523)', 'MPER Abs Baseline (N = 2)', 'MPER Abs Post-Immunization (N = 87)'))
vh_table_subset$VGene_Heavy <- factor(vh_table_subset$VGene_Heavy, levels = )


d <- data.frame(expand(vh_table_subset, VGene_Heavy, PTID2))

d <- merge(d, vh_table_subset, by = c('PTID2', 'VGene_Heavy'), all = T)
d$freq <- ifelse(is.na(d$freq), 0, d$freq)

d2 <- subset(vh_table_subset, PTID != 'MPER_Pre')
d2$PTID2 <- factor(d2$PTID2, levels = c('Pooled Baseline (N = 209,523)', 'MPER Abs Post-Immunization (N = 87)'))
d2 <- data.frame(expand(d2, VGene_Heavy, PTID2))

d2 <- merge(d2, vh_table_subset, by = c('PTID2', 'VGene_Heavy'), all.x = T)
d2$freq <- ifelse(is.na(d2$freq), 0, d2$freq)

h_sub <- ggplot(d2, aes(x = VGene_Heavy, y = freq, fill = PTID2)) + geom_bar(stat = 'identity', position = 'dodge')  + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = .5)) + 
  xlab('V-Gene') + ylab('Frequency') + labs(fill = '') + ggtitle('Heavy') + scale_x_discrete(drop = F)
h_sub

# VK 
vkgenes <- c('IGKV1-6', 'IGKV1-12', 'IGKV1-16', 'IGKV1-27', 'IGKV1-39', 'IGKV2-28', 'IGKV3-20', 'IGKV4-1')
vkgenes <- gsub('IGKV', 'VK', vkgenes)

vk_table_subset <- subset(vk_table, VGene_Light %in% vkgenes)

vk_table_subset <- pivot_longer(vk_table_subset, cols = `133.23`:Pooled, names_to = 'PTID', values_to = 'Freq' )
vk_table_subset$PTID <- gsub('\\.', '-', vk_table_subset$PTID)

vk_table_subset$gene <- gsub('VK', '', vk_table_subset$VGene_Light)
vk_table_subset <- cbind(vk_table_subset, str_split_fixed(vk_table_subset$gene, pattern = '-', n = Inf))
vk_table_subset <- vk_table_subset %>% mutate_at(c('1', '2'), as.numeric)
vk_table_subset <- vk_table_subset[order(vk_table_subset$`1`, vk_table_subset$`2`, decreasing = F),]
vk_table_subset$VGene_Light <- factor(vk_table_subset$VGene_Light, levels = unique(vk_table_subset$VGene_Light))

mper_pre$L_VGENE <- gsub('\\*[0-9]{1,}', '', mper_pre$L_VGENE) %>% gsub('IGKV', 'VK', .)
mper_post$L_VGENE <- gsub('\\*[0-9]{1,}', '', mper_post$L_VGENE) %>% gsub('IGKV', 'VK', .) %>% gsub('IGLV', 'VL', .)

mper_pre_vgene_counts <- mper_pre %>% group_by(L_VGENE) %>% summarize(n=n()) %>% mutate(freq=n/sum(n)) %>% data.frame()
mper_post_vgene_counts <- mper_post %>% group_by(L_VGENE) %>% summarize(n=n()) %>% mutate(freq=n/sum(n)) %>% data.frame()

mper_pre_vgene_counts$PTID <- 'MPER_Pre'
mper_post_vgene_counts$PTID <- 'MPER_Post'

colnames(mper_pre_vgene_counts)[1] <- c('VGene_Light')
colnames(mper_post_vgene_counts)[1] <- c('VGene_Light')

all_list <- list(pooled_vklgene_counts,mper_post_vgene_counts)

all_mper_vgene_counts <- smartbind(list = all_list)
all_mper_vgene_counts_kappa <- subset(all_mper_vgene_counts, grepl(pattern = 'VK', VGene_Light ))

vk_table_subset <- subset(all_mper_vgene_counts_kappa, VGene_Light %in% vkgenes)
vk_table_subset$PTID2 <- case_when(vk_table_subset$PTID == 'Pooled' ~ 'Pooled Baseline (N = 209,523)', 
                                   vk_table_subset$PTID == 'MPER_Post' ~ 'MPER Abs Post-Immunization (N = 87)')

vk_table_subset$PTID2 <- factor(vk_table_subset$PTID2, levels = c('Pooled Baseline (N = 209,523)', 'MPER Abs Post-Immunization (N = 87)'))
vk_table_subset$VGene_Light <- factor(vk_table_subset$VGene_Light, levels = vkgenes)

d <- data.frame(expand(vk_table_subset, VGene_Light, PTID2))

d <- merge(d, vk_table_subset, by = c('PTID2', 'VGene_Light'), all = T)
d$freq <- ifelse(is.na(d$freq), 0, d$freq)

k <- ggplot(d, aes(x = VGene_Light, y = freq, fill = PTID2)) + geom_bar(stat = 'identity', position = 'dodge')  + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = .5)) + 
  xlab('V-Gene') + ylab('Frequency') + labs(fill = '') + ggtitle('Kappa') 
k

#VL 
vlgenes <- c('IGLV2-18', 'IGLV3-1', 'IGLV3-19', 'IGLV8-61')
vlgenes <- gsub('IGLV', 'VL', vlgenes)
all_mper_vgene_counts_lambda <- subset(all_mper_vgene_counts, grepl(pattern = 'VL', VGene_Light ))

vl_table_subset <- subset(all_mper_vgene_counts_lambda, VGene_Light %in% vlgenes)
vl_table_subset$PTID2 <- case_when(vl_table_subset$PTID == 'Pooled' ~ 'Pooled Baseline (N = 209,523)', 
                                   vl_table_subset$PTID == 'MPER_Post' ~ 'MPER Abs Post-Immunization (N = 87)')

vl_table_subset$PTID2 <- factor(vl_table_subset$PTID2, levels = c('Pooled Baseline (N = 209,523)', 'MPER Abs Post-Immunization (N = 87)'))
vl_table_subset$VGene_Light <- factor(vl_table_subset$VGene_Light, levels = vlgenes)

d <- data.frame(expand(vl_table_subset, VGene_Light, PTID2))

d <- merge(d, vl_table_subset, by = c('PTID2', 'VGene_Light'), all = T)
d$freq <- ifelse(is.na(d$freq), 0, d$freq)

l <- ggplot(d, aes(x = VGene_Light, y = freq, fill = PTID2)) + geom_bar(stat = 'identity', position = 'dodge')  + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = .5)) + 
  xlab('V-Gene') + ylab('Frequency') + labs(fill = '') + ggtitle('Lambda') 
l

h2 <- h_sub + theme_bw() + theme(plot.title = element_text(hjust = .5 ), axis.text.x = element_text(angle = 45, hjust = 1)) + 
  geom_vline(xintercept = (0:10)+0.5, linetype = 'dashed', color = 'black')
h2

k2 <- k + theme_bw() + theme(plot.title = element_text(hjust = .5 ), axis.text.x = element_text(angle = 45, hjust = 1)) + 
  geom_vline(xintercept = (0:10)+0.5, linetype = 'dashed', color = 'black')
k2

l2 <- l + theme_bw() + theme(plot.title = element_text(hjust = .5 ), axis.text.x = element_text(angle = 45, hjust = 1)) + 
  geom_vline(xintercept = (0:4)+0.5, linetype = 'dashed', color = 'black')
l2

png('VGene_MPER_Baseline_Frequencies.png', height = 10, width = 8, units = 'in', res = 600)
ggarrange(h2, k2, l2, ncol=1, nrow=3, common.legend = TRUE, legend="right")
dev.off()

# Mutation Frequency
all_mu_freq_heavy <- select(all, PTID, MuFreq_Heavy)
all_mu_freq_heavy$Chain <- 'Heavy'
colnames(all_mu_freq_heavy) <- c('PTID', 'MuFreq', 'Chain')

all_mu_freq_light <- select(all, PTID, MuFreq_Light, Chain)
colnames(all_mu_freq_light) <- c('PTID', 'MuFreq', 'Chain')

all_mu_freq <- smartbind(all_mu_freq_heavy, all_mu_freq_light)

all_mu_freq$PTID2 <- 'Pooled Baseline (N = 209,523)'

#mper_pre$H_MUT_FREQ <- gsub('\\%', '', mper_pre$H_MUT_FREQ) %>% as.numeric()
#mper_pre$L_MUT_FREQ <- gsub('\\%', '', mper_pre$L_MUT_FREQ) %>% as.numeric()

# get new mufreq from immunogenetics

mpers <- read.csv('/datacommons/dhvi/10X_Genomics_data/HVTN133/Analysis/MPER_Clonal_Analysis/HVTN133_MPER_Abs.merged_clonal_analysis.csv')
#remove pre-immunization obs
mpers <- subset(mpers, !(UID %in% c('H027739+K024786', 'H028315+K025152')))

mper_post_hmufreq <- data.frame(Chain = 'Heavy', MuFreq = mpers$MuFreq_Heavy)
mper_post_klmufreq <- data.frame(Chain = mpers$Chain, MuFreq = mpers$MuFreq_Light)

mper_post_mufreq <- smartbind(list = list(mper_post_hmufreq, mper_post_klmufreq))
mper_post_mufreq$PTID2 <- 'MPER Abs Post-Immunization (N = 87)'

all_mper_mu_freq <- smartbind(list = list(all_mu_freq, mper_post_mufreq))
all_mper_mu_freq$PTID2 <- factor(all_mper_mu_freq$PTID2, levels = c('Pooled Baseline (N = 209,523)',
                                                                    'MPER Abs Post-Immunization (N = 87)'))
all_mper_mu_freq$Chain <- factor(all_mper_mu_freq$Chain, levels = c('Heavy', 'Kappa', 'Lambda'))
#all_mper_mu_freq <- rbind(c('MPER Abs Baseline (N = 2)', 4000, 'Lambda', 'MPER Abs Baseline (N = 2)'), all_mper_mu_freq)
#all_mper_mu_freq$MuFreq <- as.numeric(all_mper_mu_freq$MuFreq)


png('All_MPER_MuFreq_ViolinPlots.png', height = 6, width = 8, units = 'in', res = 600)
ggplot(all_mper_mu_freq, aes(x = Chain, fill = PTID2, y = MuFreq)) + geom_violin(draw_quantiles = c(0.25, .5, 0.75), scale = 'width') +
  coord_cartesian(ylim = c(0,.31)) + ggtitle('Mutation Frequency') + ylab('Mutation Frequency') + 
  theme_bw() + theme(plot.title = element_text(hjust = .5 ), axis.text.x = element_text(angle = 45, hjust = 1)) + 
  labs(fill = '') + 
  geom_vline(xintercept = (0:3)+0.5, linetype = 'dashed', color = 'black')
dev.off()

# CDR3Length
all_cdr3length_heavy <- select(all, PTID, CDR3Length_Heavy)
all_cdr3length_heavy$Chain <- 'Heavy'
colnames(all_cdr3length_heavy) <- c('PTID', 'CDR3Length', 'Chain')

all_cdr3length_light <- select(all, PTID, CDR3Length_Light, Chain)
colnames(all_cdr3length_light) <- c('PTID', 'CDR3Length', 'Chain')

all_cdr3length <- smartbind(all_cdr3length_heavy, all_cdr3length_light)

all_cdr3length$PTID2 <- 'Pooled Baseline (N = 209,523)'
all_cdr3length$CDR3LengthAA <- all_cdr3length$CDR3Length/3

mper_pre_hcdr3length <- data.frame(Chain = 'Heavy', CDR3LengthAA = mper_pre$H_CDR3_LEN_AA)
mper_pre_kcdr3length <- data.frame(Chain = 'Kappa', CDR3LengthAA = mper_pre$L_CDR3_LEN_AA)

mper_pre_cdr3length <- smartbind(mper_pre_hcdr3length, mper_pre_kcdr3length)
mper_pre_cdr3length$PTID2 <- 'MPER Abs Baseline (N = 2)'

mper_post_hcdr3length <- data.frame(Chain = 'Heavy', CDR3LengthAA = mper_post$H_CDR3_LEN_AA)
mper_post_kcdr3length <- data.frame(Chain = 'Kappa', CDR3LengthAA = subset(mper_post, grepl(pattern = 'K', mper_post$L_VGENE))[,'L_CDR3_LEN_AA'])
mper_post_lcdr3length <- data.frame(Chain = 'Lambda', CDR3LengthAA = subset(mper_post, grepl(pattern = 'L', mper_post$L_VGENE))[,'L_CDR3_LEN_AA'])

mper_post_cdr3length <- smartbind(list = list(mper_post_hcdr3length, mper_post_kcdr3length, mper_post_lcdr3length))
mper_post_cdr3length$PTID2 <- 'MPER Abs Post-Immunization (N = 87)'

all_mper_cdr3length <- smartbind(list = list(all_cdr3length, mper_pre_cdr3length, mper_post_cdr3length))
all_mper_cdr3length$PTID2 <- factor(all_mper_cdr3length$PTID2, levels = c('Pooled Baseline (N = 209,523)',
                                                                          'MPER Abs Baseline (N = 2)',
                                                                          'MPER Abs Post-Immunization (N = 87)'))
all_mper_cdr3length$Chain <- factor(all_mper_cdr3length$Chain, levels = c('Heavy', 'Kappa', 'Lambda'))
#all_mper_cdr3length <- rbind(c('MPER Abs Baseline (N = 2)', 4000, 'Lambda', 'MPER Abs Baseline (N = 2)', 40000), all_mper_cdr3length)
#all_mper_cdr3length$CDR3LengthAA <- as.numeric(all_mper_cdr3length$CDR3LengthAA)

heavy_cdr3length <- subset(all_mper_cdr3length, Chain == 'Heavy')
light_cdr3length <- subset(all_mper_cdr3length, Chain != 'Heavy')

heavy_mper_cdr3length_counts <- heavy_cdr3length %>% group_by(PTID2, CDR3LengthAA) %>% summarize(n=n())%>% mutate(freq = n /sum(n)) %>% data.frame()
heavy_mper_cdr3length_counts$Chain <- 'Heavy'
light_mper_cdr3length_counts <- light_cdr3length %>% group_by(PTID2, Chain, CDR3LengthAA) %>% summarize(n=n()) %>% mutate(freq = n/sum(n)) %>% data.frame()

light_mper_cdr3length_counts$N <- case_when(light_mper_cdr3length_counts$PTID2 == 'Pooled Baseline (N = 209,523)' ~ 209523,
                                            light_mper_cdr3length_counts$PTID2 == 'MPER Abs Baseline (N = 2)' ~ 2,
                                            light_mper_cdr3length_counts$PTID2 == 'MPER Abs Post-Immunization (N = 87)' ~ 87) 

light_mper_cdr3length_counts$freq <- light_mper_cdr3length_counts$n/ light_mper_cdr3length_counts$N

all_mper_cdr3length_counts <- smartbind(heavy_mper_cdr3length_counts, light_mper_cdr3length_counts)
all_mper_cdr3length_counts$Chain <- factor(all_mper_cdr3length_counts$Chain, levels = c('Heavy', 'Kappa', 'Lambda'))

all_mper_cdr3length_counts$PTID <- case_when(all_mper_cdr3length_counts$PTID2 == 'Pooled Baseline (N = 209,523)' ~ 'Pooled Baseline',
                                             all_mper_cdr3length_counts$PTID2 == 'MPER Abs Baseline (N = 2)' ~ 'MPER Abs Baseline',
                                             all_mper_cdr3length_counts$PTID2 == 'MPER Abs Post-Immunization (N = 87)' ~ 'MPER Abs Post-Immunization') 

all_mper_cdr3length_counts$PTID <- factor(all_mper_cdr3length_counts$PTID, levels = c('Pooled Baseline',
                                                                                      'MPER Abs Baseline',
                                                                                      'MPER Abs Post-Immunization'))

all_mper_cdr3length_counts <- subset(all_mper_cdr3length_counts, PTID != 'MPER Abs Baseline')
all_mper_cdr3length_counts_sub <- subset(all_mper_cdr3length_counts, Chain == 'Heavy' & CDR3LengthAA >=5 & CDR3LengthAA <=26 | 
                                           Chain %in% c('Kappa', 'Lambda') & CDR3LengthAA >=9 & CDR3LengthAA<=13)


count <- 0
breaks_fun <- function(x) {
  count <<- count + 1L
  switch(
    count,
    c(5,10,15,20,25),
    c(9,10,11,12,13),
    c(9,10,11,12,13),
    c(5,10,15,20,25),
    c(9,10,11,12,13),
    c(9,10,11,12,13)
  )
}

breaks_fun <- function(x){
  if (max(x) > 15){
    c(5,10,15,20,25)}
  else{
    c(9,10,11,12,13)}
}

png('All_MPER_CDR3Length_Barplots_Subset.png', height = 6, width = 9, units = 'in', res = 600)
ggplot(all_mper_cdr3length_counts_sub, aes(x = CDR3LengthAA, y = freq, fill = PTID2)) + geom_bar(position = 'dodge', stat = 'identity') + 
  ggtitle('CDR3 Length') + ylab('Frequency') + xlab('CDR3 Length (AA)') +
  theme_bw() + theme(plot.title = element_text(hjust = .5 )) + 
  facet_grid(PTID~Chain, scales = 'free') + labs(fill = '') + scale_x_continuous(breaks = breaks_fun)
dev.off()

