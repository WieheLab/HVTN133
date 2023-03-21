library(gtools)
library(dplyr)
library(ggplot2)
library(stringr)
library(ggridges)
library(ggpubr)

setwd('/datacommons/dhvi/10X_Genomics_data/HVTN133/Analysis/VDJ')

dat <- read.csv('HVTN133_all_VH7-4-1_pairings_by_PTID_TimePoint.csv')
View(subset(dat, VGene_Heavy == 'VH7-4-1'))

### pre and post frequencies of VH7-4-1 for all PTIDs outputted to a table

dat <- subset(dat, PTID != 'V5 - 778192591')

dat$PTID2 <- case_when(dat$PTID == 'V1 - 726212200' ~ '133-23',
                       dat$PTID == 'V2 - 821167321' ~ '133-39',
                       dat$PTID == 'V3 - 726273758' ~ '133-35',
                       dat$PTID == 'V4 - 726888509' ~ '133-33',
                       dat$PTID == 'V6 - 778232955' ~ '133-30',
                       dat$PTID == 'P1 - 726794953' ~ '133-17 (Placebo)')

dat$TimePoint2 <- ifelse(dat$TimePoint == 'Visit 2', 'Pre-Vaccine', 'Post-Vaccine')
dat$TimePoint2 <- factor(dat$TimePoint2, levels = c('Pre-Vaccine', 'Post-Vaccine'))

agg <- aggregate(freq_per_100k ~ PTID2 + TimePoint2 , data = dat, FUN = sum, drop = F)
str(agg)
agg$PTID2 <- factor(agg$PTID2, levels = c('133-23', '133-39', '133-35', '133-33', '133-30', '133-17 (Placebo)'))


png('HVTN133_VH7-4-1_freq_per_100k_by_PTID_Timepoint_New.png', height = 6, width = 8, units = 'in', res = 600)
ggplot(agg, aes(x = PTID2, y = freq_per_100k, fill = TimePoint2)) + 
  geom_bar(stat = 'identity', position = 'dodge') + 
  theme(plot.title = element_text(hjust = .5), axis.text.x = element_text(angle=45, hjust =1)) + ggtitle('VH7-4-1 Frequency by Subject & Timepoint') + 
  ylab('Frequency per 100k') + xlab('Subject') + labs(fill = 'Timepoint') + 
  scale_fill_manual(values = c('cadetblue2', 'indianred2'), drop = F)
dev.off()

pdf('HVTN133_VH7-4-1_freq_per_100k_by_PTID_Timepoint_New.pdf', height = 6, width = 8)
ggplot(agg, aes(x = PTID2, y = freq_per_100k, fill = TimePoint2)) + 
  geom_bar(stat = 'identity', position = 'dodge') + 
  theme(plot.title = element_text(hjust = .5), axis.text.x = element_text(angle=45, hjust =1)) + ggtitle('VH7-4-1 Frequency by Subject & Timepoint') + 
  ylab('Frequency per 100k') + xlab('Subject') + labs(fill = 'Timepoint') + 
  scale_fill_manual(values = c('cadetblue2', 'indianred2'), drop = F)
dev.off()

vh7_4_1 <- subset(all_sub, VGene_Heavy == 'VH7-4-1')
vh7_4_1 <- factor(vh7_4_1$PTID)
vh7_4_1_counts <- vh7_4_1 %>% group_by(PTID, TimePoint2, .drop = F) %>% summarize(quant90 = quantile(VMuFreq_Heavy, probs =  .9, na.rm = T))  
