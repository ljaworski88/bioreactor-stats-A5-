###############################################################################
# Script to handle all PCR data from a bioreactor experiment.
# The bioreactors are made according to the plans found here: [link]
# The loading pattern for the experiment was 8 hours of cyclic loading at
# between 0.4-0.8Mpa and then 16 hours of static loading at 0.2Mpa.
# Two reactors were kept at 21% oxygen tension while the other two were kept at
# 5% oxygen tension. The IVDs were collected at days 7 and 14. From previous 
# studies the reference genes that will be used here are GAPDH, RPL4, and YWHAZ.
# Genes included in the study are: MMP3, MMP13, PACE4, TIMP1, TIMP2, TIMP3,
# ADAMTS1, ADMAMTS2, HIF-1a, and T.

# Author <- Lukas Jaworski, l.jaworski@umiami.edu
# Institution <- University of Miami
# PI <- Alicia Jackson, a.jackson2@miami.edu

# Copyright University of Miami 2019
###############################################################################

###############################################################################
## Library imports
###############################################################################

library(ReadqPCR)
library(NormqPCR)
library(ggplot2)
library(gridExtra)
library(cowplot)

###############################################################################
## Functions
###############################################################################

aov.results <- function(genes, data, formula){
  # A function to get the statistically significant results without any results,
  # which are not significant
  
  
  anova.results <- list()
  for (x in genes){
    anova.results[[x]] <- aov(as.formula(gsub('{}', x, formula, fixed = T)),
                              data = data)
  }
  # ANOVA results and perform post hoc tests
  aov.summary<-list()
  for (x in genes){
    temp.summary <- summary.aov(anova.results[[x]])[[1]]
    aov.summary[[x]][[1]] <- temp.summary[temp.summary['Pr(>F)'] <= 0.1 &
                                            !is.na(temp.summary['Pr(>F)']),
                                          'Pr(>F)',
                                          drop = F]
    if (nrow(aov.summary[[x]][[1]])>0){
      temp.tukey <- TukeyHSD(anova.results[[x]])
      aov.summary[[x]][[2]] <- lapply(temp.tukey,
                                      function(x) x[x[, 'p adj'] <= 0.1, 'p adj'])
    }
  }
  return(aov.summary)
}

###############################################################################
## Script
###############################################################################
# Set working Directory
if (Sys.info()[['sysname']]=='Linux'){
  working.directory <- paste0('/home/',Sys.info()[['user']],
                              '/Dropbox/Lab/Experiments/PCR/',
                              'a5-bioreactor-ox_lvl/eds_files')
} else {
  working.directory <- paste0('C:/Users/',Sys.info()[['user']],
                              '/Dropbox/Lab/Experiments/PCR/',
                              'a5-bioreactor-ox_lvl/eds_files')
}
working.directory <- paste0('P:/Dropbox/Lab/Experiments/PCR/',
                            'a5-bioreactor-ox_lvl/eds_files')
setwd(working.directory)

# Assemble the pheno data for the experiment
pheno.data <- read.csv("sample_id.csv")
colnames(pheno.data)[colnames(pheno.data)=='id'] <- 'Sample'
colnames(pheno.data)[1] <- 'Sample'
pheno.data <- pheno.data[with(pheno.data, order(Sample)), ]
pheno.data$day <- factor(pheno.data$day)
pheno.data$pig <- factor(pheno.data$pig)
rownames(pheno.data) <- pheno.data$Sample
# Metadata for our annotated dataframe
meta <- data.frame(c('Sample','pig', 'day', 'oxygen', 'loading', 'tissue'))

# Gotta match that eSet format so boom
names(meta) <- 'labelDescription'

# need to clean up some of the sample names to have leading zeros, remove an extraneous column, and rename a gene
sample.cleaning <- read.csv('result.csv')
sample.cleaning$Contaminated <- NULL
sample.cleaning$Sample <- gsub("b([[:digit:]])$","b0\\1", sample.cleaning$Sample)
sample.cleaning$Detector <- gsub("18s","x18s", sample.cleaning$Detector)

genes.of.interest <- c("ADAMTS4", 
                       "ADAMTS5", 
                       "HIF1a", 
                       "MMP13",
                       "MMP3",
                       "PACE4", 
                       "RPL4",
                       "T", 
                       "TIMP1",
                       "TIMP2",
                       "TIMP3",
                       "x18s" )

# adding in NA values for missing values

sample.cleaning <- sample.cleaning[sample.cleaning$Detector %in% genes.of.interest,]
sample.cleaning[c('Sample', 'Detector')] <- lapply(sample.cleaning[c('Sample', 'Detector')], factor)
samples <- levels(sample.cleaning$Sample)
detectors <- levels(sample.cleaning$Detector)
allDetectors <- sample.cleaning$Detector
for (sample in samples) { # for each sample
  total.detectors <- length(allDetectors[sample.cleaning$Sample == sample])
  individual.detectors <- length(levels(allDetectors[sample.cleaning$Sample == sample]))
  tech.reps <- total.detectors/individual.detectors
  if ((tech.reps %% 1) != 0) { # if total number of replicates not a multiple of number of individual detectors
    missing.genes <- setdiff(genes.of.interest, sample.cleaning[sample.cleaning$Sample == sample, "Detector"])
    for (missing.gene in missing.genes){
      df <- data.frame(Sample=sample, Detector = missing.gene, Cq = NA)
      df[c('Sample', 'Detector')] <- lapply(df[c('Sample', 'Detector')], factor)
      sample.cleaning <- rbind(sample.cleaning, df)
    }
  }
  
}

# sample.cleaning$Cq[sample.cleaning$Cq > 39.9] <- NA
write.table(sample.cleaning, file = 'result.txt', row.names = FALSE, sep = ' ')
pheno.data <- pheno.data[pheno.data$Sample %in% levels(sample.cleaning$Sample),]
#Finally making the annotated dataframe
annotated.pheno <- AnnotatedDataFrame(pheno.data,
                                      meta, 
                                      c('sampleNames','sampleColumns'))
Cq.vals <- read.qPCR('result.txt', annotated.pheno)

# The control genes were verified using another script
hk.gene.gluc <- c('x18s', 
                  'RPL4')
dCq.vals <- deltaCq(qPCRBatch = Cq.vals, 
                    hkgs = hk.gene.gluc, 
                    calc="geom")

gene.expression <- data.frame(t(exprs(dCq.vals)))
gene.expression$day <- dCq.vals@phenoData@data$day
gene.expression$Sample <- dCq.vals@phenoData@data$Sample
gene.expression$pig <- dCq.vals@phenoData@data$pig
gene.expression$oxygen <- dCq.vals@phenoData@data$oxygen
gene.expression$loading <- dCq.vals@phenoData@data$loading
gene.expression$tissue <- dCq.vals@phenoData@data$tissue
gene.expression <- gene.expression[!(gene.expression$Sample %in% c('b33', 'b45', 'b47', 'b59', 'b91', 'b93', 'b95')),]

# Filter out NP data only

gene.expression <- gene.expression[gene.expression$tissue == 'NP', ]
gene.expression.loaded <- subset(gene.expression, loading == '1' | (loading == 'x' & day == '0'))
gene.expression.unloaded <- subset(gene.expression, loading == 'x')

# run the ANOVAs
ADAMTS4.anova <- aov(ADAMTS4 ~ oxygen * day * loading + pig,
                     data = gene.expression) 
summary.aov(ADAMTS4.anova)

ADAMTS5.anova <- aov(ADAMTS5 ~ oxygen * day * loading + pig,
                     data = gene.expression)
summary.aov(ADAMTS5.anova)

HIF1a.anova <- aov(HIF1a ~ oxygen * day * loading + pig,
                   data = gene.expression) 
summary.aov(HIF1a.anova)

MMP13.anova <- aov(MMP13 ~ oxygen * day * loading + pig,
                   data = gene.expression)
summary.aov(MMP13.anova)

MMP3.anova <- aov(MMP3 ~ oxygen * day * loading + pig,
                   data = gene.expression)
summary.aov(MMP3.anova)

PACE4.anova <- aov(PACE4 ~ oxygen * day * loading + pig,
                   data = gene.expression)
summary.aov(PACE4.anova)

T.anova <- aov(T ~ oxygen * day * loading + pig,
                   data = gene.expression)
summary.aov(T.anova)

TIMP1.anova <- aov(TIMP1 ~ oxygen * day * loading + pig,
                   data = gene.expression)
summary.aov(TIMP1.anova)

TIMP2.anova <- aov(TIMP2 ~ oxygen * day * loading + pig,
                   data = gene.expression)
summary.aov(TIMP2.anova)

TIMP3.anova <- aov(TIMP3 ~ oxygen * day * loading + pig,
                   data = gene.expression)
summary.aov(TIMP3.anova)

# Now split the based on oxygen tension as well for graphing

gene.expression.loaded.plot.ADAMTS4 <- ggplot(gene.expression.loaded, aes(day, -ADAMTS4, fill=oxygen)) +
  ggtitle('ADAMTS 4') +
  labs(x = 'Time (days)',
       y = expression(paste('Gene Expression (-', Delta,'Cq)'))) +
  geom_boxplot(position = position_dodge()) +
  scale_fill_manual(name = 'Oxygen\nTension', labels = c('21%', '5%', 'Harvest'), values = c('#005030', '#F47321', '#500c00'))

gene.expression.loaded.plot.ADAMTS5 <- ggplot(gene.expression.loaded, aes(day, -ADAMTS5, fill=oxygen)) +
  ggtitle('ADAMTS 5') +
  labs(x = 'Time (days)',
       y = expression(paste('Gene Expression (-', Delta,'Cq)'))) +
  geom_boxplot(position = position_dodge()) +
  scale_fill_manual(name = 'Oxygen\nTension', labels = c('21%', '5%', 'Harvest'), values = c('#005030', '#F47321', '#500c00'))

gene.expression.loaded.plot.HIF1a <- ggplot(gene.expression.loaded, aes(day, -HIF1a, fill=oxygen)) +
  ggtitle('HIF1a') +
  labs(x = 'Time (days)',
       y = expression(paste('Gene Expression (-', Delta,'Cq)'))) +
  geom_boxplot(position = position_dodge()) +
  scale_fill_manual(name = 'Oxygen\nTension', labels = c('21%', '5%', 'Harvest'), values = c('#005030', '#F47321', '#500c00'))

gene.expression.loaded.plot.MMP13 <- ggplot(gene.expression.loaded, aes(day, -MMP13, fill=oxygen)) +
  ggtitle('MMP 13') +
  labs(x = 'Time (days)',
       y = expression(paste('Gene Expression (-', Delta,'Cq)'))) +
  geom_boxplot(position = position_dodge()) +
  scale_fill_manual(name = 'Oxygen\nTension', labels = c('21%', '5%', 'Harvest'), values = c('#005030', '#F47321', '#500c00'))

gene.expression.loaded.plot.MMP3 <- ggplot(gene.expression.loaded, aes(day, -MMP3, fill=oxygen)) +
  ggtitle('MMP 3') +
  labs(x = 'Time (days)',
       y = expression(paste('Gene Expression (-', Delta,'Cq)'))) +
  geom_boxplot(position = position_dodge()) +
  scale_fill_manual(name = 'Oxygen\nTension', labels = c('21%', '5%', 'Harvest'), values = c('#005030', '#F47321', '#500c00'))

gene.expression.loaded.plot.PACE4 <- ggplot(gene.expression.loaded, aes(day, -PACE4, fill=oxygen)) +
  ggtitle('PACE 4') +
  labs(x = 'Time (days)',
       y = expression(paste('Gene Expression (-', Delta,'Cq)'))) +
  geom_boxplot(position = position_dodge()) +
  scale_fill_manual(name = 'Oxygen\nTension', labels = c('21%', '5%', 'Harvest'), values = c('#005030', '#F47321', '#500c00'))

gene.expression.loaded.plot.TIMP1 <- ggplot(gene.expression.loaded, aes(day, -TIMP1, fill=oxygen)) +
  ggtitle('TIMP 1') +
  labs(x = 'Time (days)',
       y = expression(paste('Gene Expression (-', Delta,'Cq)'))) +
  geom_boxplot(position = position_dodge()) +
  scale_fill_manual(name = 'Oxygen\nTension', labels = c('21%', '5%', 'Harvest'), values = c('#005030', '#F47321', '#500c00'))

gene.expression.loaded.plot.TIMP2 <- ggplot(gene.expression.loaded, aes(day, -TIMP2, fill=oxygen)) +
  ggtitle('TIMP 2') +
  labs(x = 'Time (days)',
       y = expression(paste('Gene Expression (-', Delta,'Cq)'))) +
  geom_boxplot(position = position_dodge()) +
  scale_fill_manual(name = 'Oxygen\nTension', labels = c('21%', '5%', 'Harvest'), values = c('#005030', '#F47321', '#500c00'))

gene.expression.loaded.plot.TIMP3 <- ggplot(gene.expression.loaded, aes(day, -TIMP3, fill=oxygen)) +
  ggtitle('TIMP 3') +
  labs(x = 'Time (days)',
       y = expression(paste('Gene Expression (-', Delta,'Cq)'))) +
  geom_boxplot(position = position_dodge()) +
  scale_fill_manual(name = 'Oxygen\nTension', labels = c('21%', '5%', 'Harvest'), values = c('#005030', '#F47321', '#500c00'))

gene.expression.loaded.plot.T <- ggplot(gene.expression.loaded, aes(day, -T, fill=oxygen)) +
  ggtitle('Brachyury') +
  labs(x = 'Time (days)',
       y = expression(paste('Gene Expression (-', Delta,'Cq)'))) +
  geom_boxplot(position = position_dodge()) +
  scale_fill_manual(name = 'Oxygen\nTension', labels = c('21%', '5%', 'Harvest'), values = c('#005030', '#F47321', '#500c00'))

gene.expression.unloaded.plot.ADAMTS4 <- ggplot(gene.expression.unloaded, aes(day, -ADAMTS4, fill=oxygen)) +
  ggtitle('ADAMTS 4') +
  labs(x = 'Time (days)',
       y = expression(paste('Gene Expression (-', Delta,'Cq)'))) +
  geom_boxplot(position = position_dodge()) +
  scale_fill_manual(name = 'Oxygen\nTension', labels = c('21%', '5%', 'Harvest'), values = c('#005030', '#F47321', '#500c00'))

gene.expression.unloaded.plot.ADAMTS5 <- ggplot(gene.expression.unloaded, aes(day, -ADAMTS5, fill=oxygen)) +
  ggtitle('ADAMTS 5') +
  labs(x = 'Time (days)',
       y = expression(paste('Gene Expression (-', Delta,'Cq)'))) +
  geom_boxplot(position = position_dodge()) +
  scale_fill_manual(name = 'Oxygen\nTension', labels = c('21%', '5%', 'Harvest'), values = c('#005030', '#F47321', '#500c00'))

gene.expression.unloaded.plot.HIF1a <- ggplot(gene.expression.unloaded, aes(day, -ADAMTS4, fill=oxygen)) +
  ggtitle('HIF1a') +
  labs(x = 'Time (days)',
       y = expression(paste('Gene Expression (-', Delta,'Cq)'))) +
  geom_boxplot(position = position_dodge()) +
  scale_fill_manual(name = 'Oxygen\nTension', labels = c('21%', '5%', 'Harvest'), values = c('#005030', '#F47321', '#500c00'))

gene.expression.unloaded.plot.MMP13 <- ggplot(gene.expression.unloaded, aes(day, -MMP13, fill=oxygen)) +
  ggtitle('MMP 13') +
  labs(x = 'Time (days)',
       y = expression(paste('Gene Expression (-', Delta,'Cq)'))) +
  geom_boxplot(position = position_dodge()) +
  scale_fill_manual(name = 'Oxygen\nTension', labels = c('21%', '5%', 'Harvest'), values = c('#005030', '#F47321', '#500c00'))

gene.expression.unloaded.plot.MMP3 <- ggplot(gene.expression.unloaded, aes(day, -MMP3, fill=oxygen)) +
  ggtitle('MMP 3') +
  labs(x = 'Time (days)',
       y = expression(paste('Gene Expression (-', Delta,'Cq)'))) +
  geom_boxplot(position = position_dodge()) +
  scale_fill_manual(name = 'Oxygen\nTension', labels = c('21%', '5%', 'Harvest'), values = c('#005030', '#F47321', '#500c00'))

gene.expression.unloaded.plot.PACE4 <- ggplot(gene.expression.unloaded, aes(day, -PACE4, fill=oxygen)) +
  ggtitle('PACE 4') +
  labs(x = 'Time (days)',
       y = expression(paste('Gene Expression (-', Delta,'Cq)'))) +
  geom_boxplot(position = position_dodge()) +
  scale_fill_manual(name = 'Oxygen\nTension', labels = c('21%', '5%', 'Harvest'), values = c('#005030', '#F47321', '#500c00'))

gene.expression.unloaded.plot.TIMP1 <- ggplot(gene.expression.unloaded, aes(day, -TIMP1, fill=oxygen)) +
  ggtitle('TIMP 1') +
  labs(x = 'Time (days)',
       y = expression(paste('Gene Expression (-', Delta,'Cq)'))) +
  geom_boxplot(position = position_dodge()) +
  scale_fill_manual(name = 'Oxygen\nTension', labels = c('21%', '5%', 'Harvest'), values = c('#005030', '#F47321', '#500c00'))

gene.expression.unloaded.plot.TIMP2 <- ggplot(gene.expression.unloaded, aes(day, -TIMP2, fill=oxygen)) +
  ggtitle('TIMP 2') +
  labs(x = 'Time (days)',
       y = expression(paste('Gene Expression (-', Delta,'Cq)'))) +
  geom_boxplot(position = position_dodge()) +
  scale_fill_manual(name = 'Oxygen\nTension', labels = c('21%', '5%', 'Harvest'), values = c('#005030', '#F47321', '#500c00'))

gene.expression.unloaded.plot.TIMP3 <- ggplot(gene.expression.unloaded, aes(day, -TIMP3, fill=oxygen)) +
  ggtitle('TIMP 3') +
  labs(x = 'Time (days)',
       y = expression(paste('Gene Expression (-', Delta,'Cq)'))) +
  geom_boxplot(position = position_dodge()) +
  scale_fill_manual(name = 'Oxygen\nTension', labels = c('21%', '5%', 'Harvest'), values = c('#005030', '#F47321', '#500c00'))

gene.expression.unloaded.plot.T <- ggplot(gene.expression.unloaded, aes(day, -T, fill=oxygen)) +
  ggtitle('Brachyury') +
  labs(x = 'Time (days)',
       y = expression(paste('Gene Expression (-', Delta,'Cq)'))) +
  geom_boxplot(position = position_dodge()) +
  scale_fill_manual(name = 'Oxygen\nTension', labels = c('21%', '5%', 'Harvest'), values = c('#005030', '#F47321', '#500c00'))

plot_grid(gene.expression.loaded.plot.T, gene.expression.unloaded.plot.T,
          gene.expression.loaded.plot.HIF1a, gene.expression.unloaded.plot.HIF1a,
          ncol = 2,
          labels = 'AUTO')

plot_grid(gene.expression.loaded.plot.ADAMTS4, gene.expression.unloaded.plot.ADAMTS4,
          gene.expression.loaded.plot.ADAMTS5, gene.expression.unloaded.plot.ADAMTS5,
          gene.expression.loaded.plot.MMP3, gene.expression.unloaded.plot.MMP3,
          gene.expression.loaded.plot.MMP13, gene.expression.unloaded.plot.MMP13,
          ncol = 2,
          labels = 'AUTO')

plot_grid(gene.expression.loaded.plot.TIMP1, gene.expression.unloaded.plot.TIMP1,
          gene.expression.loaded.plot.TIMP2, gene.expression.unloaded.plot.TIMP2,
          gene.expression.loaded.plot.TIMP3, gene.expression.unloaded.plot.TIMP3,
          ncol = 2,
          labels = 'AUTO')

# gene.expression.loaded.5ox <- subset(gene.expression.loaded,
#                               oxygen == '5' 
#                               | (oxygen == 'x' & day == '0'))
# gene.expression.loaded.21ox <- subset(gene.expression.loaded,
#                               oxygen == '21' 
#                               | (oxygen == 'x' & day == '0'))
# gene.expression.unloaded.5ox <- subset(gene.expression.unloaded,
#                               oxygen == '5' 
#                               | (oxygen == 'x' & day == '0'))
# gene.expression.unloaded.21ox <- subset(gene.expression.unloaded,
#                               oxygen == '21' 
#                               | (oxygen == 'x' & day == '0'))
# 
# brachyury.loaded.5ox.plot <- ggplot(gene.expression.loaded.5ox, aes(x = day,y = -T, group = pig, color = pig)) + 
#   labs(x = 'Time (days)',
#        y = expression(paste('Gene Expression (-', Delta,'Cq)'))) +
#   geom_line() +
#   geom_point() 
# plot(brachyury.loaded.5ox.plot)
# 
# brachyury.loaded.21ox.plot <- ggplot(gene.expression.loaded.21ox, aes(x = day,y = -T, group = pig, color = pig)) + 
#   labs(x = 'Time (days)',
#        y = expression(paste('Gene Expression (-', Delta,'Cq)'))) +
#   geom_line() +
#   geom_point() 
# plot(brachyury.loaded.21ox.plot)
# 
# brachyury.unloaded.5ox.plot <- ggplot(gene.expression.unloaded.5ox, aes(x = day,y = -T, group = pig, color = pig)) + 
#   labs(x = 'Time (days)',
#        y = expression(paste('Gene Expression (-', Delta,'Cq)'))) +
#   geom_line() +
#   geom_point() 
# plot(brachyury.unloaded.5ox.plot)
# 
# brachyury.unloaded.21ox.plot <- ggplot(gene.expression.unloaded.21ox, aes(x = day,y = -T, group = pig, color = pig)) + 
#   labs(x = 'Time (days)',
#        y = expression(paste('Gene Expression (-', Delta,'Cq)'))) +
#   geom_line() +
#   geom_point() 
# plot(brachyury.unloaded.21ox.plot)
# 
# plot_grid(brachyury.loaded.5ox.plot,
#           brachyury.loaded.21ox.plot,
#           brachyury.unloaded.5ox.plot,
#           brachyury.unloaded.21ox.plot,
#           ncol = 2,
#           labels = 'AUTO')
# 
# hif1a.loaded.5ox.plot <- ggplot(gene.expression.loaded.5ox, aes(x = day,y = -HIF1a, group = pig, color = pig)) + 
#   labs(x = 'Time (days)',
#        y = expression(paste('Gene Expression (-', Delta,'Cq)'))) +
#   geom_line() +
#   geom_point() 
# plot(hif1a.loaded.5ox.plot)
# 
# hif1a.loaded.21ox.plot <- ggplot(gene.expression.loaded.21ox, aes(x = day,y = -HIF1a, group = pig, color = pig)) + 
#   labs(x = 'Time (days)',
#        y = expression(paste('Gene Expression (-', Delta,'Cq)'))) +
#   geom_line() +
#   geom_point() 
# plot(hif1a.loaded.21ox.plot)
# 
# hif1a.unloaded.5ox.plot <- ggplot(gene.expression.unloaded.5ox, aes(x = day,y = -HIF1a, group = pig, color = pig)) + 
#   labs(x = 'Time (days)',
#        y = expression(paste('Gene Expression (-', Delta,'Cq)'))) +
#   geom_line() +
#   geom_point() 
# plot(hif1a.unloaded.5ox.plot)
# 
# hif1a.unloaded.21ox.plot <- ggplot(gene.expression.unloaded.21ox, aes(x = day,y = -HIF1a, group = pig, color = pig)) + 
#   labs(x = 'Time (days)',
#        y = expression(paste('Gene Expression (-', Delta,'Cq)'))) +
#   geom_line() +
#   geom_point() 
# plot(hif1a.unloaded.21ox.plot)
# 
# plot_grid(hif1a.loaded.5ox.plot,
#           hif1a.loaded.21ox.plot,
#           hif1a.unloaded.5ox.plot,
#           hif1a.unloaded.21ox.plot,
#           ncol = 2,
#           labels = 'AUTO')
# 
# timp1.loaded.5ox.plot <- ggplot(gene.expression.loaded.5ox, aes(x = day,y = -TIMP1, group = pig, color = pig)) + 
#   labs(x = 'Time (days)',
#        y = expression(paste('Gene Expression (-', Delta,'Cq)'))) +
#   geom_line() +
#   geom_point() 
# plot(timp1.loaded.5ox.plot)
# 
# timp1.loaded.21ox.plot <- ggplot(gene.expression.loaded.21ox, aes(x = day,y = -TIMP1, group = pig, color = pig)) + 
#   labs(x = 'Time (days)',
#        y = expression(paste('Gene Expression (-', Delta,'Cq)'))) +
#   geom_line() +
#   geom_point() 
# plot(timp1.loaded.21ox.plot)
# 
# timp1.unloaded.5ox.plot <- ggplot(gene.expression.unloaded.5ox, aes(x = day,y = -TIMP1, group = pig, color = pig)) + 
#   labs(x = 'Time (days)',
#        y = expression(paste('Gene Expression (-', Delta,'Cq)'))) +
#   geom_line() +
#   geom_point() 
# plot(timp1.unloaded.5ox.plot)
# 
# timp1.unloaded.21ox.plot <- ggplot(gene.expression.unloaded.21ox, aes(x = day,y = -TIMP1, group = pig, color = pig)) + 
#   labs(x = 'Time (days)',
#        y = expression(paste('Gene Expression (-', Delta,'Cq)'))) +
#   geom_line() +
#   geom_point() 
# plot(timp1.unloaded.21ox.plot)
# 
# plot_grid(timp1.loaded.5ox.plot,
#           timp1.loaded.21ox.plot,
#           timp1.unloaded.5ox.plot,
#           timp1.unloaded.21ox.plot,
#           ncol = 2,
#           labels = 'AUTO')
# 
# timp2.loaded.5ox.plot <- ggplot(gene.expression.loaded.5ox, aes(x = day,y = -TIMP2, group = pig, color = pig)) + 
#   labs(x = 'Time (days)',
#        y = expression(paste('Gene Expression (-', Delta,'Cq)'))) +
#   geom_line() +
#   geom_point() 
# plot(timp2.loaded.5ox.plot)
# 
# timp2.loaded.21ox.plot <- ggplot(gene.expression.loaded.21ox, aes(x = day,y = -TIMP2, group = pig, color = pig)) + 
#   labs(x = 'Time (days)',
#        y = expression(paste('Gene Expression (-', Delta,'Cq)'))) +
#   geom_line() +
#   geom_point() 
# plot(timp2.loaded.21ox.plot)
# 
# timp2.unloaded.5ox.plot <- ggplot(gene.expression.unloaded.5ox, aes(x = day,y = -TIMP2, group = pig, color = pig)) + 
#   labs(x = 'Time (days)',
#        y = expression(paste('Gene Expression (-', Delta,'Cq)'))) +
#   geom_line() +
#   geom_point() 
# plot(timp2.unloaded.5ox.plot)
# 
# timp2.unloaded.21ox.plot <- ggplot(gene.expression.unloaded.21ox, aes(x = day,y = -TIMP2, group = pig, color = pig)) + 
#   labs(x = 'Time (days)',
#        y = expression(paste('Gene Expression (-', Delta,'Cq)'))) +
#   geom_line() +
#   geom_point() 
# plot(timp2.unloaded.21ox.plot)
# 
# plot_grid(timp2.loaded.5ox.plot,
#           timp2.loaded.21ox.plot,
#           timp2.unloaded.5ox.plot,
#           timp2.unloaded.21ox.plot,
#           ncol = 2,
#           labels = 'AUTO')
# 
# timp3.loaded.5ox.plot <- ggplot(gene.expression.loaded.5ox, aes(x = day,y = -TIMP3, group = pig, color = pig)) + 
#   labs(x = 'Time (days)',
#        y = expression(paste('Gene Expression (-', Delta,'Cq)'))) +
#   geom_line() +
#   geom_point() 
# plot(timp3.loaded.5ox.plot)
# 
# timp3.loaded.21ox.plot <- ggplot(gene.expression.loaded.21ox, aes(x = day,y = -TIMP3, group = pig, color = pig)) + 
#   labs(x = 'Time (days)',
#        y = expression(paste('Gene Expression (-', Delta,'Cq)'))) +
#   geom_line() +
#   geom_point() 
# plot(timp3.loaded.21ox.plot)
# 
# timp3.unloaded.5ox.plot <- ggplot(gene.expression.unloaded.5ox, aes(x = day,y = -TIMP3, group = pig, color = pig)) + 
#   labs(x = 'Time (days)',
#        y = expression(paste('Gene Expression (-', Delta,'Cq)'))) +
#   geom_line() +
#   geom_point() 
# plot(timp3.unloaded.5ox.plot)
# 
# timp3.unloaded.21ox.plot <- ggplot(gene.expression.unloaded.21ox, aes(x = day,y = -TIMP3, group = pig, color = pig)) + 
#   labs(x = 'Time (days)',
#        y = expression(paste('Gene Expression (-', Delta,'Cq)'))) +
#   geom_line() +
#   geom_point() 
# plot(timp3.unloaded.21ox.plot)
# 
# plot_grid(timp3.loaded.5ox.plot,
#           timp3.loaded.21ox.plot,
#           timp3.unloaded.5ox.plot,
#           timp3.unloaded.21ox.plot,
#           ncol = 2,
#           labels = 'AUTO')
# 
# mmp3.loaded.5ox.plot <- ggplot(gene.expression.loaded.5ox, aes(x = day,y = -MMP3, group = pig, color = pig)) + 
#   labs(x = 'Time (days)',
#        y = expression(paste('Gene Expression (-', Delta,'Cq)'))) +
#   geom_line() +
#   geom_point() 
# plot(mmp3.loaded.5ox.plot)
# 
# mmp3.loaded.21ox.plot <- ggplot(gene.expression.loaded.21ox, aes(x = day,y = -MMP3, group = pig, color = pig)) + 
#   labs(x = 'Time (days)',
#        y = expression(paste('Gene Expression (-', Delta,'Cq)'))) +
#   geom_line() +
#   geom_point() 
# plot(mmp3.loaded.21ox.plot)
# 
# mmp3.unloaded.5ox.plot <- ggplot(gene.expression.unloaded.5ox, aes(x = day,y = -MMP3, group = pig, color = pig)) + 
#   labs(x = 'Time (days)',
#        y = expression(paste('Gene Expression (-', Delta,'Cq)'))) +
#   geom_line() +
#   geom_point() 
# plot(mmp3.unloaded.5ox.plot)
# 
# mmp3.unloaded.21ox.plot <- ggplot(gene.expression.unloaded.21ox, aes(x = day,y = -MMP3, group = pig, color = pig)) + 
#   labs(x = 'Time (days)',
#        y = expression(paste('Gene Expression (-', Delta,'Cq)'))) +
#   geom_line() +
#   geom_point() 
# plot(mmp3.unloaded.21ox.plot)
# 
# plot_grid(mmp3.loaded.5ox.plot,
#           mmp3.loaded.21ox.plot,
#           mmp3.unloaded.5ox.plot,
#           mmp3.unloaded.21ox.plot,
#           ncol = 2,
#           labels = 'AUTO')
# 
# mmp13.loaded.5ox.plot <- ggplot(gene.expression.loaded.5ox, aes(x = day,y = -MMP13, group = pig, color = pig)) + 
#   labs(x = 'Time (days)',
#        y = expression(paste('Gene Expression (-', Delta,'Cq)'))) +
#   geom_line() +
#   geom_point() 
# plot(mmp13.loaded.5ox.plot)
# 
# mmp13.loaded.21ox.plot <- ggplot(gene.expression.loaded.21ox, aes(x = day,y = -MMP13, group = pig, color = pig)) + 
#   labs(x = 'Time (days)',
#        y = expression(paste('Gene Expression (-', Delta,'Cq)'))) +
#   geom_line() +
#   geom_point() 
# plot(mmp13.loaded.21ox.plot)
# 
# mmp13.unloaded.5ox.plot <- ggplot(gene.expression.unloaded.5ox, aes(x = day,y = -MMP13, group = pig, color = pig)) + 
#   labs(x = 'Time (days)',
#        y = expression(paste('Gene Expression (-', Delta,'Cq)'))) +
#   geom_line() +
#   geom_point() 
# plot(mmp13.unloaded.5ox.plot)
# 
# mmp13.unloaded.21ox.plot <- ggplot(gene.expression.unloaded.21ox, aes(x = day,y = -MMP13, group = pig, color = pig)) + 
#   labs(x = 'Time (days)',
#        y = expression(paste('Gene Expression (-', Delta,'Cq)'))) +
#   geom_line() +
#   geom_point() 
# plot(mmp13.unloaded.21ox.plot)
# 
# plot_grid(mmp13.loaded.5ox.plot,
#           mmp13.loaded.21ox.plot,
#           mmp13.unloaded.5ox.plot,
#           mmp13.unloaded.21ox.plot,
#           ncol = 2,
#           labels = 'AUTO')
# 
# adamts4.loaded.5ox.plot <- ggplot(gene.expression.loaded.5ox, aes(x = day,y = -ADAMTS4, group = pig, color = pig)) + 
#   labs(x = 'Time (days)',
#        y = expression(paste('Gene Expression (-', Delta,'Cq)'))) +
#   geom_line() +
#   geom_point() 
# plot(adamts4.loaded.5ox.plot)
# 
# adamts4.loaded.21ox.plot <- ggplot(gene.expression.loaded.21ox, aes(x = day,y = -ADAMTS4, group = pig, color = pig)) + 
#   labs(x = 'Time (days)',
#        y = expression(paste('Gene Expression (-', Delta,'Cq)'))) +
#   geom_line() +
#   geom_point() 
# plot(adamts4.loaded.21ox.plot)
# 
# adamts4.unloaded.5ox.plot <- ggplot(gene.expression.unloaded.5ox, aes(x = day,y = -ADAMTS4, group = pig, color = pig)) + 
#   labs(x = 'Time (days)',
#        y = expression(paste('Gene Expression (-', Delta,'Cq)'))) +
#   geom_line() +
#   geom_point() 
# plot(adamts4.unloaded.5ox.plot)
# 
# adamts4.unloaded.21ox.plot <- ggplot(gene.expression.unloaded.21ox, aes(x = day,y = -ADAMTS4, group = pig, color = pig)) + 
#   labs(x = 'Time (days)',
#        y = expression(paste('Gene Expression (-', Delta,'Cq)'))) +
#   geom_line() +
#   geom_point() 
# plot(adamts4.unloaded.21ox.plot)
# 
# plot_grid(adamts4.loaded.5ox.plot,
#           adamts4.loaded.21ox.plot,
#           adamts4.unloaded.5ox.plot,
#           adamts4.unloaded.21ox.plot,
#           ncol = 2,
#           labels = 'AUTO')
# 
# adamts5.loaded.5ox.plot <- ggplot(gene.expression.loaded.5ox, aes(x = day, y = -ADAMTS5, group = pig, color = pig)) + 
#   labs(x = 'Time (days)',
#        y = expression(paste('Gene Expression (-', Delta,'Cq)'))) +
#   geom_line() +
#   geom_point() 
# plot(adamts5.loaded.5ox.plot)
# 
# adamts5.loaded.21ox.plot <- ggplot(gene.expression.loaded.21ox, aes(x = day, y = -ADAMTS5, group = pig, color = pig)) + 
#   labs(x = 'Time (days)',
#        y = expression(paste('Gene Expression (-', Delta,'Cq)'))) +
#   geom_line() +
#   geom_point() 
# plot(adamts5.loaded.21ox.plot)
# 
# adamts5.unloaded.5ox.plot <- ggplot(gene.expression.unloaded.5ox, aes(x = day, y = -ADAMTS5, group = pig, color = pig)) + 
#   labs(x = 'Time (days)',
#        y = expression(paste('Gene Expression (-', Delta,'Cq)'))) +
#   geom_line() +
#   geom_point() 
# plot(adamts5.unloaded.5ox.plot)
# 
# adamts5.unloaded.21ox.plot <- ggplot(gene.expression.unloaded.21ox, aes(x = day, y = -ADAMTS5, group = pig, color = pig)) + 
#   labs(x = 'Time (days)',
#        y = expression(paste('Gene Expression (-', Delta,'Cq)'))) +
#   geom_line() +
#   geom_point() 
# plot(adamts5.unloaded.21ox.plot)
# 
# plot_grid(adamts5.loaded.5ox.plot,
#           adamts5.loaded.21ox.plot,
#           adamts5.unloaded.5ox.plot,
#           adamts5.unloaded.21ox.plot,
#           ncol = 2,
#           labels = 'AUTO')
# 
# pace4.loaded.5ox.plot <- ggplot(gene.expression.loaded.5ox, aes(x = day,y = -PACE4, group = pig, color = pig)) + 
#   labs(x = 'Time (days)',
#        y = expression(paste('Gene Expression (-', Delta,'Cq)'))) +
#   geom_line() +
#   geom_point() 
# plot(pace4.loaded.5ox.plot)
# 
# pace4.loaded.21ox.plot <- ggplot(gene.expression.loaded.21ox, aes(x = day,y = -PACE4, group = pig, color = pig)) + 
#   labs(x = 'Time (days)',
#        y = expression(paste('Gene Expression (-', Delta,'Cq)'))) +
#   geom_line() +
#   geom_point() 
# plot(pace4.loaded.21ox.plot)
# 
# pace4.unloaded.5ox.plot <- ggplot(gene.expression.unloaded.5ox, aes(x = day,y = -PACE4, group = pig, color = pig)) + 
#   labs(x = 'Time (days)',
#        y = expression(paste('Gene Expression (-', Delta,'Cq)'))) +
#   geom_line() +
#   geom_point() 
# plot(pace4.unloaded.5ox.plot)
# 
# pace4.unloaded.21ox.plot <- ggplot(gene.expression.unloaded.21ox, aes(x = day,y = -PACE4, group = pig, color = pig)) + 
#   labs(x = 'Time (days)',
#        y = expression(paste('Gene Expression (-', Delta,'Cq)'))) +
#   geom_line() +
#   geom_point() 
# plot(pace4.unloaded.21ox.plot)
# 
# plot_grid(pace4.loaded.5ox.plot,
#           pace4.loaded.21ox.plot,
#           pace4.unloaded.5ox.plot,
#           pace4.unloaded.21ox.plot,
#           ncol = 2,
#           labels = 'AUTO')
# 
