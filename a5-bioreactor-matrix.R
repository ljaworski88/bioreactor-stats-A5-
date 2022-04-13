# A script used to calculate the differences between actual values of
# extracellular matrix components, specfically GAG content, using the DMMB
# assay, and collagen, using the hydroxyproline.

# Author <- Lukas Jaworski, l.jaworski@umiami.edu
# Institution <- University of Miami
# PI <- Alicia Jackson, a.jackson2@miami.edu

# Copyright University of Miami 2019
###############################################################################
## Library imports
###############################################################################

library(ggplot2)
library(ggsignif)
library(cowplot)
library(dplyr)

###############################################################################
## Functions
###############################################################################

facet_signif <- function(plot, mapping){
  #add in significance bars on a faceted plot.
  #plot: a ggplot2 object
  #mapping: which of the facets each of the significance bars goes to
  #         it is a list containing a vector for each facet with the
  #         position of the annotations in ggsignif
  
  # example: 
  # my.plot <- ggplot(test.data, aes(time, dependant.variable))+
  #   geom_boxplot()+
  #   facet_grid(. ~ grp) +
  #   ## it takes these significance annotations ##
  #   geom_signif(annotations = c('*','*', '**', '*'),
  #               y_position = c(-4.5, -1.5, -2, -1), 
  #               ## ie y_position[1] = -4.5, y_position[2] = -1.5, etc.
  #               xmin = c(1, 1, 1, 1),
  #               xmax = c(3, 3, 2, 3), 
  #               vjust = 0.5) +
  #
  ## the significance mapping shows where each annotation will go ##
  #
  # my.plot <- facet_signif(my.plot, list(c(1),c(2),c(3,4)))
  #
  ## The plot has three facets, the first on recieves annotation 1,
  ## the second facet receives the annotation 2 and the thrid facet
  ## receives annotations 3 and 4. 
  
  
  #pull out the plot data
  plot.info <- ggplot_build(plot)
  signif.drawing <- plot.info$data[[2]]
  new.drawing <- data.frame()
  new.drawing <- merge(signif.drawing, new.drawing)
  # took too long to comment this so there's some black magic below
  # it adds to each facet the number position of the annotations added 
  # for ggsignif
  for (i in 1:length(mapping)){
    pv <- signif.drawing[((signif.drawing$group %in% mapping[[i]]) &
                            (signif.drawing$PANEL == i)), ]
    
    new.drawing <- merge(new.drawing, 
                         pv, 
                         all.x = T, 
                         all.y = T)
  }
  plot.info$data[[2]] <- new.drawing
  return(ggplot_gtable(plot.info))
}

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
    aov.summary[[x]][[1]] <- temp.summary[temp.summary['Pr(>F)'] <= 0.05 &
                                            !is.na(temp.summary['Pr(>F)']),
                                          'Pr(>F)',
                                          drop = F]
    if (nrow(aov.summary[[x]][[1]])>0){
      temp.tukey <- TukeyHSD(anova.results[[x]])
      aov.summary[[x]][[2]] <- lapply(temp.tukey,
                                      function(x) x[x[, 'p adj'] <= 0.05, 'p adj'])
    }
  }
  return(aov.summary)
}
###############################################################################
## Script
###############################################################################

###############################################################################
# Admnistrative Setup
###############################################################################

# Set the working directory to where the data is kept, 
# make it workstation agnostic
if (Sys.info()[['sysname']]== 'Linux'){
  working.directory <- paste0('/home/',Sys.info()[['user']],
                              '/Dropbox/Lab/Experiments/Matrix_Characterization/',
                              'a5-bioreactor-ox_lvl')
} else {
  working.directory <- paste0('C:/Users/',Sys.info()[['user']],
                              '/Dropbox/Lab/Experiments/Matrix_Characterization/',
                              'a5-bioreactor-ox_lvl')
}
working.directory <- paste0('P:/Dropbox/Lab/Experiments/Matrix_Characterization/',
                            'a5-bioreactor-ox_lvl')
setwd(working.directory)

###############################################################################
# Global Data imports
###############################################################################

# Import the the file which contains the sample ids and matches them to
# experiemntal conditions
sample.ids <- read.csv('tissue_weight.csv')

###############################################################################
# DMMB Data import and cleaning
###############################################################################

# Retrieve the dmmb data processed by the python script and make it one script

dmmb.set1 <- read.csv('DMMB/dmmb_b01-b23.csv')
dmmb.set2 <- read.csv('DMMB/dmmb_b24-b46.csv')
dmmb.set3 <- read.csv('DMMB/dmmb_b47-b71.csv')
dmmb.set4 <- read.csv('DMMB/dmmb_b72-b96.csv')
dmmb.set5 <- read.csv('DMMB/dmmb_b97-b100.csv')
dmmb.data <- rbind(dmmb.set1, dmmb.set2, dmmb.set3, dmmb.set4, dmmb.set5)

# Replacing the data values from the original which made no physical sense
# such as: negative concentration values

dmmb.repeats <- read.csv('DMMB/dmmb_repeats.csv')
repeats.index <- match(dmmb.repeats$sample, dmmb.data$sample)
dmmb.data[repeats.index,] <- dmmb.repeats

# Merge the DMMB data with the id dataframe

dmmb.data.w.id <- merge(sample.ids, dmmb.data,
                        by.x = 'id', by.y = 'sample')

# The following equation was used to convert the concentration of GAGs found
# using the DMMB assay to the liturature concentration units of
# ug GAG/mg tissue dry weight
# assay_concentration (ug/mL) * volume_of_Papain_used (mL) / [ tissue_dry_weight (g) * 1000 (mg/g) ]

dmmb.data.w.id$gag_fraction <- (dmmb.data.w.id$concentration * dmmb.data.w.id$papain_volume) / (dmmb.data.w.id$dry_wt * 1000)

# the variables 'pig' and 'day' need to be changed from integers to factors

dmmb.data.w.id[c('pig', 'day')] <- lapply(dmmb.data.w.id[c('pig', 'day')], factor)

# Label day 0 loading as its own factor category

levels(dmmb.data.w.id$loading) <- c(levels(dmmb.data.w.id$loading), 's')
dmmb.data.w.id[dmmb.data.w.id$loading == 'x' & dmmb.data.w.id$day == '0', "loading"] <- 's'

# Dropping data columns that are not useful for analysis

dmmb.data.w.id <- dmmb.data.w.id[ , !(names(dmmb.data.w.id) %in% c('disc',
                                                                   'empty_epi',
                                                                   'wet_epi',
                                                                   'dry_epi',
                                                                   'wet_wt',
                                                                   'dry_fraction',
                                                                   'units',
                                                                   'avg_reading',
                                                                   'reading'))]

###############################################################################
# DMMB Analysis
###############################################################################

# Pull the NP data out
np.dmmb.data <- subset(dmmb.data.w.id, tissue == 'NP')

# Remove sample 33 where the tissue weight is negative
np.dmmb.data <- np.dmmb.data[!(np.dmmb.data$id == 'b33'),]

# Perform 3-way ANOVA
np.dmmb.anova <- aov(gag_fraction ~ oxygen * day * loading + pig,
                      data = np.dmmb.data)

summary.aov(np.dmmb.anova)

# Separate the loaded and unloaded sample datasets
# This is being done because there might be an artifact by comparing both sets
# to the day 0 disc
np.dmmb.loaded.data <- subset(np.dmmb.data,
                              loading == '1' 
                              | loading == 's')
np.dmmb.unloaded.data <- subset(np.dmmb.data,
                                loading == 'x'
                                | loading == 's')

# # Now split the data based on oxygen tension as well for graphing
# np.dmmb.loaded.5ox.data <- subset(np.dmmb.loaded.data,
#                               oxygen == '5'
#                               | (oxygen == 'x' & day == '0'))
# np.dmmb.loaded.21ox.data <- subset(np.dmmb.loaded.data,
#                               oxygen == '21'
#                               | (oxygen == 'x' & day == '0'))
# np.dmmb.unloaded.5ox.data <- subset(np.dmmb.unloaded.data,
#                               oxygen == '5'
#                               | (oxygen == 'x' & day == '0'))
# np.dmmb.unloaded.21ox.data <- subset(np.dmmb.unloaded.data,
#                               oxygen == '21'
#                               | (oxygen == 'x' & day == '0'))

# # Plot analysis results as Line charts due to low n
# np.dmmb.plot.loaded.5ox.days <- ggplot(np.dmmb.loaded.5ox.data, aes(x = day,y = gag_fraction, group = pig, color = pig)) + 
#   labs(x = 'Time (days)',
#        y = expression('Proteoglycan Content ('*frac(mu*'g'['GAG'], 'mg'['dry weight'])*')')) +
#   geom_line() +
#   geom_point() 
# plot(np.dmmb.plot.loaded.5ox.days)
# 
# np.dmmb.plot.loaded.21ox.days <- ggplot(np.dmmb.loaded.21ox.data, aes(x = day,y = gag_fraction, group = pig, color = pig)) + 
#   labs(x = 'Time (days)',
#        y = expression('Proteoglycan Content ('*frac(mu*'g'['GAG'], 'mg'['dry weight'])*')')) +
#   geom_line() +
#   geom_point() 
# plot(np.dmmb.plot.loaded.21ox.days)
# 
# np.dmmb.plot.unloaded.5ox.days <- ggplot(np.dmmb.unloaded.5ox.data, aes(x = day,y = gag_fraction, group = pig, color = pig)) + 
#   labs(x = 'Time (days)',
#        y = expression('Proteoglycan Content ('*frac(mu*'g'['GAG'], 'mg'['dry weight'])*')')) +
#   geom_line() +
#   geom_point() 
# plot(np.dmmb.plot.unloaded.5ox.days)
# 
# np.dmmb.plot.unloaded.21ox.days <- ggplot(np.dmmb.unloaded.21ox.data, aes(x = day,y = gag_fraction, group = pig, color = pig)) + 
#   labs(x = 'Time (days)',
#        y = expression('Proteoglycan Content ('*frac(mu*'g'['GAG'], 'mg'['dry weight'])*')')) +
#   geom_line() +
#   geom_point() 
# plot(np.dmmb.plot.unloaded.21ox.days)
# 
# plot_grid(np.dmmb.plot.loaded.5ox.days, 
#           np.dmmb.plot.loaded.21ox.days, 
#           np.dmmb.plot.unloaded.5ox.days, 
#           np.dmmb.plot.unloaded.21ox.days,
#           ncol = 2,
#           labels = 'AUTO')

## Plot NP DMMB graphs

np.dmmb.loaded.summary <- np.dmmb.loaded.data %>%
  group_by(day, oxygen) %>%
  summarise(mean_dmmb = mean(gag_fraction),
             se_dmmb = sd(gag_fraction) / sqrt(n()))

np.dmmb.unloaded.summary <- np.dmmb.unloaded.data %>%
  group_by(day, oxygen) %>%
  summarise(mean_dmmb = mean(gag_fraction),
            se_dmmb = sd(gag_fraction) / sqrt(n()))

np.dmmb.plot.loaded.days <- ggplot(np.dmmb.loaded.summary, aes(day, mean_dmmb, fill=oxygen)) +
  ggtitle('NP') +
  labs(x = 'Time (days)',
       y = expression('Proteoglycan Content ('*frac(mu*'g'['GAG'], 'mg'['dry weight'])*')')) +
  geom_col(position = position_dodge(),
           color = 'black') +
  geom_errorbar(aes(ymin=mean_dmmb, ymax=mean_dmmb+se_dmmb),
                position = position_dodge()) +
  scale_fill_manual(name = 'Oxygen\nTension', labels = c('21%', '5%', 'Harvest'), values = c('#005030', '#F47321', '#500c00'))


np.dmmb.plot.unloaded.days <- ggplot(np.dmmb.unloaded.summary, aes(day, mean_dmmb, fill=oxygen)) +
  ggtitle('NP') +
  labs(x = 'Time (days)',
       y = expression('Proteoglycan Content ('*frac(mu*'g'['GAG'], 'mg'['dry weight'])*')')) +
  geom_col(position = position_dodge(),
           color = 'black') +
  geom_errorbar(aes(ymin=mean_dmmb, ymax=mean_dmmb+se_dmmb),
                position = position_dodge()) +
  scale_fill_manual(name = 'Oxygen\nTension', labels = c('21%', '5%', 'Harvest'), values = c('#005030', '#F47321', '#500c00'))


plot_grid(np.dmmb.plot.loaded.days,
          np.dmmb.plot.unloaded.days,
          ncol = 1,
          labels = 'AUTO')

# np.dmmb.plot.loaded.oxygen <- ggplot(np.dmmb.loaded.data, aes(oxygen, gag_fraction)) +
#   labs(x = 'Oxygen Tension (%)',
#        y = expression('Proteoglycan Content ('*frac(mu*'g'['GAG'], 'mg'['dry weight'])*')')) +
#   geom_col()
# 
# np.dmmb.plot.unloaded.days <- ggplot(np.dmmb.unloaded.data, aes(day, gag_fraction)) +
#   labs(x = 'Time (days)',
#        y = expression('Proteoglycan Content ('*frac(mu*'g'['GAG'], 'mg'['dry weight'])*')')) +
#   geom_col()
# 
# np.dmmb.plot.unloaded.oxygen <- ggplot(np.dmmb.unloaded.data, aes(oxygen, gag_fraction)) +
#   labs(x = 'Oxygen Tension (%)',
#        y = expression('Proteoglycan Content ('*frac(mu*'g'['GAG'], 'mg'['dry weight'])*')')) +
#   geom_col()

### Start the analysis on the AF tissue ###
af.dmmb.data <- subset(dmmb.data.w.id,
                       tissue == 'AF')

# Remove sample 80 where the tissue weight is negative
af.dmmb.data <- af.dmmb.data[!(af.dmmb.data$id == 'b80'),]

# Perform three way ANOVA
af.dmmb.anova <- aov(gag_fraction ~ oxygen * day * loading + pig,
                     data = af.dmmb.data)

summary.aov(af.dmmb.anova)
# Separate the loaded and unloaded sample datasets
# This is being done because there might be an artifact by comparing both sets
# to the day 0 disc
af.dmmb.loaded.data <- subset(af.dmmb.data,
                              loading == '1' 
                              | loading == 's')
af.dmmb.unloaded.data <- subset(af.dmmb.data,
                                loading == 'x'
                                | loading == 's')

# # Now split the data based on oxygen tension as well for graphing
# af.dmmb.loaded.5ox.data <- subset(af.dmmb.loaded.data,
#                               oxygen == '5' 
#                               | (oxygen == 'x' & day == '0'))
# af.dmmb.loaded.21ox.data <- subset(af.dmmb.loaded.data,
#                               oxygen == '21' 
#                               | (oxygen == 'x' & day == '0'))
# af.dmmb.unloaded.5ox.data <- subset(af.dmmb.unloaded.data,
#                               oxygen == '5' 
#                               | (oxygen == 'x' & day == '0'))
# af.dmmb.unloaded.21ox.data <- subset(af.dmmb.unloaded.data,
#                               oxygen == '21' 
#                               | (oxygen == 'x' & day == '0'))

# # plotting the results
# af.dmmb.plot.loaded.5ox.days <- ggplot(af.dmmb.loaded.5ox.data, aes(x = day,y = gag_fraction, group = pig, color = pig)) + 
#   labs(x = 'Time (days)',
#        y = expression('Proteoglycan Content ('*frac(mu*'g'['GAG'], 'mg'['dry weight'])*')')) +
#   geom_line() +
#   geom_point() 
# plot(af.dmmb.plot.loaded.5ox.days)
# 
# af.dmmb.plot.loaded.21ox.days <- ggplot(af.dmmb.loaded.21ox.data, aes(x = day,y = gag_fraction, group = pig, color = pig)) + 
#   labs(x = 'Time (days)',
#        y = expression('Proteoglycan Content ('*frac(mu*'g'['GAG'], 'mg'['dry weight'])*')')) +
#   geom_line() +
#   geom_point() 
# plot(af.dmmb.plot.loaded.21ox.days)
# 
# af.dmmb.plot.unloaded.5ox.days <- ggplot(af.dmmb.unloaded.5ox.data, aes(x = day,y = gag_fraction, group = pig, color = pig)) + 
#   labs(x = 'Time (days)',
#        y = expression('Proteoglycan Content ('*frac(mu*'g'['GAG'], 'mg'['dry weight'])*')')) +
#   geom_line() +
#   geom_point() 
# plot(af.dmmb.plot.unloaded.5ox.days)
# 
# af.dmmb.plot.unloaded.21ox.days <- ggplot(af.dmmb.unloaded.21ox.data, aes(x = day,y = gag_fraction, group = pig, color = pig)) + 
#   labs(x = 'Time (days)',
#        y = expression('Proteoglycan Content ('*frac(mu*'g'['GAG'], 'mg'['dry weight'])*')')) +
#   geom_line() +
#   geom_point() 
# plot(af.dmmb.plot.unloaded.21ox.days)
# 
# plot_grid(af.dmmb.plot.loaded.5ox.days, 
#           af.dmmb.plot.loaded.21ox.days, 
#           af.dmmb.plot.unloaded.5ox.days, 
#           af.dmmb.plot.unloaded.21ox.days,
#           ncol = 2,
#           labels = 'AUTO')

## Plot AF DMMB Graphs ##
af.dmmb.loaded.summary <- af.dmmb.loaded.data %>%
  group_by(day, oxygen) %>%
  summarise(mean_dmmb = mean(gag_fraction),
            se_dmmb = sd(gag_fraction) / sqrt(n()))

af.dmmb.unloaded.summary <- af.dmmb.unloaded.data %>%
  group_by(day, oxygen) %>%
  summarise(mean_dmmb = mean(gag_fraction),
            se_dmmb = sd(gag_fraction) / sqrt(n()))

af.dmmb.plot.loaded.days <- ggplot(af.dmmb.loaded.summary, aes(day, mean_dmmb, fill=oxygen)) +
  ggtitle('AF') +
  labs(x = 'Time (days)',
       y = expression('Proteoglycan Content ('*frac(mu*'g'['GAG'], 'mg'['dry weight'])*')')) +
  geom_col(position = position_dodge(),
           color = 'black') +
  geom_errorbar(aes(ymin=mean_dmmb, ymax=mean_dmmb+se_dmmb),
                position = position_dodge()) +
  scale_fill_manual(name = 'Oxygen\nTension', labels = c('21%', '5%', 'Harvest'), values = c('#005030', '#F47321', '#500c00'))


af.dmmb.plot.unloaded.days <- ggplot(af.dmmb.unloaded.summary, aes(day, mean_dmmb, fill=oxygen)) +
  ggtitle('AF') +
  labs(x = 'Time (days)',
       y = expression('Proteoglycan Content ('*frac(mu*'g'['GAG'], 'mg'['dry weight'])*')')) +
  geom_col(position = position_dodge(),
           color = 'black') +
  geom_errorbar(aes(ymin=mean_dmmb, ymax=mean_dmmb+se_dmmb),
                position = position_dodge()) +
  scale_fill_manual(name = 'Oxygen\nTension', labels = c('21%', '5%', 'Harvest'), values = c('#005030', '#F47321', '#500c00'))


plot_grid(np.dmmb.plot.loaded.days, af.dmmb.plot.loaded.days,
          np.dmmb.plot.unloaded.days, af.dmmb.plot.unloaded.days,
          ncol = 2,
          labels = 'AUTO')
# 
# af.dmmb.plot.loaded.days <- ggplot(af.dmmb.loaded.data, aes(day, gag_fraction)) +
#   labs(x = 'Time (days)',
#        y = expression('Proteoglycan Content ('*frac(mu*'g'['GAG'], 'mg'['dry weight'])*')')) +
#   geom_col()
# 
# af.dmmb.plot.loaded.oxygen <- ggplot(af.dmmb.loaded.data, aes(oxygen, gag_fraction)) +
#   labs(x = 'Oxygen Tension (%)',
#        y = expression('Proteoglycan Content ('*frac(mu*'g'['GAG'], 'mg'['dry weight'])*')')) +
#   geom_col()
# 
# af.dmmb.plot.unloaded.days <- ggplot(af.dmmb.unloaded.data, aes(day, gag_fraction)) +
#   labs(x = 'Time (days)',
#        y = expression('Proteoglycan Content ('*frac(mu*'g'['GAG'], 'mg'['dry weight'])*')')) +
#   geom_col()
# 
# af.dmmb.plot.unloaded.oxygen <- ggplot(af.dmmb.unloaded.data, aes(oxygen, gag_fraction)) +
#   labs(x = 'Oxygen Tension (%)',
#        y = expression('Proteoglycan Content ('*frac(mu*'g'['GAG'], 'mg'['dry weight'])*')')) +
#   geom_col()
# 
# plot_grid(af.dmmb.plot.loaded.days,
#           af.dmmb.plot.loaded.days,
#           af.dmmb.plot.unloaded.days,
#           af.dmmb.plot.unloaded.days,
#           ncol = 2,
#           labels = 'AUTO')

###############################################################################
# Hydroxyproline Data import and cleaning
###############################################################################

# Import the hydroxyproline data
proline.set1 <- read.csv('Hydroxyproline/hydroxyproline_b01-b26.csv')
proline.set2 <- read.csv('Hydroxyproline/hydroxyproline_b27-b49.csv')
proline.set3 <- read.csv('Hydroxyproline/hydroxyproline_b53-b78.csv')
proline.set4 <- read.csv('Hydroxyproline/hydroxyproline_b79-b100.csv')
proline.data <- rbind(proline.set1, proline.set2, proline.set3, proline.set4)

# Replacing the data values from the original which made no physical sense
# such as: negative concentration values
proline.repeats <- read.csv('Hydroxyproline/hydroxyproline_repeats.csv')
repeats.index <- match(proline.repeats$sample, proline.data$sample)
proline.data[repeats.index,] <- proline.repeats

# Add the sample id data, previously imported for the DMMB assay
proline.data.w.id <- merge(sample.ids, proline.data,
                           by.x = 'id', by.y = 'sample')

# Transforming the data to to get proline content on a per mass dry basis (ug proline/mg dry_tissue)
proline.data.w.id$proline_content <- (proline.data.w.id$concentration * proline.data.w.id$papain_volume) / (proline.data.w.id$dry_wt * 1000)

# The pig number and day variable need to be changed to factors
proline.data.w.id[c('pig', 'day')] <- lapply(proline.data.w.id[c('pig', 'day')], factor)

# Label day 0 loading as its own factor category

levels(proline.data.w.id$loading) <- c(levels(proline.data.w.id$loading), 's')
proline.data.w.id[proline.data.w.id$loading == 'x' & proline.data.w.id$day == '0', "loading"] <- 's'

# Dropping data columns that are not useful for analysis
proline.data.w.id <- proline.data.w.id[ , !(names(proline.data.w.id) %in% c('disc',
                                                                            'empty_epi',
                                                                            'wet_epi',
                                                                            'dry_epi',
                                                                            'wet_wt',
                                                                            'dry_fraction',
                                                                            'units',
                                                                            'avg_reading',
                                                                            'reading'))]

###############################################################################
# Hydroxyproline Analysis
###############################################################################

# Separate the NP data
np.proline.data <- subset(proline.data.w.id, tissue == 'NP')
# Remove sample 33 where the tissue weight is negative
np.proline.data <- np.proline.data[!(np.proline.data$id == 'b33'),]

# Perform 3 way ANOVA
np.proline.anova <- aov(proline_content ~ oxygen * day * loading + pig,
                     data = np.proline.data)

summary.aov(np.proline.anova)

# Separate loaded and unloaded datasets to remove arificial inclusion of day 0
# in unloaded group
# Then perform the analysis and report the results

np.proline.loaded.data <- subset(np.proline.data, loading == '1' | loading == 's')
np.proline.unloaded.data <- subset(np.proline.data, loading == 'x' | loading == 's')

np.proline.loaded.summary <- np.proline.loaded.data %>%
  group_by(day, oxygen) %>%
  summarise(mean_proline = mean(proline_content),
            se_proline = sd(proline_content) / sqrt(n()))

np.proline.unloaded.summary <- np.proline.unloaded.data %>%
  group_by(day, oxygen) %>%
  summarise(mean_proline = mean(proline_content),
            se_proline = sd(proline_content) / sqrt(n()))

np.proline.plot.loaded.days <- ggplot(np.proline.loaded.summary, aes(day, mean_proline, fill=oxygen)) +
  ggtitle('NP') +
  labs(x = 'Time (days)',
       y = expression('Hydroxyproline Content ('*frac(mu*'g'['Hydroxyproline'], 'mg'['dry weight'])*')')) +
  geom_col(position = position_dodge(),
           color = 'black') +
  geom_errorbar(aes(ymin=mean_proline, ymax=mean_proline+se_proline),
                position = position_dodge()) +
  scale_fill_manual(name = 'Oxygen\nTension', labels = c('21%', '5%', 'Harvest'), values = c('#005030', '#F47321', '#500c00'))


np.proline.plot.unloaded.days <- ggplot(np.proline.unloaded.summary, aes(day, mean_proline, fill=oxygen)) +
  ggtitle('NP') +
  labs(x = 'Time (days)',
       y = expression('Hydroxyproline Content ('*frac(mu*'g'['Hydroxyproline'], 'mg'['dry weight'])*')')) +
  geom_col(position = position_dodge(),
           color = 'black') +
  geom_errorbar(aes(ymin=mean_proline, ymax=mean_proline+se_proline),
                position = position_dodge()) +
  scale_fill_manual(name = 'Oxygen\nTension', labels = c('21%', '5%', 'Harvest'), values = c('#005030', '#F47321', '#500c00'))


plot_grid(np.proline.plot.loaded.days,
          np.proline.plot.unloaded.days,
          ncol = 1,
          labels = 'AUTO')

# Now split the data based on oxygen tension as well for graphing
# np.proline.loaded.5ox.data <- subset(np.proline.loaded.data,
#                               oxygen == '5' 
#                               | (oxygen == 'x' & day == '0'))
# np.proline.loaded.21ox.data <- subset(np.proline.loaded.data,
#                               oxygen == '21' 
#                               | (oxygen == 'x' & day == '0'))
# np.proline.unloaded.5ox.data <- subset(np.proline.unloaded.data,
#                               oxygen == '5' 
#                               | (oxygen == 'x' & day == '0'))
# np.proline.unloaded.21ox.data <- subset(np.proline.unloaded.data,
#                               oxygen == '21' 
#                               | (oxygen == 'x' & day == '0'))
# 
# np.proline.plot.loaded.5ox.days <- ggplot(np.proline.loaded.5ox.data, aes(x = day,y = proline_content, group = pig, color = pig)) + 
#   labs(x = 'Time (days)',
#        y = expression('Hydroxyproline Content ('*frac(mu*'g'['GAG'], 'mg'['dry weight'])*')')) +
#   geom_line() +
#   geom_point() 
# plot(np.proline.plot.loaded.5ox.days)
# 
# np.proline.plot.loaded.21ox.days <- ggplot(np.proline.loaded.21ox.data, aes(x = day,y = proline_content, group = pig, color = pig)) + 
#   labs(x = 'Time (days)',
#        y = expression('Hydroxyproline Content ('*frac(mu*'g'['GAG'], 'mg'['dry weight'])*')')) +
#   geom_line() +
#   geom_point() 
# plot(np.proline.plot.loaded.21ox.days)
# 
# np.proline.plot.unloaded.5ox.days <- ggplot(np.proline.unloaded.5ox.data, aes(x = day,y = proline_content, group = pig, color = pig)) + 
#   labs(x = 'Time (days)',
#        y = expression('Hydroxyproline Content ('*frac(mu*'g'['GAG'], 'mg'['dry weight'])*')')) +
#   geom_line() +
#   geom_point() 
# plot(np.proline.plot.unloaded.5ox.days)
# 
# np.proline.plot.unloaded.21ox.days <- ggplot(np.proline.unloaded.21ox.data, aes(x = day,y = proline_content, group = pig, color = pig)) + 
#   labs(x = 'Time (days)',
#        y = expression('Hydroxyproline Content ('*frac(mu*'g'['GAG'], 'mg'['dry weight'])*')')) +
#   geom_line() +
#   geom_point() 
# plot(np.proline.plot.unloaded.21ox.days)
# 
# plot_grid(np.proline.plot.loaded.5ox.days, 
#           np.proline.plot.loaded.21ox.days, 
#           np.proline.plot.unloaded.5ox.days, 
#           np.proline.plot.unloaded.21ox.days,
#           ncol = 2,
#           labels = 'AUTO')

# np.proline.plot.loaded.days <- ggplot(np.proline.loaded.data, aes(day, proline_content)) + 
#   labs(x = 'Time (days)',
#        y = expression('Hydroxyproline Content ('*frac(mu*'g'['Hydroxyproline'], 'mg'['dry weight'])*')')) +
#   geom_boxplot()
# 
# np.proline.plot.loaded.oxygen <- ggplot(np.proline.loaded.data, aes(oxygen, proline_content)) + 
#   labs(x = 'Oxygen Tesnion (%)',
#        y = expression('Hydroxyproline Content ('*frac(mu*'g'['Hydroxyproline'], 'mg'['dry weight'])*')')) +
#   geom_boxplot()

# np.proline.plot.unloaded.days <- ggplot(np.proline.unloaded.data, aes(day, proline_content)) + 
#   labs(x = 'Time (days)',
#        y = expression('Hydroxyproline Content ('*frac(mu*'g'['Hydroxyproline'], 'mg'['dry weight'])*')')) +
#   geom_boxplot() 
# 
# np.proline.plot.unloaded.oxygen <- ggplot(np.proline.unloaded.data, aes(oxygen, proline_content)) + 
#   labs(x = 'Oxygen Tension (%)',
#        y = expression('Hydroxyproline Content ('*frac(mu*'g'['Hydroxyproline'], 'mg'['dry weight'])*')')) +
#   geom_boxplot() 


# Pull out the AF data
af.proline.data <- subset(proline.data.w.id, tissue == 'AF')
# Remove sample 80 where the tissue weight is negative
af.proline.data <- af.proline.data[!(af.proline.data$id == 'b80'),]

# Doing the 3-way ANOVA
af.proline.anova <- aov(proline_content ~ oxygen * day * loading + pig,
                        data = af.proline.data)

summary.aov(af.proline.anova)
# Separate loaded and unloaded datasets to remove arificial inclusion of day 0
# in unloaded group
# Then perform the analysis and report the results

af.proline.loaded.data <- subset(af.proline.data, loading == '1' | loading == 's')
af.proline.unloaded.data <- subset(af.proline.data, loading == 'x' | loading == 's')

af.proline.loaded.summary <- af.proline.loaded.data %>%
  group_by(day, oxygen) %>%
  summarise(mean_proline = mean(proline_content),
            se_proline = sd(proline_content) / sqrt(n()))

af.proline.unloaded.summary <- af.proline.unloaded.data %>%
  group_by(day, oxygen) %>%
  summarise(mean_proline = mean(proline_content),
            se_proline = sd(proline_content) / sqrt(n()))

af.proline.plot.loaded.days <- ggplot(af.proline.loaded.summary, aes(day, mean_proline, fill=oxygen)) +
  ggtitle('AF') +
  labs(x = 'Time (days)',
       y = expression('Hydroxyproline Content ('*frac(mu*'g'['Hydroxyproline'], 'mg'['dry weight'])*')')) +
  geom_col(position = position_dodge(),
           color = 'black') +
  geom_errorbar(aes(ymin=mean_proline, ymax=mean_proline+se_proline),
                position = position_dodge()) +
  scale_fill_manual(name = 'Oxygen\nTension', labels = c('21%', '5%', 'Harvest'), values = c('#005030', '#F47321', '#500c00'))


af.proline.plot.unloaded.days <- ggplot(af.proline.unloaded.summary, aes(day, mean_proline, fill=oxygen)) +
  ggtitle('AF') +
  labs(x = 'Time (days)',
       y = expression('Hydroxyproline Content ('*frac(mu*'g'['Hydroxyproline'], 'mg'['dry weight'])*')')) +
  geom_col(position = position_dodge(),
           color = 'black') +
  geom_errorbar(aes(ymin=mean_proline, ymax=mean_proline+se_proline),
                position = position_dodge()) +
  scale_fill_manual(name = 'Oxygen\nTension', labels = c('21%', '5%', 'Harvest'), values = c('#005030', '#F47321', '#500c00'))


plot_grid(np.proline.plot.loaded.days, af.proline.plot.loaded.days,
          np.proline.plot.unloaded.days, af.proline.plot.unloaded.days,
          ncol = 2,
          labels = 'AUTO')


# # Now split the data based on oxygen tension as well for graphing
# af.proline.loaded.5ox.data <- subset(af.proline.loaded.data,
#                               oxygen == '5' 
#                               | (oxygen == 'x' & day == '0'))
# af.proline.loaded.21ox.data <- subset(af.proline.loaded.data,
#                               oxygen == '21' 
#                               | (oxygen == 'x' & day == '0'))
# af.proline.unloaded.5ox.data <- subset(af.proline.unloaded.data,
#                               oxygen == '5' 
#                               | (oxygen == 'x' & day == '0'))
# af.proline.unloaded.21ox.data <- subset(af.proline.unloaded.data,
#                               oxygen == '21' 
#                               | (oxygen == 'x' & day == '0'))
# 
# af.proline.plot.loaded.5ox.days <- ggplot(af.proline.loaded.5ox.data, aes(x = day,y = proline_content, group = pig, color = pig)) + 
#   labs(x = 'Time (days)',
#        y = expression('Hydroxyproline Content ('*frac(mu*'g'['GAG'], 'mg'['dry weight'])*')')) +
#   geom_line() +
#   geom_point() 
# plot(af.proline.plot.loaded.5ox.days)
# 
# af.proline.plot.loaded.21ox.days <- ggplot(af.proline.loaded.21ox.data, aes(x = day,y = proline_content, group = pig, color = pig)) + 
#   labs(x = 'Time (days)',
#        y = expression('Hydroxyproline Content ('*frac(mu*'g'['GAG'], 'mg'['dry weight'])*')')) +
#   geom_line() +
#   geom_point() 
# plot(af.proline.plot.loaded.21ox.days)
# 
# af.proline.plot.unloaded.5ox.days <- ggplot(af.proline.unloaded.5ox.data, aes(x = day,y = proline_content, group = pig, color = pig)) + 
#   labs(x = 'Time (days)',
#        y = expression('Hydroxyproline Content ('*frac(mu*'g'['GAG'], 'mg'['dry weight'])*')')) +
#   geom_line() +
#   geom_point() 
# plot(af.proline.plot.unloaded.5ox.days)
# 
# af.proline.plot.unloaded.21ox.days <- ggplot(af.proline.unloaded.21ox.data, aes(x = day,y = proline_content, group = pig, color = pig)) + 
#   labs(x = 'Time (days)',
#        y = expression('Hydroxyproline Content ('*frac(mu*'g'['GAG'], 'mg'['dry weight'])*')')) +
#   geom_line() +
#   geom_point()  
# plot(af.proline.plot.unloaded.21ox.days)
# 
# plot_grid(af.proline.plot.loaded.5ox.days, 
#           af.proline.plot.loaded.21ox.days, 
#           af.proline.plot.unloaded.5ox.days, 
#           af.proline.plot.unloaded.21ox.days,
#           ncol = 2,
#           labels = 'AUTO')
# 
# ## Graphs for after statistically significant quantity of pigs is done
# 
# af.proline.plot.loaded.days <- ggplot(af.proline.loaded.data, aes(day, proline_content)) +
#   labs(x = 'Time (days)',
#        y = expression('Hydroxyproline Content ('*frac(mu*'g'['Hydroxyproline'], 'mg'['dry weight'])*')')) +
#   geom_boxplot()
# 
# af.proline.plot.loaded.oxygen <- ggplot(af.proline.loaded.data, aes(oxygen, proline_content)) +
#   labs(x = 'Oxygen Tension (%)',
#        y = expression('Hydroxyproline Content ('*frac(mu*'g'['Hydroxyproline'], 'mg'['dry weight'])*')')) +
#   geom_boxplot()
# 
# af.proline.plot.unloaded.days <- ggplot(af.proline.unloaded.data, aes(day, proline_content)) +
#   labs(x = 'Time (days)',
#        y = expression('Hydroxyproline Content ('*frac(mu*'g'['Hydroxyproline'], 'mg'['dry weight'])*')')) +
#   geom_boxplot()
# 
# af.proline.plot.unloaded.oxygen <- ggplot(af.proline.unloaded.data, aes(oxygen, proline_content)) +
#   labs(x = 'Oxygen Tension (%)',
#        y = expression('Hydroxyproline Content ('*frac(mu*'g'['Hydroxyproline'], 'mg'['dry weight'])*')')) +
#   geom_boxplot()

###############################################################################
# H2O content analysis
###############################################################################

# Perform 3-way ANOVA
np.h2o.anova <- aov(water_percent ~ oxygen * day * loading + pig,
                     data = np.dmmb.data)

summary.aov(np.h2o.anova)

af.h2o.anova <- aov(water_percent ~ oxygen * day * loading + pig,
                    data = af.dmmb.data)

summary.aov(af.h2o.anova)

np.h2o.loaded.summary <- np.dmmb.loaded.data %>%
  group_by(day, oxygen) %>%
  summarise(mean_h2o = mean(water_percent),
             se_h2o = sd(water_percent) / sqrt(n()))

np.h2o.unloaded.summary <- np.dmmb.unloaded.data %>%
  group_by(day, oxygen) %>%
  summarise(mean_h2o = mean(water_percent),
            se_h2o = sd(water_percent) / sqrt(n()))

np.h2o.plot.loaded.days <- ggplot(np.h2o.loaded.summary, aes(day, mean_h2o, fill=oxygen)) +
  ggtitle('NP') +
  labs(x = 'Time (days)',
       y = 'Water Content (%)') +
  geom_col(position = position_dodge(),
           color = 'black') +
  geom_errorbar(aes(ymin=mean_h2o, ymax=mean_h2o+se_h2o),
                position = position_dodge()) +
  scale_fill_manual(name = 'Oxygen\nTension', labels = c('21%', '5%', 'Harvest'), values = c('#005030', '#F47321', '#500c00'))


np.h2o.plot.unloaded.days <- ggplot(np.h2o.unloaded.summary, aes(day, mean_h2o, fill=oxygen)) +
  ggtitle('NP') +
  labs(x = 'Time (days)',
       y = 'Water Content (%)') +
  geom_col(position = position_dodge(),
           color = 'black') +
  geom_errorbar(aes(ymin=mean_h2o, ymax=mean_h2o+se_h2o),
                position = position_dodge()) +
  scale_fill_manual(name = 'Oxygen\nTension', labels = c('21%', '5%', 'Harvest'), values = c('#005030', '#F47321', '#500c00'))


# plot_grid(np.h2o.plot.loaded.days,
#           np.h2o.plot.unloaded.days,
#           ncol = 1,
#           labels = 'AUTO')

af.h2o.loaded.summary <- af.dmmb.loaded.data %>%
  group_by(day, oxygen) %>%
  summarise(mean_h2o = mean(water_percent),
             se_h2o = sd(water_percent) / sqrt(n()))

af.h2o.unloaded.summary <- af.dmmb.unloaded.data %>%
  group_by(day, oxygen) %>%
  summarise(mean_h2o = mean(water_percent),
            se_h2o = sd(water_percent) / sqrt(n()))

af.h2o.plot.loaded.days <- ggplot(af.h2o.loaded.summary, aes(day, mean_h2o, fill=oxygen)) +
  ggtitle('AF') +
  labs(x = 'Time (days)',
       y = 'Water Content (%)') +
  geom_col(position = position_dodge(),
           color = 'black') +
  geom_errorbar(aes(ymin=mean_h2o, ymax=mean_h2o+se_h2o),
                position = position_dodge()) +
  scale_fill_manual(name = 'Oxygen\nTension', labels = c('21%', '5%', 'Harvest'), values = c('#005030', '#F47321', '#500c00'))


af.h2o.plot.unloaded.days <- ggplot(af.h2o.unloaded.summary, aes(day, mean_h2o, fill=oxygen)) +
  ggtitle('AF') +
  labs(x = 'Time (days)',
       y = 'Water Content (%)') +
  geom_col(position = position_dodge(),
           color = 'black') +
  geom_errorbar(aes(ymin=mean_h2o, ymax=mean_h2o+se_h2o),
                position = position_dodge()) +
  scale_fill_manual(name = 'Oxygen\nTension', labels = c('21%', '5%', 'Harvest'), values = c('#005030', '#F47321', '#500c00'))


plot_grid(np.h2o.plot.loaded.days, af.h2o.plot.loaded.days,
          np.h2o.plot.unloaded.days, af.h2o.plot.unloaded.days,
          ncol = 2,
          labels = 'AUTO')

# np.h2o.plot.loaded.5ox.days <- ggplot(np.dmmb.loaded.5ox.data, aes(x = day,y = water_percent, group = pig, color = pig)) + 
#   labs(x = 'Time (days)',
#        y = 'Water Content (%)') +
#   geom_line() +
#   geom_point() 
# plot(np.h2o.plot.loaded.5ox.days)
# 
# np.h2o.plot.loaded.21ox.days <- ggplot(np.dmmb.loaded.21ox.data, aes(x = day,y = water_percent, group = pig, color = pig)) + 
#   labs(x = 'Time (days)',
#        y = 'Water Content (%)') +
#   geom_line() +
#   geom_point() 
# plot(np.h2o.plot.loaded.21ox.days)
# 
# np.h2o.plot.unloaded.5ox.days <- ggplot(np.dmmb.unloaded.5ox.data, aes(x = day,y = water_percent, group = pig, color = pig)) + 
#   labs(x = 'Time (days)',
#        y = 'Water Content (%)') +
#   geom_line() +
#   geom_point() 
# plot(np.h2o.plot.unloaded.5ox.days)
# 
# np.h2o.plot.unloaded.21ox.days <- ggplot(np.dmmb.unloaded.21ox.data, aes(x = day,y = water_percent, group = pig, color = pig)) + 
#   labs(x = 'Time (days)',
#        y = 'Water Content (%)') +
#   geom_line() +
#   geom_point() 
# plot(np.h2o.plot.unloaded.21ox.days)
# 
# plot_grid(np.h2o.plot.loaded.5ox.days, 
#           np.h2o.plot.loaded.21ox.days, 
#           np.h2o.plot.unloaded.5ox.days, 
#           np.h2o.plot.unloaded.21ox.days,
#           ncol = 2,
#           labels = 'AUTO')
# 
# af.h2o.plot.loaded.5ox.days <- ggplot(af.dmmb.loaded.5ox.data, aes(x = day,y = water_percent, group = pig, color = pig)) + 
#   labs(x = 'Time (days)',
#        y = 'Water Content (%)') +
#   geom_line() +
#   geom_point() 
# plot(af.h2o.plot.loaded.5ox.days)
# 
# af.h2o.plot.loaded.21ox.days <- ggplot(af.dmmb.loaded.21ox.data, aes(x = day,y = water_percent, group = pig, color = pig)) + 
#   labs(x = 'Time (days)',
#        y = 'Water Content (%)') +
#   geom_line() +
#   geom_point() 
# plot(af.h2o.plot.loaded.21ox.days)
# 
# af.h2o.plot.unloaded.5ox.days <- ggplot(af.dmmb.unloaded.5ox.data, aes(x = day,y = water_percent, group = pig, color = pig)) + 
#   labs(x = 'Time (days)',
#        y = 'Water Content (%)') +
#   geom_line() +
#   geom_point() 
# plot(af.h2o.plot.unloaded.5ox.days)
# 
# af.h2o.plot.unloaded.21ox.days <- ggplot(af.dmmb.unloaded.21ox.data, aes(x = day,y = water_percent, group = pig, color = pig)) + 
#   labs(x = 'Time (days)',
#        y = 'Water Content (%)') +
#   geom_line() +
#   geom_point() 
# plot(af.h2o.plot.unloaded.21ox.days)
# 
# plot_grid(af.h2o.plot.loaded.5ox.days, 
#           af.h2o.plot.loaded.21ox.days, 
#           af.h2o.plot.unloaded.5ox.days, 
#           af.h2o.plot.unloaded.21ox.days,
#           ncol = 2,
#           labels = 'AUTO')
