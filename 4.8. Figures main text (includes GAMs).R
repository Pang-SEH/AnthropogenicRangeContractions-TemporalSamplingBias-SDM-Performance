library(ggplot2)
library(raster)
library(gridExtra)
library(grid)
library(ggpubr)
library(viridis)
library(reshape2)
library(stringr)
library(dplyr)
library(lme4)
library(ggpmisc)
library(ggtext)
library(doParallel)
library(foreach)
setwd('D:/YourDirectory/SupplementaryFiles/DataFiles')

# mean and sd function #
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

#### Fig. 2 ################################################################################
## temporal sampling patterns (observed and simulated) ##
eg.sp <- read.csv("TSP example.csv")
sim.freq <- read.csv("sim.freq.csv")
eg.sp <- eg.sp[,c('year', 'species', 'Freq')]
eg.sp$variable[eg.sp$species == 'Cinnamomum verum'] <- "CP"
eg.sp$variable[eg.sp$species == 'Daphniphyllum glaucescens'] <- "CI"
eg.sp$variable[eg.sp$species == 'Ficus benguetensis'] <- "CR"
eg.sp$variable[eg.sp$species == 'Cinnamomum parthenoxylon'] <- "SP"
eg.sp$variable[eg.sp$species == 'Croton tiglium'] <- "SI"
eg.sp$variable[eg.sp$species == 'Knema latifolia'] <- "SR"
sim.freq <- merge.data.frame(sim.freq, eg.sp[,c('year', 'variable', 'species')])
colnames(sim.freq) <- c('year', 'variable', 'Freq', 'species')
sim.freq$type <- 'Simulated'
eg.sp$type <- 'Actual'
temppat <- rbind(sim.freq, eg.sp)
temppat$species <- factor(temppat$species, levels = c('Cinnamomum verum', 'Daphniphyllum glaucescens',
                                                      'Ficus benguetensis', 'Cinnamomum parthenoxylon',
                                                      'Croton tiglium', 'Knema latifolia'))
colnames(temppat) <- c('Year', 'variable', 'Freq', 'species', 'Pattern')

title.nm <- c('**Clustered Past (1920)**<br>*Cinnamomum verum*', '**Clustered Intermediate (1965)**<br>*Daphniphyllum glaucescens*',
              '**Clustered Recent (1991)**<br>*Ficus benguetensis*', '**Spread Past (1933)**<br>*Cinnamomum parthenoxylon*',
              '**Spread Intermediate (1950)**<br>*Croton tiglium*', '**Spread Recent (1967)**<br>*Knema latifolia*')
names(title.nm) <- c('Cinnamomum verum', 'Daphniphyllum glaucescens',
                     'Ficus benguetensis', 'Cinnamomum parthenoxylon',
                     'Croton tiglium', 'Knema latifolia')
act.sim <- ggplot(temppat, aes(x=Year, y=Freq, group = Pattern, colour = Pattern, size = Pattern, alpha = Pattern)) +
  geom_line() +
  scale_size_manual(values = c(0.5, 0.5), labels = c('Real Species', 'Simulated')) +
  scale_alpha_manual(values = c(0.7, 0.5), labels = c('Real Species', 'Simulated')) +
  scale_color_manual(values = c('black', 'red'), labels = c('Real Species', 'Simulated')) +
  scale_x_continuous(breaks = seq(1900, 2000, 20), limits = c(1900,2000)) +
  facet_wrap(~species, scales = 'free_y', labeller = labeller(species = title.nm)) +
  labs(x = 'Year', y = 'Sampling\nFrequency / Probability') +
  theme_bw() +
  theme(strip.text = ggtext::element_markdown(size = 12),
        strip.background = element_blank(),
        panel.border = element_rect(colour = 'grey50'),
        legend.position = 'bottom',
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10, colour = 'grey30'),
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        panel.grid.major = element_line(colour = 'grey98'),
        panel.grid.minor = element_blank())
act.sim

#### Fig. 3 ################################################################################
## manifestation of data misrepresentations among temporal sampling patterns ##
## temporal sampling distributions were aggregated ##
er <- read.csv('Error Rates in Sampling Data (focc) (all).csv')#[,-4]
colnames(er) <- c('species', colnames(er)[2:length(er)])

er <- melt(er, id.vars = c('species'))
er <- cbind(er, str_split_fixed(er$species, '[.]', 4))
er <- er[, c(4,5,6,7,2,3)]
colnames(er) <- c('Cniche', 'Prev', 'Bias', 'Sample', 'Error', 'Degree')
er$Bias <- factor(er$Bias, levels = c('HT', 'CA', 'SA', 'SR', 'CR', 'ST',  'CT', 'MT'))

sp.stat <- read.csv('species statistics.csv')[,-1]
colnames(sp.stat) <- c('Cniche', colnames(sp.stat)[-1])
sp.stat$ARC <- (sp.stat$D.1900 - sp.stat$D.2000)/sp.stat$D.1900

### correcting niche truncation against negative controls ###
# er.dif <- er[er$Error == 'Truncation',]
# 
# cl <- makeCluster(detectCores()-4)
# registerDoParallel(cl)
# cor.trunc <- foreach (i = 1:nrow(er.dif)) %dopar% {
#   chg <- mean(er[er$Cniche == er.dif$Cniche[i] &
#                    er$Prev == er.dif$Prev[i] &
#                    er$Bias == 'HT' &
#                    er$Error == 'Truncation', 'Degree'])
#   er.dif$Degree[i] - chg
# }
# stopCluster(cl)
# er.dif$Degree <- unlist(cor.trunc)
# er.dif <- rbind(er.dif, er[er$Error == 'Discrepancy',])
# er.dif <- merge.data.frame(er.dif, sp.stat)
# 
# save(er.dif, file = 'error.difference.corrected (all)')
load('error.difference.corrected (all)')

# plot #
er.df <- rbind(er.dif[er.dif$Error == 'Discrepancy' & !er.dif$Bias == 'MT',],
               er.dif[er.dif$Error == 'Truncation' & !er.dif$Bias == 'HT',])
er.df$weight[grepl('CA|SA', er.df$Bias)] <- 'Past'
er.df$weight[grepl('CR|SR', er.df$Bias)] <- 'Intermediate'
er.df$weight[grepl('CT|ST', er.df$Bias)] <- 'Recent'
er.df$weight[grepl('HT|MT', er.df$Bias)] <- 'Control'
er.add <- er.df[grepl('HT|MT', er.df$Bias),]
er.add$pattern <- 'Clustered'
er.df$pattern[grepl('CA|CR|CT', er.df$Bias)] <- 'Clustered'
er.df$pattern[grepl('SA|SR|ST|HT|MT', er.df$Bias)] <- 'Spread'
er.df <- rbind(er.df, er.add); rm(er.add)
head(er.df)
er.df$weight <- factor(er.df$weight, levels = c('Control', 'Past', 'Intermediate', 'Recent'))
er.df$Error <- factor(er.df$Error, levels = c('Discrepancy', 'Truncation'))
levels(er.df$Error) <- c("Discrepancy" = "Mean difference in\nhabitat suitability",
                         "Truncation" = "Corrected\n1 - Warren's I")
er.var <- ggplot(data = er.df[er.df$Prev == 'P25', ], aes(x = ARC, y = Degree, colour = weight, linetype = weight)) +
  geom_smooth(method = 'lm', se = F) +
  geom_point(alpha = 0.1, shape = 18, size = 0.3) +
  geom_hline(yintercept = 0, size = 0.5, colour = 'grey60') +
  scale_linetype_manual(values = c('dashed', 'solid', 'solid', 'solid'),
                        labels = c('Positive Control', 'Past', 'Intermediate', 'Recent')) +
  scale_colour_manual(values = c('black', "#00AFBB", "#E7B800", "#FC4E07"),
                      labels = c('Positive Control', 'Past', 'Intermediate', 'Recent')) +
  stat_poly_eq(formula = y ~ x,
               aes(label = paste(..eq.label.., ..rr.label.., sep = "*plain(\",\")~")),
               parse = TRUE, size = 3) +
  labs(y = 'Misrepresentation',
       x = 'Anthropogenic range contractions<br>(proportion)') +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.title.x = ggtext::element_markdown(),
        axis.text = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, size = 12, face = 'bold'),
        legend.position = 'bottom',
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank()) +
  facet_grid(Error~NA, switch = 'y', scales = 'free_y') +
  theme(strip.background = element_blank(),
        strip.placement = 'outside', strip.text = element_text(size = 12), strip.text.x = element_blank())
er.var

#### Fig. 4 ################################################################################
## Evaluating model performance ##
## Generalised additive model ##
ma <- read.csv('Temporal Bias scores (focc) (all) cont.csv', row.names = 1)
colnames(ma) <- c('species', colnames(ma)[2:ncol(ma)])

## data prep ##
samp.er <- read.csv('Error Rates in Sampling Data (focc) (all).csv') [,c(1,2,3,5,7)]
colnames(samp.er) <- c('species', colnames(samp.er)[2:length(samp.er)])

load('error.difference.corrected (all)')
er.dif$species <- paste(er.dif$Cniche, er.dif$Prev, er.dif$Bias, er.dif$Sample, sep = '.')
for (i in samp.er$species) {
  samp.er$Truncation[samp.er$species == i] <- er.dif$Degree[er.dif$Error == 'Truncation' & er.dif$species == i]
}

ma <- melt(ma, id.vars = 'species')
ma <- merge.data.frame(ma, samp.er, by = 'species')

ma <- cbind(ma, str_split_fixed(ma$variable, '[.]', 2))
ma <- cbind(ma, str_split_fixed(ma$species, '[.]', 4))
head(ma)
ma <- ma[, c(10,11,12,13,8,9,3,5,6)]
colnames(ma) <- c('Cniche', 'Prev', 'Bias', 'Sample',
                  'Model', 'Metric', 'Score',
                  'Discrepancy', 'Truncation')
ma$Metric <- factor(ma$Metric, levels = c('mAUC', 'mTSS', 'mPOD', 'mSR', 'mPOFD',
                                          'aFAP', 'aFAS', 'aNOI',
                                          'hFAP', 'hFAS', 'hNOI',
                                          'aAUC', 'aTSS',
                                          'hAUC', 'hTSS',
                                          'aPOD', 'aSR', 'aPOFD',
                                          'hPOD', 'hSR', 'hPOFD',
                                          'aFBI', 'hFBI'))
ma$Bias <- factor(ma$Bias, levels = c('HT', 'CA', 'SA', 'SR', 'CR', 'ST',  'CT', 'MT'))
head(ma)

## model fit ##
## model relationship with errors and parameters
## projection ##
# library(nlme)
# library(lme4)
# library(mgcv)
# library(gam)
# library(tidyverse)

## GAM to model independent effects ##
# mgam <- as.data.frame(matrix(NA,0,0))
# head(ma)
# for (m in unique(ma$Metric)) {
#   for (t in c('wLand', 'EnvO')) {
#     for (p in unique(ma$Prev)) {
#       dat <- ma[grepl(m, ma$Metric) & ma$Model == t,]
#       dat <- dat[grepl(p, dat$Prev),]
#       smod <- mgcv::gam(Score ~ s(Discrepancy) + s(Truncation),
#                         data = dat)
#       
#       Disc.pred <- tibble(Discrepancy = dat$Discrepancy,
#                           Truncation = 0,
#                           Prev = p)
#       Disc.pred <- predict(smod, newdata = Disc.pred, se.fit = T, terms = c("s(Discrepancy)", "s(Truncation)")) %>%
#         as_tibble() %>%
#         cbind(Disc.pred)
#       colnames(Disc.pred) <- c('Score', 'se.fit', 'Degree', 'Truncation', 'Prev')
#       Disc.pred$Error <- 'Discrepancy'
#       Disc.pred$Metric <- m
#       Disc.pred$Model <- t
#       Disc.pred <- Disc.pred[, c(1,2,3,6,5,7,8)]
#       
#       Trun.pred <- tibble(Discrepancy = 0,
#                           Truncation = dat$Truncation,
#                           Prev = p)
#       Trun.pred <- predict(smod, newdata = Trun.pred, se.fit = T, terms = c("s(Discrepancy)", "s(Truncation)")) %>%
#         as_tibble() %>%
#         cbind(Trun.pred)
#       colnames(Trun.pred) <- c('Score', 'se.fit', 'Discrepancy', 'Degree', 'Prev')
#       Trun.pred$Error <- 'Truncation'
#       Trun.pred$Metric <- m
#       Trun.pred$Model <- t
#       Trun.pred <- Trun.pred[, c(1,2,4,6,5,7,8)]
#       
#       mgam <- rbind(mgam, Disc.pred, Trun.pred)
#     }
#   }
# }
# setwd(res.dir)
# save(mgam, file = 'Model Performance GAM (focc)(negctrl) (all) v1')
# rm(mgam)

## independent effects on model performance ##
load('Model Performance GAM (focc)(negctrl) (all) v1')

## plotting ##

mgam <- mgam[mgam$Prev == 'P25',]
mgam.prev <- mgam[grepl('aFAP|aAUC|aPOD|aSR', mgam$Metric),]
mgam.prev$Metric <- sub('a', '', mgam.prev$Metric)
mgam.prev$Projection <- 'Contemporary'
mgam.prev2 <- mgam[grepl('hFAP|hAUC|hPOD|hSR', mgam$Metric),]
mgam.prev2$Metric <- sub('h', '', mgam.prev2$Metric)
mgam.prev2$Projection <- 'Historical'
mgam.prev <- rbind(mgam.prev, mgam.prev2)
rm(mgam.prev2)

mgam.prev$Score[grepl('POD', mgam.prev$Metric)] <- 1 - mgam.prev$Score[grepl('POD', mgam.prev$Metric)]
mgam.prev$Score[grepl('SR', mgam.prev$Metric)] <- 1 - mgam.prev$Score[grepl('SR', mgam.prev$Metric)]
mgam.prev$Metric <- as.character(mgam.prev$Metric)
mgam.prev$Metric[mgam.prev$Metric == 'POD'] <- 'OR'
mgam.prev$Metric[mgam.prev$Metric == 'SR'] <- 'CR'

mgam.prev$Metric <- factor(mgam.prev$Metric, levels = c('FAP', 'AUC', 'OR', 'CR'))
mgam.prev$Model <- factor(mgam.prev$Model, levels = c('wLand', 'EnvO'))
mgam.prev$Projection <- factor(mgam.prev$Projection, levels = c('Historical', 'Contemporary'))

mrd.plot1 <- ggplot(mgam.prev[grepl('Discrepancy', mgam.prev$Error) &
                                mgam.prev$Prev == 'P25' &
                                grepl('FAP|AUC', mgam.prev$Metric),], 
                    aes(x = Degree, y = Score, fill = Model, colour = Model)) +
  geom_ribbon(alpha = 0.1, colour = NA,
              aes(ymin = Score - 1.96 * se.fit,
                  ymax = Score + 1.96 * se.fit)) +
  geom_line(size = 1) +
  scale_fill_manual(values = c('#0072B2', '#D55E00')) +
  scale_colour_manual(values = c('#0072B2', '#D55E00')) +
  ggtitle('Probabilistic Performance') +
  ylim(0.4,1) +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 10),
        legend.position = 'none',
        plot.title = element_text(hjust = 0.5, size = 12, face = 'bold'),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank()) +
  labs(x = 'Mean difference in habitat suitability', y = NULL, fill = 'Model:', colour = 'Model:') +
  facet_grid(Metric ~ Projection, switch = 'y') +
  theme(strip.background = element_blank(),
        strip.placement = 'outside', strip.text = element_text(size = 12), 
        strip.text.y = element_text(face = 'bold'))
mrd.plot1

mrd.plot2 <- ggplot(mgam.prev[grepl('Discrepancy', mgam.prev$Error) &
                                mgam.prev$Prev == 'P25' &
                                grepl('OR|CR', mgam.prev$Metric),], 
                    aes(x = Degree, y = Score, fill = Model, colour = Model)) +
  geom_ribbon(alpha = 0.1, colour = NA,
              aes(ymin = Score - 1.96 * se.fit,
                  ymax = Score + 1.96 * se.fit)) +
  geom_line(size = 1) +
  scale_fill_manual(values = c('#0072B2', '#D55E00'), labels = c('With land use', 'Without land use')) +
  scale_colour_manual(values = c('#0072B2', '#D55E00'), labels = c('With land use', 'Without land use')) +
  ggtitle('Binary Error Rates') +
  ylim(0,0.8) +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 10),
        legend.position = 'bottom',
        plot.title = element_text(hjust = 0.5, size = 12, face = 'bold'),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank()) +
  labs(x = 'Mean difference in habitat suitability', y = NULL, fill = 'Models:', colour = 'Models:') +
  facet_grid(Metric ~ Projection, switch = 'y') +
  theme(strip.background = element_blank(),
        strip.placement = 'outside', strip.text = element_text(size = 12), 
        strip.text.y = element_text(face = 'bold'))
mrd.plot2

mrd.leg <- get_legend(mrd.plot2)
mrd.plot2 <- mrd.plot2 + theme(legend.position = 'none')
mrd.plot1 <- mrd.plot1 + theme(legend.position = 'none')
mrd.plot <-  ggarrange(ggarrange(mrd.plot1, mrd.plot2), 
                       nrow = 2, heights = c(1,0.05))

mrt.plot1 <- ggplot(mgam.prev[grepl('Truncation', mgam.prev$Error) &
                                mgam.prev$Prev == 'P25' &
                                grepl('FAP|AUC', mgam.prev$Metric),], 
                    aes(x = Degree, y = Score, fill = Model, colour = Model)) +
  geom_ribbon(alpha = 0.1, colour = NA,
              aes(ymin = Score - 1.96 * se.fit,
                  ymax = Score + 1.96 * se.fit)) +
  geom_line(size = 1) +
  scale_fill_manual(values = c('#0072B2', '#D55E00')) +
  scale_colour_manual(values = c('#0072B2', '#D55E00')) +
  ggtitle('Probabilistic Performance') +
  ylim(0.4,1) +
  scale_x_continuous(breaks = c(0, 0.1, 0.2), limits = c(0, 0.23)) +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 10),
        legend.position = 'none',
        plot.title = element_text(hjust = 0.5, size = 12, face = 'bold'),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank()) +
  labs(x = expression(paste("Corrected 1 - Warren's ", italic("I"))), y = NULL) +
  facet_grid(Metric ~ Projection, switch = 'y') +
  theme(strip.background = element_blank(),
        strip.placement = 'outside', strip.text.y = element_text(face = 'bold'),
        strip.text = element_text(size = 12))

mrt.plot1

mrt.plot2 <- ggplot(mgam.prev[grepl('Truncation', mgam.prev$Error) &
                                mgam.prev$Prev == 'P25' &
                                grepl('OR|CR', mgam.prev$Metric),], 
                    aes(x = Degree, y = Score, fill = Model, colour = Model)) +
  geom_ribbon(alpha = 0.1, colour = NA,
              aes(ymin = Score - 1.96 * se.fit,
                  ymax = Score + 1.96 * se.fit)) +
  geom_line(size = 1) +
  scale_fill_manual(values = c('#0072B2', '#D55E00'), labels = c('With land use', 'Without land use')) +
  scale_colour_manual(values = c('#0072B2', '#D55E00'), labels = c('With land use', 'Without land use')) +
  ggtitle('Binary Error Rates') +
  ylim(0,0.8) +
  scale_x_continuous(breaks = c(0, 0.1, 0.2), limits = c(0, 0.23)) +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 10),
        legend.position = 'bottom',
        plot.title = element_text(hjust = 0.5, size = 12, face = 'bold'),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank()) +
  labs(x = expression(paste("Corrected 1 - Warren's ", italic("I"))), y = NULL, fill = 'Models:', colour = 'Models:') +
  facet_grid(Metric ~ Projection, switch = 'y') +
  theme(strip.background = element_blank(),
        strip.placement = 'outside', strip.text = element_text(size = 12), 
        strip.text.y = element_text(face = 'bold'))
mrt.plot2

mrt.leg <- get_legend(mrt.plot2)
mrt.plot2 <- mrt.plot2 + theme(legend.position = 'none')
mrt.plot <- ggarrange(ggarrange(mrt.plot1, mrt.plot2), 
                      nrow = 2, heights = c(1,0.05))

mpp <- ggarrange(mrd.plot, mrt.plot, mrt.leg, 
                 labels = c("a)","b)", ""), 
                 # vjust = 0,
                 nrow = 3, heights = c(1,1,0.1))
mpp

#### Fig. 5 ################################################################################
setwd('D:/YourDirectory/SupplementaryFiles/Maps')
eg.map <- stack(list.files()[c(4,3,2,1,8,7,6,5)])
names(eg.map) <- c(1:8)

eg.df <- as.data.frame(eg.map, xy = T)
eg.df <- melt(eg.df, id.vars = c('x', 'y'))
eg.df$time[grepl('X1|X3|X5|X7', eg.df$variable)] <- 'Contemporary'
eg.df$time[grepl('X2|X4|X6|X8', eg.df$variable)] <- 'Historical'
eg.df$lulc[grepl('X1|X2|X5|X6', eg.df$variable)] <- 'With land use'
eg.df$lulc[grepl('X3|X4|X7|X8', eg.df$variable)] <- 'Without land use'
eg.df$time <- factor(eg.df$time, levels = c('Historical', 'Contemporary'))
eg.df$lulc <- factor(eg.df$lulc, levels = c('Without land use', 'With land use'))

m.plot1 <- ggplot(eg.df[grepl('X1|X2|X3|X4', eg.df$variable),],
                  aes(x = x, y = y, fill = factor(value))) +
  geom_tile() +
  theme_bw() +
  coord_equal() +
  scale_fill_manual(na.value = 'transparent',
                    values = c('#D55E00', 'grey90', '#009E73', '#CC79A7'),
                    labels = c('False Presence', 'Correct Absence', 'Correct Presence', 'False Absence', '')) +
  theme(axis.title = element_blank(),
        plot.title = ggtext::element_markdown(angle = 0, hjust = 0.5, size = 11),
        axis.text = element_blank(),
        legend.position = 'none',
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 12)) +
  facet_grid(lulc ~ time, switch = 'y')
m.plot1
m.plot2 <- ggplot(eg.df[grepl('X5|X6|X7|X8', eg.df$variable),],
                  aes(x = x, y = y, fill = factor(value))) +
  geom_tile() +
  theme_bw() +
  coord_equal() +
  scale_fill_manual(na.value = 'transparent',
                    values = c('#D55E00', 'grey90', '#009E73', '#CC79A7'),
                    labels = c('False Presence', 'Correct Absence', 'Correct Presence', 'False Absence', '')) +
  theme(axis.title = element_blank(),
        plot.title = ggtext::element_markdown(angle = 0, hjust = 0.5, size = 11),
        axis.text = element_blank(),
        legend.position = 'none',
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 12)) +
  facet_grid(lulc ~ time, switch = 'y')

m.legh <- m.plot2 + 
  scale_fill_manual(breaks = c('-1', '0', '1', '2'),
                    name = element_blank(),
                    values = c('#009E73', 'grey90', '#D55E00', '#CC79A7'),
                    labels = c('Correct Presence', 'Correct Absence', 'False Presence', 'False Absence')) +
  theme(legend.position = 'bottom')
m.legh <- get_legend(m.legh)

bi.map <- ggarrange(ggarrange(m.plot1, m.plot2, ncol = 2, #hjust = -0.2,
                              labels = c('a)', 'b)')),
                    m.legh, nrow = 2, 
                    heights = c(1,0.05))
bi.map