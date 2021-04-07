library(ggplot2)
library(gridExtra)
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

#### Fig. S1.5 ################################################################################
## Temporal sampling patterns observed (total/real species/virtual species) ##
## a)
abs.bias <- read.csv('abs.bias.csv', row.names = 1)
gbif.or <- ggplot(abs.bias, aes(x=Var1, y=Freq, group = 1)) +
  geom_line(alpha = 0.7) +
  scale_x_continuous(breaks = seq(1700, 2020, 10), limits = c(1850,2010)) +
  ylim(0,NA) +
  labs(x = 'Year', y = 'Sampling Frequency\n(Magnoliopsida)') +
  geom_smooth(method='gam', formula= y ~ s(x, bs = "cs"), se = F, size = 0.7, alpha = 0.5) +
  theme_bw() +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 10, colour = 'grey30'))
gbif.or 

## b)
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

## c)
simsamp <- read.csv("Simulated Temporal Sampling.csv")
simsamp <- as.data.frame(table(simsamp))
for (b in unique(simsamp$pattern)) {
  simsamp[simsamp$pattern == b,]$Freq <- simsamp[simsamp$pattern == b,]$Freq / sum(simsamp[simsamp$pattern == b,]$Freq)
}
simsamp$year <- as.numeric(as.character(simsamp$year))

vir.samppat <- c('**Clustered Past**', '**Clustered Intermediate**',
                 '**Clustered Recent**', '**Spread Past**',
                 '**Spread Intermediate**', '**Spread Recent**')
names(vir.samppat) <- c('CA', 'CR',
                        'CT', 'SA',
                        'SR', 'ST')

sim.samp.plot <- ggplot(simsamp, aes(x=year, y=Freq, group = 1)) +
  geom_line(alpha = 0.7) +
  scale_x_continuous(breaks = seq(1900, 2000, 20), limits = c(1900,2000)) +
  facet_wrap(~pattern, scales = 'free_y', labeller = labeller(pattern = vir.samppat)) +
  labs(x = 'Year', y = 'Virtual Species\nSampling Frequency') +
  theme_bw() +
  theme(strip.text = ggtext::element_markdown(size = 12),
        strip.background = element_blank(),
        panel.border = element_rect(colour = 'grey50'),
        legend.position = 'bottom',
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10, colour = 'grey30'),
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        panel.grid.major = element_line(colour = 'grey98'),
        panel.grid.minor = element_blank())
sim.samp.plot
## final figure
ggarrange(gbif.or, act.sim, sim.samp.plot, nrow = 3, labels = c('a)', 'b)', 'c)'), heights = c(0.8, 1.1, 0.9))

#### Fig. S1.6 ################################################################################
## randomly sampled species temporal sampling pattern ##
eg.sps <- read.csv('eg.sps.csv', row.names = 1)
# rand.sp <- sample(1:length(unique(eg.sps$species)), 24)
# # randomly selected list of species to show the temporal sampling patterns observed
# save(rand.sp, file = 'rand.sp')
load('rand.sp')
rand.samp <- eg.sps[eg.sps$species %in% unique(eg.sps$species)[rand.sp],]

# labelling table of temporal sampling pattern heuristically designated to each species
tsp.lab <- data.frame(year = 1905, Freq = aggregate(Freq~species, data = rand.samp, max)[,2] * 0.95, 
                      species = unique(rand.samp$species), 
                      tp = c('CR', 'SR', 'SR', 'CR',
                             'SI', 'SI', 'SI', 'CI',
                             'SR', 'CP', 'CR', 'SI',
                             'CR', 'CI', 'CR', 'CR',
                             'SR', 'SR', 'CI', 'SI',
                             'CI', 'SI', 'CP', 'SP'))

rsp.tsp <- ggplot(rand.samp, aes(x=year, y=Freq, group = 1)) +
  geom_line(colour = 'grey20') +
  geom_smooth(method='lm', formula= y~x, se = F, col = 'red') +
  geom_smooth(method='gam', formula= y ~ s(x, bs = "cs"), se = F) +
  scale_x_continuous(breaks = seq(1900, 2000, 50)) +
  ylim(0,NA) +
  facet_wrap(~species, scales = 'free_y', ncol = 4) +
  labs(x = 'Year', y = 'Sampling Frequency') +
  theme_bw() +
  theme(strip.text = element_text(size = 8, face = 'italic'),
        strip.background = element_blank(),
        panel.border = element_rect(colour = 'grey50'),
        legend.position = 'none',
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 8, colour = 'grey30'),
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        panel.grid.major = element_line(colour = 'grey98'),
        panel.grid.minor = element_blank()) +
  geom_text(aes(label = tp), data = tsp.lab, size = 4, colour = 'black')
rsp.tsp

#### Fig. S1.7 ################################################################################
## Raw indices of data misrepresentation and the corrected version of 1 - Warren's I ##

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
er.ucd <- merge.data.frame(er, sp.stat) # error, for uncorrected difference
er.ucd <- er.ucd[er.ucd$Error %in% c('Truncation', 'Discrepancy'),]
sub.er.dif <- er.dif[er.dif$Error == 'Truncation', ]
sub.er.dif$Error <- 'Cor.Truncation'
er.ucd <- rbind(er.ucd, sub.er.dif)
er.ucd <- er.ucd[er.ucd$Bias %in% c('HT', 'MT'),]
er.ucd$Error <- factor(er.ucd$Error, levels = c('Discrepancy', 'Truncation', 'Cor.Truncation'))
levels(er.ucd$Error) <- c("Discrepancy" = "Mean difference in\nhabitat suitability",
                          "Truncation" = "Uncorrected\n1 - Warren's I",
                          'Cor.Truncation' = "Corrected\n1 - Warren's I")
ucd.lab <- data.frame(ARC = 0.06, Degree = 0.39, 
                      Error = c('Discrepancy', 'Truncation', 'Cor.Truncation'), 
                      label = c('a)', 'b)', 'c)'))
ucd.lab$Error <- factor(ucd.lab$Error, levels = c('Discrepancy', 'Truncation', 'Cor.Truncation'))
levels(ucd.lab$Error) <- c("Discrepancy" = "Mean difference in\nhabitat suitability",
                           "Truncation" = "Uncorrected\n1 - Warren's I",
                           'Cor.Truncation' = "Corrected\n1 - Warren's I")
ucd.plot <- ggplot(er.ucd[er.ucd$Prev == 'P25', ], aes(x = ARC, y = Degree, colour = Bias)) +
  geom_smooth(method = 'lm', se = F) +
  geom_point(alpha = 0.1, shape = 18, size = 1) +
  facet_grid(.~Error) +
  theme_bw() +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, size = 12, face = 'bold'),
        legend.position = 'bottom',
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank()) +
  scale_colour_discrete(name = 'Pos/Neg Controls', labels = c('Only Start', 'Only End')) +
  labs(x = 'Anthropogenic Range Contractions (proportion)', y = NULL) +
  theme(strip.background = element_blank(),
        strip.placement = 'outside', strip.text = element_text(size = 12), 
        strip.text.y = element_text(face = 'bold')) +
  geom_text(data = ucd.lab, aes(x = ARC, y = Degree, colour = NA, label = label),
            colour = 'black', size = 6)
ucd.plot

#### Fig. S2.1 ################################################################################
## Manifestation of data misrepresentations among temporal sampling patterns (disaggregated) ##

setwd(res.dir)
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

## after correcting niche truncation indice ##
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
       x = 'Anthropogenic range contractions (proportion)') +
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
  facet_grid(Error~pattern, switch = 'y', scales = 'free_y') +
  theme(strip.background = element_blank(),
        strip.placement = 'outside', strip.text = element_text(size = 12))
er.var

#### Fig. S2.2 ################################################################################
## Manifestation of data misrepresentation among prevalence levels ##

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

## correcting for niche truncation against negative controls ##
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

ecv.plot1 <- ggplot(data = er.df[er.df$pattern == 'Clustered' & er.df$weight == 'Control',], 
                    aes(x = ARC, y = Degree, colour = Prev)) +
  geom_smooth(method = 'lm', se = F, linetype = 'dashed') +
  geom_point(alpha = 0.1, shape = 18, size = 1) +
  geom_hline(yintercept = 0, size = 0.5, colour = 'grey60') +
  stat_poly_eq(formula = y ~ x,
               aes(label = paste(..eq.label.., ..rr.label.., sep = "*plain(\",\")~")),
               parse = TRUE, size = 3) +
  labs(y = 'Misrepresentation',
       x = 'Anthropogenic range contractions<br>(proportion)') +
  ggtitle('Positive Controls') +
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

ecv.plot2 <- ggplot(data = er.df[er.df$pattern == 'Clustered' & !er.df$weight == 'Control',], 
                    aes(x = ARC, y = Degree, colour = Prev)) +
  geom_smooth(method = 'lm', se = F) +
  geom_point(alpha = 0.1, shape = 18, size = 0.3) +
  geom_hline(yintercept = 0, size = 0.5, colour = 'grey60') +
  stat_poly_eq(formula = y ~ x,
               aes(label = paste(..eq.label.., ..rr.label.., sep = "*plain(\",\")~")),
               parse = TRUE, size = 2) +
  labs(y = 'Misrepresentation',
       x = 'Anthropogenic range contractions (proportion)') +
  ggtitle('Clustered') +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.title.x = ggtext::element_markdown(),
        axis.text = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, size = 12, face = 'bold'),
        legend.position = 'none',
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank()) +
  facet_grid(Error~weight, switch = 'y', scales = 'free_y') +
  theme(strip.background = element_blank(),
        strip.placement = 'outside', strip.text = element_text(size = 12))

ecv.plot3 <- ggplot(data = er.df[er.df$pattern == 'Spread' & !er.df$weight == 'Control',], 
                    aes(x = ARC, y = Degree, colour = Prev)) +
  geom_smooth(method = 'lm', se = F) +
  geom_point(alpha = 0.1, shape = 18, size = 0.3) +
  geom_hline(yintercept = 0, size = 0.5, colour = 'grey60') +
  stat_poly_eq(formula = y ~ x,
               aes(label = paste(..eq.label.., ..rr.label.., sep = "*plain(\",\")~")),
               parse = TRUE, size = 2) +
  labs(y = 'Misrepresentation',
       x = 'Anthropogenic range contractions (proportion)') +
  ggtitle('Spread') +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.title.x = ggtext::element_markdown(),
        axis.text = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, size = 12, face = 'bold'),
        legend.position = 'none',
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank()) +
  facet_grid(Error~weight, switch = 'y', scales = 'free_y') +
  theme(strip.background = element_blank(),
        strip.placement = 'outside', strip.text = element_text(size = 12))

ecv.plot <- ggarrange(ecv.plot1, 
                      ggarrange(ecv.plot2, ecv.plot3, ncol = 1, labels = c('b)', 'c)')),
                      widths = c(0.7, 1), labels = c('a)', ''))
ecv.plot

#### Fig. S2.3 ################################################################################
## Model performance across prevalence levels ##
# Occurrence-habitat mismatching / mean difference in habitat suitability #
load('Model Performance GAM (focc)(negctrl) (all) v1')

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

# a)
mrd.plot1 <- ggplot(mgam.prev[grepl('Discrepancy', mgam.prev$Error) &
                                grepl('FAP|AUC', mgam.prev$Metric),], 
                    aes(x = Degree, y = Score, fill = Prev, colour = Prev, linetype = Model)) +
  geom_ribbon(alpha = 0.1, colour = NA,
              aes(ymin = Score - 1.96 * se.fit,
                  ymax = Score + 1.96 * se.fit)) +
  geom_line(size = 1) +
  scale_linetype(labels = c('With land use', 'Without land use')) +
  ggtitle('Probabilistic Performance') +
  ylim(0.31,1) +
  scale_x_continuous(breaks = c(0, 0.1, 0.2, 0.3, 0.4), limits = c(0, 0.41)) +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.title.x = element_text(size = 11),
        axis.text = element_text(size = 10),
        legend.position = 'none',
        plot.title = element_text(hjust = 0.5, size = 12, face = 'bold'),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank()) +
  labs(x = 'Mean difference in habitat suitability', y = NULL, fill = 'Prevalence', colour = 'Prevalence') +
  facet_grid(Metric ~ Projection, switch = 'y') +
  theme(strip.background = element_blank(),
        strip.placement = 'outside', strip.text = element_text(size = 12), 
        strip.text.y = element_text(face = 'bold'))
mrd.plot1

# b)
mrd.plot2 <- ggplot(mgam.prev[grepl('Discrepancy', mgam.prev$Error) &
                                grepl('OR|CR', mgam.prev$Metric),], 
                    aes(x = Degree, y = Score, fill = Prev, colour = Prev, linetype = Model)) +
  geom_ribbon(alpha = 0.1, colour = NA,
              aes(ymin = Score - 1.96 * se.fit,
                  ymax = Score + 1.96 * se.fit)) +
  geom_line(size = 1) +
  scale_linetype(labels = c('With land use', 'Without land use')) +
  ggtitle('Binary Error Rates') +
  ylim(0,0.8) +
  scale_x_continuous(breaks = c(0, 0.1, 0.2, 0.3, 0.4), limits = c(0, 0.41)) +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.title.x = element_text(size = 11),
        axis.text = element_text(size = 10),
        legend.position = 'right',
        plot.title = element_text(hjust = 0.5, size = 12, face = 'bold'),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank()) +
  labs(x = 'Mean difference in habitat suitability', y = NULL, fill = 'Prevalence', colour = 'Prevalence',
       linetype = 'Model') +
  facet_grid(Metric ~ Projection, switch = 'y') +
  theme(strip.background = element_blank(),
        strip.placement = 'outside', strip.text = element_text(size = 12), 
        strip.text.y = element_text(face = 'bold'))
mrd.plot2
mrd.leg <- get_legend(mrd.plot2)
mrd.plot2 <- mrd.plot2 + theme(legend.position = 'none')
mrd.plot1 <- mrd.plot1 + theme(legend.position = 'none')
mrd.plot <-  ggarrange(ggarrange(mrd.plot1, mrd.plot2, ncol = 1, labels = c('a)', 'b)')), mrd.leg, 
                       ncol = 2, widths = c(1,0.3))
mrd.plot

#### Fig. S2.4 ####
## Model performance across prevalence levels ##
# Niche truncation / corrected 1 - Warren's I #
# continuation data from part Fig. S2.3
# a)
mrt.plot1 <- ggplot(mgam.prev[grepl('Truncation', mgam.prev$Error) &
                                grepl('FAP|AUC', mgam.prev$Metric),], 
                    aes(x = Degree, y = Score, fill = Prev, colour = Prev, linetype = Model)) +
  geom_ribbon(alpha = 0.1, colour = NA,
              aes(ymin = Score - 1.96 * se.fit,
                  ymax = Score + 1.96 * se.fit)) +
  geom_line(size = 1) +
  scale_linetype(labels = c('With land use', 'Without land use')) +
  ggtitle('Probabilistic Performance') +
  ylim(0.31,1) +
  scale_x_continuous(breaks = c(0, 0.1, 0.2, 0.3, 0.4), limits = c(0, 0.42)) +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.title.x = element_text(size = 11),
        axis.text = element_text(size = 10),
        legend.position = 'right',
        plot.title = element_text(hjust = 0.5, size = 12, face = 'bold'),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank()) +
  labs(x = expression(paste("Corrected 1 - Warren's ", italic("I"))), y = NULL, 
       fill = 'Prevalence', colour = 'Prevalence', linetype = 'Model') +
  facet_grid(Metric ~ Projection, switch = 'y') +
  theme(strip.background = element_blank(),
        strip.placement = 'outside', strip.text.y = element_text(face = 'bold'),
        strip.text = element_text(size = 12))

mrt.plot1

# b)
mrt.plot2 <- ggplot(mgam.prev[grepl('Truncation', mgam.prev$Error) &
                                grepl('OR|CR', mgam.prev$Metric),], 
                    aes(x = Degree, y = Score, fill = Prev, colour = Prev, linetype = Model)) +
  geom_ribbon(alpha = 0.1, colour = NA,
              aes(ymin = Score - 1.96 * se.fit,
                  ymax = Score + 1.96 * se.fit)) +
  geom_line(size = 1) +
  scale_linetype(labels = c('With land use', 'Without land use')) +
  ggtitle('Binary Error Rates') +
  ylim(0,0.8) +
  scale_x_continuous(breaks = c(0, 0.1, 0.2, 0.3, 0.4), limits = c(0, 0.42)) +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.title.x = element_text(size = 11),
        axis.text = element_text(size = 10),
        legend.position = 'right',
        plot.title = element_text(hjust = 0.5, size = 12, face = 'bold'),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank()) +
  labs(x = expression(paste("Corrected 1 - Warren's ", italic("I"))), y = NULL, 
       fill = 'Prevalence', colour = 'Prevalence', linetype = 'Model') +
  facet_grid(Metric ~ Projection, switch = 'y') +
  theme(strip.background = element_blank(),
        strip.placement = 'outside', strip.text = element_text(size = 12), 
        strip.text.y = element_text(face = 'bold'))
mrt.plot2
mrt.leg <- get_legend(mrt.plot2)
mrt.plot2 <- mrt.plot2 + theme(legend.position = 'none')
mrt.plot1 <- mrt.plot1 + theme(legend.position = 'none')
mrt.plot <-  ggarrange(ggarrange(mrt.plot1, mrt.plot2, ncol = 1, labels = c('a)', 'b)')), mrt.leg, 
                       ncol = 2, widths = c(1,0.3))
mrt.plot
#### Fig. S2.5 ################################################################################
## sensitivity analysis for land use variable effect size (permutation importance) ##
er <- read.csv('Error Rates in Sampling Data (focc) (all).csv')#[,-4]
colnames(er) <- c('species', colnames(er)[2:length(er)])

er <- melt(er, id.vars = c('species'))
er <- cbind(er, str_split_fixed(er$species, '[.]', 4))
er <- er[, c(4,5,6,7,2,3)]
colnames(er) <- c('Cniche', 'Prev', 'Bias', 'Sample', 'Error', 'Degree')
er$Bias <- factor(er$Bias, levels = c('HT', 'CA', 'SA', 'SR', 'CR', 'ST',  'CT', 'MT'))

## Partial correlation between errors ##
load('error.difference.corrected (all)')
library(ppcor)
er.dit <- er.dif[er.dif$Error == 'Truncation', c(1,2,3,4,6,10)]
er.did <- er.dif[er.dif$Error == 'Discrepancy', c(1,2,3,4,6,10)]
names(er.dit)[5] <- "Truncation"
names(er.did)[5] <- "Discrepancy"
er.cor <- merge.data.frame(er.dit, er.did)
pcor(er.cor[er.cor$Prev == 'P25', c('ARC', 'Truncation', 'Discrepancy')])

## land use effect size ##
luvarb <- read.csv('Variable Importance (focc_cont).csv')
luvarb <- cbind(luvarb, str_split_fixed(luvarb$X, '[.]', 4))
luvarb <- luvarb[, c(4,5,6,7,2)]
names(luvarb) <- c('Cniche', 'Prev', 'Bias', 'Sample', 'PI')
luvarb <- merge.data.frame(er.cor, luvarb)
head(luvarb)
pcor(luvarb[, c('ARC', 'Discrepancy', 'PI')], method = 'spearman')
luvarb <- luvarb[order(luvarb$ARC, decreasing = F),]

# a)
pi.plot1 <- ggplot(luvarb, aes(x = Discrepancy, y = PI, colour = ARC)) +
  geom_point(size = 3, alpha = 0.3) +
  scale_color_continuous(low = 'blue', high = 'orange1', 
                         name = 'Anthropogenic range contractions',
                         limits = c(0,NA),
                         breaks = c(0, 0.2, 0.4, 0.6),
                         guide = guide_colorbar(title.position = 'top', title.hjust = 0.5,
                                                direction = 'horizontal', 
                                                barwidth = 12, barheight = 0.2)) +
  theme_bw() +
  labs(y = 'Permutation Importance (%)',
       x = 'Mean difference in habitat suitability') +
  theme(legend.position = 'none',
        legend.text = element_text(size = 10),
        panel.grid.minor = element_blank())
pi.plot1

# b)
pi.plot2 <- ggplot(luvarb, aes(x = Truncation, y = PI, colour = ARC)) +
  geom_point(size = 3, alpha = 0.3) +
  scale_color_continuous(low = 'blue', high = 'orange1', 
                         name = 'Anthropogenic range contractions',
                         limits = c(0,NA),
                         breaks = c(0, 0.2, 0.4, 0.6),
                         guide = guide_colorbar(title.position = 'top', title.hjust = 0.5,
                                                direction = 'horizontal', 
                                                barwidth = 12, barheight = 0.2)) +
  theme_bw() +
  labs(y = 'Permutation Importance (%)',
       x = "Corrected 1 - Warren's I") +
  theme(legend.position = 'bottom',
        legend.text = element_text(size = 10),
        panel.grid.minor = element_blank())
pi.plot2
pi.leg <- get_legend(pi.plot2)
pi.plot2 <- pi.plot2 + theme(legend.position = 'none')

pi <- ggarrange(ggarrange(pi.plot1, pi.plot2, labels = c('a)', 'b)')), pi.leg, ncol = 1, heights = c(1,0.2))
pi
#### Fig. S2.6 ################################################################################
## AUC scores for other prev and for model fit ##
mcon <- read.csv('Temporal Bias scores (focc) (all) cont.csv', row.names = 1)
colnames(mcon) <- c('species', colnames(mcon)[2:ncol(mcon)])
auc.df <- mcon[, c(1, grep("AUC", colnames(mcon)))]

auc.df <- melt(auc.df, id.vars = 'species')
auc.df <- cbind(auc.df, str_split_fixed(auc.df$variable, '[.]', 2))
auc.df <- cbind(auc.df, str_split_fixed(auc.df$species, '[.]', 4))
auc.df <- auc.df[, c(6,7,8,9,4,5,3)]
colnames(auc.df) <- c('Cniche', 'Prev', 'Bias', 'Sample',
                      'Model', 'Metric', 'Score')
auc.df$Metric <- factor(auc.df$Metric, levels = c('mAUC', 'aAUC','hAUC'))
auc.df$Bias <- factor(auc.df$Bias, levels = c('HT', 'CA', 'SA', 'SR', 'CR', 'ST',  'CT', 'MT'))
head(auc.df)
amt <- c('Environment Only', 'With Land Use')
names(amt) <- c('EnvO', 'wLand')

auc.plot <- ggplot(auc.df, aes(x = Metric, y = Score, colour = Prev, fill = Prev)) +
  geom_violin(scale = 'width', position = position_dodge(1)) +
  geom_hline(yintercept = 0.7, colour = 'black', size = 1, linetype = 'dashed') +
  facet_grid(~Model, labeller = labeller(Model = amt)) +
  theme_bw() +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.position = 'right',
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        panel.grid.minor.y = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 10, face = 'bold')) +
  scale_x_discrete(labels = c('Model Fit', 'Contemporary', 'Historical'), name = element_blank()) +
  ylim(0.68, 1) +
  labs(colour = 'Prevalence', fill = 'Prevalence')
auc.plot