#### examining the temporal sampling patterns within real species occurrence datasets ####
library(ggplot2)

main <- 'D:\\Work Sch\\Phd\\Databases Occurrence\\GBIF Plant\\SEA (all angio)'
setwd(main)

raw <- read.csv('GBIF_SEA_allAngio.csv')
colnames(raw)
raw$species <- as.character(raw$species)
raw$basisOfRecord <- as.character(raw$basisOfRecord)
raw.sub <- raw[!raw$species == "",]
unique(raw.sub$basisOfRecord)

# only preserved specimen
raw.sub <- raw.sub[raw.sub$basisOfRecord == "PRESERVED_SPECIMEN",]
colnames(raw.sub)
obs.bias <- raw.sub[, c(10,16, 22,23, 33)]
write.csv(obs.bias, file = 'Yearly Distribution.csv')

obs.bias <- read.csv('Yearly Distribution.csv', row.names = 1)
obs.bias$species <- as.character(obs.bias$species)
unique(obs.bias$species)
sp.freq <- as.data.frame(table(obs.bias$species))
sp.freq <- sp.freq[sp.freq$Freq >= 100 ,]
table(sp.freq$Freq)
mean(sp.freq$Freq)

species <- as.character(sp.freq$Var1)

obs.bias.s <- obs.bias[obs.bias$species %in% species,]
abs.freqbias <- as.data.frame(obs.bias.s$year[!is.na(obs.bias.s$year)])
colnames(abs.freqbias) <- 'Var1'
abs.freqbias$Var1 <- as.numeric(as.character(abs.freqbias$Var1))
colnames(abs.freqbias) <- 'Var1'
abs.bias <- as.data.frame(table(obs.bias.s$year))
obs.bias.s <- as.data.frame(table(obs.bias.s[, c(1,5)]))
obs.bias.studyperiod <- obs.bias.s[obs.bias.s$year %in% 1900:2000,]
summarySE <- function(data = NULL, measurevar, groupvars = NULL, na.rm=FALSE,
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
                     mean = median (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- plyr::rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

agg.bias <- summarySE(obs.bias.s, 'Freq', 'year')

for (sp in species) {
  obs.bias.s[ obs.bias.s$species == sp,]$Freq <- obs.bias.s[obs.bias.s$species == sp,]$Freq / sum(obs.bias.s[ obs.bias.s$species == sp,]$Freq)
}
for (sp in species) {
  obs.bias.studyperiod[obs.bias.studyperiod$species == sp, ]$Freq <- obs.bias.studyperiod[obs.bias.studyperiod$species == sp,]$Freq / sum(obs.bias.studyperiod[obs.bias.studyperiod$species == sp, ]$Freq)
}
scl.bias <- summarySE(obs.bias.s, 'Freq', 'year') # scaled frequencies
raw.bias <- as.data.frame(table(obs.bias[,c(1,5)])) # raw frequencies
scl.raw.bias <- summarySE(raw.bias, 'Freq', 'year') # scaled raw frequencies
raw.bias <- as.data.frame(table(obs.bias$year))
sp.scl <- summarySE(obs.bias.s, 'Freq', c('year', 'species')) # scale by species
eg.sps <- summarySE(obs.bias.studyperiod, 'Freq', c('year', 'species')) # example species to examine for study period
library(ggplot2)
ggplot(agg.bias, aes(x=year, y=Freq, group = 1)) +
  geom_line() +
  scale_x_discrete(breaks = seq(1700, 2020, 10))

ggplot(scl.bias, aes(x=year, y=Freq, group = 1)) +
  geom_line() +
  scale_x_discrete(breaks = seq(1700, 2020, 10))

################ absolute temporal sampling bias to be examined ####
abs.bias$Var1 <- as.numeric(as.character(abs.bias$Var1))
write.csv(abs.bias, "abs.bias.csv")

ggplot(abs.freqbias, aes(x=Var1)) +
  geom_freqpoly(bins = 66) +
  scale_x_continuous(breaks = seq(1770, 2010, 10)) +
  geom_vline(xintercept = c(1900, 2000), col = 'red')

ggplot(raw.bias, aes(x=Var1, y=(Freq), group = 1)) +
  geom_line() +
  scale_x_discrete(breaks = seq(1700, 2020, 10))

sp.scl$year <- as.numeric(as.character(sp.scl$year))
eg.sps$year <- as.numeric(as.character(eg.sps$year))
write.csv(eg.sps, file = 'eg.sps.csv')

#### example species with clear temporal sampling patterns ###############################################
## examples that were used in the main text Fig. 1
eg.sp <- eg.sps[eg.sps$species %in% c("Cinnamomum parthenoxylon", "Croton tiglium", "Knema latifolia",
                                      "Ficus benguetensis", "Daphniphyllum glaucescens", "Cinnamomum verum"),]
write.csv(eg.sp, "TSP example.csv", row.names = F)
ggplot(eg.sp, aes(x=year, y=Freq, group = 1)) +
  geom_line() +
  scale_x_continuous(breaks = seq(1900, 2000, 20), limits = c(1900,2000)) +
  facet_wrap(~species, scales = 'free_y')

setwd('D:/OneDrive - National University of Singapore/0 TempBias/Temporal Sampling')
write.csv(scl.bias[,c(1,3)], file = 'Obs_bias.csv')

#### simulating the temporal sampling patterns to be used for the observer model #########################
## temporal sampling function ##
temporal.sampling <- function(freq.table = NULL, n = 100) {
  sort(sample(freq.table$year, replace = T, prob = freq.table$Freq, n))
}

## Creating Bias distribution ############################################################################
temp.dir <- ('D:/OneDrive - National University of Singapore/0 TempBias/Temporal Sampling')
setwd(temp.dir)
## Observed temporal distribution of samples
yr <- 1900:2000
obser <- read.csv('Obs_bias.csv', row.names = 1)
obser <- obser[obser$year %in% yr,]
#### Clustered Away sampling on 1919 ## now clustered past
nonlin <- function(t, a, b, c) { (a/b) * exp( -(t-c)/b - exp(-(t-c)/b) ) }
ClusA <- obser
ClusA$Freq <- nonlin(yr, 1, 1, 1919)
plot(ClusA$Freq)
write.csv(ClusA, file = 'Bias_ClusA.csv', row.names = F)
#### Clustered Towards sampling on 1990 ## now clustered recent
ClusT <- obser
ClusT$Freq <- nonlin(yr, 1, 1, 1990)
plot(ClusT$Freq)
write.csv(ClusT, file = 'Bias_ClusT.csv', row.names = F)
#### Clustered Regular sampling on 1964 ## now clustered intermediate
ClusR <- obser
ClusR$Freq <- nonlin(yr, 1, 1, 1964)
plot(ClusR$Freq)
write.csv(ClusR, file = 'Bias_ClusR.csv', row.names = F)

#### Sloped Away sampling ## now spread past
SlopA <- obser
SlopA$Freq <- seq(0.2, 0, length.out = nrow(SlopA))
plot(SlopA$Freq)
write.csv(SlopA, file = 'Bias_SlopA.csv', row.names = F)
#### Sloped Towards sampling ## now spread recent
SlopT <- obser
SlopT$Freq <- seq(0, 0.2, length.out = nrow(SlopT))
plot(SlopT$Freq)
write.csv(SlopT, file = 'Bias_SlopT.csv', row.names = F)
#### Sloped Regular sampling ## now spread intermediate
SlopR <- obser
SlopR$Freq <- 0.2
plot(SlopR$Freq)
write.csv(SlopR, file = 'Bias_SlopR.csv', row.names = F)
#### Match sampling ## now Only End
Mat <- obser
Mat$Freq <- 0
Mat$Freq[Mat$year == 2000] <- 1
plot(Mat$Freq)
write.csv(Mat, file = 'Bias_Match.csv', row.names = F)
#### Historical sampling ## now Only Start
Hist <- obser
Hist$Freq <- 0
Hist <- rbind(c(1900, 1),Hist)
plot(Hist$Freq)
write.csv(Hist, file = 'Bias_Hist.csv', row.names = F)

#### sampling matrix ##
bias.list <- c('ClusA', 'ClusT', 'ClusR', 'SlopA', 'SlopT', 'SlopR', 'Match', 'Hist')
# Note: new names are in sequence, clustered past, clustered recent, clustered intermediate
# spread past, spread recent, spread intermediate, only end, and only start
bias.mat <- as.data.frame(matrix(NA, 100, 8))
colnames(bias.mat) <- bias.list
for (B in bias.list) {
  bias <- read.csv(paste('Bias_',B, '.csv', sep = ""))
  bias.mat[,B] <- temporal.sampling(bias, 100)
}
write.csv(bias.mat, file = 'bias_matrix.csv', row.names = F)

#### making sample example ####
nonlin <- function(t, a, b, c) { (a/b) * exp( -(t-c)/b - exp(-(t-c)/b) ) }
ClusA <- obser
ClusA$Freq <- nonlin(yr, 1, 1, 1919)
#### Clustered Towards sampling on 1990 ##
ClusT <- obser
ClusT$Freq <- nonlin(yr, 1, 1, 1990)
#### Clustered Regular sampling on 1964 ##
ClusR <- obser
ClusR$Freq <- nonlin(yr, 1, 1, 1964)

#### Sloped Away sampling ##
SlopA <- obser
SlopA$Freq <- seq(0.04, 0, length.out = nrow(SlopA))
#### Sloped Towards sampling ##
SlopT <- obser
SlopT$Freq <- seq(0, 0.04, length.out = nrow(SlopT))
#### Sloped Regular sampling ##
SlopR <- obser
SlopR$Freq <- 0.02

bias.plot <- ClusA
colnames(bias.plot) <- c('year', 'ClusA')
bias.plot$ClusT <- ClusT$Freq
bias.plot$ClusR <- ClusR$Freq
bias.plot$SlopA <- SlopA$Freq
bias.plot$SlopT <- SlopT$Freq
bias.plot$SlopR <- SlopR$Freq
colnames(bias.plot) <- c('year', 'CP', 'CR', 'CI', 'SP', 'SR', 'SI')

bias.plot <- melt(bias.plot, id.vars = 'year')
bias.plot$variable <- factor(bias.plot$variable, levels = c('SP', 'CP', 'SI', 'CI', 'SR', 'CR'))
write.csv(bias.plot, "sim.freq.csv", row.names = F)
