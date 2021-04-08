#### sampling occurrence points ####
# simulating the entire sampling process #
# sampling for each species set and for each bias
# species statistics and sampling error traits are also 
# being generated here

library(raster)
library(virtualspecies)
library(stringr)
main <- 'D:/OneDrive - National University of Singapore/0 TempBias/'

land.n.dir <- paste(main, 'LandNiche/cont', sep = '') # Land Suitability
lc.dir <- paste(main, 'Env/LUH_sim', sep = '') # Land Cover Predictors (to get true values)
env.n.dir <- paste(main, 'EnvNiche', sep = '') # Env Suitability
sampling.dir <- 'D:/S.Pang/0 TempBias/Samples' # Sampling directory (set local due to number of files)

setwd(sampling.dir)
sp.stat <- as.data.frame(list.files()) # using previous "creating sampling sheet"
colnames(sp.stat) <- 'Name'
sp.stat <- as.data.frame(str_split_fixed(sp.stat$Name, '[.]', 4)[,1:2]) # splitting for EnvNiche n Prev
sp.stat <- unique(sp.stat)
sp.stat[, 3:5] <- NA
colnames(sp.stat) <- c('Env.N', 'Prev', 'Thrs', 'D.1900', 'D.2000')
head(sp.stat)

bias <- c('CA', 'CT', 'CR', 'SA', 'ST', 'SR', 'MT', 'HT')
# Note: new names are, in order above, clustered past, clustered recent, clustered intermediate
# spread past, spread recent, spread intermediate, only end, and only start
env.niche <- paste('N', seq(1, 100, 1), sep = '')
prev.lvl <- c('P05', 'P10', 'P15', 'P25', 'P35', 'P50', 'P65')
prev <- c(0.05, 0.10,0.15,0.25,0.35,0.50,0.65)
# prevalence separated for inclusion in virtual species and for labelling

obj.prevalence <- function(obj){
  obj <- freq(obj)
  if(nrow(obj) == 3) {
    return(obj[2,2]/(obj[1,2] + obj[2,2]))
  } else {
    return(0)
  }
} # function to determine prevalence (binary distribution)

for (N in env.niche) {
  for (P in 1:length(prev)) {
    p <- prev[P] # number - small p
    P <- prev.lvl[P] # label - big P
    
    # getting species statistics, based on 1900 #
    rn <- sp.stat$Env.N == N & sp.stat$Prev == P # row number
    setwd(env.n.dir)
    Env.N <- raster(paste(N, '.tif', sep = '')) # Env Niche
    setwd(land.n.dir)
    ld <- raster(paste('Lsuit_cont_',1900,".tif", sep = "")) # historical distribution
    PA <- convertToPA(Env.N*ld,
                      PA.method = 'probability', plot = F,
                      alpha = -0.01, species.prevalence = p)
    # Prevalence is used to determine threshold which is kept constant for the species
    # Prevalence is set at historical 1900
    sp.beta <- as.numeric(PA$PA.conversion['beta']) # threshold
    sp.stat$Thrs[rn] <- sp.beta # storing threshold value
    PA <- convertToPA(Env.N*ld,
                      PA.method = 'threshold', plot = F,
                      beta = sp.stat$Thrs[rn])
    # creating species presence absence based on threshold (consistency)
    sp.stat$D.1900[rn] <- obj.prevalence(PA$pa.raster) # prevalence of species for 1900
    
    
    # for 2000 #
    ld <- raster(paste('Lsuit_cont_',2000,".tif", sep = "")) # 2000 distribution
    PA <- convertToPA(Env.N*ld,
                      PA.method = 'threshold', plot = F,
                      beta = sp.stat$Thrs[rn])
    sp.stat$D.2000[rn] <- obj.prevalence(PA$pa.raster)
  }
}
setwd(main)
setwd('Results')
write.csv(sp.stat, file = 'species statistics.csv')


#### simulating observe model and sampling of occurrence data ####
library(foreach)
library(doParallel)

main <- 'D:/OneDrive - National University of Singapore/0 TempBias'
main.l <- 'D:/S.Pang/0 TempBias'

setwd(main)
setwd('LandNiche')
setwd('cont')
alc <- stack(paste('Lsuit_cont_', 1900:2000, '.tif', sep = ''))
names(alc) <- sub('Lsuit_cont_', 'L', names(alc))

setwd(main)
setwd('EnvNiche')
cn <- stack(paste('N', 1:100, '.tif', sep = ''))
names(cn)

setwd(main)
setwd('Results')
sp.stats <- read.csv('species statistics.csv')
head(sp.stats)

# different temporal sampling patterns to be used #
bias <- c('CA', 'CT', 'CR', 'SA', 'ST', 'SR', 'MT', 'HT')
bias.list <- c('ClusA', 'ClusT', 'ClusR', 'SlopA', 'SlopT', 'SlopR', 'Match', 'Hist')
# Note: new names are in sequence, clustered past, clustered recent, clustered intermediate
# spread past, spread recent, spread intermediate, only end, and only start
temporal.sampling <- function(freq.table = NULL, n = 100) {
  sort(sample(freq.table$year, replace = T, prob = freq.table$Freq, n))
}

Prev <- paste('P', c('10', '15', '25', '35', '50'), sep = '') # labelling purpose only

cores=detectCores()
cl <- makeCluster(10)
registerDoParallel(cl)

foreach(S = 1:15, .packages = 'raster') %dopar% { # 15 replicates as some issues was observed for maxnet
  #### Random sampling of points with probabilities ##
  ## random sampling of points across region ##
  # also extracting the env suitability only #
  set.seed(2102 + S)
  S <- paste('S', S, sep = '')
  prob.samp <- as.data.frame(sampleRandom(cn, 2000, xy = T)) # cn is the environmental suitability
  
  for (B in 1:length(bias)) {
    ## setting different sampling years for each point ##
    B.mod <- bias.list[B]
    B <- bias[B]
    setwd(main)
    setwd('Temporal Sampling')
    mod.bias <- read.csv(paste('Bias_',B.mod, '.csv', sep = "")) # temporal sampling pattern file
    samp.B <- prob.samp
    samp.B$Year <- sample(temporal.sampling(mod.bias, nrow(samp.B)))
    samp.B <- samp.B[,c(1,2,103, 3:102)]
    
    ## extracting the land cover suit for each point ##
    alc.val <- cbind(samp.B[, c('x', 'y')], extract(alc, samp.B[, c('x', 'y')])) # alc is the actual land use map
    
    ## adjusting env suitability with lc suit for that year ##
    for (pt in 1:nrow(samp.B)) {
      samp.B[pt, paste('N', 1:100, sep = '')] <- samp.B[pt, paste('N', 1:100, sep = '')]*alc.val[pt, paste('L', samp.B[pt,'Year'], sep = '')]
    }
    
    samp.P <- samp.B[,c('x','y','Year')]
    samp.PA <- samp.B[,c('x','y','Year')]
    for (N in paste('N', 1:100, sep = '')) {
      for (P in Prev) {
        TH <- sp.stats[sp.stats$Env.N == N & sp.stats$Prev == P, 'Thrs'] # threshold for binary data
        samp.P[, paste(N, P, sep = '.')] <- samp.B[, N]
        samp.P[, paste(N, P, sep = '.')][samp.P[, paste(N, P, sep = '.')] <= TH] <- 0
        samp.PA[, paste(N, P, sep = '.')] <- rbinom(nrow(samp.P), 1, samp.P[, paste(N, P, sep = '.')])
      } # prevalence end loop
    } # env niche end loop
    
    setwd(main.l)
    setwd('Samples Prep')
    write.csv(samp.B, file = paste('RawProb', B, S, 'csv', sep = '.'), row.names = F)
    write.csv(samp.P, file = paste('Prob', B, S, 'csv', sep = '.'), row.names = F)
    write.csv(samp.PA, file = paste('PA', B, S, 'csv', sep = '.'), row.names = F)
  } # bias end loop
} # replicates end loop
stopCluster(cl)


#### converting to sampling datasets ####
Prev <- c('P10', 'P15', 'P25', 'P35', 'P50')

for (B in bias) {
  for (S in paste('S', 1:15, sep = '')) {
    setwd(main.l)
    setwd('Samples Prep')
    
    PA <- read.csv(paste('PA', B, S, 'csv', sep = '.'))
    prob <- read.csv(paste('Prob', B, S, 'csv', sep = '.'))
    
    for (N in paste('N', 1:100, sep = '')) {
      for (P in Prev) {
        
        samp <- prob[PA[,paste(N, P, sep = '.')] == 1 ,c('Year', 'x', 'y', paste(N, P, sep = '.'))]
        colnames(samp) <- c('Year', 'Lon', 'Lat', 'Sprob')
        
        setwd(main.l)
        setwd('Samples (focc)') # focc is fixed sampling effort of occurrence data
        write.csv(samp, file = paste(N, P, B, S, 'csv', sep = '.'))
      } # prevalence end loop
    } # environmental niche end loop
  } # sampling replicate end loop
} # bias end loop

