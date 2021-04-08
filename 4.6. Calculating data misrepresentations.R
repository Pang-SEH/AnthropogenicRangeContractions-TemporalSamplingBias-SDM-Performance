#### Analysis Sampling Set Error Rates ####
library(raster)
library(virtualspecies)
library(ecospat)
library(stringr)
main <- 'D:/OneDrive - National University of Singapore/0 TempBias/'

land.n <- paste(main, 'LandNiche/cont', sep = '')
env.n <- paste(main, 'EnvNiche', sep = '')
setwd(paste(main, 'Env/pca', sep = ''))
env <- stack('PCs_all.grd')

lc.dir <- paste(main, 'Env/LUH_sim', sep = '')
sampling.dir <- 'D:/S.Pang/0 TempBias/Samples (focc)'
sampling.prep <- 'D:/S.Pang/0 TempBias/Samples Prep'
setwd(main)
setwd('Results')
res.dir <- getwd()
sp.dist <- read.csv('Distribution Range Changes.csv', row.names = 1) # to measure truncation
colnames(sp.dist) <- sub('X', '', colnames(sp.dist)) # saving colnames as numbers will add X to front
sp.stat <- read.csv('species statistics.csv', row.names = 1)
sp.stat <- sp.stat[!grepl('P05|P65', sp.stat$Prev),]

### this function calculates 2D environmental niche overlaps for all
### possible paired combination of environmental variables (5 PC-axes)
Env.Niche.Overlap <- function(env, sp.range, sp.pts) {
  niche.d <- as.data.frame(stack(env, sp.range), na.rm = T)
  pt.d <- as.data.frame(extract(env, sp.pts[, c('Lon', 'Lat')]))
  pt.d$layer.1 <- 1
  
  require(combinat)
  comb <- combn(1:nlayers(env), 2)
  NOI.mat <- as.data.frame(matrix(NA, 2, ncol(comb)))
  row.names(NOI.mat) <- c('D', 'I')
  
  for (i in 1:ncol(comb)) {
    grid.clim.range <- ecospat.grid.clim.dyn(glob=niche.d[,comb[,i]],
                                             glob1=niche.d[which(niche.d[,6]==1), comb[,i]],
                                             sp=niche.d[which(niche.d[,6]==1), comb[,i]], R=100,
                                             th.sp=0)
    # gridding the invasive niche
    grid.clim.pts <- ecospat.grid.clim.dyn(glob=niche.d[, comb[,i]],
                                           glob1=niche.d[which(niche.d[,6]==1), comb[,i]],
                                           sp=pt.d[, comb[,i]], R=100,
                                           th.sp=0)
    
    # Compute Schoener's D, index of niche overlap
    D.overlap <- ecospat.niche.overlap (grid.clim.range, grid.clim.pts, cor=T)
    NOI.mat['D', i] <- D.overlap$D
    NOI.mat['I', i] <- D.overlap$I
  }
  return(NOI.mat)
}

library(foreach)
library(doParallel)
cores=detectCores()
cl <- makeCluster(cores[1]-2)
registerDoParallel(cl)

Sys.time()
foreach (i = 1:nrow(sp.stat), .packages = c('virtualspecies', 'ecospat', 'stringr')) %dopar% {
  TH <- sp.stat$Thrs[i] # specific threshold calculated before
  N <- as.character(sp.stat$Env.N[i]) # Env Niche
  P <- as.character(sp.stat$Prev[i]) # Prevalence level set
  setwd(env.n)
  Env.N <- raster(paste(N, '.tif', sep = '')) # pulling out Env Niche
  setwd(land.n)
  Land.N <- raster(paste("Lsuit_cont_",1900,".tif", sep = ""))
  historical <- convertToPA(Env.N*Land.N,
                            PA.method = 'threshold', plot = F,
                            beta = TH)
  
  # extracting the distribution changes over time for this species
  dist.shf <- sp.dist[sp.dist$Env.N == N & sp.dist$Prev == P, ]
  
  setwd(sampling.dir) # listing out all samples of this species
  samp.list <- list.files(pattern = paste(N,P, sep = '.'))
  samp.list <- sub('.csv', '', samp.list)
  
  er <- as.data.frame(matrix(NA, length(samp.list), 6))
  rownames(er) <- samp.list
  colnames(er) <- c('T2ER', 'Discrepancy', 'Yearly', 'Truncation', 'Trunc.D', 'Contraction')
  # list of actual and proxies to decern errors
  for (samp in samp.list) {
    setwd(sampling.dir)
    df <- read.csv(paste(samp, '.csv', sep = ''), row.names = 1)
    
    S <- str_split_fixed(samp, '[.]', 4)[4]
    setwd(sampling.prep)
    ht.df <- read.csv(paste('PA.MT', S, 'csv', sep = '.'))[, c('x', 'y', paste(N, P, sep = '.'))]
    colnames(ht.df) <- c('Lon', 'Lat', paste(N, P, 'PA', sep = '.'))
    df <- merge.data.frame(ht.df, df)
    ht.df <- read.csv(paste('RawProb.MT', S, 'csv', sep = '.'))[, c('x', 'y', N)]
    colnames(ht.df) <- c('Lon', 'Lat', paste(N, 'Prob', sep = '.'))
    df <- merge.data.frame(ht.df, df)
    rm(ht.df)
    
    ### several metrics that were initial considered to measure
    ### occurrence-habitat mismatching and niche truncation
    # number of false positives, presence points at the time of sampling but should be absence now
    er[samp, 'T2ER'] <- sum(df[,paste(N, P, 'PA', sep = '.')] == 0) / nrow(df)
    # the average misrepresentation of an occurrence's suitability value (or sampling probability)
    er[samp, 'Discrepancy'] <- mean(df$Sprob - df[, paste(N, 'Prob', sep = '.')])
    # the average Year for which occurrences come from
    er[samp, 'Yearly'] <- round(mean(df$Year))
    # average proportion of species range that could no longer be sampled at time of sampling
    DT.val <- 0
    for (y in df$Year) {
      DT.val <- DT.val + dist.shf[1, as.character(y)]
    }
    er[samp, 'Contraction'] <- DT.val / 100
    eno.mat <- Env.Niche.Overlap(env, historical$pa.raster, df)
    er[samp, 'Trunc.D'] <- 1 - rowMeans(eno.mat['D',])
    er[samp, 'Truncation'] <- 1 - rowMeans(eno.mat['I',])
    
  }
  setwd(res.dir)
  setwd('Error Rates sub res')
  write.csv(er, file = paste(N,P,'ER','csv', sep = '.'))
} # sp set is stored and save separately first, foreach efficiency and RAM space
stopCluster(cl)
Sys.time()

setwd(res.dir)
setwd('Error Rates sub res')
er.list <- list.files(pattern = 'ER') # loading each matrix of sampling error per sp
sp.er <- as.data.frame(matrix(NA, 0, 6))
colnames(sp.er) <- c('T2ER', 'Discrepancy', 'Yearly', 'Truncation', 'Trunc.D', 'Contraction')
for (i in 1:length(er.list)) {
  df <- read.csv(er.list[i], row.names = 1)
  sp.er <- rbind(sp.er, df)
}

setwd(res.dir)
write.csv(sp.er, file = 'Error Rates in Sampling Data (focc) (all).csv', row.names = T)
Sys.time()
