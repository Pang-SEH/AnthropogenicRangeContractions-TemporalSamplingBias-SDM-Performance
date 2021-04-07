library(virtualspecies)
library(biomod2)
library(PresenceAbsence)
main <- 'D:/OneDrive - National University of Singapore/0 TempBias/'
main.l <- 'D:/S.Pang/0 TempBias/'

clim.n <- paste(main, 'ClimNiche', sep = '') # environmental niche
land.n <- paste(main, 'LandNiche/cont', sep = '') # land use suitability
sp.stat <- read.csv(paste(main, 'Results/species statistics.csv', sep = '')) # thresholds for sp
sp.stat <- sp.stat[!grepl('P05|P65', sp.stat$Prev),] # thresholds to examine removing previously included prevalences

samp.dir <- paste(main.l, 'Samples (focc)', sep = '') # local storage of sample sheet (too many files)
op <- paste(main, 'Results n/', 'Temporal Bias sub res', sep = '') # output
ip <- paste(main.l, 'Models (focc)', sep = '') # input
res.dir <- paste(main, 'Results n', sep = '') # results directory

#### functions for score assessments ####
mod.true.performance <- function(VirtSpOutputFile, Pred.D, Pred.th, prefix = NULL) {
  require(virtualspecies)
  require(PresenceAbsence)
  # Virtual Species Output File must have both the pa.raster (binary distribution)
  # and the suitab.raster (probabilitic distribution)
  # formating the bi = binary distribution for PresenceAbsence
  bi <- VirtSpOutputFile$pa.raster; bi <- as.data.frame(bi); bi$pixel <- rownames(bi)
  bi <- bi[!is.na(bi$layer),]; bi <- bi[,c(2,1)]
  
  P.df <- as.data.frame(Pred.D); P.df$pixel <- rownames(P.df); P.df <- P.df[!is.na(P.df[,1]),]
  P.df <- P.df[,c(2,1)]; colnames(P.df) <- c('pixel', 'proj')
  
  tp.df <- merge.data.frame(bi, P.df)
  colnames(tp.df) <- c('plotID', 'Observed', 'Predicted1')
  
  s <- as.data.frame(matrix(NA, 1, 9))
  # the scores can be changed and is flexible
  colnames(s) <- c('AUC', 'TSS', 'POD', 'SR', 'POFD', 'FBI', 'FAP', 'FAS', 'NOI')
  df.score <- presence.absence.accuracy(tp.df, 
                                       threshold = Pred.th, 
                                       st.dev = F)
  # a = hits, b = false alarm, c = miss, d = correct negatives
  a <- sum(tp.df$Observed == 1 & tp.df$Predicted1 >= Pred.th)
  c <- sum(tp.df$Observed == 1 & tp.df$Predicted1 < Pred.th)
  b <- sum(tp.df$Observed == 0 & tp.df$Predicted1 >= Pred.th)
  d <- sum(tp.df$Observed == 0 & tp.df$Predicted1 < Pred.th)
  
  
  s[1, 'AUC'] <- df.score$AUC
  s[1, 'TSS'] <- a/(a+c) + d/(b+d) - 1
  s[1, 'POD'] <- a/(a+c)
  s[1, 'SR'] <- a/(a+b)
  s[1, 'POFD'] <- b/(b+d)
  s[1, 'FBI'] <- (a+b)/(a+c)
  
  ## made for comparisons using probability only ##
  # heavier focus on the probabilities overall (decernibility across env gradient)
  # may inflate performance for more range restricted species (more negatives to match)
  f.df <- as.data.frame(VirtSpOutputFile$suitab.raster * VirtSpOutputFile$pa.raster)
  f.df$pixel <- rownames(f.df); f.df <- f.df[!is.na(f.df[,1]),]; f.df <- f.df[, c(2,1)]
  cor.df <- merge.data.frame(f.df, P.df)
  # scores for continuous measures of probabilities (Functional Accuracy Pearsons)
  # and for (Niche Overlap Index) using default warren's adjusted Hellinger's Distance (I)
  s[1, 'FAP'] <- cor(cor.df[,2:3], use = 'pairwise.complete.obs', method = 'pearson')[1,2]
  s[1, 'FAS'] <- cor(cor.df[,2:3], use = 'pairwise.complete.obs', method = 'spearman')[1,2]
  s[1, 'NOI'] <- dismo::nicheOverlap(VirtSpOutputFile$suitab.raster * VirtSpOutputFile$pa.raster, Pred.D, checkNegatives = F)
  colnames(s) <- paste(prefix, colnames(s), sep = '')
  
  # ## made for comparisons using binary X probability ##
  # # heavier focus on the probabilities within and not outside
  # # (inflate errors of over prediction and vice versa)
  # f.df <- as.data.frame(VirtSpOutputFile$suitab.raster * VirtSpOutputFile$pa.raster)
  # f.df$pixel <- rownames(f.df); f.df <- f.df[!is.na(f.df[,1]),]; f.df <- f.df[, c(2,1)]
  # cor.df <- merge.data.frame(f.df, P.df)
  # s[1, 'FAP'] <- cor(cor.df[,2:3], use = 'pairwise.complete.obs', method = 'pearson')[1,2]
  # s[1, 'NOI'] <- dismo::nicheOverlap(VirtSpOutputFile$suitab.raster * VirtSpOutputFile$pa.raster, Pred.D, checkNegatives = F)
  # colnames(s) <- paste(prefix, colnames(s), sep = '')
  
  return(s)
}

## species list ##
setwd(paste(main, 'Env/pca', sep = ''))
bg.data <- read.csv('bg.data.csv')
bias <- c('CA', 'CT', 'CR', 'SA', 'ST', 'SR', 'MT', 'HT')
sims <- c('wLand', 'EnvO')
land.format <- c('cont') # require for previous iterations of the analysis left as it is for convinence
#### starting loop ####
library(foreach)
library(doParallel)
cl <- makeCluster(detectCores()-2)
registerDoParallel(cl)

num.col <- 2*(4 + 2*9) # num of col per sampling score sheet with all 2 simulations

for (L in land.format) {
  # for each land cover map format as either cont or cate
  foreach (i = 1:nrow(sp.stat), .packages = c('virtualspecies', 'biomod2', 'PresenceAbsence')) %dopar% {
    TH <- sp.stat$Thrs[i] # binary range threshold
    N <- sp.stat$Clim.N[i] # virtual species environmental niche
    P <- sp.stat$Prev[i] # prevalence
    setwd(clim.n)
    Clim.N <- raster(paste(N, '.tif', sep = ''))
    setwd(land.n) # land cover suitability, NOT FORMAT
    Land.N <- raster(paste("Lsuit_cont_",2000,".tif", sep = ""))
    # actual is land use for 2000, contemporary/end
    actual <- convertToPA(Clim.N*Land.N,
                          PA.method = 'threshold', plot = F,
                          beta = TH)
    Land.N <- raster(paste("Lsuit_cont_",1900,".tif", sep = ""))
    # historical is land use for 1900, historical/start
    historical <- convertToPA(Clim.N*Land.N,
                              PA.method = 'threshold', plot = F,
                              beta = TH)
    
    sp.res <- as.data.frame(matrix(NA, 0, num.col)) # species results for easier compilation later
    for (B in bias) {
      for (S in paste('S', 1:15, sep = "")) {
        sp <- paste(N, P, B, S, sep = ".")
        setwd(samp.dir) # for each sampling dataset
        if(file.exists(paste(sp, 'csv', sep = '.'))) {
          occ.pt <- read.csv(paste(sp, 'csv', sep = '.'))[,c('Lon', 'Lat')]
          
          # obj df to store scores (model fit) (prediction) (projection) for each sp
          # and for each of the simulation (modelling procedure including controls)
          sp.row <- as.data.frame(matrix(NA, 1, 0))
          for (sim in sims) {
            setwd(ip) # input folder
            setwd(L) # land cover format used
            setwd(sp) # the species folder
            setwd(sim) # the simulations folder (procedure)
            mod.th <- read.csv('Thresholds.csv') # threshold as determined by the model
            mod.th <- mod.th[mod.th$Threshold == 'Maximum training sensitivity plus specificity',
                             'Cloglog.value']
            mod.proj.a <- NA # model prediction a = actual (for evaluation)
            mod.proj.h <- NA # model projection h = historical (for evaluation)
            if (sim == 'wLand' | sim == 'nDisc' | sim == 'rPres'){
              mod.proj.a <- raster('Pred_2000_cloglog.tif')
              mod.proj.h <- raster('Pred_1900_cloglog.tif')
            }
            if (sim == 'EnvO' | sim == 'rLand') { # single outputs (landcover "unchanged")
              mod.proj.a <- raster('Pred_cloglog.tif')
              mod.proj.h <- raster('Pred_cloglog.tif')
            }
            
            confusion.matrix <- read.csv('confusion_matrix.csv', row.names = 1)['Maximum training sensitivity plus specificity',]
            colnames(confusion.matrix) <- c('th', 'a', 'b', 'c', 'd')
            
            mod.scores <- read.csv('train_score.csv')
            mod.scores$mPOD <- confusion.matrix$a / (confusion.matrix$a + confusion.matrix$c)
            mod.scores$mSR <- confusion.matrix$a / (confusion.matrix$a + confusion.matrix$b)
            mod.scores$mPOFD <- confusion.matrix$b / (confusion.matrix$b + confusion.matrix$d)
            colnames(mod.scores) <- paste(sim, colnames(mod.scores), sep = '.')
            
            mod.perf <- mod.true.performance(actual, mod.proj.a, mod.th, 'a') # prediction
            colnames(mod.perf) <- paste(sim, colnames(mod.perf), sep = '.')
            
            mod.hist <- mod.true.performance(historical, mod.proj.h, mod.th, 'h') # projection
            colnames(mod.hist) <- paste(sim, colnames(mod.hist), sep = '.')
            
            sp.row <- cbind(sp.row, mod.scores, mod.perf, mod.hist) # scores + performance * 2 = per sim
          }
          row.names(sp.row) <- sp
          sp.res <- rbind(sp.res, sp.row)
        }
      }
    }
    setwd(op)
    write.csv(sp.res, file = paste(N,P,L,'csv', sep = '.')) # sub outputs for compile later
  }
  setwd(op)
  sp.mod.res.list <- list.files(pattern = L) # list to compile by lc format
  mod.res <- as.data.frame(matrix(NA, 0, num.col))
  for (i in 1:length(sp.mod.res.list)) {
    df <- read.csv(sp.mod.res.list[i])
    mod.res <- rbind(mod.res, df)
  }
  setwd(res.dir)
  # model scores for this land cover format
  write.csv(mod.res, file = paste('Temporal Bias scores (focc) (all) ', L, '.csv', sep = ''), row.names = T)
}
stopCluster(cl)
