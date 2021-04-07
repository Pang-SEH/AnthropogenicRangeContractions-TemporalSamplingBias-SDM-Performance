#### Creating specific Niche types
library(raster)
library(virtualspecies)
library(data.table)
library(foreach)
library(doParallel)

base.1 <- 'D:/OneDrive - National University of Singapore/0 TempBias/'
L.dir <- 'D:/S.Pang/0 TempBias'
main <- base.1

### pc axes for environmental variables ##
pca.dir <- paste(main, 'Env/pca', sep = '')
setwd(pca.dir)
PCs <- stack('PCs_all.grd')

#### Creating land use niches ####
landcat <- paste(main, 'Env', 'LUH_sim', sep = '/')
land.n.dir <- paste(main, 'LandNiche', sep = '')
## land suitability functions ##
setwd(landcat)
setwd('cont')
luh_cont <- raster('LUH_cont_2000.tif')
names(luh_cont) <- 'luh_cont'

formatFunctions(x = luh_cont,
                luh_cont = c(fun = 'logisticFun',
                             alpha = 0.1, beta = 0.5))
suit_lg <- function(x) {
  1 / (1 + (exp((x - 0.5)/0.1)))
}

plot(suit_lg(luh_cont))

## land use suitability map for each year

years <- c(seq(1500, 1880, 20), 1900:2000)
for (yr in years) {
  setwd(landcat)
  setwd('cont')
  r.cont <- raster(paste('LUH_cont_', yr, '.tif', sep = ''))
  r.cont <- suit_lg(r.cont)
  setwd(land.n.dir)
  setwd('cont')
  writeRaster(r.cont, paste('Lsuit_cont_', yr, '.tif', sep = ''), overwrite = T)
}

#### Creating Climatic Niches ####
setwd(L.dir)
setwd('EnvNicheObjects')
obj.dir <- getwd()
setwd(L.dir)
setwd('EnvNichePrediction')
en.dir <- getwd()

cores=detectCores()
cl <- makeCluster(cores[1]-2)
registerDoParallel(cl)

# important for realistic.sp = T, which sequentially creates specie response curves to prevent unrealistic species
# i.e., no overlap between variables, see ?generateRandomSp for more detail
foreach (num = 1:10000, .packages = 'virtualspecies') %dopar% {
  Clim.N <- generateRandomSp(PCs, approach = 'response', 
                             rescale = T, rescale.each.response = T,
                             realistic.sp = T, species.type = 'multiplicative',
                             niche.breadth = 'wide',
                             sample.points = F, nb.points = 10000,
                             convert.to.PA = F, PA.method = 'probability', alpha = -0.1, adjust.alpha = F,
                             beta = 'probability', species.prevalence = NULL,
                             plot = T)
  convertToPA(Clim.N, PA.method = 'threshold', beta = "random", alpha = -0.1)
  # setting limits to exclude extreme niches where almost all of SEA is suitable or unsuitable
  # this prevents an excessive amount of extreme niches which slows down later analyses
  if (freq(Clim.N$suitab.raster >= 0.7)[2,2] > 1000 & freq(Clim.N$suitab.raster <= 0.3)[2,2] > 500) {
    setwd(obj.dir)
    save(Clim.N, file = paste('niche', num, sep = ""))
    setwd(en.dir)
    writeRaster(Clim.N$suitab.raster, filename = paste('niche', num, '.tif', sep = ''))
    print(paste("Saved", "Niche", num))
  }
}
stopCluster(cl)

#### Selecting Climatic Niches ####
setwd(en.dir)
niches <- stack(list.files(pattern = 'niche'))
library(usdm)
cn_cor <- vifcor(niches, th = 0.85)
cn_cor
save(cn_cor, file = 'cn_cor')
load('cn_cor')
cn_cor <- as.character(cn_cor@results$Variables)
length(cn_cor)
plot(niches[[cn_cor[1:12]]])
plot(niches[[cn_cor[13:24]]])
plot(niches[[cn_cor[25:36]]])
plot(niches[[cn_cor[37:48]]])
plot(niches[[cn_cor[49:length(cn_cor)]]])
freq(niches[[cn_cor[46]]])
#### removing too narrow niches ####
# second iteration of filtering out extreme niches (some had suitability bands of only a few pixels still)
nl <- NA
for (i in cn_cor) {
  if(freq(niches[[i]] > 0.7)[2,2] > 1000 & freq(niches[[i]] < 0.3)[2,2] > 1000) {
    nl <- c(nl, i)
  } 
}
nl <- NA
for (i in cn_cor) {
  if(freq(niches[[i]] < 0.001)[2,2] < 15000 & freq(niches[[i]] > 0.999)[2,2] < 2000) {
    nl <- c(nl, i)
  } 
} 
# empty space less than 35% prevalence for P65 and full 1 suit less than 5% prevalence for P05
# this filterings were based on the initial set of prevalence values to be tested
# this ensures that realistically narrow niches (0.05 prevalence) were still considered

nl <- nl[2:length(nl)]
save(nl, file = 'nl')
load('nl')
#### selected climatic niches #####
env.niches <- sample(nl, 100, replace = F) # randomly selected 100 environmental niches
setwd(L.dir)
setwd('EnvNiche')
save(env.niches, file = 'env.niches')
load('env.niches')

for (i in 1:length(env.niches)) {
  r <- niches[[env.niches[i]]]
  writeRaster(r, filename = paste('N', i, '.tif', sep = ''), overwrite = T)
} # Final name of niches

### species distributions were left deconstructed because creating and saving
### all species distributions would take up too much space
### 100 species 7 prevalence 101 years binary/probabilistic = 141,400 rasters