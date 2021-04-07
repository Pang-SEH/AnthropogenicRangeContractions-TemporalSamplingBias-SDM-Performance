library(rgdal)
library(raster)
library(plyr)
### Data preparation v5.0 ###

##### Climatic Variables #####
base.1 <- 'D:/OneDrive - National University of Singapore/0 TempBias'

main <- base.1
setwd(main)
setwd('Env')
shp <- readOGR('SEA_cropped_rms.shp') # shape file of Southeast Asia
setwd('SEA (clim)')
b.list <- list.files(pattern = 'bio')
g.env <- getwd()
s.env <- paste(main, '/Env/SEA (clim)', sep = '')

for (b in b.list){ # cropping bio to shape file
  setwd(g.env)
  bio <- mask(crop(raster(b), shp) , shp)
  setwd(s.env)
  writeRaster(bio, b)
}

ref_lay <- raster('bio01.tif')

soil.dir <- 'D:/OneDrive - National University of Singapore/Soil Predictors'
setwd(soil.dir)
b.list <- list.files(pattern = '.tif') # see Table S1.2 for list of soil predictors

for (b in b.list){ # Cropping soil variables, resampling and recrop final
  setwd(soil.dir)
  bio <- mask(crop(raster(b), shp) , shp)
  bio <- resample(bio, ref_lay, 'bilinear')
  bio <- mask(crop(bio, ref_lay), ref_lay)
  
  setwd(s.env)
  writeRaster(bio, b)
}
ref_lay <- bio

setwd(main)
setwd('Env')
setwd('SEA (clim)')
b.list <- list.files(pattern = 'bio')
g.env <- getwd()

for (b in b.list){ # re-cropping bio to new raster extent
  setwd(g.env)
  bio <- mask(crop(raster(b), ref_lay) , ref_lay)
  setwd(s.env)
  writeRaster(bio, b, overwrite = T)
}

## PCA for cimate and soil ##
setwd(s.env)
env <- stack(list.files()) # stack of climate and soil variables
env <- as.data.frame(env, na.rm = T)

setwd(main)
setwd('Env')
setwd('pca')
pca <- prcomp(env, scale=T) # important scale = T
save(pca, file = "PCA analysis") # PCA saved
pca.sum <- summary(pca)

write.csv(pca.sum$rotation, file = "PCA_loadings.csv") # contributions towards each PC axis
write.csv(pca.sum$importance, file = "PCA_Variances.csv") # Cumulative variance accounted for

## prediction and PC axis ##
num.axis <- as.numeric(sum(pca.sum$importance[3,] <= 0.85)) + 1 # conservative 85%
setwd(s.env)
env <- stack(list.files())
PCs <- predict(env, pca, index=1:num.axis) # new PC maps for all
names(PCs) <- paste('PC', 1:num.axis, sep = "") # simple rename

setwd(main)
setwd('Env')
setwd('pca')
writeRaster(PCs, filename = 'PCs_all.grd') # .grd saves name

#### land use variable ####
base <- 'D:/OneDrive - National University of Singapore/'
main <- paste(base, '0 TempBias', sep = '')

setwd(main)
setwd('Env')
shp <- readOGR('SEA_cropped_rms.shp') # shape file of Southeast Asia west of Wallace's line

setwd('pca')
env <- list.files(pattern = 'grd')
env <- stack(env)
ref_lay <- env$PC1
##### Land Use variables #####
library(ncdf4) # nc file
d.land <- paste(base, 'LUH2_v2h_Historic', sep = '')
setwd(d.land)
st.list <- c("primf", "primn" ,"secdf", "secdn", "urban" ,"c3ann", "c4ann", "c3per", "c4per" ,"c3nfx" ,"pastr" ,"range")
luh.urban <- brick('states.nc', varname='urban')
years <- c(seq(1500, 1880, 20), 1900:2000)
yearsC <- years - 849
library(foreach)
library(doParallel)
cl <- makeCluster(12)
registerDoParallel(cl)
## resampling of land use raster first ##
foreach (st = st.list, .packages = c('raster', 'ncdf4')) %dopar% {
  setwd(d.land)
  luh <- brick('states.nc', varname=st)
  luh <- luh[[yearsC]]
  luh <- resample(luh, ref_lay, 'bilinear')
  luh <- mask(crop(luh, ref_lay), ref_lay)
  names(luh) <- paste(st, years, sep = '_')
  setwd(paste(main, 'Env', 'SEA (LUH)', sep = '/'))
  writeRaster(luh, filename = paste(st, '.grd', sep = ""))
}
stopCluster(cl)
closeAllConnections()

# summing up the anthropogenic land use classes #
luh.dir <- paste(main, 'Env', 'SEA (LUH)', sep = '/')
setwd(luh.dir)

anthro <- stack('urban.grd') + stack('c3ann.grd') + stack('c4ann.grd') + stack('c3per.grd') +
  stack ('c4per.grd') + stack('c3nfx.grd') + stack('pastr.grd') + stack('range.grd')
names(anthro) <- paste('anthr', years, sep = '_')
writeRaster(anthro, 'anthr.grd')

## creating individual maps for easier access ##
op <- paste(main, 'Env', 'LUH_sim', sep = '/')
setwd(op)
years <- c(seq(1500, 1880, 20), 1900:2000)
for (i in years) {
  r.cont <- anthro[[paste('anthr', i, sep = '_')]]
  setwd(op)
  setwd('cont')
  writeRaster(r.cont, paste('LUH_cont_', i, '.tif', sep = ''))
}


