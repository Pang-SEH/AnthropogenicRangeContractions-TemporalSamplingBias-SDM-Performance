library(SDMtune)
library(ENMeval)
main <- 'D:/OneDrive - National University of Singapore/0 TempBias/'
main.l <- 'D:/S.Pang/0 TempBias/'

#### environmental predictors ####
pca.dir <- paste(main, 'Env/pca', sep = '')
setwd(pca.dir)
PCs <- stack('PCs_all.grd')
bg.data <- read.csv('bg.data_5.csv')

samp.dir <- paste(main.l, 'Samples (focc)', sep = '')
opb <- paste(main.l, 'Models (focc)', sep = '')

setwd(samp.dir)
sp.list <- list.files()
sp.list <- sub('.csv', '', sp.list)

#### modeling for loop ####
library(foreach)
library(doParallel)
cl <- makeCluster(detectCores() - 2)
registerDoParallel(cl)

foreach (sp = sp.list, .packages = c('SDMtune', 'ENMeval', 'rJava')) %dopar% {
  for (L in c('cont')) {
    ## setting land cover format ##
    landcat <- paste(main, 'Env/LUH_sim', sep = '')
    setwd(landcat)
    setwd(L) # land cover format to be used
    landcat <- raster(paste('LUH_', L, '_2000.tif', sep = ''))
    names(landcat) <- 'luh'
    hlandcat <- raster(paste('LUH_', L, '_1900.tif', sep = ''))
    names(hlandcat) <- 'luh'
    
    cate.cont <- NULL # assumes non-categorical first
    if (L == 'cate') {
      landcat <- as.factor(landcat)
      hlandcat <- as.factor(hlandcat)
      cate.cont <- 'luh'
    } # is L is categorical, we set the luh (land cover) variables as a categorical one
    
    setwd(opb)
    if (!dir.exists(L)) {dir.create(L)}
    setwd(L) # outputs set for each type of landcover format
    op <- getwd()
    
    #### Starting Extraction and Modelling for each sample ####
    ## data extraction ##
    # st <- proc.time()
    setwd(samp.dir)
    # sampling occurrence points
    samp <- read.csv(paste(sp, '.csv', sep = ''))[, c('Lon', 'Lat')]
    # true land cover values for each occurrence point (land cover format specific)
    # samp.land <- read.csv(paste(sp, '.csv', sep = ''))[, paste('S', L, sep = '')]
    # creating species dir #
    setwd(op)
    if (!dir.exists(sp)) {dir.create(sp)}
    setwd(sp)
    
    #### modelling with land cover predictor (with Land) ####
    Proj.nm <- 'wLand'
    setwd(op); setwd(sp)
    if (!dir.exists(Proj.nm)) {dir.create(Proj.nm)}
    setwd(Proj.nm)
    # # #
    swd <- NA; enm <- NA
    # extracting data for species
    # categorical follows cate.cont (null or luh)
    swd <- prepareSWD(species = sp, p = samp, a = bg.data,
                      env = stack(PCs, landcat), categorical = cate.cont)
    enm <- train(method = 'Maxnet', data = swd)
    
    s <- as.data.frame(matrix(NA, 1, 2))
    colnames(s) <- c('mAUC', 'mTSS') # model fit scores
    s[1, 'mAUC'] <- SDMtune::auc(enm)
    s[1, 'mTSS'] <- SDMtune::tss(enm)
    
    confusion.matrix <- SDMtune::confMatrix(enm, th = thresholds(enm, type = 'cloglog')[,2], type = 'cloglog')
    rownames(confusion.matrix) <- thresholds(enm, type = 'cloglog')[,1]
    
    # training model scores
    write.csv(s, file = 'train_score.csv', row.names = F)
    # Variable Importance (reliance on LandCover)
    write.csv(varImp(enm, 10), file = 'VarImp.csv', row.names = F)
    # Different Threshold (mainly Max TSS) for binary conversion
    write.csv(thresholds(enm, type = 'cloglog'), file = 'Thresholds.csv', row.names = F)
    # Confusion matrix for the different thresholds
    write.csv(confusion.matrix, file = 'confusion_matrix.csv')
    # Prediction (Outputs using contemporary maps 2000)
    writeRaster(predict(enm, stack(PCs, landcat), type = 'cloglog', clamp = F),
                filename = 'Pred_2000_cloglog.tif', overwrite = T)
    # Projection (Outputs using historical maps 1900)
    writeRaster(predict(enm, stack(PCs, hlandcat), type = 'cloglog', clamp = F),
                filename = 'Pred_1900_cloglog.tif', overwrite = T)
    
    ### model with environments only ##
    Proj.nm <- 'EnvO'
    setwd(op); setwd(sp)
    if (!dir.exists(Proj.nm)) {dir.create(Proj.nm)}
    setwd(Proj.nm)
    # # #
    swd <- NA; enm <- NA
    # full reset, data prep without ANY landcover
    swd <- prepareSWD(species = sp, p = samp, a = bg.data,
                      env = PCs)
    enm <- train(method = 'Maxnet', data = swd)

    s <- as.data.frame(matrix(NA, 1, 2))
    colnames(s) <- c('mAUC', 'mTSS')
    s[1, 'mAUC'] <- SDMtune::auc(enm)
    s[1, 'mTSS'] <- SDMtune::tss(enm)
    
    confusion.matrix <- NA
    confusion.matrix <- SDMtune::confMatrix(enm, th = thresholds(enm, type = 'cloglog')[,2], type = 'cloglog')
    rownames(confusion.matrix) <- thresholds(enm, type = 'cloglog')[,1]
    
    write.csv(s, file = 'train_score.csv', row.names = F)
    write.csv(varImp(enm, 10), file = 'VarImp.csv', row.names = F)
    write.csv(confusion.matrix, file = 'confusion_matrix.csv')
    write.csv(thresholds(enm, type = 'cloglog'), file = 'Thresholds.csv', row.names = F)
    # again, only single projection as Environmental Variables were assumed unchanged
    writeRaster(predict(enm, PCs, type = 'cloglog', clamp = F),
                filename = 'Pred_cloglog.tif', overwrite = T)
    
    setwd(op)
  }
}
stopCluster(cl)
closeAllConnections()

### note: there are some issues with running on maxnet
# a convergence error consistently pops up that prevents the development of
# some models.. An example of the error is below

# Error in { : task 8 failed - "index larger than maximal 168"

# The developers are still trying to resolve this error, but further testing concludes
# that there is no effect on the final model
# No distinct difference was noted between replicates that had this error and those that did not
# when  maxent.jar was used instead.
# In summary, some models might fail, but excluding them out would have no effect on the results
