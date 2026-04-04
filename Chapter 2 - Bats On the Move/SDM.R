####R packages needed for analyses
library(randomForest)
library(raster)
library(rgeos)
library(maptools)
library(dismo)
library(kernlab)
library(sp)
library(ecospat)
library(ENMeval)
library(maps)
library(pROC)
library(sdm)

####Change to folder where all occurrences are
setwd("")

####Reads in all tiff files (environmental variables)
files1 <- list.files(path="./",pattern = 'sum.tif$', full.names = TRUE)
pred=stack(files1)
ext=extent(pred)
r=raster(ext=ext, res=c(1,1))

setwd("")
files2 <- list.files(path="./",pattern = '.csv$', full.names = TRUE)


##########Creating data matrix to collect data
a1=matrix(nrow=length(files2),ncol=16)
a2=matrix(nrow=length(files2),ncol=36)
colnames(a1)=c('Season','Species','AUC_rf','AUC_gl','AUC_mx','AUC_bc','thr_AUC','AUC_AUCw_en','AUC_tssw_en','tss_rf','tss_gl','tss_mx','tss_bc','thr_tss','tss_AUCw_en','tss_tssw_en')
colnames(a2)=c("Season","Species","rf_forest","rf_karst","rf_Mines","rf_landcover","rf_ndvi","rf_elev","rf_prec","rf_srad","rf_tavg","rf_vapr","rf_wind","gl_intercept","gl_forest","gl_karst","gl_Mines","gl_landcover","gl_ndvi","gl_elev","gl_prec","gl_srad","gl_tavg","gl_vapr","gl_wind","mx_forest","mx_karst","mx_Mines","mx_landcover","mx_ndvi","mx_elev","mx_prec","mx_srad","mx_tavg","mx_vapr","mx_wind")


for(i in 331:length(files2)){

  ###Sets the java path so rJava runs properly (only run this if rJava doesn't run without)
  Sys.setenv(JAVA_HOME='C:/Program Files/Java/jre1.8.0_431')
  library(rJava)


  ####Reads in location data, and extracts values from each environmental layer
  loc=read.csv(files2[i])
  loc=loc[,-1]
  pres=raster::extract(pred, loc)
  set.seed(0)
  backgr=randomPoints(pred, 500)
    
  ###psuedo absence points
  absvals=raster::extract(pred, backgr)
  pb=c(rep(1, nrow(pres)), rep(0, nrow(absvals)))
  sdmdata=data.frame(cbind(pb, rbind(pres, absvals)))
    
  group <- kfold(loc, 5)
  pres_train=loc[group != 1, ]
  pres_test=loc[group == 1, ]
  names(pred)=c("forest","karst","global_mines","landcover","ndvi","elev","prec","srad","tavg","vapr","wind")
  
  #################Starting actual predictions	
  ###BioClim Model
  bc=bioclim(pred, pres_train)
  backg=randomPoints(pred, n=1000, ext=ext, extf = 1.25)
  colnames(backg) = c('Longitude', 'Latitude')
  group=kfold(backg, 5)
  backg_train=backg[group != 1, ]
  backg_test=backg[group == 1, ]
    
  ####Actual prediction from BC
  pb=predict(pred, bc, ext=ext, progress='')
    
    
  ####maxent
  ###Need Java installed
  ###Need MaxENT in the dismo folder
  ###MaxENT will need to be downloaded separately and name where it is saved
  jar=paste(system.file(package="dismo"), "./java/maxent.jar", sep='')
  xm=dismo::maxent(pred, pres_train)
  px=dismo::predict(pred, xm, ext=ext, progress='')
  
   
  ####RAndom Forest
  #colnames(backg_train)=c("lon","lat")
  names(pred)=c("forest","karst","global_mines","landcover","ndvi","elev","prec","srad","tavg","vapr","wind")
  train=rbind(pres_train, backg_train)
  pb_train=c(rep(1, nrow(pres_train)), rep(0, nrow(backg_train)))
  envtrain=raster::extract(pred, train)
  envtrain=data.frame( cbind(pa=pb_train, envtrain) )
  envtrain=na.omit(envtrain)
  model=factor(pa)~forest+karst+global_mines+landcover+ndvi+elev+prec+srad+tavg+vapr+wind
  rf1=randomForest(model, data=envtrain)
  pr=dismo::predict(pred, rf1, ext=ext)
    
  ####GLM
  glm1=glm(pa~forest+karst+global_mines+landcover+srad+tavg+vapr+wind+elev+ndvi+prec, data = envtrain)
  pg=dismo::predict(pred,glm1,ext=ext)
  
  
  ###Weighted Mean (for ensemble model [AUC])
  emx=evaluate(pres_test, backg_test, xm, pred)

  test=rbind(pres_test,backg_test)
  test2=raster::extract(pr,test)
  test2[is.na(test2)]=0
  test_k=c(rep(1, nrow(pres_test)), rep(0, nrow(backg_test)))
  erf=pROC::auc(test_k,test2)
  
  egl=evaluate(pres_test, backg_test, glm1, pred)
  ebc=evaluate(pres_test, backg_test, bc, pred)
    
  auc <- sapply(list(egl, emx, ebc), function(x) x@auc)
  auc=c(erf,auc)
  w <- (auc-0.5)^2
  models=stack(pr,pg,px,pb)
  m1 <- weighted.mean(models, w)
    
  ###Weighted Mean (for ensemble model [TSS])
  colnames(backgr)=c("Longitude","Latitude")
  tss=rbind(loc,backgr)
  tss_pts=SpatialPoints(cbind(tss[,1],tss[,2]))
  tss2=cbind(tss,sdmdata$pb)
    
  tss.mx_values=raster::extract(px,tss_pts,method='simple')
  tss.glm_values=raster::extract(pg,tss_pts,method='simple')
  tss.rf_values=raster::extract(pr,tss_pts,method='simple')
  tss.bc_values=raster::extract(pb,tss_pts,method='simple')
    
  tss2=cbind(tss,sdmdata$pb,tss.mx_values,tss.bc_values,tss.glm_values,tss.rf_values)
  tss_final=na.omit(tss2)
    
  tss.mx=ecospat.max.tss(tss_final[,4],tss_final[,3])
  tss.bc=ecospat.max.tss(tss_final[,5],tss_final[,3])
  tss.glm=ecospat.max.tss(tss_final[,6],tss_final[,3])
  tss.rf=ecospat.max.tss(tss_final[,7],tss_final[,3])
    
  tss_values=c(tss.rf$max.TSS,tss.glm$max.TSS,tss.mx$max.TSS,tss.bc$max.TSS)
  tss_values2=as.numeric(tss_values[1:4])
  models=stack(pr,pg,px,pb)
  m2 <- weighted.mean(models, tss_values2)
    
  #####Evaluating ensemble models
  tss.en_auc_values=raster::extract(m1,tss_pts,method='simple')
  tss.en_tss_values=raster::extract(m2,tss_pts,method='simple')
  tss.en_p1=cbind(tss,sdmdata$pb,tss.en_auc_values,tss.en_tss_values)
  tss.en_p2=na.omit(tss.en_p1)
  tss.en.auc=ecospat.max.tss(tss.en_p2[,4],tss.en_p2[,3])
  tss.en.tss=ecospat.max.tss(tss.en_p2[,5],tss.en_p2[,3])
    
    
  tss_en_no_neg=tss.en_p2
  tss_en_no_neg[,3]=ifelse(tss_en_no_neg[,3]<0,0,tss_en_no_neg[,3])
  tss_en_no_neg[,4]=ifelse(tss_en_no_neg[,4]<0,0,tss_en_no_neg[,4])
  tss_en_no_neg[,5]=ifelse(tss_en_no_neg[,5]<0,0,tss_en_no_neg[,5])
    
  auc.en.auc=pROC::auc(tss_en_no_neg[,3],tss_en_no_neg[,4])
  auc.en.tss=pROC::auc(tss_en_no_neg[,3],tss_en_no_neg[,5])
    
    
  ###Generate threshold
  tr1=ecospat.mpa(m1,pres_test,perc = 0.9)
  tr2=ecospat.mpa(m2,pres_test,perc = 0.9)
    
  ###Generate p/np raster
  m1_tr=m1 > tr1
  m2_tr=m2 > tr2
    
  species=files2[i]
  pattern <- "./\\s*(.*?)\\s*summer"
  result <- regmatches(species, regexec(pattern, species))
  species=result[[1]][2]
  ###Add data to matrix
  ###Tells us which models performed best
  ###We are looking for largest values in each one starting with "AUC" [3:9] and "TSS" [10:16]
  a1[i,1]="summer"
  a1[i,2]=species
  a1[i,3:6]=w
  a1[i,7]=tr1
  a1[i,8]=auc.en.auc-0.5
  a1[i,9]=auc.en.tss-0.5
  a1[i,10:13]=tss_values2
  a1[i,14]=tr2
  a1[i,15]=tss.en.auc$max.TSS
  a1[i,16]=tss.en.tss$max.TSS
    

  mx_importance=ecospat.maxentvarimport (model=xm, dfvar=envtrain[,-1], nperm=5)
    
  ###variable importance matrix
  a2[i,1]="summer"
  a2[i,2]=species
  a2[i,3:13]=importance(rf1)
  a2[i,14:25]=glm1$coefficients
  a2[i,26:36]=mx_importance
  
  print(i)
}

