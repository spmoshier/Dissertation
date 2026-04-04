#set working directory
setwd("")

#load in libraries
library(ggplot2) #plotting
library(terra) #load in tif raster files, xyCellFRom & ext functions
library(ENMTools) #raster.overlap function
library(geosphere) #distGeo function

#load in SDMs (TIF files)
sum_files<-list.files(path="./Summer_SDMs/",pattern = '.tif$', full.names = T)
summer<-lapply(sum_files, raster) #builds a raster stack object
win_files<-list.files(path="./Winter_SDMs/",pattern = '.tif$', full.names = T)
winter<-lapply(win_files, raster) #builds a raster stack object
species<-tools::file_path_sans_ext(basename(win_files))

#create new data frame for data
df<-data.frame()
df[length(summer),1:6]=NA
colnames(df)=c("species","summer_tif","winter_tif","D","I","rank.corr")
df[,1]<-species
df[,2]<-sum_files
df[,3]<-win_files

##redo extents for sdms that did not clip properly
#set working directory
setwd("")

#load in SDMs (TIF files) - make sure to use terra rather than raster package
sum_files<-list.files(path="./Summer_SDMs/",pattern = '.tif$', full.names = T)
summer<-lapply(sum_files, terra::rast) #builds a raster stack object
win_files<-list.files(path="./Winter_SDMs/",pattern = '.tif$', full.names = T)
winter<-lapply(win_files, terra::rast) #builds a raster stack object
species<-tools::file_path_sans_ext(basename(win_files))

#create new data frame for data
df<-data.frame()
df[length(summer),1:6]=NA
colnames(df)=c("species","summer_tif","winter_tif","D","I","rank.corr")
df[,1]<-species
df[,2]<-sum_files
df[,3]<-win_files

#new extent for Artibeus literatus
ext(summer[[1]])
plot(summer[[1]])
summer[[1]]<-crop(summer[[1]], ext(-112.5, -25, -32.9166666666667, 30.5833333333333))
winter[[1]]<-crop(winter[[1]], ext(-112.5, -25, -32.9166666666667, 30.5833333333333))
plot(summer[[1]])
#new extent for Eptesicus serotinus
ext(summer[[2]])
plot(summer[[2]])
summer[[2]]<-crop(summer[[2]], ext(-123.166666666667,-40, 16.9166666666667, 64.75))
winter[[2]]<-crop(winter[[2]], ext(-123.166666666667,-40, 16.9166666666667, 64.75))
plot(summer[[2]])
#new extent for Eumops auripendulus
ext(summer[[3]])
plot(summer[[3]])
summer[[3]]<-crop(summer[[3]], ext(-97.9166666666667, -30, -33, 26))
winter[[3]]<-crop(winter[[3]], ext(-97.9166666666667, -30, -33, 26))
plot(summer[[3]])
#new extent for Lasiurus borealis
ext(summer[[4]])
plot(summer[[4]])
summer[[4]]<-crop(summer[[4]], ext(-127.666666666667, -25, -47, 57.8333333333333))
winter[[4]]<-crop(winter[[4]], ext(-127.666666666667, -25, -47, 57.8333333333333))
plot(summer[[4]])
#new extent for Myotis lucifugus
ext(summer[[5]])
plot(summer[[5]])
summer[[5]]<-crop(summer[[5]], ext(-160.75, -45, 14.4166666666667, 72.4166666666667))
winter[[5]]<-crop(winter[[5]], ext(-160.75, -45, 14.4166666666667, 72.4166666666667))
plot(summer[[5]])
#new extent for Myotis myotis
ext(summer[[6]])
plot(summer[[6]])
summer[[6]]<-crop(summer[[6]], ext(-20, 48.4166666666667, 10.5, 63.5))
winter[[6]]<-crop(winter[[6]], ext(-20, 48.4166666666667, 10.5, 63.5))
plot(summer[[6]])
#new extent for Myotis mystacinus
ext(summer[[7]])
plot(summer[[7]])
summer[[7]]<-crop(summer[[7]], ext(-30, 124.083333333333, -22.0833333333333, 68.9166666666667))
winter[[7]]<-crop(winter[[7]], ext(-30, 124.083333333333, -22.0833333333333, 68.9166666666667))
plot(summer[[7]])
#new extent for Nyctalus leisleri
ext(summer[[8]])
plot(summer[[8]])
summer[[8]]<-crop(summer[[8]], ext(-20, 76.1666666666667, 29.0833333333333, 71.5))
winter[[8]]<-crop(winter[[8]], ext(-20, 76.1666666666667, 29.0833333333333, 71.5))
plot(winter[[8]])
#new extent for Rhinolophus ferrumequinum
ext(summer[[9]])
plot(summer[[9]])
summer[[9]]<-crop(summer[[9]], ext(-20, 146.666666666667, 21.25, 58.25))
winter[[9]]<-crop(winter[[9]], ext(-20, 146.666666666667, 21.25, 58.25))
plot(summer[[9]])
#new extent for Rousettus aegyptiacus
ext(summer[[10]])
plot(summer[[10]])
plot(winter[[10]])
summer[[10]]<-crop(summer[[10]], ext(-20, 69.0833333333333, -38.9166666666667, 45))
winter[[10]]<-crop(winter[[10]], ext(-20, 69.0833333333333, -38.9166666666667, 45))
plot(winter[[10]])
#new extent for Tadarida brasiliensis
ext(summer[[11]])
plot(summer[[11]])
summer[[11]]<-crop(summer[[11]], ext(-128.166666666667, -25, -47, 47.8333333333333))
winter[[11]]<-crop(winter[[11]], ext(-128.166666666667, -25, -47, 47.8333333333333))
plot(summer[[11]])

# Export rasters as GeoTIFF files (recommended for broad compatibility)
setwd("C:/Users/shelb/OneDrive/Desktop/Chapter 1 Stuff/sdm_redos/Summer_SDMs_redone")
for (i in 1:length(species)){
  writeRaster(summer[[i]], filename = paste0(species[[i]],"_summer.tif"), overwrite=T)
}

setwd("C:/Users/shelb/OneDrive/Desktop/Chapter 1 Stuff/sdm_redos/Winter_SDMs_redone")
for (i in 1:length(species)){
  writeRaster(winter[[i]], filename = paste0(species[[i]],".tif"), overwrite=T)
}

#redo centroid and niche overlap calculations
setwd("")

#Load in rasters!
sum_files<-list.files(path="./Summer_SDMs_redone/",pattern = '.tif$', full.names = T)
summer<-lapply(sum_files, terra::rast) #builds a raster stack object
win_files<-list.files(path="./Winter_SDMs_redone/",pattern = '.tif$', full.names = T)
winter<-lapply(win_files, terra::rast) #builds a raster stack object
species<-tools::file_path_sans_ext(basename(win_files))

#Create dataframe to fill
df<-data.frame()
df[length(species),1:12]=NA
colnames(df)=c("Latin.Name","summer_tif","winter_tif","D","I","rank.corr","sum_lon","sum_lat","win_lon","win_lat","distance","bearing")
df[,1]<-species
df[,2]<-sum_files
df[,3]<-win_files

#Loop through rasters to calculate niche overlap/centroids/centroid distance & fill into dataset
for (i in 1:length(species)){
  #Calculates niche overlap (3 metrics) & writes them to df
  over<-ENMTools::raster.overlap(summer[[i]],winter[[i]])
  df$D[[i]]<-over$D
  df$I[[i]]<-over$I
  df$rank.corr[[i]]<-over$rank.cor
  
  #Calculates centroids of the presence area & writes the coordinates to df
  sum_cen=colMeans(terra::xyFromCell(summer[[i]], which(terra::values(summer[[i]])==1)))
  win_cen=colMeans(terra::xyFromCell(winter[[i]], which(terra::values(winter[[i]])==1)))
  
  df$sum_lon[[i]]<-sum_cen[[1]]
  df$sum_lat[[i]]<-sum_cen[[2]]
  df$win_lon[[i]]<-win_cen[[1]]
  df$win_lat[[i]]<-win_cen[[2]]
  
  #Calculates distance between centroids and direction & writes them to df
  dis=distGeo(win_cen,sum_cen)
  dir=bearing(win_cen,sum_cen)
  df$distance[[i]]<-dis
  df$bearing[[i]]<-dir
  
  #Calculates latitudinal extent
  sum_lat=colMeans(terra::ext(summer[[i]], which(terra::values(summer[[i]])==1)))
  win_lat=colMeans(terra::ext(winter[[i]], which(terra::values(winter[[i]])==1)))
}
