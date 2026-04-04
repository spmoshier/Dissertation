library(raster)
library(geosphere)

###Creates list of rasters (with and without full names)
folders=list.dirs(path='.')

####Data frame to write calculations into
dat=data.frame()
dat[1,1:7]=NA
colnames(dat)=c("species","sum_lon","sum_lat","win_lon","win_lat","distance","bearing")

for(j in 2:length(folders)){
  
  
  ####Reads in locations to determine cropped extent
  files=list.files(path=folders[j], pattern='thin1.csv$', full.names=TRUE)
  
  dat1=read.csv(files[1])
  dat2=read.csv(files[2])
  
  
  ####Creates lists of max and min of each direction
  lat=c(max(dat1$decimalLatitude),max(dat2$decimalLatitude),min(dat1$decimalLatitude),min(dat2$decimalLatitude))
  lon=c(max(dat1$decimalLongitude),max(dat2$decimalLongitude),min(dat1$decimalLongitude),min(dat2$decimalLongitude))
  
  ####Creates extent for cropping (has 5deg buffer of sampled locations)
  e=extent(min(lon)-5,max(lon)+5,min(lat)-5,max(lat)+5)
  
  
  ####Reads in suitable/not suitable rasters for centroid calculation
  files2=list.files(path=folders[j],pattern='np.tif$',full.names=TRUE)
  
  r_sum1=raster(files2[1])
  r_win1=raster(files2[2])
  
  ####Crops suitable/ns rasters to range of samples with 5deg buffer
  r_sum2=crop(r_sum1, e)
  r_win2=crop(r_win1, e)
  
  ####Calculates centroids
  sum_cen=colMeans(xyFromCell(r_sum2, which(r_sum2[]==1)))
  win_cen=colMeans(xyFromCell(r_win2, which(r_win2[]==1)))
  
  ####Calculates distance between centroids and direction
  dis=distGeo(win_cen,sum_cen)
  dir=bearing(win_cen,sum_cen)
  
  ####Writes data to data frame
  dat[j,1]=sub('./', '', folders[j])
  dat[j,2:3]=sum_cen
  dat[j,4:5]=win_cen
  dat[j,6]=dis
  dat[j,7]=dir
  
  
  ####Reads in 0-1 rasters to put in the correct folder setup
  files3=list.files(path=folders[j],pattern='r.tif$',full.names=TRUE)
  r_sum3=raster(files2[1])
  r_win3=raster(files2[2])
  
  ####Crops rasters to range of samples with 5deg buffer
  r_sum4=crop(r_sum1, e)
  r_win4=crop(r_win1, e)
  
  
  ####Writes out initial rasters
  fn1_sum=paste("./z_Input1/Summer_SDMs",sub('.', '', folders[j]),'_summer.tif',sep='')
  fn1_win=paste("./z_Input1/Winter_SDMs",sub('.', '', folders[j]),'.tif',sep='')
  
  writeRaster(r_sum3,fn1_sum)
  writeRaster(r_win3,fn1_win)
  
  ####Writes out cropped rasters
  fn2_sum=paste("./z_Input2/Summer_SDMs/",sub('.', '', folders[j]),'_summer.tif',sep='')
  fn2_win=paste("./z_Input2/Winter_SDMs/",sub('.', '', folders[j]),'.tif',sep='')
  
  writeRaster(r_sum4,fn2_sum)
  writeRaster(r_win4,fn2_win)
  
}
