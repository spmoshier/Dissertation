#load in iucn shapefiles - MAMMALS_TERRESTRIAL_ONLY
setwd("")
#load in libraries for the entire set of operations
library(sf) #needed for spatial union, geometry calculations
library(dplyr)
library(lwgeom)
library(sp) #needed for range characteristic calculations - bbox function

#filter iucn dataset to only our species
species<-MO_filter$Latin_Name
write.csv(species, "./species.csv")
setdiff(species, bats$sci_name)

#update taxonomy in the base file and run some checks to ensure that our species match with the iucn species
ourspecies<-read.csv("./species_iucnmatch.csv")
mammals<-st_read("MAMMALS_TERRESTRIAL_ONLY.shp") #reads in the iucn multipolygon shapefile
bats<-mammals[mammals$order_=="CHIROPTERA",]
#figure out who is in our final species list that doesn't match the iucn list
setdiff(ourspecies$new, bats$sci_name)

ourbats<-mammals %>%
  filter(sci_name %in% ourspecies$new)
#check - it returned only the expected differences so it's good to go
length(unique(ourbats$sci_name))
names<-unique(ourbats$sci_name)
setdiff(names, ourspecies$new)
setdiff(ourspecies$new, names)

#merge the polygons by species
sf_use_s2(FALSE) #keeps it from throwing an error for dissolving some of the large polygons
merged_bats <- ourbats |>
  group_by(sci_name) |>
  summarise(
    # You can also summarize other attributes if needed (e.g., total area, count)
    total_area = sum(SHAPE_Area, na.rm = TRUE),
    total_leng = sum(SHAPE_Leng, na.rm = TRUE),
    n_polygons = n(),
 )

#now summarize a bunch of data for the polygons to attribute to the bats later on
df<-data.frame()
df[length(merged_bats$sci_name),1:5]=NA
colnames(df)=c("sci_name","total_area_m2","lat_extent_deg","mid_lat", "mid_lon")
df$sci_name<-merged_bats$sci_name

#Loop through rasters to calculate niche overlap/centroids/centroid distance & fill into dataset
for (i in 1:length(merged_bats$sci_name)){
  #Calculates polygon area, latitudinal extent & mid lat coordinate
  area<-st_area(merged_bats[i,1])
  df$total_area_m2[[i]]<-area
  
  bbox<-st_bbox(merged_bats[i,1])
  min_lat<-bbox["ymin"]
  max_lat<-bbox["ymax"]
  latext<-max_lat-min_lat
  df$lat_extent_deg[[i]]<-latext
  
  cent<-st_centroid(merged_bats[i,1])
  coord<-st_coordinates(cent)
  df$mid_lon[[i]]<-coord[1,"X"]
  df$mid_lat[[i]]<-coord[1,"Y"]
}

#merge back with our species so we can match it to the niche overlap dataframe with all the other stats
colnames(ourspecies)[3]<-"sci_name"
colnames(ourspecies)[2]<-"Latin_Name"
bat_range_stats<-merge(df,ourspecies, by="sci_name",all=T) #full join of the two datasets
write.csv(bat_range_stats,"")
