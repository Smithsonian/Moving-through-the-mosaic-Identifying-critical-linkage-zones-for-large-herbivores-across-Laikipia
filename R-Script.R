############################################################################################################################
# Moving through the mosaic: Identifying critical linkage zones for large herbivores across a multiple-use African landscape
# Ramiro D. Crego, Harry. B.M. Wells, Kimani S. Ndung’u, Lauren Evans, Redempta Njeri Nduguta, Muthiuru A. Chege, 
# Michael B. Brown, Joseph O. Ogutu, Gordon O. Ojwang, Julian Fennessy, David O'Connor, Jenna Stacy-Dawes, Daniel I. Rubenstein,
# Dino J. Martins, Peter Leimgruber, and Jared A. Stabach

# This code generates the resistance surfaces and runs the Circuitscape analyses.

# Libraries
library(sp)
library(rgdal)
library(raster)
library(spatialEco)
library(maptools)
library(rgeos)
library(move)
library(dismo)
library(plyr)

############################################ 
## Define study area and coordinate system #
############################################

CRS<-("+proj=utm +zone=37 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")

# Define extent
extL<-extent(186885.9, 321651.3, -32707.34, 96041.66) #Laikipia extent

# Create a 100 m resolution raster
a <- raster(ext=extL, res=100, vals=0)

# Load Study Area
StudyArea<-readOGR(dsn = "./GIS_Layers", layer = "StArea")
StudyAreaRast <- rasterize(StudyArea, a)
plot(StudyAreaRast)

############################
### Determine core areas ###
############################

centroids<-readOGR(dsn = "./GIS_Layers", layer = "Centroids_UTM37N")
Centroids <- rasterize(x = centroids, y = StudyAreaRast, field="cons")
projection(Centroids)<-CRS
plot(Centroids, add=T)
writeRaster(Centroids,"Centroids_rast.asc",overwrite=TRUE)

# Eliminate Pyramids conservancy for elephant analysis (Elephants are excluded from this property)
Elcentroid<-centroids[-15,]
ElCentroid <- rasterize(x = Elcentroid, y = StudyAreaRast, field="cons")
projection(ElCentroid)<-CRS
plot(ElCentroid, add=T)
writeRaster(ElCentroid,"Ele_Centroids_rast.asc",overwrite=TRUE)

############
# Elephant #
############

Eleph<-brick("./GIS_Layers/Elephant.tif")
# Resample to 100x100 m resolution
Eleph_100m<-resample(Eleph,a)
projection(Eleph_100m)<-CRS
# Estimate mean occupancy
Eleph <- (Eleph_100m$Eleph.1 + Eleph_100m$Eleph.2 + Eleph_100m$Eleph.3 + 
            Eleph_100m$Eleph.4 + Eleph_100m$Eleph.5 + Eleph_100m$Eleph.6 +
            Eleph_100m$Eleph.7 + Eleph_100m$Eleph.8)/8

#Invert raster
ElephResistance<-((Eleph - 1)*-1)
ElephResistance <- ElephResistance * 100
ElephResistance <- round(ElephResistance,0)
plot(ElephResistance)

#Import anthropogenic areas derived from Jacobson raster
Antr<-readOGR(dsn = "./GIS_Layers", layer = "AnthropogenicAreas_UTM37N")
Antr<-spTransform(Antr, CRS)
AntrRaster<-rasterize(Antr, StudyAreaRast, Antr$DN, background=NA)
#Assign a resistance value 
AntrRaster<-AntrRaster*95

#Combine rasters
ElephResistance2 <- cover(AntrRaster, ElephResistance)
plot(ElephResistance2)

#Add fences
# Reclassify the fece data to assign the resistance value
#Assign resistance values
# 1 = Porcupine fence
# 2 = Tall fence
# 3 = Tall netted fence
# 4 = Cattle fence
# 5 = Cattle fence owner repeals wildlife
# 6 = unknown fence
# 7 = ditch

fencesele<-readOGR(dsn = "./GIS_Layers", layer = "FencesUTM37NEle") #Fence with fence gaps into Pyramid removed
fencesele$Class<-as.numeric(fencesele$Class)
Fencesele<-rasterize(fencesele, StudyAreaRast, fencesele$Class, background=NA)

fromToArray <- c(1,60,2,-9999,3,-9999,4,50,5,-9999,6,50,7,72)
fromToMatrix <- matrix(fromToArray, ncol=2, byrow=TRUE)

FencesEleph<-reclassify(Fencesele,fromToMatrix)

#Add fences data keeping pixels with max resistance value
ElephResistance3<- mosaic(ElephResistance2,FencesEleph, fun=max)

#Add fences with infinite value
ElephResistance3[FencesEleph == -9999] <- -9999
plot(ElephResistance3)

#Clip and save
ElephResistanceFin<-mask(ElephResistance3,StudyAreaRast)
plot(ElephResistanceFin)
writeRaster(ElephResistanceFin,"./ElephResist_rast.asc",overwrite=TRUE)

########################
### Run circuitscape  ##
########################

# Make a place holder for the cs_run.exe path
CS_exe <- 'C:/"Program Files"/Circuitscape/cs_run.exe' 

# Make an .ini file
CS_ini <- c("[circuitscape options]",            
            "data_type = raster",
            "scenario = pairwise",
            "set_focal_node_currents_to_zero = True",
            "write_cur_maps = True",
            "write_cum_cur_map_only = True",
            "log_transform_maps = False",
            paste(c("point_file =",
                    "habitat_file =",
                    "output_file ="),
                  paste(getwd(),c("Ele_Centroids_rast.asc",
                                  "ElephResist_rast.asc",
                                  "Eleph_CS.out"),
                        sep="/")))

# Write it to your working directory
writeLines(CS_ini,"myini.ini")

# Make the CS run cmd
CS_run <- paste(CS_exe, paste(getwd(),"myini.ini",sep="/")) # Make the cmd

# Run the command
system(CS_run)


###########
# Giraffe #
###########

Giraffe<-brick("./GIS_Layers/Giraffe.tif")

# Resample to 100x100 m resolution
Giraffe_100m<-resample(Giraffe,a)
projection(Giraffe_100m)<-CRS
Giraffe<- (Giraffe_100m$Giraffe.1 + Giraffe_100m$Giraffe.2 + Giraffe_100m$Giraffe.3 + 
             Giraffe_100m$Giraffe.4 + Giraffe_100m$Giraffe.5 + Giraffe_100m$Giraffe.6 +
             Giraffe_100m$Giraffe.7 + Giraffe_100m$Giraffe.8)/8

GiraffeResistance <- ((Giraffe - 1)* -1)
GiraffeResistance <- GiraffeResistance * 100
GiraffeResistance <- round(GiraffeResistance,0)
plot(GiraffeResistance)

#Add Anthropogenic areas 
GiraffeResistance2 <- cover(AntrRaster, GiraffeResistance)
plot(GiraffeResistance2)

#Add fences
# Reclassify the fece data to assign the resistance value
#Assign resistance values
# 1 = Porcupine fence
# 2 = Tall fence
# 3 = Tall netted fence
# 4 = Cattle fence
# 5 = Cattle fence owner repeals wildlife
# 6 = unknown fence
# 7 = ditch
fences<-readOGR(dsn = "./GIS_Layers", layer = "FencesUTM37N")
fences$Class<-as.numeric(fences$Class)
Fences<-rasterize(fences, StudyAreaRast, fences$Class, background=NA)

fromToArrayG <- c(1,65,2,-9999,3,-9999,4,61,5,-9999,6,65,7,82)
fromToArrayG <- matrix(fromToArrayG, ncol=2, byrow=TRUE)
fromToMatrixG

FencesGir<-reclassify(Fences,fromToMatrix2)

#Add fences data keeping pixels with max resistance value
GiraffeResistance3<- mosaic(GiraffeResistance2,FencesGir, fun=max)

#Add fences with infinite value
GiraffeResistance3[FencesGir ==  -9999] <- -9999
plot(GiraffeResistance3)

#Clip and save
GiraffeResistanceFin<-mask(GiraffeResistance3,StudyAreaRast)
plot(GiraffeResistanceFin)
writeRaster(GiraffeResistanceFin,"GirResist_rast.asc",overwrite=TRUE)

########################
### Run Circuitscape  ##
########################

# Make a place holder for the cs_run.exe path
CS_exe <- 'C:/"Program Files"/Circuitscape/cs_run.exe' # Don't forget the "Program Files" problem

# Make an .ini file
CS_ini <- c("[circuitscape options]",            
            "data_type = raster",
            "scenario = pairwise",
            "set_focal_node_currents_to_zero = True",
            "write_cur_maps = True",
            "write_cum_cur_map_only = True",
            "log_transform_maps = False",
            paste(c("point_file =",
                    "habitat_file =",
                    "output_file ="),
                  paste(getwd(),c("Centroids_rast.asc",
                                  "GirResist_rast.asc",
                                  "Giraffe_CS.out"),
                        sep="/")))

# Write it to your working directory
writeLines(CS_ini,"myini.ini")

# Make the CS run cmd
CS_run <- paste(CS_exe, paste(getwd(),"myini.ini",sep="/")) # Make the cmd

# Run the command
system(CS_run)


################
# Plains zebra #
################

PlZebra<-brick("./GIS_Layers/PlainsZebra.tif")

# Resample to 100x100 m resolution
PlainsZebra_100m<-resample(PlZebra,a)
projection(PlainsZebra_100m)<-CRS
plot(PlainsZebra_100m)
PlainsZebra<- (PlainsZebra_100m$ComZebra.1 + PlainsZebra_100m$ComZebra.2 + PlainsZebra_100m$ComZebra.3 + 
              PlainsZebra_100m$ComZebra.4 + PlainsZebra_100m$ComZebra.5 + PlainsZebra_100m$ComZebra.6 +
              PlainsZebra_100m$ComZebra.7 + PlainsZebra_100m$ComZebra.8)/8

PlainsZebraResistance <- ((PlainsZebra-1)*-1)
PlainsZebraResistance <- PlainsZebraResistance * 100
PlainsZebraResistance <- round(PlainsZebraResistance,0)
plot(PlainsZebraResistance)

#Add Anthropogenic areas 
PlainsZebraResistance2 <- cover(AntrRaster, PlainsZebraResistance)
plot(PlainsZebraResistance2)

#Add fences
# Reclassify the fece data to assign the resistance value
#Assign resistance values
# 1 = Porcupine fence
# 2 = Tall fence
# 3 = Tall netted fence
# 4 = Cattle fence
# 5 = Cattle fence owner repeals wildlife
# 6 = unknown fence
# 7 = ditch

fromToArrayPZ <- c(1,40,2,-9999,3,-9999,4,50,5,-9999,6,40,7,53)
fromToMatrixPZ <- matrix(fromToArrayPZ, ncol=2, byrow=TRUE)
fromToMatrixPZ

FencesPlainsZebra<-reclassify(Fences,fromToMatrixPZ)

#Add fences data keeping pixels with max resistance value
PlainsZebraResistance3<- mosaic(PlainsZebraResistance2,FencesPlainsZebra, fun=max)

#Add fences with infinite value
PlainsZebraResistance3[FencesPlainsZebra ==  -9999] <- -9999
plot(PlainsZebraResistance3)

PlainsZebraResistanceFin<-mask(PlainsZebraResistance3,StudyAreaRast)
plot(PlainsZebraResistanceFin)
writeRaster(PlainsZebraResistanceFin,"PlainsZebraResist_rast.asc",overwrite=TRUE)


########################
### Run Circuitscape  ##
########################

#With fences
# Make a place holder for the cs_run.exe path
CS_exe <- 'C:/"Program Files"/Circuitscape/cs_run.exe' # Don't forget the "Program Files" problem

# Make an .ini file
CS_ini <- c("[circuitscape options]",            
            "data_type = raster",
            "scenario = pairwise",
            "set_focal_node_currents_to_zero = True",
            "write_cur_maps = True",
            "write_cum_cur_map_only = True",
            "log_transform_maps = False",
            paste(c("point_file =",
                    "habitat_file =",
                    "output_file ="),
                  paste(getwd(),c("Centroids_rast.asc",
                                  "PlainsZebraResist_rast.asc",
                                  "PlainsZebra_CS.out"),
                        sep="/")))

# Write it to your working directory
writeLines(CS_ini,"myini.ini")

# Make the CS run cmd
CS_run <- paste(CS_exe, paste(getwd(),"myini.ini",sep="/")) # Make the cmd

# Run the command
system(CS_run)


################
# Grevys zebra #
################

GrevyZebra<-brick("./GIS_Layers/GrevyZebra.tif")
# Resample to 100x100 m resolution
GrevyZebra_100m<-resample(GrevyZebra,a)
projection(GrevyZebra_100m)<-CRS
plot(GrevyZebra)
GrevyZebra<- (GrevyZebra_100m$GrevyZebra.1 + GrevyZebra_100m$GrevyZebra.2 + GrevyZebra_100m$GrevyZebra.3 + 
                GrevyZebra_100m$GrevyZebra.4 + GrevyZebra_100m$GrevyZebra.5 + GrevyZebra_100m$GrevyZebra.6 +
                GrevyZebra_100m$GrevyZebra.7 + GrevyZebra_100m$GrevyZebra.8)/8

GrevyZebraResistance <- ((GrevyZebra-1)*-1)
GrevyZebraResistance <- GrevyZebraResistance * 100
GrevyZebraResistance <- round(GrevyZebraResistance,0)
plot(GrevyZebraResistance)

#Add Anthropogenic areas 
GrevyZebraResistance2 <- cover(AntrRaster, GrevyZebraResistance)
plot(GrevyZebraResistance2)

#Add fences
# Reclassify the fece data to assign the resistance value
#Assign resistance values
# 1 = Porcupine fence
# 2 = Tall fence
# 3 = Tall netted fence
# 4 = Cattle fence
# 5 = Cattle fence owner repeals wildlife
# 6 = unknown fence
# 7 = ditch

fromToArrayGZ <- c(1,40,2,-9999,3,-9999,4,65,5,-9999,6,40,7,53)
fromToMatrixGZ <- matrix(fromToArrayGZ, ncol=2, byrow=TRUE)
fromToMatrixGZ

FencesGrevyZebra<-reclassify(Fences,fromToMatrixGZ)

#Add fences data keeping pixels with max resistance value
GrevyZebraResistance3<- mosaic(GrevyZebraResistance2,FencesGrevyZebra, fun=max)

#Add fences with infinite value
GrevyZebraResistance3[FencesGrevyZebra ==  -9999] <- -9999
plot(GrevyZebraResistance3)

GrevyZebraResistanceFin<-mask(GrevyZebraResistance3,StudyAreaRast)
plot(GrevyZebraResistanceFin)
writeRaster(GrevyZebraResistanceFin,"GrevyZebraResist_rast.asc",overwrite=TRUE)

########################
### Run Circuitscape  ##
########################

#With fences
# Make a place holder for the cs_run.exe path
CS_exe <- 'C:/"Program Files"/Circuitscape/cs_run.exe' # Don't forget the "Program Files" problem

# Make an .ini file
CS_ini <- c("[circuitscape options]",            
            "data_type = raster",
            "scenario = pairwise",
            "set_focal_node_currents_to_zero = True",
            "write_cur_maps = True",
            "write_cum_cur_map_only = True",
            "log_transform_maps = False",
            paste(c("point_file =",
                    "habitat_file =",
                    "output_file ="),
                  paste(getwd(),c("Centroids_rast.asc",
                                  "GrevyZebraResist_rast.asc",
                                  "GrevyZebra_CS.out"),
                        sep="/")))

# Write it to your working directory
writeLines(CS_ini,"myini.ini")

# Make the CS run cmd
CS_run <- paste(CS_exe, paste(getwd(),"myini.ini",sep="/")) # Make the cmd

# Run the command
system(CS_run)


################################
# Multi-species level analysis #
################################

# Load raster of estimated richness from occupancy models
richness<-brick("./GIS_Layers/SpRichness.tif")

# Resample to 100x100 m resolution
richness_100m<-resample(richness,a)
projection(richness_100m)<-CRS

rich <- (richness_100m$SpRichness.1 + richness_100m$SpRichness.2 + richness_100m$SpRichness.3 + 
           richness_100m$SpRichness.4 + richness_100m$SpRichness.5 + richness_100m$SpRichness.6 +
           richness_100m$SpRichness.7 + richness_100m$SpRichness.8)/8

#### Reclassify richness
Res<- (rich - 3.663573)/(10.20257-3.663573)*100
Res[Res<=1]<-1
Resistance<-(Res-100)*-1
Resistance<-round(Resistance,0)
plot(Resistance)


## Add anthropogenic areas
Resistance1 <- cover(AntrRaster, Resistance)
plot(Resistance1)

#Add fences

#Assign resistance values
# 1 = Porcupine fence
# 2 = Tall fence
# 3 = Tall netted fence
# 4 = Cattle fence
# 5 = Cattle fence owner repeals wildlife
# 6 = unknown fence
# 7 = ditch

fromToArraySR <- c(1,39,2,-9999,3,-9999,4,42,5,-9999,6,39,7,62)
fromToMatrixSR <- matrix(fromToArraySR, ncol=2, byrow=TRUE)

Fences2<-reclassify(Fences,fromToMatrixSR)
Fences2

#Add fences data keeping pixels with max resistance value
Resistance2<- mosaic(Resistance1,Fences2, fun=max)

#Add fences with infinite value
Resistance2[Fences2 ==  -9999] <- -9999
ResistanceFin<-round(Resistance2)
ResistanceFin<-mask(ResistanceFin,StudyAreaRast)
plot(ResistanceFin)

writeRaster(ResistanceFin,"comm_resist_rast.asc",overwrite=TRUE)

########################
### Run circuitscape  ##
########################

# Make a place holder for the cs_run.exe path
CS_exe <- 'C:/"Program Files"/Circuitscape/cs_run.exe' # Don't forget the "Program Files" problem

# Make an .ini file
CS_ini <- c("[circuitscape options]",            
            "data_type = raster",
            "scenario = pairwise",
            "set_focal_node_currents_to_zero = True",
            "write_cur_maps = True",
            "write_cum_cur_map_only = True",
            "log_transform_maps = False",
            paste(c("point_file =",
                    "habitat_file =",
                    "output_file ="),
                  paste(getwd(),c("Centroids_rast.asc",
                                  "comm_resist_rast.asc",
                                  "Community_CS.out"),
                        sep="/")))

# Write it to your working directory
writeLines(CS_ini,"myini.ini")

# Make the CS run cmd
CS_run <- paste(CS_exe, paste(getwd(),"myini.ini",sep="/")) # Make the cmd

# Run the command
system(CS_run)