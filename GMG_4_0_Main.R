#Procedural World Generator#
####Initialize####
setwd(dirname(rstudioapi::callFun("getActiveDocumentContext")$path)) #Working Directory
#Read in all required data
source("GMG_4_0_Functions.R")
source("GMG_4_0_Constants.R")

gRinstall() #BioConductor Packages to speed processing
if (!require("pacman")) {install.packages("pacman")} #Package Manager
pacman::p_load(beepr, bitops, colorspace, data.table, dplyr, ggplot2, ggrepel, gRbase, png, stringr, tidyverse) #Packages needed for script



####User Defined Variables####

#Custom Image: Uncomment the line below, and enter in the PNG filename of the image you would like to use instead of the milky way. This requires the image to be 576x576 to have the same scale.

#readGalaxy(filename = "YourFileHere.png")


####Packages####
# pacman::p_load(,data.table,ggrepel,magick,plotrix, cluster, factoextra, ggalt, beepr, ggforce, stringr)

#####Running Script####
#Wallart time.
#map <- generateMap()

#As separate functions, the jumpvectors can be calculated separately, which may be safer for now.
wall <- snapgaldata(14, planets=T)
wallc<- (snapgaldata(14, planets=T, garden=T))
jumpdata <- jumpvector(wall, test=T)

jmpp <- nationmapfullmap(wall, gardenworlddata = wallc, ps = jumpdata, size = 24,save = T, filename = "WallMap_9_Nations_Borderless.png", savecsv=F)
