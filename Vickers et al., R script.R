#### Scripts for Vickers et al. Sensitivity of migratory connectivity metrics to spatial sampling design ####

# This script file simulates breeding and non-breeding locations of a migratory popualiton with three different strengths of migratory connectiivty
# This is done for a single contiguous populaiton scenario and a patchy populaiton scenario
# Several spatial sampling designs are then applied to assess potential bias
# Code for manuscript figures is also included

# Certain elements of this script were conducted on a high permance cluster. 
# As such, these elements may take considerable time to run.

# These scripts were devloped on:
# R version 3.6.1 (2019-07-05) -- "Action of the Toes"
# Copyright (C) 2019 The R Foundation for Statistical Computing
# Platform: x86_64-apple-darwin15.6.0 (64-bit)

# Packages ----

if (!require('sp')) install.packages('sp'); library(sp)
if (!require('raster')) install.packages('raster'); library(raster)
if (!require('geosphere')) install.packages('geosphere'); library(geosphere)
if (!require('circular')) install.packages('circular'); library(circular)
if (!require('ade4')) install.packages('ade4'); library(ade4)
if (!require('reshape2')) install.packages('reshape2'); library(reshape2)
if (!require('plyr')) install.packages('plyr'); library(plyr)
if (!require('tidyr')) install.packages('tidyr'); library(tidyr)
if (!require('MigConnectivity')) install.packages('MigConnectivity'); library(MigConnectivity)
if (!require('devtools')) install.packages('devtools'); library(devtools)
if (!require('rgdal')) install.packages('rgdal'); library(rgdal)
if (!require('readxl')) install.packages('readxl'); library(readxl)
if (!require('tidyverse')) install.packages('tidyverse'); library(tidyverse)
if (!require('tidyr')) install.packages('tidyr'); library(tidyr)
if (!require('ggplot2')) install.packages('ggplot2'); library(ggplot2)
if (!require('doBy')) install.packages('doBy'); library(doBy)
if (!require('cowplot')) install.packages('cowplot'); library(cowplot)
if (!require('data.table')) install.packages('data.table'); library(data.table)
if (!require('ggpubr')) install.packages('ggpubr'); library(ggpubr)
options(scipen=99999999)

# Section 1: Simulating migratory populations (single breeding and non-breeding population) ----

# Each time this section is run, the resultant populaitons will vary slightly.
# The specific populations used for downstream analysis are availbale within the document repository.

# 1.1: Set up data and objects for a for loop #

# Set up the environment:
xdim <- 100 # x dimension of space
ydim <- xdim*2 # y dimension of space
grid <- raster(ncols=ydim,nrows=ydim) # create a grid of cells in raster format
extent(grid) <- c(0,1,0,ydim/xdim) # set up extent of the grid object
grid[]<-1 
xy <- data.frame(xyFromCell(grid,1:ncell(grid))) #extract coords of cell centres

n <- 11000 # set a total population size
n2 <- 10000 # set population size after restricting population

# breeding range 
brd_latrange <- c(quantile(xy$y,0.67),quantile(xy$y,1)) # Set latitudinal range of breeding population
brd_lonrange <- c(quantile(xy$x,0),quantile(xy$x,1)) # Set longitudinal range of breeding population

# non-breeding range
win_latrange <- c(quantile(xy$y,0),quantile(xy$y,0.33))  # Set latitudinal range of non-breeding population
win_lonrange  <- c(quantile(xy$x,0),quantile(xy$x,1))   # Set longitudinal range of non-breeding population

# Set scenario:
con.levs <- c(1,1.75,2.5) # Set values that will determine different strength levels of migratory connectivity in the population

# Set direction of travel for individuals. Setting this outside of the for loop means direction of movement is maintained across migratory connectivity strengths
direction <- runif(n,-180,180)

# Create total pop:

brd_pops <- list()
win_pops <- list()

# 1.2: for loop #

for(j in 1:3){
  #breeding locations:
  brd_locs <- data.frame(x=runif(n,brd_lonrange[1],brd_lonrange[2]),
                         y=runif(n,brd_latrange[1],brd_latrange[2]))

  #Simulate non-breeding locations:
  midbrd <- c(mean(brd_lonrange),mean(brd_latrange)) #midpoint of breeding range
  midwin <- c(mean(win_lonrange),mean(win_latrange)) #midpoint of non-breeding range
  mid_dist <- distm(midbrd,midwin)
  win_start <- destPointRhumb(brd_locs,b=180,d=mid_dist) #starting points
  win_dispersal <- rlnorm(n,8-con.levs[j],1)*25 # add non-breeding dispersal
  win_locs<- destPointRhumb(win_start, direction, d=win_dispersal)
  
  # Extreme outliers within the non-breeding zone are removed
  win_locs2 <- as.data.frame(win_locs)
  win_locs <- subset(win_locs2, subset = lon>(-1) & lon<2 & lat>(-1) & lat<1.67)
  # Removed those individuals from breeding dataframe
  brd_locs2 <- subset(brd_locs, rownames(brd_locs) %in% rownames(win_locs))
  # limit population size to 10,000
  brd_locs2<-brd_locs2[1:n2,]
  win_locs<-win_locs[1:n2,]
  
  brd_pops[[j]] <- brd_locs2
  win_pops[[j]] <- win_locs
}    

# Breeding and non-breeding populations were compiled into a single .csv file available in the file repository

# Section 2: Simulating migratory populations (patchy breeding and non-breeding populations) ----

# Each time this section is run, the resultant populaitons will vary slightly.
# The specific populations used for downstream analysis are availbale within the document repository.

rm(list = ls())

# 2.1: Set up data and objects for a for loop #

# Set up the environment:
xdim <- 100 # x dimension of space
ydim <- xdim*2 # y dimension of space
grid <- raster(ncols=ydim,nrows=ydim) # create a grid of cells in raster format
extent(grid) <- c(0,1,0,ydim/xdim) # set up extent of the grid object
grid[]<-1 
xy <- data.frame(xyFromCell(grid,1:ncell(grid))) #extract coords of cell centres

n <- 40000 # set a total population size
n2 <- 10000 # set population size after restricting population

# breeding range 
brd_latrange <- c(quantile(xy$y,0.67),quantile(xy$y,1)) # Set latitudinal range of breeding population
brd_lonrange <- c(quantile(xy$x,0),quantile(xy$x,1)) # Set longitudinal range of breeding population

# non-breeding range
win_latrange <- c(quantile(xy$y,0),quantile(xy$y,0.33))  # Set latitudinal range of non-breeding population
win_lonrange  <- c(quantile(xy$x,0),quantile(xy$x,1))   # Set longitudinal range of non-breeding population

# Set scenario:
con.levs <- c(1,1.75,2.5) # Set values that will determine different strength levels of migratory connectivity in the population

# Set direction of travel for individuals. Setting this outside of the for loop means direction of movement is maintained across migratory connectivity strengths
direction <- runif(n,-180,180)

# Create total pop:

brd_pops <- list()
win_pops <- list()

# 2.2: for loop #

for(j in 1:3){
  
  # simulatebreeding locations
  x1 <- runif(n/4,brd_lonrange[1],(brd_lonrange[1]+((brd_lonrange[2]-brd_lonrange[1])*.33)))
  x2 <- runif(n/4,brd_lonrange[1],(brd_lonrange[1]+((brd_lonrange[2]-brd_lonrange[1])*.33)))
  x3 <- runif(n/4,(brd_lonrange[1]+((brd_lonrange[2]-brd_lonrange[1])*.66)),brd_lonrange[2])
  x4 <- runif(n/4,(brd_lonrange[1]+((brd_lonrange[2]-brd_lonrange[1])*.66)),brd_lonrange[2])
  x <- c(x1,x2,x3,x4)
  
  y1 <- runif(n/4,brd_latrange[1],(brd_latrange[1]+((brd_latrange[2]-brd_latrange[1])*.33)))
  y2 <- runif(n/4,brd_latrange[1],(brd_latrange[1]+((brd_latrange[2]-brd_latrange[1])*.33)))
  y3 <- runif(n/4,(brd_latrange[1]+((brd_latrange[2]-brd_latrange[1])*.66)),brd_latrange[2])
  y4 <- runif(n/4,(brd_latrange[1]+((brd_latrange[2]-brd_latrange[1])*.66)),brd_latrange[2])
  y <- c(y1,y3,y2,y4)
  
  brd_locs <- data.frame(x=x,
                         y=y)
  
  #Simulate non-breeding locations
  midbrd <- c(mean(brd_lonrange),mean(brd_latrange)) # midpoint of breeding range
  midwin <- c(mean(win_lonrange),mean(win_latrange)) # midpoint of non-breeding range
  mid_dist <- distm(midbrd,midwin)
  win_start <- destPointRhumb(brd_locs,b=180,d=mid_dist)#starting points
  win_dispersal <- rlnorm(n,8-con.levs[j],1)*25 #add winter dispersal
  win_locs<- destPointRhumb(win_start, direction, d=win_dispersal)

  # Restrict non-breeding indivduals to patchy distribution
  win_locs2 <- as.data.frame(win_locs)
  
  win_locs_1 <- subset(win_locs2, lon > 0 & lon < .33 & lat > 0.4411333 & lat < 0.6617)
  win_locs_2 <- subset(win_locs2, lon > .66 & lon < 1 & lat > 0.4411333 & lat < 0.6617)
  win_locs_3 <- subset(win_locs2, lon > 0 & lon < .33 & lat > 0 & lat < 0.2205667)
  win_locs_4 <- subset(win_locs2, lon > .66 & lon < 1 & lat > 0 & lat < 0.2205667)
  
  win_locs <- rbind(win_locs_1,
                    win_locs_2,
                    win_locs_3,
                    win_locs_4)
  
  
  win_locs<-win_locs[sample((rownames(win_locs)), n2),]
  
  brd_locs2<-brd_locs[as.numeric(rownames(win_locs)),]
  
  brd_pops[[j]] <- brd_locs2
  win_pops[[j]] <- win_locs
}    

# Breeding and non-breeding populations were compiled into a single .csv file available in the file repository


# Section 3: Applying sampling designs (Single breeding and non-breeding popualtion) ----

# This section was iterated 100 times on a high performance cluster. 
# Mantel correaltions are highly memory intensive due to the exponential increase in matrix size with larger populations.
# As a result, a single iteration may take several hours to complete.

rm(list = ls())

## Set up data and objects for a for loop ##

# Set up the environment:
xdim <- 100 # x dimension of space
ydim <- xdim*2 # y dimension of space
grid <- raster(ncols=ydim,nrows=ydim) # create a grid of cells in raster format
extent(grid) <- c(0,1,0,ydim/xdim) # set up extent of the grid object
grid[]<-1 
xy <- data.frame(xyFromCell(grid,1:ncell(grid))) #extract coords of cell centres

# breeding range 
brd_latrange <- c(quantile(xy$y,0.67),quantile(xy$y,1)) 
brd_lonrange <- c(quantile(xy$x,0),quantile(xy$x,1)) 
# winter range
win_latrange <- c(quantile(xy$y,0),quantile(xy$y,0.33)) 
win_lonrange  <- c(quantile(xy$x,0),quantile(xy$x,1)) 

# Set scenario:
n <- 11000 # set a total population size
n2 <- 10000 # set population size after restricting population
con.levs <- seq(1,5,length.out=3) # connectivity levels - 1= weak, 5= strong (dictates mean (on a log scale) of lognormal dist)
nreps <- 100 # set the number of repeats 
samp_size <- 200
bigboxsize <- seq(0.5,1.5,length.out=3)
smallbox <- 0.12

# Create empty lists, data frames, and set object starting numbers. All needs to be done before each run of the loop
winloclist <- list() 
count <- 0 # resets count object back to 0
pls <- list()
counter <- 0 # resets counter object back to 0
ps <- list()
area <- list()
slp <- list()
q <- 5 # Counter used to fill data_sheet correctly 
f <- 8 # Counter used to fill data_sheet correctly 

# Create empty matrix for data
data_sheet <- data.frame(matrix(nrow=10,ncol=6)) # produces a data frame for the calculated mantel stats to be stored in. 4= number of variations in sampling spread/size/number
names(data_sheet) <- c("Mantel","Cohen","Samp_Size","Spread", "Area", "Scenario")
data_sheet$Scenario <- c('All', 'Allsmall', 'Allmedium', 'Alllarge', 'Spread','Spread','Spread','Area','Area','Area')
data_sheet_list <- list()

# Spread and area scenarios to allow descriptor filling in data_sheet
Spread <- c("Low","Medium", "High")
Area <- c("Small","Medium", "Large")

# Setting up for Cohen
nBreeding <- 3 #number of breeding regions
nNonBreeding <- 3 #number of non-breeding regions

# Set up extent of the boxes that delimit the spread of sampling.
# big box 1 (xmin, xmax, ymin, ymax)
bb1 <-  c(0.315925, 0.684075, 1.545161, 1.78814)
# big box 2 (xmin, xmax, ymin, ymax)
bb2 <-  c(0.19155, 0.80845, 1.463073, 1.870227)
# big box 3 (xmin, xmax, ymin, ymax)
bb3 <-  c(0.067175, 0.932825, 1.380986, 1.952315)
bb <- list(bb1,bb2,bb3)

setwd("") # Set working directory to folder with population breeding and non-breeding location files

brd_MC1 <- read.csv("brd_MC1.csv")
win_MC1 <- read.csv("win_MC1.csv")
brd_MC2 <- read.csv("brd_MC2.csv")
win_MC2 <- read.csv("win_MC2.csv")
brd_MC3 <- read.csv("brd_MC3.csv")
win_MC3 <- read.csv("win_MC3.csv")

brd <- list(brd_MC1, brd_MC2, brd_MC3)
win <- list(win_MC1, win_MC2, win_MC3)

for(g in 1:3){ # loops thorugh levels of migratory connectivity
  
  brd_locs <- brd[[g]]
  win_locs <- win[[g]]  
  
  #### Mantel Methodology for entire population #### - computationally intensive step -
  brddists <- dist(data.frame(brd_locs[1:2]))  
  windists <- dist(data.frame(win_locs[1:2]))
  dfbrd <- melt(as.matrix(brddists), varnames = c("row", "col"))
  dfwin <- melt(as.matrix(windists), varnames = c("row", "col"))
  data_sheet[1,1] <-  mantel.rtest(brddists,windists,nrepet=1)
  
  #### Cohen Methodology for entire population ####
  brdcohen <- data.frame(brd_locs)
  wincohen <- data.frame(win_locs)
  names(brdcohen) <- c('x','y','zone1')
  names(wincohen) <- c('x','y','zone2')
  brdcohen$zone1 <- 1
  brdcohen$zone1[brdcohen$x>quantile(brdcohen$x,0.33)] <- 2
  brdcohen$zone1[brdcohen$x>quantile(brdcohen$x,0.67)] <- 3
  wincohen$zone2 <- 1
  wincohen$zone2[wincohen$x>quantile(wincohen$x,0.33)] <- 2
  wincohen$zone2[wincohen$x>quantile(wincohen$x,0.67)] <- 3
  brdwincohen <- cbind(brdcohen,wincohen)
  brdwincohen$zone1 <- as.factor(brdwincohen$zone1)
  brdwincohen$zone2 <- as.factor(brdwincohen$zone2)
  bwcz1 <- subset(brdwincohen, zone1 == 1)
  bwcz2 <- subset(brdwincohen, zone1 == 2)
  bwcz3 <- subset(brdwincohen, zone1 == 3)
  bwcz1summary <- summary(bwcz1$zone2)
  bwcz1summary <- as.data.frame(bwcz1summary)
  step1 <- bwcz1summary[1,]/sum(summary(bwcz1$zone2))
  step2 <- bwcz1summary[2,]/sum(summary(bwcz1$zone2))
  step3 <- bwcz1summary[3,]/sum(summary(bwcz1$zone2))
  bwcz2summary <- summary(bwcz2$zone2)
  bwcz2summary <- as.data.frame(bwcz2summary)
  step4 <- bwcz2summary[1,]/sum(summary(bwcz2$zone2))
  step5 <- bwcz2summary[2,]/sum(summary(bwcz2$zone2))
  step6 <- bwcz2summary[3,]/sum(summary(bwcz2$zone2))
  bwcz3summary <- summary(bwcz3$zone2)
  bwcz3summary <- as.data.frame(bwcz3summary)
  step7 <- bwcz3summary[1,]/sum(summary(bwcz3$zone2))
  step8 <- bwcz3summary[2,]/sum(summary(bwcz3$zone2))
  step9 <- bwcz3summary[3,]/sum(summary(bwcz3$zone2))
  
  #transition probabilities form each breeding to each non-breeding region
  Cohenspsi <- matrix(c(step1,step2,step3,step4,step5,step6,step7,step8,step9), nBreeding, nNonBreeding, byrow = TRUE)
  Cohenspsi[is.na(Cohenspsi)] <- 0
  
  # Define the relative abundance within the three breeding regions
  relN <- c(sum(summary(bwcz1$zone2)/10000),sum(summary(bwcz2$zone2)/10000),sum(summary(bwcz3$zone2)/10000))
  
  # set distance matrices for cohen
  breedDist <- matrix(c((quantile(brdcohen$x,0.165)-quantile(brdcohen$x,0.165)), (quantile(brdcohen$x,0.495)-quantile(brdcohen$x,0.165)), (quantile(brdcohen$x,0.835)-quantile(brdcohen$x,0.165)),
                        (quantile(brdcohen$x,0.495)-quantile(brdcohen$x,0.165)), (quantile(brdcohen$x,0.495)-quantile(brdcohen$x,0.495)), (quantile(brdcohen$x,0.835)-quantile(brdcohen$x,0.495)),
                        (quantile(brdcohen$x,0.835)-quantile(brdcohen$x,0.165)), (quantile(brdcohen$x,0.835)-quantile(brdcohen$x,0.495)), (quantile(brdcohen$x,0.835)-quantile(brdcohen$x,0.835))), nBreeding, nBreeding)
  nonBreedDist <- matrix(c((quantile(wincohen$x,0.165)-quantile(wincohen$x,0.165)), (quantile(wincohen$x,0.495)-quantile(wincohen$x,0.165)), (quantile(wincohen$x,0.835)-quantile(wincohen$x,0.165)),
                           (quantile(wincohen$x,0.495)-quantile(wincohen$x,0.165)), (quantile(wincohen$x,0.495)-quantile(wincohen$x,0.495)), (quantile(wincohen$x,0.835)-quantile(wincohen$x,0.495)),
                           (quantile(wincohen$x,0.835)-quantile(wincohen$x,0.165)), (quantile(wincohen$x,0.835)-quantile(wincohen$x,0.495)), (quantile(wincohen$x,0.835)-quantile(wincohen$x,0.835))), nNonBreeding, nNonBreeding)
  
  # breeding polygons
  
  x1 <- c(quantile(brdcohen$x,0),  quantile(brdcohen$x,0),  quantile(brdcohen$x,0.33), quantile(brdcohen$x,0.33), quantile(brdcohen$x,0))
  y <- c(quantile(brdcohen$y,0),  quantile(brdcohen$y,1),  quantile(brdcohen$y,1), quantile(brdcohen$y,0), quantile(brdcohen$y,0))
  xy1 <- cbind(x1, y)
  
  x2 <- c(quantile(brdcohen$x,0.33),  quantile(brdcohen$x,0.33),  quantile(brdcohen$x,0.66), quantile(brdcohen$x,0.66), quantile(brdcohen$x,0.33))
  xy2 <- cbind(x2, y)
  
  x3 <- c(quantile(brdcohen$x,0.66),  quantile(brdcohen$x,0.66),  quantile(brdcohen$x,1), quantile(brdcohen$x,1), quantile(brdcohen$x,0.66))
  xy3 <- cbind(x3, y)
  
  p1 <- Polygon(xy1)
  p2 <- Polygon(xy2)
  p3 <- Polygon(xy3)
  ps1 <- Polygons(list(p1), ID = 1)
  ps2 <- Polygons(list(p2), ID = 2)
  ps3 <- Polygons(list(p3), ID = 3)
  
  originsites <- SpatialPolygons(list(ps1,ps2,ps3))
  
  x1 <- c(quantile(wincohen$x,0),  quantile(wincohen$x,0),  quantile(wincohen$x,0.33), quantile(wincohen$x,0.33), quantile(wincohen$x,0))
  y <- c(quantile(wincohen$y,0),  quantile(wincohen$y,1),  quantile(wincohen$y,1), quantile(wincohen$y,0), quantile(wincohen$y,0))
  xy1 <- cbind(x1, y)
  
  x2 <- c(quantile(wincohen$x,0.33),  quantile(wincohen$x,0.33),  quantile(wincohen$x,0.66), quantile(wincohen$x,0.66), quantile(wincohen$x,0.33))
  xy2 <- cbind(x2, y)
  
  x3 <- c(quantile(wincohen$x,0.66),  quantile(wincohen$x,0.66),  quantile(wincohen$x,1), quantile(wincohen$x,1), quantile(wincohen$x,0.66))
  xy3 <- cbind(x3, y)
  
  p1 <- Polygon(xy1)
  p2 <- Polygon(xy2)
  p3 <- Polygon(xy3)
  ps1 <- Polygons(list(p1), ID = 1)
  ps2 <- Polygons(list(p2), ID = 2)
  ps3 <- Polygons(list(p3), ID = 3) 
  
  targetsites <- SpatialPolygons(list(ps1,ps2,ps3))
  
  proj4string(originsites) <- CRS("+proj=eqc +lat_ts=0 +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +a=6371007 +b=6371007 +units=m +no_defs ")
  raster::crs(targetsites)<-raster::crs(originsites)
  
  spatbreeding <- SpatialPoints(brdcohen[,1:2])
  raster::crs(spatbreeding)<-raster::crs(originsites)
  spatwinter <- SpatialPoints(wincohen[,1:2])
  raster::crs(spatwinter)<-raster::crs(originsites)
  
  GPS_mc<-estMC(isGL=FALSE, # Logical vector: light-level geolocator(T)/GPS(F)
                targetDist = nonBreedDist, # targetSites distance matrix
                originDist = breedDist, # originSites distance matrix
                originRelAbund = relN, #Origin relative abund.
                targetSites = targetsites, # Non-breeding target sites
                originSites = originsites, # Breeding origin sites
                originPoints = spatbreeding, # Capture Locations
                targetPoints = spatwinter, # Device target locations
                verbose = 0,   # output options
                nSamples = 1000) # This is set low for example
  
  data_sheet[1,2] <- GPS_mc$meanMC
  
  data_sheet[1,3] <- n2
  
  ####################################################
  # Sampling the population --------------
  
  # Delimit sampling areas (breeding range):
  
  for(t in 1:3){ # loop through spread scenarios
    
    # Create Big box - whole study area (no sampling outside this)
    bigbox <- bigboxsize[t] # sensible range 0.5 to 1.5. This would be added and small box fixed to allow big box to change but fix small box size.
    bigboxdims <- c(bigbox*diff(brd_lonrange),bigbox*diff(brd_latrange))
    midbrd <- c(mean(brd_lonrange),mean(brd_latrange))#midpoint of breeding range
    midwin <- c(mean(win_lonrange),mean(win_latrange))
    area1 <- c(midbrd[1]-(0.5*bigboxdims[1]),midbrd[1]+(0.5*bigboxdims[1]),
               midbrd[2]-(0.5*bigboxdims[2]),midbrd[2]+(0.5*bigboxdims[2]))
    
    # #### Mantel methodlogy for inds within big box extent ####
    bbbox <- bb[[t]]
    bigxcoords <- c(bbbox[1], bbbox[1], bbbox[2], bbbox[2], bbbox[1])
    bigycoords <- c(bbbox[3],bbbox[4],bbbox[4],bbbox[3],bbbox[3])
    bigboxdata <- data.frame(x=bigxcoords,y=bigycoords)
    bigboxpoly <- Polygon(bigboxdata)
    bigboxlist <- Polygons(list(bigboxpoly),1)
    bigboxplot <- SpatialPolygons(list(bigboxlist))
    try(coordinates(brd_locs)<-c("x","y"), silent=T)
    biginds <- as.numeric(which(over(brd_locs,bigboxplot)==1))
    big_brd <- data.frame(brd_locs[biginds,])
    big_win <- data.frame(win_locs[biginds,])
    bigboxbreed <- dist(big_brd[1:2]) 
    bigboxwint <- dist(big_win[1:2])
    bigboxmantel <-  mantel.rtest(bigboxbreed,bigboxwint,nrepet=1)
    data_sheet[(1+t),1] <- bigboxmantel$obs 
    
    #### Cohen methodology for inds within big box extent ####
    
    brdBBcohen <- data.frame(brd_locs[biginds,])
    winBBcohen <- data.frame(win_locs[biginds,])
    names(brdBBcohen) <- c('x','y','zone1')
    names(winBBcohen) <- c('x','y','zone2')
    brdBBcohen$zone1 <- 1
    brdBBcohen$zone1[brdBBcohen$x>quantile(brdBBcohen$x,0.33)] <- 2
    brdBBcohen$zone1[brdBBcohen$x>quantile(brdBBcohen$x,0.67)] <- 3
    winBBcohen$zone2 <- 1
    winBBcohen$zone2[winBBcohen$x>quantile(winBBcohen$x,0.33)] <- 2
    winBBcohen$zone2[winBBcohen$x>quantile(winBBcohen$x,0.67)] <- 3
    brdwinBBcohen <- cbind(brdBBcohen,winBBcohen)
    brdwinBBcohen$zone1 <- as.factor(brdwinBBcohen$zone1)
    brdwinBBcohen$zone2 <- as.factor(brdwinBBcohen$zone2)
    bwBBcz1 <- subset(brdwinBBcohen, zone1 == 1)
    bwBBcz2 <- subset(brdwinBBcohen, zone1 == 2)
    bwBBcz3 <- subset(brdwinBBcohen, zone1 == 3)
    bwBBcz1summary <- summary(bwBBcz1$zone2)
    bwBBcz1summary <- as.data.frame(bwBBcz1summary)
    steps1 <- bwBBcz1summary[1,]/sum(summary(bwBBcz1$zone2))
    steps2 <- bwBBcz1summary[2,]/sum(summary(bwBBcz1$zone2))
    steps3 <- bwBBcz1summary[3,]/sum(summary(bwBBcz1$zone2))
    bwBBcz2summary <- summary(bwBBcz2$zone2)
    bwBBcz2summary <- as.data.frame(bwBBcz2summary)
    steps4 <- bwBBcz2summary[1,]/sum(summary(bwBBcz2$zone2))
    steps5 <- bwBBcz2summary[2,]/sum(summary(bwBBcz2$zone2))
    steps6 <- bwBBcz2summary[3,]/sum(summary(bwBBcz2$zone2))
    bwBBcz3summary <- summary(bwBBcz3$zone2)
    bwBBcz3summary <- as.data.frame(bwBBcz3summary)
    steps7 <- bwBBcz3summary[1,]/sum(summary(bwBBcz3$zone2))
    steps8 <- bwBBcz3summary[2,]/sum(summary(bwBBcz3$zone2))
    steps9 <- bwBBcz3summary[3,]/sum(summary(bwBBcz3$zone2))
    
    #transition probabilities form each breeding to each non-breeding region
    BBspsi <- matrix(c(steps1,steps2,steps3,steps4,steps5,steps6,steps7,steps8,steps9), nBreeding, nNonBreeding, byrow = TRUE)
    BBspsi[is.na(BBspsi)] <- 0
    
    # Define the relative abundance within the three breeding regions
    BBsrelN <- c(sum(summary(bwBBcz1$zone2)/length(biginds)),sum(summary(bwBBcz2$zone2)/length(biginds)),sum(summary(bwBBcz3$zone2)/length(biginds)))
    
    # set distance matrices for cohen
    breedDist <- matrix(c((quantile(brdBBcohen$x,0.165)-quantile(brdBBcohen$x,0.165)), (quantile(brdBBcohen$x,0.495)-quantile(brdBBcohen$x,0.165)), (quantile(brdBBcohen$x,0.835)-quantile(brdBBcohen$x,0.165)),
                          (quantile(brdBBcohen$x,0.495)-quantile(brdBBcohen$x,0.165)), (quantile(brdBBcohen$x,0.495)-quantile(brdBBcohen$x,0.495)), (quantile(brdBBcohen$x,0.835)-quantile(brdBBcohen$x,0.495)),
                          (quantile(brdBBcohen$x,0.835)-quantile(brdBBcohen$x,0.165)), (quantile(brdBBcohen$x,0.835)-quantile(brdBBcohen$x,0.495)), (quantile(brdBBcohen$x,0.835)-quantile(brdBBcohen$x,0.835))), nBreeding, nBreeding)
    nonBreedDist <- matrix(c((quantile(winBBcohen$x,0.165)-quantile(winBBcohen$x,0.165)), (quantile(winBBcohen$x,0.495)-quantile(winBBcohen$x,0.165)), (quantile(winBBcohen$x,0.835)-quantile(winBBcohen$x,0.165)),
                             (quantile(winBBcohen$x,0.495)-quantile(winBBcohen$x,0.165)), (quantile(winBBcohen$x,0.495)-quantile(winBBcohen$x,0.495)), (quantile(winBBcohen$x,0.835)-quantile(winBBcohen$x,0.495)),
                             (quantile(winBBcohen$x,0.835)-quantile(winBBcohen$x,0.165)), (quantile(winBBcohen$x,0.835)-quantile(winBBcohen$x,0.495)), (quantile(winBBcohen$x,0.835)-quantile(winBBcohen$x,0.835))), nNonBreeding, nNonBreeding)
    
    x1 <- c(quantile(brdBBcohen$x,0),  quantile(brdBBcohen$x,0),  quantile(brdBBcohen$x,0.33), quantile(brdBBcohen$x,0.33), quantile(brdBBcohen$x,0))
    y <- c(quantile(brdBBcohen$y,0),  quantile(brdBBcohen$y,1),  quantile(brdBBcohen$y,1), quantile(brdBBcohen$y,0), quantile(brdBBcohen$y,0))
    xy1 <- cbind(x1, y)
    
    x2 <- c(quantile(brdBBcohen$x,0.33),  quantile(brdBBcohen$x,0.33),  quantile(brdBBcohen$x,0.66), quantile(brdBBcohen$x,0.66), quantile(brdBBcohen$x,0.33))
    xy2 <- cbind(x2, y)
    
    x3 <- c(quantile(brdBBcohen$x,0.66),  quantile(brdBBcohen$x,0.66),  quantile(brdBBcohen$x,1), quantile(brdBBcohen$x,1), quantile(brdBBcohen$x,0.66))
    xy3 <- cbind(x3, y)
    
    p1 <- Polygon(xy1)
    p2 <- Polygon(xy2)
    p3 <- Polygon(xy3)
    ps1 <- Polygons(list(p1), ID = 1)
    ps2 <- Polygons(list(p2), ID = 2)
    ps3 <- Polygons(list(p3), ID = 3)
    
    originsites <- SpatialPolygons(list(ps1,ps2,ps3))
    
    x1 <- c(quantile(winBBcohen$x,0),  quantile(winBBcohen$x,0),  quantile(winBBcohen$x,0.33), quantile(winBBcohen$x,0.33), quantile(winBBcohen$x,0))
    y <- c(quantile(winBBcohen$y,0),  quantile(winBBcohen$y,1),  quantile(winBBcohen$y,1), quantile(winBBcohen$y,0), quantile(winBBcohen$y,0))
    xy1 <- cbind(x1, y)
    
    x2 <- c(quantile(winBBcohen$x,0.33),  quantile(winBBcohen$x,0.33),  quantile(winBBcohen$x,0.66), quantile(winBBcohen$x,0.66), quantile(winBBcohen$x,0.33))
    xy2 <- cbind(x2, y)
    
    x3 <- c(quantile(winBBcohen$x,0.66),  quantile(winBBcohen$x,0.66),  quantile(winBBcohen$x,1), quantile(winBBcohen$x,1), quantile(winBBcohen$x,0.66))
    xy3 <- cbind(x3, y)
    
    p1 <- Polygon(xy1)
    p2 <- Polygon(xy2)
    p3 <- Polygon(xy3)
    ps1 <- Polygons(list(p1), ID = 1)
    ps2 <- Polygons(list(p2), ID = 2)
    ps3 <- Polygons(list(p3), ID = 3) 
    
    targetsites <- SpatialPolygons(list(ps1,ps2,ps3))
    
    proj4string(originsites) <- CRS("+proj=eqc +lat_ts=0 +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +a=6371007 +b=6371007 +units=m +no_defs ")
    raster::crs(targetsites)<-raster::crs(originsites)
    
    spatbreeding <- SpatialPoints(brdBBcohen[,1:2])
    raster::crs(spatbreeding)<-raster::crs(originsites)
    spatwinter <- SpatialPoints(winBBcohen[,1:2])
    raster::crs(spatwinter)<-raster::crs(originsites)
    
    GPS_mc<-estMC(isGL=FALSE, # Logical vector: light-level geolocator(T)/GPS(F)
                  targetDist = nonBreedDist, # targetSites distance matrix
                  originDist = breedDist, # originSites distance matrix
                  originRelAbund = BBsrelN, #Origin relative abund.
                  targetSites = targetsites, # Non-breeding target sites
                  originSites = originsites, # Breeding origin sites
                  originPoints = spatbreeding, # Capture Locations
                  targetPoints = spatwinter, # Device target locations
                  verbose = 0,   # output options
                  nSamples = 1000) # This is set low for example
    
    data_sheet[(1+t),2] <- GPS_mc$meanMC
    data_sheet[(1+t),3] <- length(biginds)
    data_sheet[(1+t),4] <- Spread[t]
    data_sheet[(1+t),5] <- Area[t]
    
    # Create smaller boxes of a given size that will be where indivudals are sample from
    nboxes <- sqrt(9) # 9 = number of sampled boxes
    smallboxdims <- c(smallbox*diff(brd_lonrange),smallbox*diff(brd_latrange))
    xcoords <- seq(area1[1],area1[2],length.out=nboxes+2)
    ycoords <- seq(area1[3],area1[4],length.out=nboxes+2)
    xcoords <- xcoords[2:(length(xcoords)-1)]
    ycoords <- ycoords[2:(length(ycoords)-1)]
    
    # create small boxes
    for(i in 1:length(xcoords)){
      for(j in 1:length(ycoords)){
        counter<-counter+1
        midpoint <- c(xcoords[i],ycoords[j])
        boxcorners <- c(midpoint[1]-(0.5*smallboxdims[1]),midpoint[1]+(0.5*smallboxdims[1]),
                        midpoint[2]-(0.5*smallboxdims[2]),midpoint[2]+(0.5*smallboxdims[2]))
        xcs <- c(boxcorners[1],boxcorners[1],boxcorners[2],boxcorners[2],boxcorners[1])
        ycs <- c(boxcorners[3],boxcorners[4],boxcorners[4],boxcorners[3],boxcorners[3])
        xym <- data.frame(x=xcs,y=ycs)
        p <- Polygon(xym)
        pls[[counter]]<-p
      }
    }
    slp <- pls
    ps <- Polygons(slp,1)
    area <- SpatialPolygons(list(ps))
    pls <- list()
    counter <- 0 # reset to zero
    count <- 0 # reset to zero
    
    # loop through sample sizes (now fixed at 200)
    sampleN <- samp_size
    
    ##### Sample populations and calculate mantel statistics
    fixN <- T 
    try(coordinates(brd_locs)<-c("x","y"), silent=T)
    
    # individuals for spread scenario ---------------------------------
    inds <- as.numeric(which(over(brd_locs,area)==1)) #pick individuals in the study area using over() from sp package
    samp_n<- ifelse(fixN==T,min(sampleN,length(inds)),floor(sampleprop*length(inds))) #number of individuals to sample
    samp_inds <- inds[sample(1:length(inds),samp_n,replace=F)] #randomly sample the pop
    try(coordinates(win_locs)<-c("lon","lat"), silent=T) 
    samp_brd <- data.frame(brd_locs[samp_inds,])
    samp_win <- data.frame(win_locs[samp_inds,])
    brddistssubset <- dist(samp_brd[1:2]) # 
    windistssubset <- dist(samp_win[1:2])
    mantelsubset <-  mantel.rtest(brddistssubset,windistssubset,nrepet=1)
    data_sheet[q,1] <-mantelsubset$obs #extract mantel statistic - needs updating if scenarios run in a loop
    data_sheet[q,3] <- length(samp_inds)
    
    #### Cohen methodology for sampled inds ####
    brdsubcohen <- data.frame(brd_locs[samp_inds,])
    winsubcohen <- data.frame(win_locs[samp_inds,])
    names(brdsubcohen) <- c('x','y','zone1')
    names(winsubcohen) <- c('x','y','zone2')
    brdsubcohen$zone1 <- 1
    brdsubcohen$zone1[brdsubcohen$x>quantile(brdsubcohen$x,0.33)] <- 2
    brdsubcohen$zone1[brdsubcohen$x>quantile(brdsubcohen$x,0.67)] <- 3
    winsubcohen$zone2 <- 1
    winsubcohen$zone2[winsubcohen$x>quantile(winsubcohen$x,0.33)] <- 2
    winsubcohen$zone2[winsubcohen$x>quantile(winsubcohen$x,0.67)] <- 3
    brdwinsubcohen <- cbind(brdsubcohen,winsubcohen)
    brdwinsubcohen$zone1 <- as.factor(brdwinsubcohen$zone1)
    brdwinsubcohen$zone2 <- as.factor(brdwinsubcohen$zone2)
    bwscz1 <- subset(brdwinsubcohen, zone1 == 1)
    bwscz2 <- subset(brdwinsubcohen, zone1 == 2)
    bwscz3 <- subset(brdwinsubcohen, zone1 == 3)
    bwscz1summary <- summary(bwscz1$zone2)
    bwscz1summary <- as.data.frame(bwscz1summary)
    stp1 <- bwscz1summary[1,]/sum(summary(bwscz1$zone2))
    stp2 <- bwscz1summary[2,]/sum(summary(bwscz1$zone2))
    stp3 <- bwscz1summary[3,]/sum(summary(bwscz1$zone2))
    bwscz2summary <- summary(bwscz2$zone2)
    bwscz2summary <- as.data.frame(bwscz2summary)
    stp4 <- bwscz2summary[1,]/sum(summary(bwscz2$zone2))
    stp5 <- bwscz2summary[2,]/sum(summary(bwscz2$zone2))
    stp6 <- bwscz2summary[3,]/sum(summary(bwscz2$zone2))
    bwscz3summary <- summary(bwscz3$zone2)
    bwscz3summary <- as.data.frame(bwscz3summary)
    stp7 <- bwscz3summary[1,]/sum(summary(bwscz3$zone2))
    stp8 <- bwscz3summary[2,]/sum(summary(bwscz3$zone2))
    stp9 <- bwscz3summary[3,]/sum(summary(bwscz3$zone2))
    
    #transition probabilities form each breeding to each non-breeding region
    spsi <- matrix(c(stp1,stp2,stp3,stp4,stp5,stp6,stp7,stp8,stp9), nBreeding, nNonBreeding, byrow = TRUE)
    spsi[is.na(spsi)] <- 0
    
    # Define the relative abundance within the three breeding regions
    srelN <- c(sum(summary(bwscz1$zone2)/samp_n),sum(summary(bwscz2$zone2)/samp_n),sum(summary(bwscz3$zone2)/samp_n))
    
    # set distance matrices for cohen
    breedDist <- matrix(c((quantile(brdsubcohen$x,0.165)-quantile(brdsubcohen$x,0.165)), (quantile(brdsubcohen$x,0.495)-quantile(brdsubcohen$x,0.165)), (quantile(brdsubcohen$x,0.835)-quantile(brdsubcohen$x,0.165)),
                          (quantile(brdsubcohen$x,0.495)-quantile(brdsubcohen$x,0.165)), (quantile(brdsubcohen$x,0.495)-quantile(brdsubcohen$x,0.495)), (quantile(brdsubcohen$x,0.835)-quantile(brdsubcohen$x,0.495)),
                          (quantile(brdsubcohen$x,0.835)-quantile(brdsubcohen$x,0.165)), (quantile(brdsubcohen$x,0.835)-quantile(brdsubcohen$x,0.495)), (quantile(brdsubcohen$x,0.835)-quantile(brdsubcohen$x,0.835))), nBreeding, nBreeding)
    nonBreedDist <- matrix(c((quantile(winsubcohen$x,0.165)-quantile(winsubcohen$x,0.165)), (quantile(winsubcohen$x,0.495)-quantile(winsubcohen$x,0.165)), (quantile(winsubcohen$x,0.835)-quantile(winsubcohen$x,0.165)),
                             (quantile(winsubcohen$x,0.495)-quantile(winsubcohen$x,0.165)), (quantile(winsubcohen$x,0.495)-quantile(winsubcohen$x,0.495)), (quantile(winsubcohen$x,0.835)-quantile(winsubcohen$x,0.495)),
                             (quantile(winsubcohen$x,0.835)-quantile(winsubcohen$x,0.165)), (quantile(winsubcohen$x,0.835)-quantile(winsubcohen$x,0.495)), (quantile(winsubcohen$x,0.835)-quantile(winsubcohen$x,0.835))), nNonBreeding, nNonBreeding)
    
    
    x1 <- c(quantile(brdsubcohen$x,0),  quantile(brdsubcohen$x,0),  quantile(brdsubcohen$x,0.33), quantile(brdsubcohen$x,0.33), quantile(brdsubcohen$x,0))
    y <- c(quantile(brdsubcohen$y,0),  quantile(brdsubcohen$y,1),  quantile(brdsubcohen$y,1), quantile(brdsubcohen$y,0), quantile(brdsubcohen$y,0))
    xy1 <- cbind(x1, y)
    
    x2 <- c(quantile(brdsubcohen$x,0.33),  quantile(brdsubcohen$x,0.33),  quantile(brdsubcohen$x,0.66), quantile(brdsubcohen$x,0.66), quantile(brdsubcohen$x,0.33))
    xy2 <- cbind(x2, y)
    
    x3 <- c(quantile(brdsubcohen$x,0.66),  quantile(brdsubcohen$x,0.66),  quantile(brdsubcohen$x,1), quantile(brdsubcohen$x,1), quantile(brdsubcohen$x,0.66))
    xy3 <- cbind(x3, y)
    
    p1 <- Polygon(xy1)
    p2 <- Polygon(xy2)
    p3 <- Polygon(xy3)
    ps1 <- Polygons(list(p1), ID = 1)
    ps2 <- Polygons(list(p2), ID = 2)
    ps3 <- Polygons(list(p3), ID = 3)
    
    originsites <- SpatialPolygons(list(ps1,ps2,ps3))
    
    x1 <- c(quantile(winsubcohen$x,0),  quantile(winsubcohen$x,0),  quantile(winsubcohen$x,0.33), quantile(winsubcohen$x,0.33), quantile(winsubcohen$x,0))
    y <- c(quantile(winsubcohen$y,0),  quantile(winsubcohen$y,1),  quantile(winsubcohen$y,1), quantile(winsubcohen$y,0), quantile(winsubcohen$y,0))
    xy1 <- cbind(x1, y)
    
    x2 <- c(quantile(winsubcohen$x,0.33),  quantile(winsubcohen$x,0.33),  quantile(winsubcohen$x,0.66), quantile(winsubcohen$x,0.66), quantile(winsubcohen$x,0.33))
    xy2 <- cbind(x2, y)
    
    x3 <- c(quantile(winsubcohen$x,0.66),  quantile(winsubcohen$x,0.66),  quantile(winsubcohen$x,1), quantile(winsubcohen$x,1), quantile(winsubcohen$x,0.66))
    xy3 <- cbind(x3, y)
    
    p1 <- Polygon(xy1)
    p2 <- Polygon(xy2)
    p3 <- Polygon(xy3)
    ps1 <- Polygons(list(p1), ID = 1)
    ps2 <- Polygons(list(p2), ID = 2)
    ps3 <- Polygons(list(p3), ID = 3) 
    
    targetsites <- SpatialPolygons(list(ps1,ps2,ps3))
    
    proj4string(originsites) <- CRS("+proj=eqc +lat_ts=0 +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +a=6371007 +b=6371007 +units=m +no_defs ")
    raster::crs(targetsites)<-raster::crs(originsites)
    
    spatbreeding <- SpatialPoints(brdsubcohen[,1:2])
    raster::crs(spatbreeding)<-raster::crs(originsites)
    spatwinter <- SpatialPoints(winsubcohen[,1:2])
    raster::crs(spatwinter)<-raster::crs(originsites)
    
    GPS_mc<-estMC(isGL=FALSE, # Logical vector: light-level geolocator(T)/GPS(F)
                  targetDist = nonBreedDist, # targetSites distance matrix
                  originDist = breedDist, # originSites distance matrix
                  originRelAbund = BBsrelN, #Origin relative abund.
                  targetSites = targetsites, # Non-breeding target sites
                  originSites = originsites, # Breeding origin sites
                  originPoints = spatbreeding, # Capture Locations
                  targetPoints = spatwinter, # Device target locations
                  verbose = 0,   # output options
                  nSamples = 1000) # This is set low for example
    
    data_sheet[q,2] <- GPS_mc$meanMC
    data_sheet[q,4] <- Spread[t]
    data_sheet[q,5] <- Area[t]
    
    # indiivudals for area scenario ---------------------------------
    bigsamp_n<- ifelse(fixN==T,min(sampleN,length(biginds)),floor(sampleprop*length(biginds))) #number of individuals to sample
    bigsamp_inds <- biginds[sample(1:length(biginds),bigsamp_n,replace=F)] #randomly sample the pop
    bigsamp_brd <- data.frame(brd_locs[bigsamp_inds,])
    bigsamp_win <- data.frame(win_locs[bigsamp_inds,]) 
    brddistssubset <- dist(bigsamp_brd[1:2]) # 
    windistssubset <- dist(bigsamp_win[1:2])
    mantelsubset <-  mantel.rtest(brddistssubset,windistssubset,nrepet=1)
    data_sheet[f,1] <- mantelsubset$obs #extract mantel statistic 
    data_sheet[f,3] <- length(bigsamp_inds)
    
    #### Cohen methodology for sampled inds ####
    brdsubcohen <- data.frame(brd_locs[bigsamp_inds,])
    winsubcohen <- data.frame(win_locs[bigsamp_inds,])
    names(brdsubcohen) <- c('x','y','zone1')
    names(winsubcohen) <- c('x','y','zone2')
    brdsubcohen$zone1 <- 1
    brdsubcohen$zone1[brdsubcohen$x>quantile(brdsubcohen$x,0.33)] <- 2
    brdsubcohen$zone1[brdsubcohen$x>quantile(brdsubcohen$x,0.67)] <- 3
    winsubcohen$zone2 <- 1
    winsubcohen$zone2[winsubcohen$x>quantile(winsubcohen$x,0.33)] <- 2
    winsubcohen$zone2[winsubcohen$x>quantile(winsubcohen$x,0.67)] <- 3
    brdwinsubcohen <- cbind(brdsubcohen,winsubcohen)
    brdwinsubcohen$zone1 <- as.factor(brdwinsubcohen$zone1)
    brdwinsubcohen$zone2 <- as.factor(brdwinsubcohen$zone2)
    bwscz1 <- subset(brdwinsubcohen, zone1 == 1)
    bwscz2 <- subset(brdwinsubcohen, zone1 == 2)
    bwscz3 <- subset(brdwinsubcohen, zone1 == 3)
    bwscz1summary <- summary(bwscz1$zone2)
    bwscz1summary <- as.data.frame(bwscz1summary)
    stp1 <- bwscz1summary[1,]/sum(summary(bwscz1$zone2))
    stp2 <- bwscz1summary[2,]/sum(summary(bwscz1$zone2))
    stp3 <- bwscz1summary[3,]/sum(summary(bwscz1$zone2))
    bwscz2summary <- summary(bwscz2$zone2)
    bwscz2summary <- as.data.frame(bwscz2summary)
    stp4 <- bwscz2summary[1,]/sum(summary(bwscz2$zone2))
    stp5 <- bwscz2summary[2,]/sum(summary(bwscz2$zone2))
    stp6 <- bwscz2summary[3,]/sum(summary(bwscz2$zone2))
    bwscz3summary <- summary(bwscz3$zone2)
    bwscz3summary <- as.data.frame(bwscz3summary)
    stp7 <- bwscz3summary[1,]/sum(summary(bwscz3$zone2))
    stp8 <- bwscz3summary[2,]/sum(summary(bwscz3$zone2))
    stp9 <- bwscz3summary[3,]/sum(summary(bwscz3$zone2))
    
    #transition probabilities form each breeding to each non-breeding region
    spsi <- matrix(c(stp1,stp2,stp3,stp4,stp5,stp6,stp7,stp8,stp9), nBreeding, nNonBreeding, byrow = TRUE)
    spsi[is.na(spsi)] <- 0
    
    # Define the relative abundance within the three breeding regions
    srelN <- c(sum(summary(bwscz1$zone2)/samp_n),sum(summary(bwscz2$zone2)/samp_n),sum(summary(bwscz3$zone2)/samp_n))
    
    # set distance matrices for cohen
    breedDist <- matrix(c((quantile(brdsubcohen$x,0.165)-quantile(brdsubcohen$x,0.165)), (quantile(brdsubcohen$x,0.495)-quantile(brdsubcohen$x,0.165)), (quantile(brdsubcohen$x,0.835)-quantile(brdsubcohen$x,0.165)),
                          (quantile(brdsubcohen$x,0.495)-quantile(brdsubcohen$x,0.165)), (quantile(brdsubcohen$x,0.495)-quantile(brdsubcohen$x,0.495)), (quantile(brdsubcohen$x,0.835)-quantile(brdsubcohen$x,0.495)),
                          (quantile(brdsubcohen$x,0.835)-quantile(brdsubcohen$x,0.165)), (quantile(brdsubcohen$x,0.835)-quantile(brdsubcohen$x,0.495)), (quantile(brdsubcohen$x,0.835)-quantile(brdsubcohen$x,0.835))), nBreeding, nBreeding)
    nonBreedDist <- matrix(c((quantile(winsubcohen$x,0.165)-quantile(winsubcohen$x,0.165)), (quantile(winsubcohen$x,0.495)-quantile(winsubcohen$x,0.165)), (quantile(winsubcohen$x,0.835)-quantile(winsubcohen$x,0.165)),
                             (quantile(winsubcohen$x,0.495)-quantile(winsubcohen$x,0.165)), (quantile(winsubcohen$x,0.495)-quantile(winsubcohen$x,0.495)), (quantile(winsubcohen$x,0.835)-quantile(winsubcohen$x,0.495)),
                             (quantile(winsubcohen$x,0.835)-quantile(winsubcohen$x,0.165)), (quantile(winsubcohen$x,0.835)-quantile(winsubcohen$x,0.495)), (quantile(winsubcohen$x,0.835)-quantile(winsubcohen$x,0.835))), nNonBreeding, nNonBreeding)
    
    
    x1 <- c(quantile(brdsubcohen$x,0),  quantile(brdsubcohen$x,0),  quantile(brdsubcohen$x,0.33), quantile(brdsubcohen$x,0.33), quantile(brdsubcohen$x,0))
    y <- c(quantile(brdsubcohen$y,0),  quantile(brdsubcohen$y,1),  quantile(brdsubcohen$y,1), quantile(brdsubcohen$y,0), quantile(brdsubcohen$y,0))
    xy1 <- cbind(x1, y)
    
    x2 <- c(quantile(brdsubcohen$x,0.33),  quantile(brdsubcohen$x,0.33),  quantile(brdsubcohen$x,0.66), quantile(brdsubcohen$x,0.66), quantile(brdsubcohen$x,0.33))
    xy2 <- cbind(x2, y)
    
    x3 <- c(quantile(brdsubcohen$x,0.66),  quantile(brdsubcohen$x,0.66),  quantile(brdsubcohen$x,1), quantile(brdsubcohen$x,1), quantile(brdsubcohen$x,0.66))
    xy3 <- cbind(x3, y)
    
    p1 <- Polygon(xy1)
    p2 <- Polygon(xy2)
    p3 <- Polygon(xy3)
    ps1 <- Polygons(list(p1), ID = 1)
    ps2 <- Polygons(list(p2), ID = 2)
    ps3 <- Polygons(list(p3), ID = 3)
    
    originsites <- SpatialPolygons(list(ps1,ps2,ps3))
    
    x1 <- c(quantile(winsubcohen$x,0),  quantile(winsubcohen$x,0),  quantile(winsubcohen$x,0.33), quantile(winsubcohen$x,0.33), quantile(winsubcohen$x,0))
    y <- c(quantile(winsubcohen$y,0),  quantile(winsubcohen$y,1),  quantile(winsubcohen$y,1), quantile(winsubcohen$y,0), quantile(winsubcohen$y,0))
    xy1 <- cbind(x1, y)
    
    x2 <- c(quantile(winsubcohen$x,0.33),  quantile(winsubcohen$x,0.33),  quantile(winsubcohen$x,0.66), quantile(winsubcohen$x,0.66), quantile(winsubcohen$x,0.33))
    xy2 <- cbind(x2, y)
    
    x3 <- c(quantile(winsubcohen$x,0.66),  quantile(winsubcohen$x,0.66),  quantile(winsubcohen$x,1), quantile(winsubcohen$x,1), quantile(winsubcohen$x,0.66))
    xy3 <- cbind(x3, y)
    
    p1 <- Polygon(xy1)
    p2 <- Polygon(xy2)
    p3 <- Polygon(xy3)
    ps1 <- Polygons(list(p1), ID = 1)
    ps2 <- Polygons(list(p2), ID = 2)
    ps3 <- Polygons(list(p3), ID = 3) 
    
    targetsites <- SpatialPolygons(list(ps1,ps2,ps3))
    
    proj4string(originsites) <- CRS("+proj=eqc +lat_ts=0 +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +a=6371007 +b=6371007 +units=m +no_defs ")
    raster::crs(targetsites)<-raster::crs(originsites)
    
    spatbreeding <- SpatialPoints(brdsubcohen[,1:2])
    raster::crs(spatbreeding)<-raster::crs(originsites)
    spatwinter <- SpatialPoints(winsubcohen[,1:2])
    raster::crs(spatwinter)<-raster::crs(originsites)
    
    GPS_mc<-estMC(isGL=FALSE, # Logical vector: light-level geolocator(T)/GPS(F)
                  targetDist = nonBreedDist, # targetSites distance matrix
                  originDist = breedDist, # originSites distance matrix
                  originRelAbund = BBsrelN, #Origin relative abund.
                  targetSites = targetsites, # Non-breeding target sites
                  originSites = originsites, # Breeding origin sites
                  originPoints = spatbreeding, # Capture Locations
                  targetPoints = spatwinter, # Device target locations
                  verbose = 0,   # output options
                  nSamples = 1000) # This is set low for example
    
    data_sheet[f,2] <- GPS_mc$meanMC
    data_sheet[f,4] <- Spread[t]
    data_sheet[f,5] <- Area[t]
    f <- f+1
    q <- q+1
    
  }
  data_sheet_list[[g]] <- data_sheet  
  q <- 5
  f <- 8
}


# 100 replicates were completed and the data compiled into a single .csv file available in the file depository


# Section 4: Sample size analysis ----

# This section was iterated 100 times on a high performance cluster to produce 100 replicates of the sampling procedure. 
# Mantel correlations are highly memory intensive due to the exponential increase in matrix size with larger populations.
# As a result, a single iteration may take several hours to complete.
# Without mantel correlations this section can be run quickly

rm(list = ls())
setwd("") # Set working directory to locartion of files

# Read in breeding and non-breeding locations created in section 1
brd_MC1 <- read.csv("brd_MC1.csv")
win_MC1 <- read.csv("win_MC1.csv")
brd_MC2 <- read.csv("brd_MC2.csv")
win_MC2 <- read.csv("win_MC2.csv")
brd_MC3 <- read.csv("brd_MC3.csv")
win_MC3 <- read.csv("win_MC3.csv")

brd <- list(brd_MC1, brd_MC2, brd_MC3)
win <- list(win_MC1, win_MC2, win_MC3)

data_sheet <- data.frame(matrix(nrow=6,ncol=6))
names(data_sheet) <- c("Mantel_MC1", "Mantel_MC2", "Mantel_MC3",
                       "Cohen_MC1", "Cohen_MC2", "Cohen_MC3")

sample_size <- c(10, 50, 100,1000,2500,5000)

nBreeding <- 3 # number of breeding regions
nNonBreeding <- 3 # number of non-breeding regions

for (a in 1:3){
  
  brd_locs <- brd[[a]]
  win_locs <- win[[a]]
  
  for (i in 1:6){
    inds <- seq(1:10000) #pick individuals in the study area using over() from sp package
    samp_inds <- inds[sample(1:10000,sample_size[i],replace=F)] #randomly sample the pop        
    try(coordinates(win_locs)<-c("lon","lat"), silent=T) 
    
    # Mantel
    
    samp_brd <- data.frame(brd_locs[samp_inds,])
    samp_win <- data.frame(win_locs[samp_inds,])
    brddistssubset <- dist(samp_brd[1:2]) 
    windistssubset <- dist(samp_win[1:2])
    mantelsubset <-  mantel.rtest(brddistssubset,windistssubset,nrepet=1)
    data_sheet[i,a] <- mantelsubset$obs
    
    # Cohen
    
    brdsubcohen <- data.frame(brd_locs[samp_inds,])
    winsubcohen <- data.frame(win_locs[samp_inds,])
    names(brdsubcohen) <- c('x','y','zone1')
    names(winsubcohen) <- c('x','y','zone2')
    
    brdsubcohen$zone1 <- 1
    brdsubcohen$zone1[brdsubcohen$x>quantile(brdsubcohen$x,0.33)] <- 2
    brdsubcohen$zone1[brdsubcohen$x>quantile(brdsubcohen$x,0.67)] <- 3
    
    winsubcohen$zone2 <- 1
    winsubcohen$zone2[winsubcohen$x>quantile(winsubcohen$x,0.33)] <- 2
    winsubcohen$zone2[winsubcohen$x>quantile(winsubcohen$x,0.67)] <- 3
    
    brdwinsubcohen <- cbind(brdsubcohen,winsubcohen)
    brdwinsubcohen$zone1 <- as.factor(brdwinsubcohen$zone1)
    brdwinsubcohen$zone2 <- as.factor(brdwinsubcohen$zone2)
    bwscz1 <- subset(brdwinsubcohen, zone1 == 1)
    bwscz2 <- subset(brdwinsubcohen, zone1 == 2)
    bwscz3 <- subset(brdwinsubcohen, zone1 == 3)
    bwscz1summary <- summary(bwscz1$zone2)
    bwscz1summary <- as.data.frame(bwscz1summary)
    stp1 <- bwscz1summary[1,]/sum(summary(bwscz1$zone2))
    stp2 <- bwscz1summary[2,]/sum(summary(bwscz1$zone2))
    stp3 <- bwscz1summary[3,]/sum(summary(bwscz1$zone2))
    bwscz2summary <- summary(bwscz2$zone2)
    bwscz2summary <- as.data.frame(bwscz2summary)
    stp4 <- bwscz2summary[1,]/sum(summary(bwscz2$zone2))
    stp5 <- bwscz2summary[2,]/sum(summary(bwscz2$zone2))
    stp6 <- bwscz2summary[3,]/sum(summary(bwscz2$zone2))
    bwscz3summary <- summary(bwscz3$zone2)
    bwscz3summary <- as.data.frame(bwscz3summary)
    stp7 <- bwscz3summary[1,]/sum(summary(bwscz3$zone2))
    stp8 <- bwscz3summary[2,]/sum(summary(bwscz3$zone2))
    stp9 <- bwscz3summary[3,]/sum(summary(bwscz3$zone2))
    
    #transition probabilities form each breeding to each non-breeding region
    spsi <- matrix(c(stp1,stp2,stp3,stp4,stp5,stp6,stp7,stp8,stp9), nBreeding, nNonBreeding, byrow = TRUE)
    spsi[is.na(spsi)] <- 0
    
    # Define the relative abundance within the three breeding regions
    srelN <- c(sum(summary(bwscz1$zone2)/sample_size[i]),sum(summary(bwscz2$zone2)/sample_size[i]),sum(summary(bwscz3$zone2)/sample_size[i]))
    
    
    
    breedDist <- matrix(c((quantile(brdsubcohen$x,0.165)-quantile(brdsubcohen$x,0.165)), (quantile(brdsubcohen$x,0.495)-quantile(brdsubcohen$x,0.165)), (quantile(brdsubcohen$x,0.835)-quantile(brdsubcohen$x,0.165)),
                          (quantile(brdsubcohen$x,0.495)-quantile(brdsubcohen$x,0.165)), (quantile(brdsubcohen$x,0.495)-quantile(brdsubcohen$x,0.495)), (quantile(brdsubcohen$x,0.835)-quantile(brdsubcohen$x,0.495)),
                          (quantile(brdsubcohen$x,0.835)-quantile(brdsubcohen$x,0.165)), (quantile(brdsubcohen$x,0.835)-quantile(brdsubcohen$x,0.495)), (quantile(brdsubcohen$x,0.835)-quantile(brdsubcohen$x,0.835))), nBreeding, nBreeding)
    nonBreedDist <- matrix(c((quantile(winsubcohen$x,0.165)-quantile(winsubcohen$x,0.165)), (quantile(winsubcohen$x,0.495)-quantile(winsubcohen$x,0.165)), (quantile(winsubcohen$x,0.835)-quantile(winsubcohen$x,0.165)),
                             (quantile(winsubcohen$x,0.495)-quantile(winsubcohen$x,0.165)), (quantile(winsubcohen$x,0.495)-quantile(winsubcohen$x,0.495)), (quantile(winsubcohen$x,0.835)-quantile(winsubcohen$x,0.495)),
                             (quantile(winsubcohen$x,0.835)-quantile(winsubcohen$x,0.165)), (quantile(winsubcohen$x,0.835)-quantile(winsubcohen$x,0.495)), (quantile(winsubcohen$x,0.835)-quantile(winsubcohen$x,0.835))), nNonBreeding, nNonBreeding)
    
    
    
    x1 <- c(quantile(brdsubcohen$x,0),  quantile(brdsubcohen$x,0),  quantile(brdsubcohen$x,0.33), quantile(brdsubcohen$x,0.33), quantile(brdsubcohen$x,0))
    y <- c(quantile(brdsubcohen$y,0),  quantile(brdsubcohen$y,1),  quantile(brdsubcohen$y,1), quantile(brdsubcohen$y,0), quantile(brdsubcohen$y,0))
    xy1 <- cbind(x1, y)
    
    x2 <- c(quantile(brdsubcohen$x,0.33),  quantile(brdsubcohen$x,0.33),  quantile(brdsubcohen$x,0.66), quantile(brdsubcohen$x,0.66), quantile(brdsubcohen$x,0.33))
    xy2 <- cbind(x2, y)
    
    x3 <- c(quantile(brdsubcohen$x,0.66),  quantile(brdsubcohen$x,0.66),  quantile(brdsubcohen$x,1), quantile(brdsubcohen$x,1), quantile(brdsubcohen$x,0.66))
    xy3 <- cbind(x3, y)
    
    p1 <- Polygon(xy1)
    p2 <- Polygon(xy2)
    p3 <- Polygon(xy3)
    ps1 <- Polygons(list(p1), ID = 1)
    ps2 <- Polygons(list(p2), ID = 2)
    ps3 <- Polygons(list(p3), ID = 3)
    
    originsites <- SpatialPolygons(list(ps1,ps2,ps3))
    
    x1 <- c(quantile(winsubcohen$x,0),  quantile(winsubcohen$x,0),  quantile(winsubcohen$x,0.33), quantile(winsubcohen$x,0.33), quantile(winsubcohen$x,0))
    y <- c(quantile(winsubcohen$y,0),  quantile(winsubcohen$y,1),  quantile(winsubcohen$y,1), quantile(winsubcohen$y,0), quantile(winsubcohen$y,0))
    xy1 <- cbind(x1, y)
    
    x2 <- c(quantile(winsubcohen$x,0.33),  quantile(winsubcohen$x,0.33),  quantile(winsubcohen$x,0.66), quantile(winsubcohen$x,0.66), quantile(winsubcohen$x,0.33))
    xy2 <- cbind(x2, y)
    
    x3 <- c(quantile(winsubcohen$x,0.66),  quantile(winsubcohen$x,0.66),  quantile(winsubcohen$x,1), quantile(winsubcohen$x,1), quantile(winsubcohen$x,0.66))
    xy3 <- cbind(x3, y)
    
    p1 <- Polygon(xy1)
    p2 <- Polygon(xy2)
    p3 <- Polygon(xy3)
    ps1 <- Polygons(list(p1), ID = 1)
    ps2 <- Polygons(list(p2), ID = 2)
    ps3 <- Polygons(list(p3), ID = 3) 
    
    targetsites <- SpatialPolygons(list(ps1,ps2,ps3))
    
    proj4string(originsites) <- CRS("+proj=eqc +lat_ts=0 +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +a=6371007 +b=6371007 +units=m +no_defs ")
    #proj4string(targetsites) <- CRS("+proj=eqc +lat_ts=0 +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +a=6371007 +b=6371007 +units=m +no_defs ")
    raster::crs(targetsites)<-raster::crs(originsites)
    
    spatbreeding <- SpatialPoints(brdsubcohen[,1:2])
    #proj4string(spatbreeding) <- CRS("+proj=eqc +lat_ts=0 +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +a=6371007 +b=6371007 +units=m +no_defs ")
    raster::crs(spatbreeding)<-raster::crs(originsites)
    spatwinter <- SpatialPoints(winsubcohen[,1:2])
    #proj4string(spatwinter) <- CRS("+proj=eqc +lat_ts=0 +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +a=6371007 +b=6371007 +units=m +no_defs ")
    raster::crs(spatwinter)<-raster::crs(originsites)
    
    GPS_mc<-estMC(isGL=FALSE, # Logical vector: light-level geolocator(T)/GPS(F)
                  targetDist = nonBreedDist, # targetSites distance matrix
                  originDist = breedDist, # originSites distance matrix
                  originRelAbund = srelN, #Origin relative abund.
                  targetSites = targetsites, # Non-breeding target sites
                  originSites = originsites, # Breeding origin sites
                  originPoints = spatbreeding, # Capture Locations
                  targetPoints = spatwinter, # Device target locations
                  verbose = 0,   # output options
                  nSamples = 1000) # This is set low for example
    
    data_sheet[i,(a+3)] <- GPS_mc$meanMC
    
  }
  data_sheet$Sampled <- sample_size
}

# 100 replicates were completed and the data compiled into a single .csv file available in the file depository

# Section 5: Applying sampling designs (patchy breeding and non-breeding popualtions) ----

# This section was iterated 100 times on a high performance cluster. 
# Mantel correaltions are highly memory intensive due to the exponential increase in matrix size with larger populations.
# As a result, a single iteration may take several hours to complete.

rm(list = ls())

#### Set up data and objects for a for loop #

n2 <- 10000 # population size

bigboxsize <- c(.25,.5,.75)

# Create empty lists, data frames, and set object starting numbers. All needs to be done before each run of the loop
winloclist <- list() 
count <- 0 # resets count object back to 0
pls <- list()
counter <- 0 # resets counter object back to 0
ps <- list()
area <- list()
slp <- list()
q <- 5 # Counter used to fill data_sheet correctly 
f <- 8 # Counter used to fill data_sheet correctly 

# Create empty matrix for data
data_sheet <- data.frame(matrix(nrow=7,ncol=6)) # produces a data frame for the calculated mantel stats to be stored in. 4= number of variations in sampling spread/size/number
names(data_sheet) <- c("Mantel","Cohen","Samp_Size","Spread", "Area", "Scenario")
data_sheet$Scenario <- c('All', 'Allsmall', 'Allmedium', 'Alllarge','Area','Area','Area')
data_sheet_list <- list()

# Spread and area scenarios to allow descriptor filling in data_sheet
Spread <- c("Low","Medium", "High")
Area <- c("Small","Medium", "Large")

# Setting up for Cohen
nBreeding <- 4 #number of breeding regions
nNonBreeding <- 4 #number of non-breeding regions

# Set up extent of the boxes that delimit the spread of sampling.
# big box 1 (xmin, xmax, ymin, ymax)
bb1 <-  c(0.0825, 0.2475, 1.835, 1.945)
# big box 2 (xmin, xmax, ymin, ymax)
bb2 <-  c(0.19155, 0.80845, 1.463073, 1.870227)
# big box 3 (xmin, xmax, ymin, ymax)
bb3 <-  c(0.067175, 0.932825, 1.380986, 1.952315)
box1 <- list(bb1,bb2,bb3)

bb1 <-  c(0.315925, 0.684075, 1.545161, 1.78814)
# big box 2 (xmin, xmax, ymin, ymax)
bb2 <-  c(0.19155, 0.80845, 1.463073, 1.870227)
# big box 3 (xmin, xmax, ymin, ymax)
bb3 <-  c(0.067175, 0.932825, 1.380986, 1.952315)
box2 <- list(bb1,bb2,bb3)

bb1 <-  c(0.315925, 0.684075, 1.545161, 1.78814)
# big box 2 (xmin, xmax, ymin, ymax)
bb2 <-  c(0.19155, 0.80845, 1.463073, 1.870227)
# big box 3 (xmin, xmax, ymin, ymax)
bb3 <-  c(0.067175, 0.932825, 1.380986, 1.952315)
box3 <- list(bb1,bb2,bb3)

bb1 <-  c(0.315925, 0.684075, 1.545161, 1.78814)
# big box 2 (xmin, xmax, ymin, ymax)
bb2 <-  c(0.19155, 0.80845, 1.463073, 1.870227)
# big box 3 (xmin, xmax, ymin, ymax)
bb3 <-  c(0.067175, 0.932825, 1.380986, 1.952315)
box4 <- list(bb1,bb2,bb3)

## The for loop ##

# Set working directory to folder containing patchy population breeding and non-breeding location files

brd_MC1 <- read.csv("brd_MC1.csv")
win_MC1 <- read.csv("win_MC1.csv")
brd_MC2 <- read.csv("brd_MC2.csv")
win_MC2 <- read.csv("win_MC2.csv")
brd_MC3 <- read.csv("brd_MC3.csv")
win_MC3 <- read.csv("win_MC3.csv")

brd <- list(brd_MC1, brd_MC2, brd_MC3)
win <- list(win_MC1, win_MC2, win_MC3)

  
for(g in 1:3){ # loops thorugh levels of migratory connectivity
  
  brd_locs <- brd[[g]]
  win_locs <- win[[g]]  
  
  #### Mantel Methodology for entire population #### - computationally intensive step -
  brddists <- dist(data.frame(brd_locs[1:2]))  
  windists <- dist(data.frame(win_locs[1:2]))

  dfbrd <- melt(as.matrix(brddists), varnames = c("row", "col"))
  dfwin <- melt(as.matrix(windists), varnames = c("row", "col"))
  data_sheet[1,1] <-  mantel.rtest(brddists,windists,nrepet=1)
  
  #### Cohen Methodology for entire population ####
  brdcohen <- data.frame(brd_locs)
  wincohen <- data.frame(win_locs)
  names(brdcohen) <- c('x','y','zone1')
  names(wincohen) <- c('x','y','zone2')
  
  brdwincohen <- cbind(brdcohen,wincohen)
  brdwincohen$zone1 <- as.factor(brdwincohen$zone1)
  brdwincohen$zone2 <- as.factor(brdwincohen$zone2)
  bwcz1 <- subset(brdwincohen, zone1 == 1)
  bwcz2 <- subset(brdwincohen, zone1 == 2)
  bwcz3 <- subset(brdwincohen, zone1 == 3)
  bwcz4 <- subset(brdwincohen, zone1 == 4)
  
  # Define the relative abundance within the three breeding regions
  relN <- c(sum(summary(bwcz1$zone2)/10000),sum(summary(bwcz2$zone2)/10000),sum(summary(bwcz3$zone2)/10000), sum(summary(bwcz4$zone2)/10000))
  
  # set distance matrices for cohen
  breedDist <- matrix(c(0,.66,.45,.8114801,
                        .66,0,.8114801,.45,
                        .45,.8114801,0,.66,
                        .8114801,.45,.66,0), nBreeding, nBreeding)
  nonBreedDist <- matrix(c(0,.66,.45,.8114801,
                           .66,0,.8114801,.45,
                           .45,.8114801,0,.66,
                           .8114801,.45,.66,0), nBreeding, nBreeding)
  
  # breeding polygons
  
  x1 <- c(0,0,.34,.34,0)
  y1 <- c(1.77,2,2,1.77,1.77)
  xy1 <- cbind(x1, y1)
  
  x2 <- c(.65,.65,1,1,.65)
  y2 <- c(1.77,2,2,1.77,1.77)
  xy2 <- cbind(x2, y2)
  
  x3 <- c(0,0,.34,.34,0)
  y3 <- c(1.33,1.56,1.56,1.33,1.33)
  xy3 <- cbind(x3, y3)
  
  x4 <- c(.65,.65,1,1,.65)
  y4 <- c(1.33,1.56,1.56,1.33,1.33)
  xy4 <- cbind(x4, y4)
  
  p1 <- Polygon(xy1)
  p2 <- Polygon(xy2)
  p3 <- Polygon(xy3)
  p4 <- Polygon(xy4)
  
  ps1 <- Polygons(list(p1), ID = 1)
  ps2 <- Polygons(list(p2), ID = 2)
  ps3 <- Polygons(list(p3), ID = 3)
  ps4 <- Polygons(list(p4), ID = 4)
  
  originsites <- SpatialPolygons(list(ps1,ps2,ps3,ps4))
  
  # target sites
  
  x1 <- c(0,0,.34,.34,0)
  y1 <- c(.44,.67,.67,.44,.44)
  xy1 <- cbind(x1, y1)
  
  x2 <- c(.65,.65,1,1,.65)
  y2 <- c(.44,.67,.67,.44,.44)
  xy2 <- cbind(x2, y2)
  
  x3 <- c(0,0,.34,.34,0)
  y3 <- c(0,.23,.23,0,0)
  xy3 <- cbind(x3, y3)
  
  x4 <- c(.65,.65,1,1,.65)
  y4 <- c(0,.23,.23,0,0)
  xy4 <- cbind(x4, y4)
  
  p1 <- Polygon(xy1)
  p2 <- Polygon(xy2)
  p3 <- Polygon(xy3)
  p4 <- Polygon(xy4)
  
  ps1 <- Polygons(list(p1), ID = 1)
  ps2 <- Polygons(list(p2), ID = 2)
  ps3 <- Polygons(list(p3), ID = 3)
  ps4 <- Polygons(list(p4), ID = 4)
  
  targetsites <- SpatialPolygons(list(ps1,ps2,ps3,ps4))
  
  proj4string(originsites) <- CRS("+proj=eqc +lat_ts=0 +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +a=6371007 +b=6371007 +units=m +no_defs ")
  raster::crs(targetsites)<-raster::crs(originsites)
  
  spatbreeding <- SpatialPoints(brdcohen[,1:2])
  raster::crs(spatbreeding)<-raster::crs(originsites)
  spatwinter <- SpatialPoints(wincohen[,1:2])
  raster::crs(spatwinter)<-raster::crs(originsites)
  
  GPS_mc<-estMC(isGL=FALSE, # Logical vector: light-level geolocator(T)/GPS(F)
                targetDist = nonBreedDist, # targetSites distance matrix
                originDist = breedDist, # originSites distance matrix
                originRelAbund = relN, #Origin relative abund.
                targetSites = targetsites, # Non-breeding target sites
                originSites = originsites, # Breeding origin sites
                originPoints = spatbreeding, # Capture Locations
                targetPoints = spatwinter, # Device target locations
                verbose = 0,   # output options
                nSamples = 1000) # This is set low for example
  
  data_sheet[1,2] <- GPS_mc$meanMC
  
  data_sheet[1,3] <- n2
  
  ####################################################
  # Sampling the population --------------
  
  # Delimit sampling areas (breeding range):
  
  for(t in 1:3){ # loop through spread scenarios
    
    # Create Big box - whole study area (no sampling outside this)
    bigbox <- bigboxsize[t] # sensible range 0.5 to 1.5. This would be added and small box fixed to allow big box to change but fix small box size.
    bigboxdims <- c(bigbox*.33,bigbox*.22)
    
    midbrd1 <- c(.165,1.89)#midpoint of breeding range
    midwin1 <- c(.165,.55)
    area1 <- c(midbrd1[1]-(0.5*bigboxdims[1]),midbrd1[1]+(0.5*bigboxdims[1]),
               midbrd1[2]-(0.5*bigboxdims[2]),midbrd1[2]+(0.5*bigboxdims[2]))
    xcs <- c(area1[1],area1[1],area1[2],area1[2],area1[1])
    ycs <- c(area1[3],area1[4],area1[4],area1[3],area1[3])
    xym <- data.frame(x=xcs,y=ycs)
    p1 <- Polygon(xym)
    pls[[1]]<-p1
    
    midbrd2 <- c(.825, 1.89)#midpoint of breeding range
    midwin2 <- c(.825,.55)
    area2 <- c(midbrd2[1]-(0.5*bigboxdims[1]),midbrd2[1]+(0.5*bigboxdims[1]),
               midbrd2[2]-(0.5*bigboxdims[2]),midbrd2[2]+(0.5*bigboxdims[2]))
    xcs <- c(area2[1],area2[1],area2[2],area2[2],area2[1])
    ycs <- c(area2[3],area2[4],area2[4],area2[3],area2[3])
    xym <- data.frame(x=xcs,y=ycs)
    p2 <- Polygon(xym)
    pls[[2]]<-p2
    
    midbrd3 <- c(.165,1.44)#midpoint of breeding range
    midwin3 <- c(.165,.11)
    area3 <- c(midbrd3[1]-(0.5*bigboxdims[1]),midbrd3[1]+(0.5*bigboxdims[1]),
               midbrd3[2]-(0.5*bigboxdims[2]),midbrd3[2]+(0.5*bigboxdims[2]))
    xcs <- c(area3[1],area3[1],area3[2],area3[2],area3[1])
    ycs <- c(area3[3],area3[4],area3[4],area3[3],area3[3])
    xym <- data.frame(x=xcs,y=ycs)
    p3 <- Polygon(xym)
    pls[[3]]<-p3
    
    midbrd4 <- c(.825,1.44)#midpoint of breeding range
    midwin4 <- c(.825,.11)
    area4 <- c(midbrd4[1]-(0.5*bigboxdims[1]),midbrd4[1]+(0.5*bigboxdims[1]),
               midbrd4[2]-(0.5*bigboxdims[2]),midbrd4[2]+(0.5*bigboxdims[2]))
    xcs <- c(area4[1],area4[1],area4[2],area4[2],area4[1])
    ycs <- c(area4[3],area4[4],area4[4],area4[3],area4[3])
    xym <- data.frame(x=xcs,y=ycs)
    p4 <- Polygon(xym)
    pls[[4]]<-p4
    
    
    slp <- pls
    ps <- Polygons(slp,1)
    area <- SpatialPolygons(list(ps))
    pls <- list()
    
    
    # loop through sample sizes (now fixed at 200)
    sampleN <- 200
    
    ##### Sample populations and calculate mantel statistics
    fixN <- T 
    try(coordinates(brd_locs)<-c("x","y"), silent=T)
    
    # individuals for spread scenario ---------------------------------
    inds <- as.numeric(which(over(brd_locs,area)==1)) #pick individuals in the study area using over() from sp package
    samp_n<- ifelse(fixN==T,min(sampleN,length(inds)),floor(sampleprop*length(inds))) #number of individuals to sample
    samp_inds <- inds[sample(1:length(inds),samp_n,replace=F)] #randomly sample the pop
    try(coordinates(win_locs)<-c("lon","lat"), silent=T) 
    samp_brd <- data.frame(brd_locs[samp_inds,])
    samp_win <- data.frame(win_locs[samp_inds,])
    brddistssubset <- dist(samp_brd[1:2]) # 
    windistssubset <- dist(samp_win[1:2])
    mantelsubset <-  mantel.rtest(brddistssubset,windistssubset,nrepet=1)
    data_sheet[q,1] <- mantelsubset$obs #extract mantel statistic - needs updating if scenarios run in a loop
    data_sheet[q,3] <- length(samp_inds)
    
    #### Cohen methodology for sampled inds ####
    brdsubcohen <- data.frame(brd_locs[samp_inds,])
    winsubcohen <- data.frame(win_locs[samp_inds,])
    names(brdsubcohen) <- c('x','y','zone1')
    names(winsubcohen) <- c('x','y','zone2')
    
    brdwinsubcohen <- cbind(brdsubcohen,winsubcohen)
    brdwinsubcohen$zone1 <- as.factor(brdwinsubcohen$zone1)
    brdwinsubcohen$zone2 <- as.factor(brdwinsubcohen$zone2)
    bwscz1 <- subset(brdwinsubcohen, zone1 == 1)
    bwscz2 <- subset(brdwinsubcohen, zone1 == 2)
    bwscz3 <- subset(brdwinsubcohen, zone1 == 3)
    bwscz4 <- subset(brdwinsubcohen, zone1 == 4)
    
    # Define the relative abundance within the three breeding regions
    srelN <- c(sum(summary(bwscz1$zone2)/samp_n),sum(summary(bwscz2$zone2)/samp_n),sum(summary(bwscz3$zone2)/samp_n), sum(summary(bwscz4$zone2)/samp_n))
    
    # set distance matrices for cohen
    breedDist <- matrix(c(0,.66,.45,.8114801,
                          .66,0,.8114801,.45,
                          .45,.8114801,0,.66,
                          .8114801,.45,.66,0), nBreeding, nBreeding)
    nonBreedDist <- matrix(c(0,.66,.45,.8114801,
                             .66,0,.8114801,.45,
                             .45,.8114801,0,.66,
                             .8114801,.45,.66,0), nBreeding, nBreeding)
    
    
    ps1 <- Polygons(list(p1), ID = 1)
    ps2 <- Polygons(list(p2), ID = 2)
    ps3 <- Polygons(list(p3), ID = 3)
    ps4 <- Polygons(list(p4), ID = 4)
    
    originsites <- SpatialPolygons(list(ps1,ps2,ps3,ps4))
    
    
    x1 <- c(0,0,.34,.34,0)
    y1 <- c(.44,.67,.67,.44,.44)
    xy1 <- cbind(x1, y1)
    
    x2 <- c(.65,.65,1,1,.65)
    y2 <- c(.44,.67,.67,.44,.44)
    xy2 <- cbind(x2, y2)
    
    x3 <- c(0,0,.34,.34,0)
    y3 <- c(0,.23,.23,0,0)
    xy3 <- cbind(x3, y3)
    
    x4 <- c(.65,.65,1,1,.65)
    y4 <- c(0,.23,.23,0,0)
    xy4 <- cbind(x4, y4)
    
    p1 <- Polygon(xy1)
    p2 <- Polygon(xy2)
    p3 <- Polygon(xy3)
    p4 <- Polygon(xy4)
    
    ps1 <- Polygons(list(p1), ID = 1)
    ps2 <- Polygons(list(p2), ID = 2)
    ps3 <- Polygons(list(p3), ID = 3)
    ps4 <- Polygons(list(p4), ID = 4)
    
    targetsites <- SpatialPolygons(list(ps1,ps2,ps3,ps4))
    
    proj4string(originsites) <- CRS("+proj=eqc +lat_ts=0 +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +a=6371007 +b=6371007 +units=m +no_defs ")
    raster::crs(targetsites)<-raster::crs(originsites)
    
    spatbreeding <- SpatialPoints(brdsubcohen[,1:2])
    raster::crs(spatbreeding)<-raster::crs(originsites)
    spatwinter <- SpatialPoints(winsubcohen[,1:2])
    raster::crs(spatwinter)<-raster::crs(originsites)
    
    GPS_mc<-estMC(isGL=FALSE, # Logical vector: light-level geolocator(T)/GPS(F)
                  targetDist = nonBreedDist, # targetSites distance matrix
                  originDist = breedDist, # originSites distance matrix
                  originRelAbund = srelN, #Origin relative abund.
                  targetSites = targetsites, # Non-breeding target sites
                  originSites = originsites, # Breeding origin sites
                  originPoints = spatbreeding, # Capture Locations
                  targetPoints = spatwinter, # Device target locations
                  verbose = 0,   # output options
                  nSamples = 1000) # This is set low for example
    
    data_sheet[q,2] <- GPS_mc$meanMC
    data_sheet[q,5] <- Area[t]
    
    # All indivudals within area scenario
    
    inds <- as.numeric(which(over(brd_locs,area)==1)) #pick individuals in the study area using over() from sp package
    try(coordinates(win_locs)<-c("lon","lat"), silent=T) 
    samp_brd <- data.frame(brd_locs[inds,])
    samp_win <- data.frame(win_locs[inds,])
    brddistssubset <- dist(samp_brd[1:2]) # 
    windistssubset <- dist(samp_win[1:2])
    mantelsubset <-  mantel.rtest(brddistssubset,windistssubset,nrepet=1)
    data_sheet[t+1,1] <- mantelsubset$obs #extract mantel statistic - needs updating if scenarios run in a loop
    data_sheet[t+1,3] <- length(inds)
    
    #### Cohen methodology for sampled inds ####
    brdsubcohen <- data.frame(brd_locs[inds,])
    winsubcohen <- data.frame(win_locs[inds,])
    names(brdsubcohen) <- c('x','y','zone1')
    names(winsubcohen) <- c('x','y','zone2')
    
    brdwinsubcohen <- cbind(brdsubcohen,winsubcohen)
    brdwinsubcohen$zone1 <- as.factor(brdwinsubcohen$zone1)
    brdwinsubcohen$zone2 <- as.factor(brdwinsubcohen$zone2)
    bwscz1 <- subset(brdwinsubcohen, zone1 == 1)
    bwscz2 <- subset(brdwinsubcohen, zone1 == 2)
    bwscz3 <- subset(brdwinsubcohen, zone1 == 3)
    bwscz4 <- subset(brdwinsubcohen, zone1 == 4)
    
    proj4string(originsites) <- CRS("+proj=eqc +lat_ts=0 +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +a=6371007 +b=6371007 +units=m +no_defs ")
    raster::crs(targetsites)<-raster::crs(originsites)
    
    spatbreeding <- SpatialPoints(brdsubcohen[,1:2])
    raster::crs(spatbreeding)<-raster::crs(originsites)
    spatwinter <- SpatialPoints(winsubcohen[,1:2])
    raster::crs(spatwinter)<-raster::crs(originsites)
    
    # Define the relative abundance within the three breeding regions
    srelN2 <- c(sum(summary(bwscz1$zone2)/length(inds)),sum(summary(bwscz2$zone2)/length(inds)),sum(summary(bwscz3$zone2)/length(inds)), sum(summary(bwscz4$zone2)/length(inds)))
    
    GPS_mc<-estMC(isGL=FALSE, # Logical vector: light-level geolocator(T)/GPS(F)
                  targetDist = nonBreedDist, # targetSites distance matrix
                  originDist = breedDist, # originSites distance matrix
                  originRelAbund = srelN2, #Origin relative abund.
                  targetSites = targetsites, # Non-breeding target sites
                  originSites = originsites, # Breeding origin sites
                  originPoints = spatbreeding, # Capture Locations
                  targetPoints = spatwinter, # Device target locations
                  verbose = 0,   # output options
                  nSamples = 1000)
    
    data_sheet[t+1,2] <- GPS_mc$meanMC
    data_sheet[t+1,5] <- Area[t]
    
    q <- q+1
    
  }
  data_sheet_list[[g]] <- data_sheet  
}

# 100 replicates were completed and the data compiled into a single .csv file available in the file depository

# Section 6: Post analysis data handling and figures ----

# Sampling processes with sections 3–5 are random and as such specific results will vary slightly.
# Specific output files (after compliling) are provided in the file repository for this section.

rm(list = ls())
gc()

## Run code chunks 1. and 2., after this each figure chunk can be ran ##

# 6.1. Spread and Area dataset production ----

All_data <- read.csv(file.choose()) # read in outputdata.csv file

All_data$Samp_Size <- as.numeric(All_data$Samp_Size)

Area <- subset(All_data, Scenario=="Area")
Area$Area <- factor(Area$Area, levels=c("Small", "Medium", "Large"))

Spread <- subset(All_data, Scenario=="Spread")
Spread$Spread <- factor(Spread$Spread, levels=c("Low", "Medium", "High"))

# Create summaries of mean Mantel scores
a <- summaryBy(Mantel ~ Area + MC + Scenario, Area, FUN = c(mean, sd))
a2 <- summaryBy(Mantel ~ Spread + MC + Scenario, Spread, FUN = c(mean, sd))

# Create summaries of mean Cohen scores
b <- summaryBy(Cohen ~ Area + MC + Scenario, Area, FUN = c(mean, sd))
b2 <- summaryBy(Cohen ~ Spread + MC + Scenario, Spread, FUN = c(mean, sd))

# Create data frame of full population Mantel and Cohen scores, Spread and Area vars are needed for plotting purposes
step1 <- subset(All_data, Samp_Size == 10000 & Class == 1)
fullscores <- rbind(step1,step1,step1)
fullscores$Area <- c('Small','Small','Small', 'Medium','Medium','Medium','Large','Large','Large')
fullscores$Spread <- c('Low', 'Low','Low', 'Medium','Medium','Medium','High','High','High')

# Create data frame of population at spatial extent of sampling Mantel and Cohen scores, Spread and Area vars are needed for plotting purposes

areascores <- subset(All_data, Samp_Size > 200 & Samp_Size < 10000 & Class == 1)

breeding <- read.csv(file.choose()) # brd_MC1 from the contiguous pop files folder
nonbreeding <- read.csv(file.choose()) # win_MC1 from the contiguous pop files folder

# Datasets without MC1 

ALT_All_data <- subset(All_data, MC != 1)
ALTArea <- subset(Area, MC != 1)
ALTSpread <- subset(Spread, MC != 1)
ALTa <- subset(a, MC != 1)
ALTa2 <- subset(a2, MC != 1)
ALTfullscores <- subset(fullscores, MC != 1)

# 6.2. Sample size dataset production ----

Data <- read.csv(file.choose()) # read in samplingdata.csv file

Meltdatabias <- reshape2::melt(Data[1:9], id='Sampled')

MantelMDB <- subset(Meltdatabias, variable %in% c("Mantel_MC1","Mantel_MC2","Mantel_MC3"))
CohenMDB <- subset(Meltdatabias, variable %in% c("Cohen_MC1","Cohen_MC2","Cohen_MC3"))

MDBa <- summaryBy(value ~ variable + Sampled, MantelMDB, FUN = c(mean, sd))
MDBb <- summaryBy(value ~ variable + Sampled, CohenMDB, FUN = c(mean, sd))

# 6.3 Figure 1 ----

rows <- sample(nrow(breeding), 1000)
breeding1 <- breeding[rows, ]
nonbreeding1<- nonbreeding[rows, ]

bdist <- dist(breeding1[1:2])
bsampdensity <- reshape2::melt(as.matrix(bdist), varnames = c("row", "col"))
bsampdensity$Season <- "Breeding"
nbdist <- dist(nonbreeding1[1:2])
nbsampdensity <- reshape2::melt(as.matrix(nbdist), varnames = c("row", "col"))
nbsampdensity$Season <- "Non-breeding"
alldensity <- rbind(bsampdensity, nbsampdensity) # Very compuationally intensive, creates large dataset
mu <- ddply(alldensity, c("Season"), summarise, Median=median(value))

f1breeding <- breeding[1:3]
rows <- sample(nrow(f1breeding), 250)
f1breeding <- f1breeding[rows, ]
f1breeding$rows <- rownames(f1breeding)
f1nonbreeding <- nonbreeding1
f1nonbreeding <- f1nonbreeding[rows, ]
f1nonbreeding$rows <- f1breeding$row

f1breeding$Season <- "Breeding"
f1nonbreeding$Season <- "Non-breeding"
colnames(f1nonbreeding) <- c('x','y', 'zone','rows','Season')
points <- rbind(f1breeding, f1nonbreeding)

x <- c(0,1)
y<- c(0, 2)
dat <- data.frame(x,y)
bb1 <-  c(0.315925, 0.684075, 1.545161, 1.78814)
x2 <- c(0.315925,0.315925,0.684075,0.684075)
y2 <- c(1.545161,1.78814,1.78814,1.545161)
datsmall <- data.frame(x2,y2)
bb2 <-  c(0.19155, 0.80845, 1.463073, 1.870227)
x3 <- c(0.19155,0.19155,0.80845,0.80845)
y3 <- c(1.463073,1.870227,1.870227,1.463073)
datmed <- data.frame(x3,y3)
bb3 <-  c(0.067175, 0.932825, 1.380986, 1.952315)
x4 <- c(0.067175,0.067175,0.932825,0.932825)
y4 <- c(1.380986,1.952315,1.952315,1.380986)
datlarge <- data.frame(x4,y4)

# Whole population 

p1 <- ggplot(dat, aes(x=x, y=y)) + 
  theme_void() +
  geom_line(data=points,aes(x=x, y=y, group=rows), colour="gray") +
  geom_point(data=points,aes(x=x, y=y, colour=Season)) +
  scale_fill_identity() +
  theme(legend.position = "none") +
  scale_color_manual(values=c("#E69F00", "#56B4E9"))

p2 <- ggplot(alldensity, aes(x=value, colour=Season)) + 
  geom_density(size=1) +
  scale_colour_manual(values=c("#E69F00", "#56B4E9"), name="Scenario") +
  geom_vline(data=mu, aes(xintercept=Median, colour=Season), linetype="dashed",show.legend=FALSE) +
  theme_pubr() +
  labs(x = "Pairwise distance", y = "Density") +
  theme(legend.position="none",axis.text = element_text(size= 18), axis.title = element_text(size=18), legend.text=element_text(size=18),legend.title=element_text(size=18)) +
  scale_x_continuous(
    breaks = seq(0, 2, 1),
    limits=c(0, 2)) +
  ylim(c(0,5))

# Small box size plot

smallbox <- Polygon(datsmall)
smallbox <- list(smallbox)
smallboxlist <- Polygons(smallbox,1)
sbp <- SpatialPolygons(list(smallboxlist))

cobreeding <- f1breeding
try(coordinates(cobreeding)<-c("x","y"), silent=T)
boxinds <- as.numeric(which(over(cobreeding,sbp)==1))
boxbreeding <- f1breeding[boxinds,]
boxnonbreeding <- f1nonbreeding[boxinds,]
boxpoints <- rbind(boxbreeding, boxnonbreeding)

bdist <- dist(boxbreeding[1:2])
bsampdensity <- reshape2::melt(as.matrix(bdist), varnames = c("row", "col"))
bsampdensity$Season <- "Breeding"
bsampdensity$Sample <- "Sampled"
nbdist <- dist(boxnonbreeding[1:2])
nbsampdensity <- reshape2::melt(as.matrix(nbdist), varnames = c("row", "col"))
nbsampdensity$Season <- "Non-breeding"
nbsampdensity$Sample <- "Sampled"

smalldensity <- rbind(bsampdensity, nbsampdensity)
smallmu <- ddply(smalldensity, c("Season"), summarise, Median=median(value))

p3 <- ggplot(dat, aes(x=x, y=y)) + 
  theme_void() +
  geom_point(data=points,aes(x=x, y=y, alpha=0.1), colour="lightgray") +
  geom_line(data=points,aes(x=x, y=y, group=rows, alpha=0)) +
  geom_line(data=boxpoints,aes(x=x, y=y, group=rows), colour="black", alpha=0.5) +
  geom_point(data=boxpoints,aes(x=x, y=y, colour=Season), size=1.5) +
  geom_polygon(data=datsmall, aes(x=x2,y=y2, fill=NA), colour="#E69F00") +
  theme(legend.position = "none") +
  scale_fill_identity() +
  scale_color_manual(values=c("#E69F00", "#56B4E9"))

p4 <- ggplot(smalldensity, aes(x=value, colour=Season)) + 
  geom_density(size=1) +
  scale_colour_manual(values=c("#E69F00", "#56B4E9"), name="Scenario") +
  geom_vline(data=smallmu, aes(xintercept=Median, colour=Season), linetype="dashed",show.legend=FALSE) +
  theme_pubr() +
  labs(x = "Pairwise distance", y = "Density") +
  theme(legend.position="none",axis.text = element_text(size= 18), axis.title = element_text(size=18), legend.text=element_text(size=18),legend.title=element_text(size=18)) +
  scale_x_continuous(
    breaks = seq(0, 2, 1),
    limits=c(0, 2)) +
  ylim(c(0,5))

# Medium box size

medbox <- Polygon(datmed)
medbox <- list(medbox)
medboxlist <- Polygons(medbox,1)
mbp <- SpatialPolygons(list(medboxlist))
cobreeding <- f1breeding
try(coordinates(cobreeding)<-c("x","y"), silent=T)
boxinds <- as.numeric(which(over(cobreeding,mbp)==1))
boxbreeding <- f1breeding[boxinds,]
boxnonbreeding <- f1nonbreeding[boxinds,]
boxpoints <- rbind(boxbreeding, boxnonbreeding)

bdist <- dist(boxbreeding[1:2])
bsampdensity <- reshape2::melt(as.matrix(bdist), varnames = c("row", "col"))
bsampdensity$Season <- "Breeding"
bsampdensity$Sample <- "Sampled"
nbdist <- dist(boxnonbreeding[1:2])
nbsampdensity <- reshape2::melt(as.matrix(nbdist), varnames = c("row", "col"))
nbsampdensity$Season <- "Non-breeding"
nbsampdensity$Sample <- "Sampled"

meddensity <- rbind(bsampdensity, nbsampdensity)
medmu <- ddply(meddensity, c("Season"), summarise, Median=median(value))

p5 <- ggplot(dat, aes(x=x, y=y)) + 
  theme_void() +
  geom_point(data=points,aes(x=x, y=y, alpha=0.1), colour="lightgray") +
  geom_line(data=points,aes(x=x, y=y, group=rows, alpha=0)) +
  geom_line(data=boxpoints,aes(x=x, y=y, group=rows), colour="black", alpha=0.5) +
  geom_point(data=boxpoints,aes(x=x, y=y, colour=Season), size=1.5) +
  geom_polygon(data=datmed, aes(x=x3,y=y3, fill=NA), colour="#E69F00") +
  theme(legend.position = "none") +
  scale_fill_identity() +
  scale_color_manual(values=c("#E69F00", "#56B4E9"))

p6 <- ggplot(meddensity, aes(x=value, colour=Season)) + 
  geom_density(size=1) +
  scale_colour_manual(values=c("#E69F00", "#56B4E9"), name="Scenario") +
  geom_vline(data=medmu, aes(xintercept=Median, colour=Season), linetype="dashed",show.legend=FALSE) +
  theme_pubr() +
  labs(x = "Pairwise distance", y = "Density") +
  theme(legend.position="none",axis.text = element_text(size= 18), axis.title = element_text(size=18), legend.text=element_text(size=18),legend.title=element_text(size=18)) +
  scale_x_continuous(
    breaks = seq(0, 2, 1),
    limits=c(0, 2)) +
  ylim(c(0,5))

# Large box size plot

larbox <- Polygon(datlarge)
larbox <- list(larbox)
larboxlist <- Polygons(larbox,1)
lbp <- SpatialPolygons(list(larboxlist))
cobreeding <- f1breeding
try(coordinates(cobreeding)<-c("x","y"), silent=T)
boxinds <- as.numeric(which(over(cobreeding,lbp)==1))
boxbreeding <- f1breeding[boxinds,]
boxnonbreeding <- f1nonbreeding[boxinds,]
boxpoints <- rbind(boxbreeding, boxnonbreeding)

bdist <- dist(boxbreeding[1:2])
bsampdensity <- reshape2::melt(as.matrix(bdist), varnames = c("row", "col"))
bsampdensity$Season <- "Breeding"
bsampdensity$Sample <- "Sampled"
nbdist <- dist(boxnonbreeding[1:2])
nbsampdensity <- reshape2::melt(as.matrix(nbdist), varnames = c("row", "col"))
nbsampdensity$Season <- "Non-breeding"
nbsampdensity$Sample <- "Sampled"

largedensity <- rbind(bsampdensity, nbsampdensity)
largemu <- ddply(largedensity, c("Season","Sample"), summarise, Median=median(value))

p7 <- ggplot(dat, aes(x=x, y=y)) + 
  theme_void() +
  geom_point(data=points,aes(x=x, y=y, alpha=0.1), colour="lightgray") +
  geom_line(data=points,aes(x=x, y=y, group=rows, alpha=0)) +
  geom_line(data=boxpoints,aes(x=x, y=y, group=rows), colour="black", alpha=0.5) +
  geom_point(data=boxpoints,aes(x=x, y=y, colour=Season), size=1.5) +
  geom_polygon(data=datlarge, aes(x=x4,y=y4, fill=NA), colour="#E69F00") +
  theme(legend.position = "none") +
  scale_fill_identity() +
  scale_color_manual(values=c("#E69F00", "#56B4E9"))

p8 <- ggplot(largedensity, aes(x=value, colour=Season)) + 
  geom_density(size=1, show.legend=FALSE) +
  scale_colour_manual(values=c("#E69F00", "#56B4E9"), name="Scenario") +
  geom_vline(data=largemu, aes(xintercept=Median, colour=Season), linetype="dashed",show.legend=FALSE) +
  theme_pubr() +
  labs(x = "Pairwise distance", y = "Density") +
  stat_density(aes(x=value, colour=Season),
               geom="line",position="identity", size = 0) + 
  guides(colour = guide_legend(override.aes=list(size=2))) +
  theme(legend.position="bottom",axis.text = element_text(size= 18), axis.title = element_text(size=18), legend.text=element_text(size=18),legend.title=element_blank()) + 
  scale_x_continuous(breaks = seq(0, 2, 1),limits=c(0, 2)) + ylim(c(0,5))


plot1 <- ggdraw() +
  draw_plot(p1, x = 0, y = .75, width = .7, height = .25) +
  draw_plot(p2, x = .7, y = .75, width = .3, height = .25) +
  draw_plot(p3, x = 0, y = .5, width = .7, height = .25) +
  draw_plot(p4, x = .7, y = .5, width = .3, height = .25) +
  draw_plot(p5, x = 0, y = .25, width = .7, height = .25) +
  draw_plot(p6, x = .7, y = .25, width = .3, height = .25) +
  draw_plot(p7, x = 0, y = 0, width = .7, height = .25) +
  draw_plot(p8, x = .7, y = 0, width = .3, height = .25) +
  draw_plot_label(label = c("A", "B", "C", "D","E","F","G","H"), size = 20, x = c(0, .7, 0, .7, 0, .7, 0, .7), y = c(1,1,.75,.75,.5,.5,.25,.25))
plot1

#ggsave("Figure_1.pdf", plot = plot1, dpi = 600, width=10)

lowspread <- read.csv(file.choose()) # lowspreadboxes.csv 
medspread <- read.csv(file.choose()) # medspreadboxes.csv
highspread <- read.csv(file.choose()) # highspreadboxes.csv

rows <- sample(nrow(breeding), 1000)
breeding <- breeding[rows, ]
nonbreeding <- nonbreeding[rows, ]

bdist <- dist(breeding[1:2])
bsampdensity <- reshape2::melt(as.matrix(bdist), varnames = c("row", "col"))
bsampdensity$Season <- "Breeding"
nbdist <- dist(nonbreeding[1:2])
nbsampdensity <- reshape2::melt(as.matrix(nbdist), varnames = c("row", "col"))
nbsampdensity$Season <- "Non-breeding"
alldensity <- rbind(bsampdensity, nbsampdensity) # Very compuationally intensive, creates large dataset
mu <- ddply(alldensity, c("Season"), summarise, Median=median(value))

f1breeding <- breeding[1:3]
rows <- sample(nrow(f1breeding), 250)
f1breeding <- f1breeding[rows, ]
f1breeding$rows <- rownames(f1breeding)
f1nonbreeding <- nonbreeding
f1nonbreeding <- f1nonbreeding[rows, ]
f1nonbreeding$rows <- f1breeding$row

f1breeding$Season <- "Breeding"
f1nonbreeding$Season <- "Non-breeding"
colnames(f1nonbreeding) <- c('x','y', 'zone','rows','Season')
points <- rbind(f1breeding, f1nonbreeding)


x <- c(0,1)
y<- c(0, 2)
dat <- data.frame(x,y)
bb1 <-  c(0.315925, 0.684075, 1.545161, 1.78814)
x2 <- c(0.315925,0.315925,0.684075,0.684075)
y2 <- c(1.545161,1.78814,1.78814,1.545161)
datsmall <- data.frame(x2,y2)
bb2 <-  c(0.19155, 0.80845, 1.463073, 1.870227)
x3 <- c(0.19155,0.19155,0.80845,0.80845)
y3 <- c(1.463073,1.870227,1.870227,1.463073)
datmed <- data.frame(x3,y3)
bb3 <-  c(0.067175, 0.932825, 1.380986, 1.952315)
x4 <- c(0.067175,0.067175,0.932825,0.932825)
y4 <- c(1.380986,1.952315,1.952315,1.380986)
datlarge <- data.frame(x4,y4)

# Whole population 

p1 <- ggplot(dat, aes(x=x, y=y)) + 
  theme_void() +
  geom_line(data=points,aes(x=x, y=y, group=rows), colour="gray") +
  geom_point(data=points,aes(x=x, y=y, colour=Season)) +
  scale_fill_identity() +
  theme(legend.position = "none") +
  scale_color_manual(values=c("#E69F00", "#56B4E9"))

p2 <- ggplot(alldensity, aes(x=value, colour=Season)) + 
  geom_density(size=1) +
  scale_colour_manual(values=c("#E69F00", "#56B4E9"), name="Scenario") +
  geom_vline(data=mu, aes(xintercept=Median, colour=Season), linetype="dashed",show.legend=FALSE) +
  theme_pubr() +
  labs(x = "Pairwise distance", y = "Density") +
  theme(legend.position="none",axis.text = element_text(size= 18), axis.title = element_text(size=18), legend.text=element_text(size=18),legend.title=element_text(size=18)) +
  scale_x_continuous(
    breaks = seq(0, 2, 1),
    limits=c(0, 2)) +
  ylim(c(0,5))

# Small box size plot

xdim <- 100 # x dimension of space
ydim <- xdim*2 # y dimension of space
grid <- raster(ncols=ydim,nrows=ydim) # create a grid of cells in raster format
extent(grid) <- c(0,1,0,ydim/xdim) # set up extent of the grid object
grid[]<-1 
xy <- data.frame(xyFromCell(grid,1:ncell(grid))) #extract coords of cell centres


brd_latrange <- c(quantile(xy$y,0.67),quantile(xy$y,1)) 
brd_lonrange <- c(quantile(xy$x,0),quantile(xy$x,1)) 
bigboxsize <- seq(0.5,1.5,length.out=3)

rows <- sample(nrow(breeding), 1000)
breeding2 <- breeding[rows, ]
breeding2$rows <- rownames(breeding2)

nonbreeding2 <- nonbreeding[rows, ]
nonbreeding2$rows <- rownames(nonbreeding2)


bb1 <-  c(0.315925, 0.684075, 1.545161, 1.78814)
# big box 2 (xmin, xmax, ymin, ymax)
bb2 <-  c(0.19155, 0.80845, 1.463073, 1.870227)
# big box 3 (xmin, xmax, ymin, ymax)
bb3 <-  c(0.067175, 0.932825, 1.380986, 1.952315)
bb <- list(bb1,bb2,bb3)


bigbox <- bigboxsize[1] # cycle through 3 options
bigboxdims <- c(bigbox*diff(brd_lonrange),bigbox*diff(brd_latrange))
midbrd <- c(mean(brd_lonrange),mean(brd_latrange))#midpoint of breeding range
area1 <- c(midbrd[1]-(0.5*bigboxdims[1]),midbrd[1]+(0.5*bigboxdims[1]),
           midbrd[2]-(0.5*bigboxdims[2]),midbrd[2]+(0.5*bigboxdims[2]))

nboxes <- sqrt(9) # 9 = number of sampled boxes
smallbox <- 0.12
smallboxdims <- c(smallbox*diff(brd_lonrange),smallbox*diff(brd_latrange))
xcoords <- seq(area1[1],area1[2],length.out=nboxes+2)
ycoords <- seq(area1[3],area1[4],length.out=nboxes+2)
xcoords <- xcoords[2:(length(xcoords)-1)]
ycoords <- ycoords[2:(length(ycoords)-1)]

counter <- 0 # reset to zero
count <- 0 
pls <- list()

for(i in 1:length(xcoords)){
  for(j in 1:length(ycoords)){
    counter<-counter+1
    midpoint <- c(xcoords[i],ycoords[j])
    boxcorners <- c(midpoint[1]-(0.5*smallboxdims[1]),midpoint[1]+(0.5*smallboxdims[1]),
                    midpoint[2]-(0.5*smallboxdims[2]),midpoint[2]+(0.5*smallboxdims[2]))
    xcs <- c(boxcorners[1],boxcorners[1],boxcorners[2],boxcorners[2],boxcorners[1])
    ycs <- c(boxcorners[3],boxcorners[4],boxcorners[4],boxcorners[3],boxcorners[3])
    xym <- data.frame(x=xcs,y=ycs)
    p <- Polygon(xym)
    pls[[counter]]<-p
  }
}
ps <- Polygons(pls,1)
area <- SpatialPolygons(list(ps))
pls <- list()
counter <- 0 # reset to zero
count <- 0 # reset to zero

fixN <- T 
try(coordinates(breeding2)<-c("x","y"), silent=T)
sampleN <- 200

# individuals for spread scenario
inds <- as.numeric(which(over(breeding2,area)==1)) #pick individuals in the study area using over() from sp package
samp_n<- ifelse(fixN==T,min(sampleN,length(inds)),floor(sampleprop*length(inds))) #number of individuals to sample
samp_inds <- inds[sample(1:length(inds),samp_n,replace=F)] #randomly sample the pop

tempbrd <- data.frame(breeding2[samp_inds,])
brddist <- dist(tempbrd[1:2])
df2brd <- reshape2::melt(as.matrix(brddist), varnames = c("row", "col"))


tempwin <- nonbreeding2[samp_inds, ]
colnames(tempwin)[1:2] <- c('x','y')
tempwin$Season <- 'Non-breeding'
tempbrd$Season <- 'Breeding'
tempbrd$rows <- tempwin$rows
temp <- rbind(tempbrd[,c(1,2,4,6)],tempwin[,c(1,2,4,5)])

windist <- dist(tempwin[1:2])
df2win <- reshape2::melt(as.matrix(windist), varnames = c("row", "col"))

df2brd <- subset(df2brd, value > 0)
df2win <- subset(df2win, value > 0)

df2brd$Season <- 'Breeding'
df2win$Season <- 'Non-breeding'

smalldensity <- rbind(df2brd, df2win)
smallmu <- ddply(smalldensity, c("Season"), summarise, Median=median(value))


p3 <- ggplot(dat, aes(x=x, y=y)) + 
  theme_void() +
  geom_point(data=points,aes(x=x, y=y, alpha=0.1), colour="lightgray") +
  geom_line(data=points,aes(x=x, y=y, group=rows, alpha=0)) +
  geom_line(data=temp,aes(x=x, y=y, group=rows), colour="black", alpha=0.5) +
  geom_point(data=temp,aes(x=x, y=y, colour=Season), size=1.5) +
  geom_polygon(data=lowspread, aes(x=x,y=y, fill=NA, group=Box), colour="#E69F00") +
  theme(legend.position = "none") +
  scale_fill_identity() +
  scale_color_manual(values=c("#E69F00", "#56B4E9"))

p4 <- ggplot(smalldensity, aes(x=value, colour=Season)) + 
  geom_density(size=1) +
  scale_colour_manual(values=c("#E69F00", "#56B4E9"), name="Scenario") +
  geom_vline(data=smallmu, aes(xintercept=Median, colour=Season), linetype="dashed",show.legend=FALSE) +
  theme_pubr() +
  labs(x = "Pairwise distance", y = "Density") +
  theme(legend.position="none",axis.text = element_text(size= 18), axis.title = element_text(size=18), legend.text=element_text(size=18),legend.title=element_text(size=18)) +
  scale_x_continuous(
    breaks = seq(0, 2, 1),
    limits=c(0, 2)) +
  ylim(c(0,5))

# Medium box size

brd_latrange <- c(quantile(xy$y,0.67),quantile(xy$y,1)) 
brd_lonrange <- c(quantile(xy$x,0),quantile(xy$x,1)) 
bigboxsize <- seq(0.5,1.5,length.out=3)

rows <- sample(nrow(breeding), 1000)
breeding2 <- breeding[rows, ]
breeding2$rows <- rownames(breeding2)

nonbreeding2 <- nonbreeding[rows, ]
nonbreeding2$rows <- rownames(nonbreeding2)


bb1 <-  c(0.315925, 0.684075, 1.545161, 1.78814)
# big box 2 (xmin, xmax, ymin, ymax)
bb2 <-  c(0.19155, 0.80845, 1.463073, 1.870227)
# big box 3 (xmin, xmax, ymin, ymax)
bb3 <-  c(0.067175, 0.932825, 1.380986, 1.952315)
bb <- list(bb1,bb2,bb3)


bigbox <- bigboxsize[2] # cycle through 3 options
bigboxdims <- c(bigbox*diff(brd_lonrange),bigbox*diff(brd_latrange))
midbrd <- c(mean(brd_lonrange),mean(brd_latrange))#midpoint of breeding range
area1 <- c(midbrd[1]-(0.5*bigboxdims[1]),midbrd[1]+(0.5*bigboxdims[1]),
           midbrd[2]-(0.5*bigboxdims[2]),midbrd[2]+(0.5*bigboxdims[2]))

nboxes <- sqrt(9) # 9 = number of sampled boxes
smallbox <- 0.12
smallboxdims <- c(smallbox*diff(brd_lonrange),smallbox*diff(brd_latrange))
xcoords <- seq(area1[1],area1[2],length.out=nboxes+2)
ycoords <- seq(area1[3],area1[4],length.out=nboxes+2)
xcoords <- xcoords[2:(length(xcoords)-1)]
ycoords <- ycoords[2:(length(ycoords)-1)]

for(i in 1:length(xcoords)){
  for(j in 1:length(ycoords)){
    counter<-counter+1
    midpoint <- c(xcoords[i],ycoords[j])
    boxcorners <- c(midpoint[1]-(0.5*smallboxdims[1]),midpoint[1]+(0.5*smallboxdims[1]),
                    midpoint[2]-(0.5*smallboxdims[2]),midpoint[2]+(0.5*smallboxdims[2]))
    xcs <- c(boxcorners[1],boxcorners[1],boxcorners[2],boxcorners[2],boxcorners[1])
    ycs <- c(boxcorners[3],boxcorners[4],boxcorners[4],boxcorners[3],boxcorners[3])
    xym <- data.frame(x=xcs,y=ycs)
    p <- Polygon(xym)
    pls[[counter]]<-p
  }
}
ps <- Polygons(pls,1)
area <- SpatialPolygons(list(ps))
pls <- list()
counter <- 0 # reset to zero
count <- 0 # reset to zero

fixN <- T 
try(coordinates(breeding2)<-c("x","y"), silent=T)
sampleN <- 200

# individuals for spread scenario
inds <- as.numeric(which(over(breeding2,area)==1)) #pick individuals in the study area using over() from sp package
samp_n<- ifelse(fixN==T,min(sampleN,length(inds)),floor(sampleprop*length(inds))) #number of individuals to sample
samp_inds <- inds[sample(1:length(inds),samp_n,replace=F)] #randomly sample the pop

tempbrd <- data.frame(breeding2[samp_inds,])
brddist <- dist(tempbrd[1:2])
df2brd <- reshape2::melt(as.matrix(brddist), varnames = c("row", "col"))


tempwin <- nonbreeding2[samp_inds, ]
colnames(tempwin)[1:2] <- c('x','y')
tempwin$Season <- 'Non-breeding'
tempbrd$Season <- 'Breeding'
tempbrd$rows <- tempwin$rows
temp <- rbind(tempbrd[,c(1,2,4,6)],tempwin[,c(1,2,4,5)])

windist <- dist(tempwin[1:2])
df2win <- reshape2::melt(as.matrix(windist), varnames = c("row", "col"))

df2brd <- subset(df2brd, value > 0)
df2win <- subset(df2win, value > 0)

df2brd$Season <- 'Breeding'
df2win$Season <- 'Non-breeding'

meddensity <- rbind(df2brd, df2win)
medmu <- ddply(smalldensity, c("Season"), summarise, Median=median(value))


p5 <- ggplot(dat, aes(x=x, y=y)) + 
  theme_void() +
  geom_point(data=points,aes(x=x, y=y, alpha=0.1), colour="lightgray") +
  geom_line(data=points,aes(x=x, y=y, group=rows, alpha=0)) +
  geom_line(data=temp,aes(x=x, y=y, group=rows), colour="black", alpha=0.5) +
  geom_point(data=temp,aes(x=x, y=y, colour=Season), size=1.5) +
  geom_polygon(data=medspread, aes(x=x,y=y, fill=NA, group=Box), colour="#E69F00") +
  theme(legend.position = "none") +
  scale_fill_identity() +
  scale_color_manual(values=c("#E69F00", "#56B4E9"))

p6 <- ggplot(meddensity, aes(x=value, colour=Season)) + 
  geom_density(size=1) +
  scale_colour_manual(values=c("#E69F00", "#56B4E9"), name="Scenario") +
  geom_vline(data=medmu, aes(xintercept=Median, colour=Season), linetype="dashed",show.legend=FALSE) +
  theme_pubr() +
  labs(x = "Pairwise distance", y = "Density") +
  theme(legend.position="none",axis.text = element_text(size= 18), axis.title = element_text(size=18), legend.text=element_text(size=18),legend.title=element_text(size=18)) +
  scale_x_continuous(
    breaks = seq(0, 2, 1),
    limits=c(0, 2)) +
  ylim(c(0,5))

# Large box size plot

brd_latrange <- c(quantile(xy$y,0.67),quantile(xy$y,1)) 
brd_lonrange <- c(quantile(xy$x,0),quantile(xy$x,1)) 
bigboxsize <- seq(0.5,1.5,length.out=3)

rows <- sample(nrow(breeding), 1000)
breeding2 <- breeding[rows, ]
breeding2$rows <- rownames(breeding2)

nonbreeding2 <- nonbreeding[rows, ]
nonbreeding2$rows <- rownames(nonbreeding2)


bb1 <-  c(0.315925, 0.684075, 1.545161, 1.78814)
# big box 2 (xmin, xmax, ymin, ymax)
bb2 <-  c(0.19155, 0.80845, 1.463073, 1.870227)
# big box 3 (xmin, xmax, ymin, ymax)
bb3 <-  c(0.067175, 0.932825, 1.380986, 1.952315)
bb <- list(bb1,bb2,bb3)


bigbox <- bigboxsize[3] # cycle through 3 options
bigboxdims <- c(bigbox*diff(brd_lonrange),bigbox*diff(brd_latrange))
midbrd <- c(mean(brd_lonrange),mean(brd_latrange))#midpoint of breeding range
area1 <- c(midbrd[1]-(0.5*bigboxdims[1]),midbrd[1]+(0.5*bigboxdims[1]),
           midbrd[2]-(0.5*bigboxdims[2]),midbrd[2]+(0.5*bigboxdims[2]))

nboxes <- sqrt(9) # 9 = number of sampled boxes
smallbox <- 0.12
smallboxdims <- c(smallbox*diff(brd_lonrange),smallbox*diff(brd_latrange))
xcoords <- seq(area1[1],area1[2],length.out=nboxes+2)
ycoords <- seq(area1[3],area1[4],length.out=nboxes+2)
xcoords <- xcoords[2:(length(xcoords)-1)]
ycoords <- ycoords[2:(length(ycoords)-1)]

for(i in 1:length(xcoords)){
  for(j in 1:length(ycoords)){
    counter<-counter+1
    midpoint <- c(xcoords[i],ycoords[j])
    boxcorners <- c(midpoint[1]-(0.5*smallboxdims[1]),midpoint[1]+(0.5*smallboxdims[1]),
                    midpoint[2]-(0.5*smallboxdims[2]),midpoint[2]+(0.5*smallboxdims[2]))
    xcs <- c(boxcorners[1],boxcorners[1],boxcorners[2],boxcorners[2],boxcorners[1])
    ycs <- c(boxcorners[3],boxcorners[4],boxcorners[4],boxcorners[3],boxcorners[3])
    xym <- data.frame(x=xcs,y=ycs)
    p <- Polygon(xym)
    pls[[counter]]<-p
  }
}
ps <- Polygons(pls,1)
area <- SpatialPolygons(list(ps))
pls <- list()
counter <- 0 # reset to zero
count <- 0 # reset to zero

fixN <- T 
try(coordinates(breeding2)<-c("x","y"), silent=T)
sampleN <- 200

# individuals for spread scenario
inds <- as.numeric(which(over(breeding2,area)==1)) #pick individuals in the study area using over() from sp package
samp_n<- ifelse(fixN==T,min(sampleN,length(inds)),floor(sampleprop*length(inds))) #number of individuals to sample
samp_inds <- inds[sample(1:length(inds),samp_n,replace=F)] #randomly sample the pop

tempbrd <- data.frame(breeding2[samp_inds,])
brddist <- dist(tempbrd[1:2])
df2brd <- reshape2::melt(as.matrix(brddist), varnames = c("row", "col"))


tempwin <- nonbreeding2[samp_inds, ]
colnames(tempwin)[1:2] <- c('x','y')
tempwin$Season <- 'Non-breeding'
tempbrd$Season <- 'Breeding'
tempbrd$rows <- tempwin$rows
temp <- rbind(tempbrd[,c(1,2,4,6)],tempwin[,c(1,2,4,5)])

windist <- dist(tempwin[1:2])
df2win <- reshape2::melt(as.matrix(windist), varnames = c("row", "col"))

df2brd <- subset(df2brd, value > 0)
df2win <- subset(df2win, value > 0)

df2brd$Season <- 'Breeding'
df2win$Season <- 'Non-breeding'

largedensity <- rbind(df2brd, df2win)
largemu <- ddply(largedensity, c("Season"), summarise, Median=median(value))

p7 <- ggplot(dat, aes(x=x, y=y)) + 
  theme_void() +
  geom_point(data=points,aes(x=x, y=y, alpha=0.1), colour="lightgray") +
  geom_line(data=points,aes(x=x, y=y, group=rows, alpha=0)) +
  geom_line(data=temp,aes(x=x, y=y, group=rows), colour="black", alpha=0.5) +
  geom_point(data=temp,aes(x=x, y=y, colour=Season), size=1.5) +
  geom_polygon(data=highspread, aes(x=x,y=y, fill=NA, group=Box), colour="#E69F00") +
  theme(legend.position = "none") +
  scale_fill_identity() +
  scale_color_manual(values=c("#E69F00", "#56B4E9"))

p8 <- ggplot(largedensity, aes(x=value, colour=Season)) + 
  geom_density(size=1, show.legend=FALSE) +
  scale_colour_manual(values=c("#E69F00", "#56B4E9"), name="Scenario") +
  geom_vline(data=largemu, aes(xintercept=Median, colour=Season), linetype="dashed",show.legend=FALSE) +
  theme_pubr() +
  labs(x = "Pairwise distance", y = "Density") +
  stat_density(aes(x=value, colour=Season),
               geom="line",position="identity", size = 0) + 
  guides(colour = guide_legend(override.aes=list(size=2))) +
  theme(legend.position="none",axis.text = element_text(size= 18), axis.title = element_text(size=18), legend.text=element_text(size=16),legend.title=element_blank()) +
  scale_x_continuous(
    breaks = seq(0, 2, 1),
    limits=c(0, 2)) +
  ylim(c(0,5))


plot <- ggdraw() +
  draw_plot(p1, x = 0, y = .75, width = .7, height = .25) +
  draw_plot(p2, x = .675, y = .75, width = .3, height = .25) +
  draw_plot(p3, x = 0, y = .5, width = .7, height = .25) +
  draw_plot(p4, x = .675, y = .5, width = .3, height = .25) +
  draw_plot(p5, x = 0, y = .25, width = .7, height = .25) +
  draw_plot(p6, x = .675, y = .25, width = .3, height = .25) +
  draw_plot(p7, x = 0, y = 0, width = .7, height = .25) +
  draw_plot(p8, x = .675, y = 0, width = .3, height = .25) +
  draw_plot_label(label = c("A", "B", "C", "D","E","F","G","H"), size = 20, x = c(0, .675, 0, .675, 0, .67, 0, .675), y = c(1,1,.75,.75,.5,.5,.25,.25))
plot
#ggsave("Figure_1.5.pdf", plot = plot, dpi = 600, width=10)

plot3 <- ggdraw() +
  draw_plot(plot1, x = 0, y = 0, width = .5, height = .95) +
  draw_plot(plot, x = .5, y = 0, width = .5, height = .95) +
  draw_plot_label(label = c("i. Single sampling site", "ii. Multiple discrete sampling sites"), size = 20, x = c(0.05, .5), y = c(1,1))
plot
ggsave("Figure_1.pdf", plot = plot3, dpi = 600, width=20, height = 12)

# 6.4 Figure 2 ----
allwhole <- rbind(fullscores[1:3,], areascores)
allwhole$Prop <- c(rep(1,3), rep(c(0.36815, 0.6169, 0.86565),3))
allwhole2 <- allwhole %>% gather(Method, Score, c(Mantel, Cohen))

allwhole2$Method <- c(rep('Mantel', 12), rep('Transition\nprobability',12))

# need to replace cohen in figure legend

p1 <- ggplot(allwhole2, aes(x=Prop, y=Score, by=factor(MC), colour=factor(MC))) +
  geom_hline(yintercept=0, linetype="dashed", 
             color = "black") +
  geom_point(size=2.5) +
  #geom_smooth(method="loess", aes(linetype=Method)) +
  geom_line(data=allwhole2, aes(linetype=Method), size=1) +
  labs(y = "Migratory Connectivity score", x = "Proportion of area") +
  theme_pubr() +
  theme(legend.position="bottom", axis.text = element_text(size= 20), axis.title = element_text(size=25), legend.text=element_text(size=20),legend.title=element_text(size=20)) +
  scale_colour_brewer(palette="Dark2", name="Migratory\nConnectivity", labels=c("Low", "Medium", "High"), guide = guide_legend(override.aes = list(linetype=NA, shape = c(16,16,16), size=4))) +
  coord_cartesian(ylim = c(-.1, 1))
p1
#ggsave("Figure_2.pdf", plot = p1, dpi = 600)
# 6.5 Figure 3 ----
p1 <- ggplot(Area, aes(x=Area, y=Mantel, by=as.factor(MC), colour=as.factor(MC))) +
  geom_hline(yintercept=0, linetype="dashed", 
             color = "black") +
  geom_jitter(show.legend = T, alpha = 0.1, position=position_jitterdodge()) +
  geom_point(data=a, aes(x = Area, y = Mantel.mean, fill="Sampled\nindividuals"), size=2.5, position=position_dodge(width=0.75)) +
  geom_errorbar(data = a, mapping = aes(x = Area, y = Mantel.mean, ymin = Mantel.mean - Mantel.sd, ymax = Mantel.mean + Mantel.sd), size=1, width=.2, position=position_dodge(width=0.75)) +
  geom_point(data=fullscores, aes(x = Area, y = Mantel, fill="Whole\npopulation"), size=5, stroke=1.5, position=position_dodge(width=0.75), shape=4) +
  labs(x="Area", y = "Mantel score") +
  theme_pubr() +
  theme(legend.position='bottom', axis.text = element_text(size= 25), axis.title = element_text(size=30), legend.text=element_text(size=25),legend.title=element_text(size=27)) +
  scale_x_discrete(labels=c("Small", "Medium", "Large")) +
  scale_color_manual(values=c('#1B9E77',"#D95F02","#7570B3"), labels=c("Low","Medium", "High")) +
  coord_cartesian(ylim = c(-0.1, 1)) +
  labs(fill='') +
  guides(colour=guide_legend(title="Migratory\nconnectivity:", override.aes = list(shape=c(NA,NA,NA), size=6)), linetype=FALSE, fill=guide_legend(title="", override.aes = list(shape=c(16,4))))


p2 <- ggplot(Area, aes(x=Area, y=Cohen, by=as.factor(MC), colour=as.factor(MC))) +
  geom_hline(yintercept=0, linetype="dashed", 
             color = "black") +
  geom_jitter(show.legend = T, alpha = 0.1, position=position_jitterdodge()) +
  geom_point(data=b, aes(x = Area, y = Cohen.mean, fill="Sampled individuals"), size=2.5, position=position_dodge(width=0.75)) +
  geom_errorbar(data = b, mapping = aes(x = Area, y = Cohen.mean, ymin = Cohen.mean - Cohen.sd, ymax = Cohen.mean + Cohen.sd), size=1, width=.2, position=position_dodge(width=0.75)) +
  geom_point(data=fullscores, aes(x = Area, y = Cohen, fill="Whole population"), size=5, stroke=1.5, position=position_dodge(width=0.75), shape=4) +
  labs(x="Area", y = "Transition probability score") +
  theme_pubr() +
  theme(legend.position="none", axis.text = element_text(size= 25), axis.title = element_text(size=30), legend.text=element_text(size=25),legend.title=element_text(size=27)) +
  scale_x_discrete(labels=c("Small", "Medium", "Large")) +
  scale_color_manual(values=c("#1B9E77", "#D95F02", "#7570B3"), labels=c("Low","Medium", "High")) +
  coord_cartesian(ylim = c(-0.1, 1)) +
  labs(fill='') +
  guides(colour=guide_legend(title="Migratory connectivity:", override.aes = list(shape=c(NA,NA,NA), size=3)), linetype=FALSE, fill=guide_legend(title="", override.aes = list(shape=c(4,16))))

plot <- ggdraw() +
  draw_plot(p2, x = .5, y = 0.095, width = .5, height = .905) +
  draw_plot(p1, x = 0, y = 0, width = .5, height = 1) +
  draw_plot_label(label = c("A", "B"), size = 20, x = c(0, .5), y = c(1, 1))
plot

#ggsave("Figure_3.pdf", plot = plot, dpi = 600, width=23)

# 6.6 Figure 4 ----

p1 <- ggplot(Spread, aes(x=Spread, y=Mantel, by=as.factor(MC), colour=as.factor(MC))) +
  geom_hline(yintercept=0, linetype="dashed", 
             color = "black") +
  geom_jitter(show.legend = T, alpha = 0.075, position=position_jitterdodge()) +
  geom_point(data=a2, aes(x = Spread, y = Mantel.mean, fill="Sampled\nindividuals"), size=2.5, position=position_dodge(width=0.75)) +
  geom_errorbar(data = a2, mapping = aes(x = Spread, y = Mantel.mean, ymin = Mantel.mean - Mantel.sd, ymax = Mantel.mean + Mantel.sd), size=1, width=.2, position=position_dodge(width=0.75)) +
  geom_point(data=areascores, aes(x = Spread, y = Mantel, fill="All individuals within\nspatial extent of\nsampling"), size=5, stroke=1.5, position=position_dodge(width=0.75), shape=4) +
  labs(x="Spread", y = "Mantel score") +
  theme_pubr() +
  theme(legend.position="bottom", axis.text = element_text(size= 25), axis.title = element_text(size=30), legend.text=element_text(size=23),legend.title=element_text(size=25)) +
  scale_x_discrete(labels=c("Low", "Medium", "High")) +
  scale_color_manual(values=c("#1B9E77", "#D95F02", "#7570B3"), labels=c("Low", "Medium", "High")) +
  coord_cartesian(ylim = c(-0.1, 1)) +
  labs(fill='') +
  guides(colour=guide_legend(title="Migratory\nconnectivity:", override.aes = list(shape=c(NA,NA,NA), size=3)), linetype=FALSE, fill=guide_legend(title="", override.aes = list(shape=c(4,16))))

p2 <- ggplot(Spread, aes(x=Spread, y=Cohen, by=as.factor(MC), colour=as.factor(MC))) +
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  geom_jitter(show.legend = T, alpha = 0.075, position=position_jitterdodge()) +
  geom_point(data=b2, aes(x = Spread, y = Cohen.mean, fill="Sampled individuals"), size=2.5, position=position_dodge(width=0.75)) +
  geom_errorbar(data = b2, mapping = aes(x = Spread, y = Cohen.mean, ymin = Cohen.mean - Cohen.sd, ymax = Cohen.mean + Cohen.sd), size=1, width=.2, position=position_dodge(width=0.75)) +
  geom_point(data=areascores, aes(x = Spread, y = Cohen, fill="All individuals within spatial\nextent of sampling"), size=5, stroke=1.5, position=position_dodge(width=0.75), shape=4) +
  labs(x="Spread", y = "Transition probability score") +
  theme_pubr() +
  theme(legend.position="none", axis.text = element_text(size= 25), axis.title = element_text(size=30), legend.text=element_text(size=12),legend.title=element_text(size=12)) +
  scale_x_discrete(labels=c("Low", "Medium", "High")) +
  scale_color_manual(values=c("#1B9E77", "#D95F02", "#7570B3"), labels=c("Low", "Medium", "High"), 
                     name=c("Migratory Connecitivity =",''), guide = guide_legend(override.aes = list(linetype=NA, shape = c(16,16,16)))) +
  coord_cartesian(ylim = c(-0.1, 1)) +
  labs(fill='') 

plot <- ggdraw() +
  draw_plot(p2, x = .5, y = 0.116, width = .5, height = .884) +
  draw_plot(p1, x = 0, y = 0, width = .5, height = 1) +
  draw_plot_label(label = c("A", "B"), size = 20, x = c(0, .5), y = c(1, 1))
plot
#ggsave("Figure_4.pdf", plot = plot, dpi = 600, width=22)


# 6.7 Figure 5 ----

p1 <- ggplot(MantelMDB, aes(x=as.factor(variable), y=value, by=as.factor(Sampled), colour=as.factor(Sampled))) +
  geom_segment(aes(x =0.6, y = 0.3256315, xend = 1.4, yend = 0.3256315), linetype="dashed", 
               color = "black") +
  geom_segment(aes(x =1.6, y = 0.5812184, xend = 2.4, yend = 0.5812184), linetype="dashed", 
               color = "black") +
  geom_segment(aes(x =2.6, y = 0.8332163, xend = 3.4, yend = 0.8332163), linetype="dashed", 
               color = "black") +
  geom_jitter(show.legend = T, alpha = 0.1, position=position_jitterdodge()) +
  geom_point(data=MDBa, aes(x = variable, y = value.mean), size=2.5, position=position_dodge(width=0.75)) +
  geom_errorbar(data = MDBa, mapping = aes(x = variable, y = value.mean, ymin = value.mean - value.sd, ymax = value.mean + value.sd), size=1, width=.2, position=position_dodge(width=0.75)) +
  labs(x="Migratory connectivity", y = "Mantel score") +
  theme_pubr() +
  theme(legend.position="none", axis.text = element_text(size= 25), axis.title = element_text(size=30), legend.text=element_text(size=25),legend.title=element_text(size=25)) +
  scale_x_discrete(labels=c("Low (0.33)", "Medium (0.58)", "High (0.83)")) +
  scale_colour_brewer(palette="Dark2", name="Sample size", 
                      labels=c("10 (0.1%)", "50 (0.5%)", "100 (1%)", "1000 (10%)", "2500 (25%)", "5000 (50%)"), guide = guide_legend(override.aes = list(linetype=NA, shape = c(16,16,16)))) +
  coord_cartesian(ylim = c(-.25, 1)) +
  scale_y_continuous(breaks=c(-.2,0,.2,.4,.6,.8,1))


p2 <- ggplot(CohenMDB, aes(x=as.factor(variable), y=value, by=as.factor(Sampled), colour=as.factor(Sampled))) +
  geom_segment(aes(x =0.6, y = 0.40800918 , xend = 1.4, yend = 0.40800918 ), linetype="dashed", 
               color = "black") +
  geom_segment(aes(x =1.6, y = 0.6378398, xend = 2.4, yend = 0.6378398), linetype="dashed", 
               color = "black") +
  geom_segment(aes(x =2.6, y = 0.7981142, xend = 3.4, yend = 0.7981142), linetype="dashed", 
               color = "black") +
  geom_jitter(show.legend = T, alpha = 0.1, position=position_jitterdodge()) +
  geom_point(data=MDBb, aes(x = variable, y = value.mean), size=2.5, position=position_dodge(width=0.75)) +
  geom_errorbar(data = MDBb, mapping = aes(x = variable, y = value.mean, ymin = value.mean - value.sd, ymax = value.mean + value.sd), size=1, width=.2, position=position_dodge(width=0.75)) +
  labs(x="Migratory Connectivity", y = "Transition probability\nscore") +
  scale_x_discrete(labels=c("Low (0.41)", "Medium (0.63)", "High (0.80)")) +
  theme_pubr() +
  theme(legend.position="bottom", axis.text = element_text(size= 25), axis.title = element_text(size=30), legend.text=element_text(size=25),legend.title=element_text(size=27)) +
  scale_colour_brewer(palette="Dark2", name="Sample size", 
                      labels=c("10 (0.1%)", "50 (0.5%)", "100 (1%)", "1000 (10%)", "2500 (25%)", "5000 (50%)"), guide = guide_legend(override.aes = list(linetype=NA, shape = c(16,16,16)))) +
  coord_cartesian(ylim = c(-.25, 1)) +
  scale_y_continuous(breaks=c(-.2,0,.2,.4,.6,.8,1))

plot <- ggdraw() +
  draw_plot(p1, x = 0.02, y = .53, width = .98, height = .47) +
  draw_plot(p2, x = 0, y = 0, width = 1, height = .53) +
  draw_plot_label(label = c("A", "B"), size = 20, x = c(0, 0), y = c(1, 0.55))
plot
#ggsave("Figure_5.pdf", plot = plot, dpi = 600, width=22)


# 6.8 Figure 6 ---- 

All_data <- read.csv(file.choose()) # metapopscores.csv

Area <- subset(All_data, Samp_Size == 200)
b <- summaryBy(Cohen ~ Area + MC, Area, FUN = c(mean, sd))
b2 <- summaryBy(Mantel ~ Area + MC, Area, FUN = c(mean, sd))

Allinds <- subset(All_data, Samp_Size == 10000)
c <- summaryBy(Cohen ~ Area + MC, Allinds, FUN = c(mean, sd))
c2 <- summaryBy(Mantel ~ Area + MC, Allinds, FUN = c(mean, sd))

Area$Area <- factor(Area$Area, levels=c("Small", "Medium", "Large"))

step1 <- subset(All_data, Samp_Size == 10000)
d <- summaryBy(Cohen ~ MC, step1, FUN = c(mean, sd))
d2 <- summaryBy(Mantel ~ MC, step1, FUN = c(mean, sd))
fullscoresCoh <- rbind(d,d,d)
fullscoresMan <- rbind(d2,d2,d2)
fullscoresCoh$Area <- c('Small','Small','Small', 'Medium','Medium','Medium','Large','Large','Large')
fullscoresMan$Area <- c('Small','Small','Small', 'Medium','Medium','Medium','Large','Large','Large')

noextrap <- subset(All_data, Samp_Size < 10000 & Samp_Size > 200)
e <- summaryBy(Cohen ~ Area + MC, noextrap, FUN = c(mean, sd))
e2 <- summaryBy(Mantel ~ Area + MC, noextrap, FUN = c(mean, sd))

p2 <- ggplot(Area, aes(x=Area, y=Cohen, by=as.factor(MC), colour=as.factor(MC))) +
  geom_jitter(show.legend = T, alpha = 0.075, position=position_jitterdodge()) +
  geom_point(data=b, aes(x = Area, y = Cohen.mean, fill="Sampled individuals"), size=2.5, position=position_dodge(width=0.75)) +
  geom_errorbar(data = b, mapping = aes(x = Area, y = Cohen.mean, ymin = Cohen.mean - Cohen.sd, ymax = Cohen.mean + Cohen.sd), size=1, width=.2, position=position_dodge(width=0.75)) +
  geom_point(data=fullscores, aes(x = Area, y = Cohen.mean, fill="All individuals"), size=5, stroke=1.5, position=position_dodge(width=0.75), shape=4) +
  labs(x="Area", y = "Transition probability\nscore") +
  theme_pubr() +
  theme(legend.position="none", axis.text = element_text(size= 23), axis.title = element_text(size=25), legend.text=element_text(size=15),legend.title=element_text(size=15)) +
  scale_x_discrete(labels=c("Small", "Medium", "Large")) +
  scale_color_manual(values=c("#1B9E77", "#D95F02", "#7570B3"), labels=c("Low", "Medium", "High")) +
  coord_cartesian(ylim = c(.6, 1)) +
  labs(fill='') +
  guides(colour=guide_legend(title="Migratory connectivity:", override.aes = list(shape=c(NA,NA,NA), size=3)), linetype=FALSE, fill=guide_legend(title="", override.aes = list(shape=c(4,16))))

p4 <- ggplot(Area, aes(x=Area, y=Cohen, by=as.factor(MC), colour=as.factor(MC))) +
  geom_jitter(show.legend = T, alpha = 0.075, position=position_jitterdodge()) +
  geom_point(data=b, aes(x = Area, y = Cohen.mean, fill="Sampled individuals"), size=2.5, position=position_dodge(width=0.75)) +
  geom_errorbar(data = b, mapping = aes(x = Area, y = Cohen.mean, ymin = Cohen.mean - Cohen.sd, ymax = Cohen.mean + Cohen.sd), size=1, width=.2, position=position_dodge(width=0.75)) +
  geom_point(data=e, aes(x = Area, y = Cohen.mean, fill="All individuals"), size=5, stroke=1.5, position=position_dodge(width=0.75), shape=4) +
  labs(x="Area", y = "Transition probability\nscore") +
  theme_pubr() +
  theme(legend.position="none", axis.text = element_text(size= 23), axis.title = element_text(size=25), legend.text=element_text(size=15),legend.title=element_text(size=15)) +
  scale_x_discrete(labels=c("Small", "Medium", "Large")) +
  scale_color_manual(values=c("#1B9E77", "#D95F02", "#7570B3"), labels=c("Low", "Medium", "High")) +
  coord_cartesian(ylim = c(.6, 1)) +
  labs(fill='') +
  guides(colour=guide_legend(title="Migratory connectivity:", override.aes = list(shape=c(NA,NA,NA), size=3)), linetype=FALSE, fill=guide_legend(title="", override.aes = list(shape=c(4,16))))

p1 <- ggplot(Area, aes(x=Area, y=Mantel, by=as.factor(MC), colour=as.factor(MC))) +
  geom_jitter(show.legend = T, alpha = 0.075, position=position_jitterdodge()) +
  geom_point(data=b2, aes(x = Area, y = Mantel.mean, fill="Sampled\nindividuals"), size=2.5, position=position_dodge(width=0.75)) +
  geom_errorbar(data = b2, mapping = aes(x = Area, y = Mantel.mean, ymin = Mantel.mean - Mantel.sd, ymax = Mantel.mean + Mantel.sd), size=1, width=.2, position=position_dodge(width=0.75)) +
  geom_point(data=fullscoresMan, aes(x = Area, y = Mantel.mean, fill="All individuals"), size=5, stroke=1.5, position=position_dodge(width=0.75), shape=4) +
  labs(x="Area", y = "Mantel score") +
  theme_pubr() +
  theme(legend.position="bottom", axis.text = element_text(size= 23), axis.title = element_text(size=25), legend.text=element_text(size=18),legend.title=element_text(size=18)) +
  scale_x_discrete(labels=c("Small", "Medium", "Large")) +
  scale_color_manual(values=c("#1B9E77", "#D95F02", "#7570B3"), labels=c("Low", "Medium", "High")) +
  coord_cartesian(ylim = c(.6, 1)) +
  labs(fill='') +
  guides(colour=guide_legend(title="Migratory\nconnectivity:", override.aes = list(shape=c(NA,NA,NA), size=3)), linetype=FALSE, fill=guide_legend(title="", override.aes = list(shape=c(4,16))))

p3 <- ggplot(Area, aes(x=Area, y=Mantel, by=as.factor(MC), colour=as.factor(MC))) +
  geom_jitter(show.legend = T, alpha = 0.075, position=position_jitterdodge()) +
  geom_point(data=b2, aes(x = Area, y = Mantel.mean, fill="Sampled\nindividuals"), size=2.5, position=position_dodge(width=0.75)) +
  geom_errorbar(data = b2, mapping = aes(x = Area, y = Mantel.mean, ymin = Mantel.mean - Mantel.sd, ymax = Mantel.mean + Mantel.sd), size=1, width=.2, position=position_dodge(width=0.75)) +
  geom_point(data=e2, aes(x = Area, y = Mantel.mean, fill="All individuals within\nspatial extent of\nsampling"), size=5, stroke=1.5, position=position_dodge(width=0.75), shape=4) +
  labs(x="Area", y = "Mantel score") +
  theme_pubr() +
  theme(legend.position="bottom", axis.text = element_text(size= 23), axis.title = element_text(size=25), legend.text=element_text(size=18),legend.title=element_text(size=18)) +
  scale_x_discrete(labels=c("Small", "Medium", "Large")) +
  scale_color_manual(values=c("#1B9E77", "#D95F02", "#7570B3"), labels=c("Low", "Medium", "High")) +
  coord_cartesian(ylim = c(.6, 1)) +
  labs(fill='') +
  guides(colour=guide_legend(title="Migratory\nconnectivity:", override.aes = list(shape=c(NA,NA,NA), size=3)), linetype=FALSE, fill=guide_legend(title="", override.aes = list(shape=c(4,16))))

plot <- ggdraw() +
  draw_plot(p2, x = .5, y = .57, width = .5, height = .43) +
  draw_plot(p4, x = .5, y = 0.095, width = .5, height = .405) +
  draw_plot(p3, x = 0, y = 0, width = .5, height = .5) +
  draw_plot(p1, x = 0, y = .5, width = .5, height = .5) +
  draw_plot_label(label = c("A", "B", "C", "D"), size = 20, x = c(0, .5, 0, .5), y = c(1,1,.5,.5))
plot
#ggsave("Figure_6.pdf", plot = plot, dpi = 600, width=18)

# 6.9 Figure SM1 ----
brd <- data.frame(breeding)
win <- data.frame(nonbreeding)

x <- c(0,0,1,1)
y <- c(1.33,2,2,1.33)
breedrange <- data.frame(x,y)

a <- brd[75,1:2]
a2 <- a
a2$y <- a2$y-2.2
b <- win[75,1:2]
colnames(b) <- c('x','y')
c <- rbind(a,b)
c$season <- c('Breeding','Non-breeding')
c[2,2] <- 0.7644167-1
c1 <- rbind(a,a2)
c1$season <- c('Breeding','Non-breeding')
c2 <- rbind(b,a2)

# movement south
p1 <- ggplot(brd, aes(x=x, y=y)) +
  geom_point(size=0.1, colour="white") +
  geom_point(data=nonbreeding, aes(x=lon, y=lat-1), colour="white", size=0.1) +
  geom_point(data=c1, aes(x=x, y=y, colour=season), size=2) +
  geom_polygon(data=breedrange, aes(x=x,y=y), fill=NA,colour="#E69F00") +
  geom_text(aes(label="Breeding area", y=2.1, x=0.12), colour="#E69F00") +
  theme_void() +
  scale_color_manual(values=c("#E69F00", "#56B4E9")) +
  theme(legend.position = "none", panel.border = element_rect(colour = "gray85", fill=NA, size=1)) +
  geom_segment(data = c1, aes(x = x[1], y = y[1], xend = x[2], yend = y[2]+0.02), colour = "black", arrow = arrow(length = unit(0.1, "cm")))

# distance and direction
p2 <- ggplot(brd, aes(x=x, y=y)) +
  geom_point(size=0.1, colour="white") +
  geom_point(data=nonbreeding, aes(x=lon, y=lat-1), colour="white", size=0.1) +
  geom_point(data=a, aes(x=x, y=y), colour="#E69F00", size=2) +
  geom_point(data=a2, aes(x=x, y=y), colour="#56B4E9", size=2, shape=1) +
  geom_point(data=b, aes(x=x, y=y), colour="#56B4E9", size=2) +
  geom_polygon(data=breedrange, aes(x=x,y=y), fill=NA,colour="#E69F00") +
  geom_text(aes(label="Breeding area", y=2.1, x=0.12), colour="#E69F00") +
  theme_void() +
  scale_color_manual(values=c("#56B4E9", "#E69F00")) +
  theme(legend.position = "none", panel.border = element_rect(colour = "gray85", fill=NA, size=1)) +
  geom_segment(data = c2, aes(x = x[2]+0.01, y = y[2]+0.01, xend = x[1]-0.015, yend = y[1]-0.02), colour = "black", arrow = arrow(length = unit(0.1, "cm")))

polar <- structure(list(degree = c(0L, 45L, 90L,180L,270L,360L), value = c(0, 1, 0, 0, 0, 0)), .Names = c("degree","value"), class = "data.frame", row.names = c(NA, -6L))
base <- ggplot(polar, aes(x=degree, y=value)) + theme_minimal() +
  scale_x_continuous(name="", breaks=seq(0,270,length.out=4), labels = c('N', 'E', 'S', 'W')) +
  scale_y_continuous(name="") +
  theme(axis.text.y = element_blank(), axis.ticks = element_blank())
p <- base + coord_polar(start=0)

compass <- p + geom_segment(aes(y=0, xend=degree, yend=value), colour="#56B4E9", size=2)

dists <- rlnorm(10000,7,1)*25
dists <- as.data.frame(dists)
dists <- subset(dists, dists<500000)

distance <- ggplot(dists, aes(x=dists)) + 
  geom_density(fill=NA, colour="black") +
  theme_minimal() +
  labs(x = "Distance", y = "Density") +
  theme(legend.position="none") +
  theme(axis.text.y = element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank()) +
  geom_segment(aes(x = 100000, y = 0.000010, xend = 100000, yend = 0.0000001), colour = "#56B4E9", size=1, arrow = arrow(length = unit(0.1, "cm")))

# Full pop
full <- ggplot(brd, aes(x=x, y=y)) +
  geom_point(size=0.1, colour="#E69F00") +
  geom_point(data=nonbreeding, aes(x=lon, y=lat-1), colour="#56B4E9", size=0.1) +
  theme_void() +
  scale_color_manual(values=c("#E69F00", "#56B4E9")) +
  theme(legend.position = "none", panel.border = element_rect(colour = "gray85", fill=NA, size=1))

plot <- ggdraw() +
  draw_plot(p1, x = 0.05, y = 0.6, width = .45, height = 0.39) +
  draw_plot(p2, x = .5, y = 0.6, width = .45, height = 0.39) +
  draw_plot(compass, x = 0, y = 0.37, width = .5, height = .24) +
  draw_plot(distance, x = .5, y = 0.4, width = .3, height = .2) +
  draw_plot(full, x = 0.275, y = 0.01, width = 0.45, height = .39) +
  draw_plot_label(label = c("A", "B", "C", "D"), size = 20, x = c(0.05, .5,.1, 0.15), y = c(0.99, 0.99,0.6,0.4))
plot
# ggsave("Figure_SM1.pdf", plot = plot, dpi = 600, width=12)


# 6.10 Figure SM2 ----

brd_locs <- read.csv(file.choose()) # brd_MC1.csv from patchy pop files
win_locs <- read.csv(file.choose()) # win_MC1.csv from patchy pop files

rows <- sample(nrow(brd_locs), 1000)
breeding1 <- brd_locs[rows, ]
nonbreeding1<- win_locs[rows, ]

bdist <- dist(breeding1[1:2])
bsampdensity <- reshape2::melt(as.matrix(bdist), varnames = c("row", "col"))
bsampdensity$Season <- "Breeding"
nbdist <- dist(nonbreeding1[1:2])
nbsampdensity <- reshape2::melt(as.matrix(nbdist), varnames = c("row", "col"))
nbsampdensity$Season <- "Non-breeding"
alldensity <- rbind(bsampdensity, nbsampdensity) # Very compuationally intensive, creates large dataset
mu <- ddply(alldensity, c("Season"), summarise, Median=median(value))

f1breeding <- breeding1[1:3]
rows <- sample(nrow(f1breeding), 250)
f1breeding <- f1breeding[rows, ]
f1breeding$rows <- rownames(f1breeding)
f1nonbreeding <- nonbreeding1
f1nonbreeding <- f1nonbreeding[rows, ]
f1nonbreeding$rows <- f1breeding$row

f1breeding$Season <- "Breeding"
f1nonbreeding$Season <- "Non-breeding"
colnames(f1nonbreeding) <- c('x','y', 'zone','rows','Season')
points <- rbind(f1breeding, f1nonbreeding)

x <- c(0,1)
y<- c(0, 2)
dat <- data.frame(x,y)

# Whole population 

plot1 <- ggplot(dat, aes(x=x, y=y)) + 
  theme_void() +
  geom_line(data=points,aes(x=x, y=y, group=rows), colour="gray") +
  geom_point(data=points,aes(x=x, y=y, colour=Season)) +
  scale_fill_identity() +
  theme(legend.position = "none") +
  scale_color_manual(values=c("#E69F00", "#56B4E9"))

plot2 <- ggplot(alldensity, aes(x=value, colour=Season)) + 
  geom_density(size=1) +
  scale_colour_manual(values=c("#E69F00", "#56B4E9"), name="Scenario") +
  geom_vline(data=mu, aes(xintercept=Median, colour=Season), linetype="dashed",show.legend=FALSE) +
  theme_pubr() +
  labs(x = "Pairwise distance", y = "Density") +
  theme(legend.position="none",axis.text = element_text(size= 18), axis.title = element_text(size=18), legend.text=element_text(size=20),legend.title=element_text(size=20)) +
  scale_x_continuous(
    breaks = seq(0, 2, 1),
    limits=c(0, 2)) +
  ylim(c(0,3))

# Small box size plot

bigbox <- bigboxsize[1] # sensible range 0.5 to 1.5. This would be added and small box fixed to allow big box to change but fix small box size.
bigboxdims <- c(bigbox*.33,bigbox*.22)

midbrd1 <- c(.165,1.89)#midpoint of breeding range
midwin1 <- c(.165,.55)
area1 <- c(midbrd1[1]-(0.5*bigboxdims[1]),midbrd1[1]+(0.5*bigboxdims[1]),
           midbrd1[2]-(0.5*bigboxdims[2]),midbrd1[2]+(0.5*bigboxdims[2]))
xcs <- c(area1[1],area1[1],area1[2],area1[2],area1[1])
ycs <- c(area1[3],area1[4],area1[4],area1[3],area1[3])
xym1 <- data.frame(x=xcs,y=ycs)
p1 <- Polygon(xym1)
pls[[1]]<-p1

midbrd2 <- c(.825, 1.89)#midpoint of breeding range
midwin2 <- c(.825,.55)
area2 <- c(midbrd2[1]-(0.5*bigboxdims[1]),midbrd2[1]+(0.5*bigboxdims[1]),
           midbrd2[2]-(0.5*bigboxdims[2]),midbrd2[2]+(0.5*bigboxdims[2]))
xcs <- c(area2[1],area2[1],area2[2],area2[2],area2[1])
ycs <- c(area2[3],area2[4],area2[4],area2[3],area2[3])
xym2 <- data.frame(x=xcs,y=ycs)
p2 <- Polygon(xym2)
pls[[2]]<-p2

midbrd3 <- c(.165,1.44)#midpoint of breeding range
midwin3 <- c(.165,.11)
area3 <- c(midbrd3[1]-(0.5*bigboxdims[1]),midbrd3[1]+(0.5*bigboxdims[1]),
           midbrd3[2]-(0.5*bigboxdims[2]),midbrd3[2]+(0.5*bigboxdims[2]))
xcs <- c(area3[1],area3[1],area3[2],area3[2],area3[1])
ycs <- c(area3[3],area3[4],area3[4],area3[3],area3[3])
xym3 <- data.frame(x=xcs,y=ycs)
p3 <- Polygon(xym3)
pls[[3]]<-p3

midbrd4 <- c(.825,1.44)#midpoint of breeding range
midwin4 <- c(.825,.11)
area4 <- c(midbrd4[1]-(0.5*bigboxdims[1]),midbrd4[1]+(0.5*bigboxdims[1]),
           midbrd4[2]-(0.5*bigboxdims[2]),midbrd4[2]+(0.5*bigboxdims[2]))
xcs <- c(area4[1],area4[1],area4[2],area4[2],area4[1])
ycs <- c(area4[3],area4[4],area4[4],area4[3],area4[3])
xym4 <- data.frame(x=xcs,y=ycs)
p4 <- Polygon(xym4)
pls[[4]]<-p4

slp <- pls
ps <- Polygons(slp,1)
area <- SpatialPolygons(list(ps))
pls <- list()

# loop through sample sizes (now fixed at 200)
sampleN <- 200

##### Sample populations and calculate mantel statistics
fixN <- T 
try(coordinates(breeding1)<-c("x","y"), silent=T)

inds <- as.numeric(which(over(breeding1,area)==1)) #pick individuals in the study area using over() from sp package
samp_n<- ifelse(fixN==T,min(sampleN,length(inds)),floor(sampleprop*length(inds))) #number of individuals to sample
samp_inds <- inds[sample(1:length(inds),samp_n,replace=F)] #randomly sample the pop
try(coordinates(nonbreeding1)<-c("lon","lat"), silent=T) 
samp_brd <- data.frame(breeding1[samp_inds,])
samp_win <- data.frame(nonbreeding1[samp_inds,])
colnames(samp_win) <- c('x','y','zone','optional')
samp_brd$rows <- rownames(samp_brd)
samp_win$rows <- rownames(samp_win)
samp_brd$Season <- "Breeding"
samp_win$Season <- "Non-breeding"

boxpoints <- rbind(samp_brd, samp_win)

bdist <- dist(samp_brd[1:2])
bsampdensity <- reshape2::melt(as.matrix(bdist), varnames = c("row", "col"))
bsampdensity$Season <- "Breeding"
bsampdensity$Sample <- "Sampled"
nbdist <- dist(samp_win[1:2])
nbsampdensity <- reshape2::melt(as.matrix(nbdist), varnames = c("row", "col"))
nbsampdensity$Season <- "Non-breeding"
nbsampdensity$Sample <- "Sampled"

smalldensity <- rbind(bsampdensity, nbsampdensity)
smallmu <- ddply(smalldensity, c("Season"), summarise, Median=median(value))


xym1$Box <- 1
xym2$Box <- 2
xym3$Box <- 3
xym4$Box <- 4
areasmall <- rbind(xym1,xym2,xym3,xym4)


plot3 <- ggplot(dat, aes(x=x, y=y)) + 
  theme_void() +
  geom_point(data=points,aes(x=x, y=y, alpha=0.1), colour="lightgray") +
  geom_line(data=points,aes(x=x, y=y, group=rows, alpha=0)) +
  geom_line(data=boxpoints,aes(x=x, y=y, group=rows), colour="black", alpha=0.5) +
  geom_point(data=boxpoints,aes(x=x, y=y, colour=Season), size=1.5) +
  geom_polygon(data=areasmall, aes(x=x,y=y, fill=NA, group=Box), colour="#E69F00") +
  theme(legend.position = "none") +
  scale_fill_identity() +
  scale_color_manual(values=c("#E69F00", "#56B4E9"))

plot4 <- ggplot(smalldensity, aes(x=value, colour=Season)) + 
  geom_density(size=1) +
  scale_colour_manual(values=c("#E69F00", "#56B4E9"), name="Scenario") +
  geom_vline(data=smallmu, aes(xintercept=Median, colour=Season), linetype="dashed",show.legend=FALSE) +
  theme_pubr() +
  labs(x = "Pairwise distance", y = "Density") +
  theme(legend.position="none",axis.text = element_text(size= 18), axis.title = element_text(size=18), legend.text=element_text(size=20),legend.title=element_text(size=20)) +
  scale_x_continuous(
    breaks = seq(0, 2, 1),
    limits=c(0, 2)) +
  ylim(c(0,3))

# Medium box size

bigbox <- bigboxsize[2] # sensible range 0.5 to 1.5. This would be added and small box fixed to allow big box to change but fix small box size.
bigboxdims <- c(bigbox*.33,bigbox*.22)

midbrd1 <- c(.165,1.89)#midpoint of breeding range
midwin1 <- c(.165,.55)
area1 <- c(midbrd1[1]-(0.5*bigboxdims[1]),midbrd1[1]+(0.5*bigboxdims[1]),
           midbrd1[2]-(0.5*bigboxdims[2]),midbrd1[2]+(0.5*bigboxdims[2]))
xcs <- c(area1[1],area1[1],area1[2],area1[2],area1[1])
ycs <- c(area1[3],area1[4],area1[4],area1[3],area1[3])
xym1 <- data.frame(x=xcs,y=ycs)
p1 <- Polygon(xym1)
pls[[1]]<-p1

midbrd2 <- c(.825, 1.89)#midpoint of breeding range
midwin2 <- c(.825,.55)
area2 <- c(midbrd2[1]-(0.5*bigboxdims[1]),midbrd2[1]+(0.5*bigboxdims[1]),
           midbrd2[2]-(0.5*bigboxdims[2]),midbrd2[2]+(0.5*bigboxdims[2]))
xcs <- c(area2[1],area2[1],area2[2],area2[2],area2[1])
ycs <- c(area2[3],area2[4],area2[4],area2[3],area2[3])
xym2 <- data.frame(x=xcs,y=ycs)
p2 <- Polygon(xym2)
pls[[2]]<-p2

midbrd3 <- c(.165,1.44)#midpoint of breeding range
midwin3 <- c(.165,.11)
area3 <- c(midbrd3[1]-(0.5*bigboxdims[1]),midbrd3[1]+(0.5*bigboxdims[1]),
           midbrd3[2]-(0.5*bigboxdims[2]),midbrd3[2]+(0.5*bigboxdims[2]))
xcs <- c(area3[1],area3[1],area3[2],area3[2],area3[1])
ycs <- c(area3[3],area3[4],area3[4],area3[3],area3[3])
xym3 <- data.frame(x=xcs,y=ycs)
p3 <- Polygon(xym3)
pls[[3]]<-p3

midbrd4 <- c(.825,1.44)#midpoint of breeding range
midwin4 <- c(.825,.11)
area4 <- c(midbrd4[1]-(0.5*bigboxdims[1]),midbrd4[1]+(0.5*bigboxdims[1]),
           midbrd4[2]-(0.5*bigboxdims[2]),midbrd4[2]+(0.5*bigboxdims[2]))
xcs <- c(area4[1],area4[1],area4[2],area4[2],area4[1])
ycs <- c(area4[3],area4[4],area4[4],area4[3],area4[3])
xym4 <- data.frame(x=xcs,y=ycs)
p4 <- Polygon(xym4)
pls[[4]]<-p4


slp <- pls
ps <- Polygons(slp,1)
area <- SpatialPolygons(list(ps))
pls <- list()


# loop through sample sizes (now fixed at 200)
sampleN <- 200

##### Sample populations and calculate mantel statistics
fixN <- T 
try(coordinates(breeding1)<-c("x","y"), silent=T)

inds <- as.numeric(which(over(breeding1,area)==1)) #pick individuals in the study area using over() from sp package
samp_n<- ifelse(fixN==T,min(sampleN,length(inds)),floor(sampleprop*length(inds))) #number of individuals to sample
samp_inds <- inds[sample(1:length(inds),samp_n,replace=F)] #randomly sample the pop
try(coordinates(nonbreeding1)<-c("lon","lat"), silent=T) 
samp_brd <- data.frame(breeding1[samp_inds,])
samp_win <- data.frame(nonbreeding1[samp_inds,])
colnames(samp_win) <- c('x','y','zone','optional')
samp_brd$rows <- rownames(samp_brd)
samp_win$rows <- rownames(samp_win)
samp_brd$Season <- "Breeding"
samp_win$Season <- "Non-breeding"

boxpoints <- rbind(samp_brd, samp_win)

bdist <- dist(samp_brd[1:2])
bsampdensity <- reshape2::melt(as.matrix(bdist), varnames = c("row", "col"))
bsampdensity$Season <- "Breeding"
bsampdensity$Sample <- "Sampled"
nbdist <- dist(samp_win[1:2])
nbsampdensity <- reshape2::melt(as.matrix(nbdist), varnames = c("row", "col"))
nbsampdensity$Season <- "Non-breeding"
nbsampdensity$Sample <- "Sampled"

meddensity <- rbind(bsampdensity, nbsampdensity)
medmu <- ddply(smalldensity, c("Season"), summarise, Median=median(value))


xym1$Box <- 1
xym2$Box <- 2
xym3$Box <- 3
xym4$Box <- 4
areamed <- rbind(xym1,xym2,xym3,xym4)

p5 <- ggplot(dat, aes(x=x, y=y)) + 
  theme_void() +
  geom_point(data=points,aes(x=x, y=y, alpha=0.1), colour="lightgray") +
  geom_line(data=points,aes(x=x, y=y, group=rows, alpha=0)) +
  geom_line(data=boxpoints,aes(x=x, y=y, group=rows), colour="black", alpha=0.5) +
  geom_point(data=boxpoints,aes(x=x, y=y, colour=Season), size=1.5) +
  geom_polygon(data=areamed, aes(x=x,y=y, fill=NA, group=Box), colour="#E69F00") +
  theme(legend.position = "none") +
  scale_fill_identity() +
  scale_color_manual(values=c("#E69F00", "#56B4E9"))

p6 <- ggplot(meddensity, aes(x=value, colour=Season)) + 
  geom_density(size=1) +
  scale_colour_manual(values=c("#E69F00", "#56B4E9"), name="Scenario") +
  geom_vline(data=medmu, aes(xintercept=Median, colour=Season), linetype="dashed",show.legend=FALSE) +
  theme_pubr() +
  labs(x = "Pairwise distance", y = "Density") +
  theme(legend.position="none",axis.text = element_text(size= 18), axis.title = element_text(size=18), legend.text=element_text(size=20),legend.title=element_text(size=20)) +
  scale_x_continuous(
    breaks = seq(0, 2, 1),
    limits=c(0, 2)) +
  ylim(c(0,3))

# Large box size plot

bigbox <- bigboxsize[3] # sensible range 0.5 to 1.5. This would be added and small box fixed to allow big box to change but fix small box size.
bigboxdims <- c(bigbox*.33,bigbox*.22)

midbrd1 <- c(.165,1.89)#midpoint of breeding range
midwin1 <- c(.165,.55)
area1 <- c(midbrd1[1]-(0.5*bigboxdims[1]),midbrd1[1]+(0.5*bigboxdims[1]),
           midbrd1[2]-(0.5*bigboxdims[2]),midbrd1[2]+(0.5*bigboxdims[2]))
xcs <- c(area1[1],area1[1],area1[2],area1[2],area1[1])
ycs <- c(area1[3],area1[4],area1[4],area1[3],area1[3])
xym1 <- data.frame(x=xcs,y=ycs)
p1 <- Polygon(xym1)
pls[[1]]<-p1

midbrd2 <- c(.825, 1.89)#midpoint of breeding range
midwin2 <- c(.825,.55)
area2 <- c(midbrd2[1]-(0.5*bigboxdims[1]),midbrd2[1]+(0.5*bigboxdims[1]),
           midbrd2[2]-(0.5*bigboxdims[2]),midbrd2[2]+(0.5*bigboxdims[2]))
xcs <- c(area2[1],area2[1],area2[2],area2[2],area2[1])
ycs <- c(area2[3],area2[4],area2[4],area2[3],area2[3])
xym2 <- data.frame(x=xcs,y=ycs)
p2 <- Polygon(xym2)
pls[[2]]<-p2

midbrd3 <- c(.165,1.44)#midpoint of breeding range
midwin3 <- c(.165,.11)
area3 <- c(midbrd3[1]-(0.5*bigboxdims[1]),midbrd3[1]+(0.5*bigboxdims[1]),
           midbrd3[2]-(0.5*bigboxdims[2]),midbrd3[2]+(0.5*bigboxdims[2]))
xcs <- c(area3[1],area3[1],area3[2],area3[2],area3[1])
ycs <- c(area3[3],area3[4],area3[4],area3[3],area3[3])
xym3 <- data.frame(x=xcs,y=ycs)
p3 <- Polygon(xym3)
pls[[3]]<-p3

midbrd4 <- c(.825,1.44)#midpoint of breeding range
midwin4 <- c(.825,.11)
area4 <- c(midbrd4[1]-(0.5*bigboxdims[1]),midbrd4[1]+(0.5*bigboxdims[1]),
           midbrd4[2]-(0.5*bigboxdims[2]),midbrd4[2]+(0.5*bigboxdims[2]))
xcs <- c(area4[1],area4[1],area4[2],area4[2],area4[1])
ycs <- c(area4[3],area4[4],area4[4],area4[3],area4[3])
xym4 <- data.frame(x=xcs,y=ycs)
p4 <- Polygon(xym4)
pls[[4]]<-p4


slp <- pls
ps <- Polygons(slp,1)
area <- SpatialPolygons(list(ps))
pls <- list()


# loop through sample sizes (now fixed at 200)
sampleN <- 200

##### Sample populations and calculate mantel statistics
fixN <- T 
try(coordinates(breeding1)<-c("x","y"), silent=T)

inds <- as.numeric(which(over(breeding1,area)==1)) #pick individuals in the study area using over() from sp package
samp_n<- ifelse(fixN==T,min(sampleN,length(inds)),floor(sampleprop*length(inds))) #number of individuals to sample
samp_inds <- inds[sample(1:length(inds),samp_n,replace=F)] #randomly sample the pop
try(coordinates(nonbreeding1)<-c("lon","lat"), silent=T) 
samp_brd <- data.frame(breeding1[samp_inds,])
samp_win <- data.frame(nonbreeding1[samp_inds,])
colnames(samp_win) <- c('x','y','zone','optional')
samp_brd$rows <- rownames(samp_brd)
samp_win$rows <- rownames(samp_win)
samp_brd$Season <- "Breeding"
samp_win$Season <- "Non-breeding"

boxpoints <- rbind(samp_brd, samp_win)

bdist <- dist(samp_brd[1:2])
bsampdensity <- reshape2::melt(as.matrix(bdist), varnames = c("row", "col"))
bsampdensity$Season <- "Breeding"
bsampdensity$Sample <- "Sampled"
nbdist <- dist(samp_win[1:2])
nbsampdensity <- reshape2::melt(as.matrix(nbdist), varnames = c("row", "col"))
nbsampdensity$Season <- "Non-breeding"
nbsampdensity$Sample <- "Sampled"

largedensity <- rbind(bsampdensity, nbsampdensity)
largemu <- ddply(smalldensity, c("Season"), summarise, Median=median(value))


xym1$Box <- 1
xym2$Box <- 2
xym3$Box <- 3
xym4$Box <- 4
arealarge <- rbind(xym1,xym2,xym3,xym4)

p7 <- ggplot(dat, aes(x=x, y=y)) + 
  theme_void() +
  geom_point(data=points,aes(x=x, y=y, alpha=0.1), colour="lightgray") +
  geom_line(data=points,aes(x=x, y=y, group=rows, alpha=0)) +
  geom_line(data=boxpoints,aes(x=x, y=y, group=rows), colour="black", alpha=0.5) +
  geom_point(data=boxpoints,aes(x=x, y=y, colour=Season), size=1.5) +
  geom_polygon(data=arealarge, aes(x=x,y=y, fill=NA, group=Box), colour="#E69F00") +
  theme(legend.position = "none") +
  scale_fill_identity() +
  scale_color_manual(values=c("#E69F00", "#56B4E9"))

p8 <- ggplot(largedensity, aes(x=value, colour=Season)) + 
  geom_density(size=1) +
  scale_colour_manual(values=c("#E69F00", "#56B4E9"), name="Scenario") +
  geom_vline(data=largemu, aes(xintercept=Median, colour=Season), linetype="dashed",show.legend=FALSE) +
  theme_pubr() +
  labs(x = "Pairwise distance", y = "Density") +
  theme(legend.position="bottom",axis.text = element_text(size= 18), axis.title = element_text(size=18), legend.text=element_text(size=10),legend.title=element_text(size=12)) +
  scale_x_continuous(
    breaks = seq(0, 2, 1),
    limits=c(0, 2)) +
  ylim(c(0,3))


allplot1 <- ggdraw() +
  draw_plot(plot1, x = .05, y = .75, width = .5, height = .25) +
  draw_plot(plot2, x = .55, y = .75, width = .44, height = .25) +
  draw_plot(plot3, x = .05, y = .5, width = .5, height = .25) +
  draw_plot(plot4, x = .55, y = .5, width = .44, height = .25) +
  draw_plot(p5, x = .05, y = .25, width = .5, height = .25) +
  draw_plot(p6, x = .55, y = .25, width = .44, height = .25) +
  draw_plot(p7, x = .05, y = 0, width = .5, height = .25) +
  draw_plot(p8, x = .55, y = 0, width = .44, height = .25) +
  draw_plot_label(label = c("A", "B", "C", "D","E","F","G","H"), size = 20, x = c(0, .55, 0, .55, 0, .55, 0, .55), y = c(1,1,.75,.75,.5,.5,.25,.25))
allplot1

#ggsave("Figure_SM2.pdf", plot = allplot1, dpi = 600, width=10)

# 6.11 Figure SM3 ----

p1 <- ggplot(Area, aes(x=Area, y=Mantel, by=as.factor(MC), colour=as.factor(MC))) +
  geom_hline(yintercept=0, linetype="dashed", 
             color = "black") +
  geom_jitter(show.legend = T, alpha = 0.1, position=position_jitterdodge()) +
  geom_point(data=a, aes(x = Area, y = Mantel.mean, fill="Sampled\nindividuals"), size=2.5, position=position_dodge(width=0.75)) +
  geom_errorbar(data = a, mapping = aes(x = Area, y = Mantel.mean, ymin = Mantel.mean - Mantel.sd, ymax = Mantel.mean + Mantel.sd), size=1, width=.2, position=position_dodge(width=0.75)) +
  geom_point(data=areascores, aes(x = Area, y = Mantel, fill="All individuals within\nspatial extent of\nsampling"), size=5, stroke=1.5, position=position_dodge(width=0.75), shape=4) +
  labs(x="Area", y = "Mantel score") +
  theme_pubr() +
  theme(legend.position="bottom", axis.text = element_text(size= 25), axis.title = element_text(size=30), legend.text=element_text(size=23),legend.title=element_text(size=25)) +
  scale_x_discrete(labels=c("Small", "Medium", "Large")) +
  scale_color_manual(values=c("#1B9E77", "#D95F02", "#7570B3"), labels=c("Low", "Medium", "High")) +
  coord_cartesian(ylim = c(-0.1, 1)) +
  labs(fill='') +
  guides(colour=guide_legend(title="Migratory\nconnectivity:", override.aes = list(shape=c(NA,NA,NA), size=3)), linetype=FALSE, fill=guide_legend(title="", override.aes = list(shape=c(4,16))))

p2 <- ggplot(Area, aes(x=Area, y=Cohen, by=as.factor(MC), colour=as.factor(MC))) +
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  geom_jitter(show.legend = T, alpha = 0.1, position=position_jitterdodge()) +
  geom_point(data=b, aes(x = Area, y = Cohen.mean, fill="Sampled individuals"), size=2.5, position=position_dodge(width=0.75)) +
  geom_errorbar(data = b, mapping = aes(x = Area, y = Cohen.mean, ymin = Cohen.mean - Cohen.sd, ymax = Cohen.mean + Cohen.sd), size=1, width=.2, position=position_dodge(width=0.75)) +
  geom_point(data=areascores, aes(x = Area, y = Cohen, fill="All individuals within spatial\nextent of sampling"), size=5, stroke=1.5, position=position_dodge(width=0.75), shape=4) +
  labs(x="Area", y = "Transition probability score") +
  theme_pubr() +
  theme(legend.position="none", axis.text = element_text(size= 25), axis.title = element_text(size=30), legend.text=element_text(size=15),legend.title=element_text(size=15)) +
  scale_x_discrete(labels=c("Small", "Medium", "Large")) +
  scale_color_manual(values=c("#1B9E77", "#D95F02", "#7570B3"), labels=c("Low", "Medium", "High")) +
  coord_cartesian(ylim = c(-0.1, 1)) +
  labs(fill='') +
  guides(colour=guide_legend(title="Migratory connectivity:", override.aes = list(shape=c(NA,NA,NA), size=3)), linetype=FALSE, fill=guide_legend(title="", override.aes = list(shape=c(4,16))))

plot <- ggdraw() +
  draw_plot(p2, x = .5, y = 0.116, width = .5, height = .884) +
  draw_plot(p1, x = 0, y = 0, width = .5, height = 1) +
  draw_plot_label(label = c("A", "B"), size = 20, x = c(0, .5), y = c(1, 1))
plot
#ggsave("Figure_SM3.pdf", plot = plot, dpi = 600, width=22)

# 6.12 Figure SM4 ----
setwd("~/OneDrive - University of East Anglia/Mig connectivity paper/Spread files")
lowspread <- read.csv('lowspreadboxes.csv')
medspread <- read.csv('medspreadboxes.csv')
highspread <- read.csv('highspreadboxes.csv')

xdim <- 100 # x dimension of space
ydim <- xdim*2 # y dimension of space
grid <- raster(ncols=ydim,nrows=ydim) # create a grid of cells in raster format
extent(grid) <- c(0,1,0,ydim/xdim) # set up extent of the grid object
grid[]<-1 
xy <- data.frame(xyFromCell(grid,1:ncell(grid))) #extract coords of cell centres

brd_latrange <- c(quantile(xy$y,0.67),quantile(xy$y,1)) 
brd_lonrange <- c(quantile(xy$x,0),quantile(xy$x,1)) 
bigboxsize <- seq(0.5,1.5,length.out=3)

rows <- sample(nrow(breeding), 1000)
breeding2 <- breeding[rows, ]
breeding2$rows <- rownames(breeding2)

nonbreeding2 <- nonbreeding[rows, ]
nonbreeding2$rows <- rownames(nonbreeding2)

x <- c(0,1)
y<- c(0, 2)
dat <- data.frame(x,y)

counter <- 0 
count <- 0
pls <- list()

listhistbrd <- list()
listhistwin <- list()
mulist <- list()
mu2list <- list()

bb1 <-  c(0.315925, 0.684075, 1.545161, 1.78814)
# big box 2 (xmin, xmax, ymin, ymax)
bb2 <-  c(0.19155, 0.80845, 1.463073, 1.870227)
# big box 3 (xmin, xmax, ymin, ymax)
bb3 <-  c(0.067175, 0.932825, 1.380986, 1.952315)
bb <- list(bb1,bb2,bb3)
big_brd_list <- list()

for (k in 1:3){
  
  bigbox <- bigboxsize[k] # cycle through 3 options
  bigboxdims <- c(bigbox*diff(brd_lonrange),bigbox*diff(brd_latrange))
  midbrd <- c(mean(brd_lonrange),mean(brd_latrange))#midpoint of breeding range
  area1 <- c(midbrd[1]-(0.5*bigboxdims[1]),midbrd[1]+(0.5*bigboxdims[1]),
             midbrd[2]-(0.5*bigboxdims[2]),midbrd[2]+(0.5*bigboxdims[2]))
  
  nboxes <- sqrt(9) # 9 = number of sampled boxes
  smallbox <- 0.12
  smallboxdims <- c(smallbox*diff(brd_lonrange),smallbox*diff(brd_latrange))
  xcoords <- seq(area1[1],area1[2],length.out=nboxes+2)
  ycoords <- seq(area1[3],area1[4],length.out=nboxes+2)
  xcoords <- xcoords[2:(length(xcoords)-1)]
  ycoords <- ycoords[2:(length(ycoords)-1)]
  
  for(i in 1:length(xcoords)){
    for(j in 1:length(ycoords)){
      counter<-counter+1
      midpoint <- c(xcoords[i],ycoords[j])
      boxcorners <- c(midpoint[1]-(0.5*smallboxdims[1]),midpoint[1]+(0.5*smallboxdims[1]),
                      midpoint[2]-(0.5*smallboxdims[2]),midpoint[2]+(0.5*smallboxdims[2]))
      xcs <- c(boxcorners[1],boxcorners[1],boxcorners[2],boxcorners[2],boxcorners[1])
      ycs <- c(boxcorners[3],boxcorners[4],boxcorners[4],boxcorners[3],boxcorners[3])
      xym <- data.frame(x=xcs,y=ycs)
      p <- Polygon(xym)
      pls[[counter]]<-p
    }
  }
  ps <- Polygons(pls,1)
  area <- SpatialPolygons(list(ps))
  pls <- list()
  counter <- 0 # reset to zero
  count <- 0 # reset to zero
  
  ##### Sample populations and calculate mantel statistics
  fixN <- T 
  try(coordinates(breeding)<-c("x","y"), silent=T)
  sampleN <- 200
  
  # individuals for spread scenario
  inds <- as.numeric(which(over(breeding,area)==1)) #pick individuals in the study area using over() from sp package
  samp_n<- ifelse(fixN==T,min(sampleN,length(inds)),floor(sampleprop*length(inds))) #number of individuals to sample
  samp_inds <- inds[sample(1:length(inds),samp_n,replace=F)] #randomly sample the pop
  
  tempbrd <- data.frame(breeding[samp_inds,])
  tempbrddist <- dist(tempbrd[1:2])
  df2brd <- melt(as.matrix(tempbrddist), varnames = c("row", "col"))
  
  df2brd$Sample <- "Sampled individuals"
  df2brd <- subset(df2brd, value > 0)
  
  bbbox <- bb[[k]]
  bigxcoords <- c(bbbox[1], bbbox[1], bbbox[2], bbbox[2], bbbox[1])
  bigycoords <- c(bbbox[3],bbbox[4],bbbox[4],bbbox[3],bbbox[3])
  bigboxdata <- data.frame(x=bigxcoords,y=bigycoords)
  bigboxpoly <- Polygon(bigboxdata)
  bigboxlist <- Polygons(list(bigboxpoly),1)
  bigboxplot <- SpatialPolygons(list(bigboxlist))
  try(coordinates(breeding)<-c("x","y"), silent=T)
  biginds <- as.numeric(which(over(breeding,bigboxplot)==1))
  big_brd <- data.frame(breeding[biginds,])
  big_brd_list[[k]] <- big_brd
  big_win <- data.frame(nonbreeding[biginds,])
  bigboxbreed <- dist(big_brd[1:2]) 
  bigboxwint <- dist(big_win[1:2])
  
  dfbrd <- melt(as.matrix(bigboxbreed), varnames = c("row", "col"))
  dfbrd$Sample <- "Whole population"
  df2brd$Sample <- "Sampled individuals"
  histbrd <- rbind(dfbrd, df2brd)
  histbrd <- subset(histbrd,value > 0)
  histbrd$Season <- "Breeding" 
  
  dfbrd <- subset(dfbrd, value > 0)
  histbrd <- rbind(dfbrd, df2brd)
  
  dfwin <- melt(as.matrix(bigboxwint), varnames = c("row", "col"))
  tempwin <- data.frame(nonbreeding[samp_inds,])
  tempwindist <- dist(tempwin[1:2])
  df2win <- melt(as.matrix(tempwindist), varnames = c("row", "col"))
  
  dfbrd$Sample <- "Whole population"
  df2brd$Sample <- "Sampled individuals"
  histbrd <- rbind(dfbrd, df2brd)
  histbrd <- subset(histbrd,value > 0)
  histbrd$Season <- "Breeding" 
  
  dfwin$Sample <- "Whole population"
  df2win$Sample <- "Sampled individuals"
  histwin <- rbind(dfwin, df2win)
  histwin <- subset(histwin,value > 0)
  histwin$Season <- "Non-breeding" 
  
  mu <- ddply(histbrd, c("Sample"), summarise, Median=median(value))
  mu2 <- ddply(histwin, c("Sample"), summarise, Median=median(value))
  
  listhistbrd[[k]] <- histbrd
  listhistwin[[k]] <- histwin
  mulist[[k]] <- mu
  mu2list[[k]] <- mu2
  
}

p1 <- ggplot(dat, aes(x=x, y=y)) + 
  geom_point(data=data.frame(big_brd_list[[1]]), size=0.1) +
  geom_polygon(data=lowspread, aes(x=x,y=y, group=Box, fill="#E69F00", alpha=0.1)) +
  theme_void() +
  theme(legend.position = "none") +
  scale_fill_identity() +
  xlim(c(0,1)) +
  ylim(c(1.33,2))

p2 <- ggplot(dat, aes(x=x, y=y)) + 
  geom_point(data=data.frame(big_brd_list[[2]]), size=0.1) +
  geom_polygon(data=medspread, aes(x=x,y=y, group=Box, fill="#E69F00", alpha=0.1)) +
  theme_void() +
  theme(legend.position = "none") +
  scale_fill_identity() +
  xlim(c(0,1)) +
  ylim(c(1.33,2))

p3 <- ggplot(dat, aes(x=x, y=y)) + 
  geom_point(data=data.frame(big_brd_list[[3]]), size=0.1) +
  geom_polygon(data=highspread, aes(x=x,y=y, group=Box, fill="#E69F00", alpha=0.1)) +
  theme_void() +
  theme(legend.position = "none") +
  scale_fill_identity() +
  xlim(c(0,1)) +
  ylim(c(1.33,2))

p4 <- ggplot(listhistbrd[[1]], aes(x=value)) + 
  geom_density(data=listhistbrd[[1]], aes(x=value, linetype=Sample, colour="#E69F00"),size=1) +
  geom_density(data=listhistwin[[1]], aes(x=value, linetype=Sample, colour="#56B4E9"), size=1) +
  scale_colour_manual(values=c("#56B4E9", "#E69F00"), name="Scenario", labels=c("Breeding", "Non-breeding")) +
  geom_vline(data=mulist[[1]], aes(xintercept=Median, linetype=Sample, colour="#E69F00")) +
  geom_vline(data=mu2list[[1]], aes(xintercept=Median, linetype=Sample, colour="#56B4E9")) +
  theme_minimal() +
  labs(x = "Pairwise distance", y = "Percentage of sample") +
  theme(legend.position="none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text = element_text(size= 20), axis.title = element_text(size=20), legend.text=element_text(size=12),legend.title=element_text(size=12)) +
  xlim(c(0,2.5))+
  ylim(c(0,5))

p5 <- ggplot(listhistbrd[[2]], aes(x=value)) + 
  geom_density(data=listhistbrd[[2]], aes(x=value, linetype=Sample, colour="#E69F00"),size=1,show_guide=FALSE) +
  geom_density(data=listhistwin[[2]], aes(x=value, linetype=Sample, colour="#56B4E9"), size=1,show_guide=FALSE) +
  stat_density(data=listhistbrd[[2]], aes(x=value, linetype=Sample, colour="#E69F00"),
               geom="line",position="identity") +
  stat_density(data=listhistwin[[2]], aes(x=value, linetype=Sample, colour="#56B4E9"),
               geom="line",position="identity") +
  scale_colour_manual(values=c("#56B4E9", "#E69F00"), name="Scenario", labels=c("Breeding", "Non-breeding")) +
  geom_vline(data=mulist[[2]], aes(xintercept=Median, linetype=Sample, colour="#E69F00")) +
  geom_vline(data=mu2list[[2]], aes(xintercept=Median, linetype=Sample, colour="#56B4E9")) +
  theme_minimal() +
  labs(x = "Pairwise distance", y = "Percentage of sample") +
  theme(legend.position="bottom", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text = element_text(size= 20), axis.title = element_text(size=20), legend.text=element_text(size=20),legend.title=element_text(size=20, face="bold")) +
  xlim(c(0,2.5))+
  ylim(c(0,5))

p6 <- ggplot(listhistbrd[[3]], aes(x=value)) + 
  geom_density(data=listhistbrd[[3]], aes(x=value, linetype=Sample, colour="#E69F00"),size=1) +
  geom_density(data=listhistwin[[3]], aes(x=value, linetype=Sample, colour="#56B4E9"), size=1) +
  scale_colour_manual(values=c("#56B4E9", "#E69F00"), name="Scenario", labels=c("Breeding", "Non-breeding")) +
  geom_vline(data=mulist[[3]], aes(xintercept=Median, linetype=Sample, colour="#E69F00")) +
  geom_vline(data=mu2list[[3]], aes(xintercept=Median, linetype=Sample, colour="#56B4E9")) +
  theme_minimal() +
  labs(x = "Pairwise distance", y = "Percentage of sample") +
  theme(legend.position="none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text = element_text(size= 20), axis.title = element_text(size=20), legend.text=element_text(size=12),legend.title=element_text(size=12)) +
  xlim(c(0,2.5)) +
  ylim(c(0,5))

spreads <- subset(Spread, MC == "1")  
spreadssum <- summaryBy(Mantel ~ Spread, spreads, FUN = c(mean, sd))

p10 <- ggplot(spreads, aes(x=Spread, y=Mantel, col=Spread)) +
  geom_hline(yintercept=0, linetype="dashed", 
             color = "black", size=1) +
  geom_jitter(aes(colour = "black"), show.legend = T, alpha = 0.1, position=position_jitterdodge()) +
  geom_point(data=spreadssum, aes(x = Spread, y = Mantel.mean, fill="Sampled individuals"), size=2.5, position=position_dodge(width=0.75)) +
  geom_errorbar(data = spreadssum, mapping = aes(x = Spread, y = Mantel.mean, ymin = Mantel.mean - Mantel.sd, ymax = Mantel.mean + Mantel.sd), size=1, width=.2, position=position_dodge(width=0.75)) +
  geom_point(data=subset(areascores,MC==1), aes(x = Spread, y = Mantel, fill="Whole population"), size=5, stroke=1.5, position=position_dodge(width=0.75), shape=4) +
  labs(x="Spread", y = "Mantel score") +
  theme_pubr() +
  theme(legend.position=c(.85,.95), axis.text = element_text(size= 22), axis.title = element_text(size=27), legend.text=element_text(size=15),legend.title=element_text(size=15)) +
  scale_x_discrete(labels=c('Low', 'Medium', 'High')) +
  scale_color_manual(values=c('black','#1B9E77',"#1B9E77","#1B9E77"), labels=c('Low', 'Medium', 'High')) +
  coord_cartesian(ylim = c(-.1, 1)) +
  guides(colour=FALSE, linetype=FALSE, fill=guide_legend(title="", override.aes = list(shape=c(16,4))))

plot <- ggdraw() +
  draw_plot(p1, x = 0.14, y = 0.77, width = .16, height = .2) +
  draw_plot(p2, x = .48, y = 0.77, width = .16, height = .2) +
  draw_plot(p3, x = .8, y = 0.77, width = .16, height = .2) +
  draw_plot(p4, x = 0, y = 0.555, width = .33, height = .445) +
  draw_plot(p5, x = .33, y = 0.5, width = .33, height = .5) +
  draw_plot(p6, x = .66, y = 0.555, width = .33, height = .445) +
  draw_plot(p10, x = 0, y = 0, width = 1, height = .5) +
  draw_plot_label(label = c("A", "B", "C", "D"), size = 20, x = c(0, .33,.66, 0), y = c(1, 1,1,0.5))
plot
#ggsave("Figure_SM4.pdf", plot = plot, dpi = 600, width=22)

# 6.13 Figure SM5 ----

lowspread <- read.csv(file.choose()) # 'lowspreadboxes.csv'
medspread <- read.csv(file.choose())  # 'medspreadboxes.csv'
highspread <- read.csv(file.choose()) # 'highspreadboxes.csv'

xdim <- 100 # x dimension of space
ydim <- xdim*2 # y dimension of space
grid <- raster(ncols=ydim,nrows=ydim) # create a grid of cells in raster format
extent(grid) <- c(0,1,0,ydim/xdim) # set up extent of the grid object
grid[]<-1 
xy <- data.frame(xyFromCell(grid,1:ncell(grid))) #extract coords of cell centres

brd_latrange <- c(quantile(xy$y,0.67),quantile(xy$y,1)) 
brd_lonrange <- c(quantile(xy$x,0),quantile(xy$x,1)) 
bigboxsize <- seq(0.5,1.5,length.out=3)

rows <- sample(nrow(brd_locs), 1000)
brd_locs <- breeding
win_locs <- nonbreeding

breeding <- brd_locs[rows, ]
breeding$rows <- rownames(breeding)

nonbreeding <- win_locs[rows, ]
nonbreeding$rows <- rownames(nonbreeding)

x <- c(0,1)
y<- c(0, 2)
dat <- data.frame(x,y)

counter <- 0 
count <- 0
pls <- list()

listhistbrd <- list()
listhistwin <- list()
muAlist <- list()
mu2Alist <- list()
muBlist <- list()
mu2Blist <- list()

bb1 <-  c(0.315925, 0.684075, 1.545161, 1.78814)
# big box 2 (xmin, xmax, ymin, ymax)
bb2 <-  c(0.19155, 0.80845, 1.463073, 1.870227)
# big box 3 (xmin, xmax, ymin, ymax)
bb3 <-  c(0.067175, 0.932825, 1.380986, 1.952315)
bb <- list(bb1,bb2,bb3)
big_brd_list <- list()
medians <- data.frame(matrix(nrow=404,ncol=4)) # produces a data frame for the calculated mantel stats to be stored in. 4= number of variations in sampling spread/size/number
names(medians) <- c("Sample","Season","Median", "Divide")
medianslist <- list()
meanmedianslist <- list()

for (k in 1:3){
  
  bigbox <- bigboxsize[k] # cycle through 3 options
  bigboxdims <- c(bigbox*diff(brd_lonrange),bigbox*diff(brd_latrange))
  midbrd <- c(mean(brd_lonrange),mean(brd_latrange))#midpoint of breeding range
  area1 <- c(midbrd[1]-(0.5*bigboxdims[1]),midbrd[1]+(0.5*bigboxdims[1]),
             midbrd[2]-(0.5*bigboxdims[2]),midbrd[2]+(0.5*bigboxdims[2]))
  
  nboxes <- sqrt(9) # 9 = number of sampled boxes
  smallbox <- 0.12
  smallboxdims <- c(smallbox*diff(brd_lonrange),smallbox*diff(brd_latrange))
  xcoords <- seq(area1[1],area1[2],length.out=nboxes+2)
  ycoords <- seq(area1[3],area1[4],length.out=nboxes+2)
  xcoords <- xcoords[2:(length(xcoords)-1)]
  ycoords <- ycoords[2:(length(ycoords)-1)]
  
  for(i in 1:length(xcoords)){
    for(j in 1:length(ycoords)){
      counter<-counter+1
      midpoint <- c(xcoords[i],ycoords[j])
      boxcorners <- c(midpoint[1]-(0.5*smallboxdims[1]),midpoint[1]+(0.5*smallboxdims[1]),
                      midpoint[2]-(0.5*smallboxdims[2]),midpoint[2]+(0.5*smallboxdims[2]))
      xcs <- c(boxcorners[1],boxcorners[1],boxcorners[2],boxcorners[2],boxcorners[1])
      ycs <- c(boxcorners[3],boxcorners[4],boxcorners[4],boxcorners[3],boxcorners[3])
      xym <- data.frame(x=xcs,y=ycs)
      p <- Polygon(xym)
      pls[[counter]]<-p
    }
  }
  ps <- Polygons(pls,1)
  area <- SpatialPolygons(list(ps))
  pls <- list()
  counter <- 0 # reset to zero
  count <- 0 # reset to zero
  
  ##### Sample populations and calculate mantel statistics
  fixN <- T 
  try(coordinates(brd_locs)<-c("x","y"), silent=T)
  sampleN <- 200
  
  # individuals for spread scenario
  inds <- as.numeric(which(over(brd_locs,area)==1)) #pick individuals in the study area using over() from sp package
  samp_n<- ifelse(fixN==T,min(sampleN,length(inds)),floor(sampleprop*length(inds))) #number of individuals to sample
  
  
  for(z in 1:100){
    samp_inds <- inds[sample(1:length(inds),samp_n,replace=F)] #randomly sample the pop
    
    tempbrd <- data.frame(brd_locs[samp_inds,])
    
    tempbrd$dist1 <- abs(tempbrd$x - quantile(tempbrd$x,0.33))
    tempbrd$dist2 <- abs(tempbrd$x - quantile(tempbrd$x,0.66))
    df2brd <- tempbrd
    df2brd$Sample <- "Sampled\nindividuals"
    
    tempwin <- data.frame(win_locs[samp_inds,])
    
    tempwin$dist1 <- abs(tempwin$lon - quantile(tempwin$lon,0.33))
    tempwin$dist2 <- abs(tempwin$lon - quantile(tempwin$lon,0.66))
    df2win <- tempwin
    df2win$Sample <- "Sampled\nindividuals"
    
    medians[z,3] <- median(df2brd$dist1)
    medians[z+100,3] <- median(df2brd$dist2)
    medians[z+200,3] <- median(df2win$dist1)
    medians[z+300,3] <- median(df2win$dist2)
    medians[1:400,1] <- "Sampled\nindividuals"
    medians[1:200,2] <- "Breeding"
    medians[201:400,2] <- "Non-Breeding"
    medians[z,4] <- '1'
    medians[z+100,4] <- '2' 
    medians[z+200,4] <- '1'
    medians[z+300,4] <- '2'
    
  }
  
  bbbox <- bb[[k]]
  bigxcoords <- c(bbbox[1], bbbox[1], bbbox[2], bbbox[2], bbbox[1])
  bigycoords <- c(bbbox[3],bbbox[4],bbbox[4],bbbox[3],bbbox[3])
  bigboxdata <- data.frame(x=bigxcoords,y=bigycoords)
  bigboxpoly <- Polygon(bigboxdata)
  bigboxlist <- Polygons(list(bigboxpoly),1)
  bigboxplot <- SpatialPolygons(list(bigboxlist))
  try(coordinates(brd_locs)<-c("x","y"), silent=T)
  biginds <- as.numeric(which(over(brd_locs,bigboxplot)==1))
  big_brd <- data.frame(brd_locs[biginds,])
  big_brd_list[[k]] <- big_brd
  
  big_brd$dist1 <- abs(big_brd$x - quantile(big_brd$x,0.33))
  big_brd$dist2 <- abs(big_brd$x - quantile(big_brd$x,0.66))
  
  big_win <- data.frame(win_locs[biginds,])
  big_win$dist1 <- abs(big_win$lon - quantile(big_win$lon,0.33))
  big_win$dist2 <- abs(big_win$lon - quantile(big_win$lon,0.66))
  
  medians[401,3] <- median(big_brd$dist1)
  medians[402,3] <- median(big_brd$dist2)
  medians[403,3] <- median(big_win$dist1)
  medians[404,3] <- median(big_win$dist2)
  medians[401:404,1] <- "Whole\npopulation"
  medians[401:402,2] <- "Breeding"
  medians[403:404,2] <- "Non-Breeding"
  medians[401,4] <- '1'
  medians[402,4] <- '2'
  medians[403,4] <- '1'
  medians[404,4] <- '2'
  
  medianslist[[k]] <- medians
  
  meanmedians <- summaryBy(Median ~ Sample + Season + Divide, medians[1:400,], FUN = c(mean, sd))
  meanmedianslist[[k]] <- meanmedians
}

a <- medianslist[[1]]
b <- medianslist[[2]]
c <- medianslist[[3]]

a2 <- meanmedianslist[[1]]
b2 <- meanmedianslist[[2]]
c2 <- meanmedianslist[[3]]

a3 <- join(a[401:404,],a2, type="full")
a3[5:8,3] <- a3[5:8,5]
b3 <- join(b[401:404,],b2, type="full")
b3[5:8,3] <- b3[5:8,5]
c3 <- join(c[401:404,],c2, type="full")
c3[5:8,3] <- c3[5:8,5]

p1 <- ggplot() +
  geom_line(a3, mapping=aes(x=Sample, y=Median,group=interaction(Divide,Season))) +
  geom_point(data=subset(a, Sample == "Whole\npopulation"), mapping=aes(x=Sample, y=Median, col=Season, shape=Divide), size=5) +
  geom_jitter(data=subset(a, Sample == "Sampled\nindividuals"), aes(x=Sample, y=Median, col=Season, shape=Divide), show.legend = T, alpha = 0.4) +
  geom_point(data=a2, mapping=aes(x=Sample, y=Median.mean, col=Season, shape=Divide), size=5) +
  labs(x="", y = "Median longitudinal distance to region divide\n('Transition distance')") +
  theme_pubr() +
  theme(legend.position='none', axis.text = element_text(size= 22), axis.title = element_text(size=30), legend.text=element_text(size=15),legend.title=element_text(size=15)) +
  scale_color_manual(values=c("#E69F00","#56B4E9")) +
  ylim(c(0,.4))

p2 <- ggplot() +
  geom_line(b3, mapping=aes(x=Sample, y=Median,group=interaction(Divide,Season))) +
  geom_point(data=subset(b, Sample == "Whole\npopulation"), mapping=aes(x=Sample, y=Median, col=Season, shape=Divide), size=5) +
  geom_jitter(data=subset(b, Sample == "Sampled\nindividuals"), aes(x=Sample, y=Median, col=Season, shape=Divide), show.legend = T, alpha = 0.4) +
  geom_point(data=b2, mapping=aes(x=Sample, y=Median.mean, col=Season, shape=Divide), size=5) +
  labs(x="", y = "") +
  theme_pubr() +
  theme(legend.position='bottom', axis.text = element_text(size= 22), axis.title = element_text(size=30), legend.text=element_text(size=22),legend.title=element_text(size=25)) +
  scale_color_manual(values=c("#E69F00", "#56B4E9")) +
  ylim(c(0,.4))

p3 <- ggplot() +
  geom_line(c3, mapping=aes(x=Sample, y=Median,group=interaction(Divide,Season))) +
  geom_point(data=subset(c, Sample == "Whole\npopulation"), mapping=aes(x=Sample, y=Median, col=Season, shape=Divide), size=5) +
  geom_jitter(data=subset(c, Sample == "Sampled\nindividuals"), aes(x=Sample, y=Median, col=Season, shape=Divide), show.legend = T, alpha = 0.4) +
  geom_point(data=c2, mapping=aes(x=Sample, y=Median.mean, col=Season, shape=Divide), size=5) +
  labs(x="", y = "") +
  theme_pubr() +
  theme(legend.position='none', axis.text = element_text(size= 22), axis.title = element_text(size=30), legend.text=element_text(size=15),legend.title=element_text(size=15)) +
  scale_color_manual(values=c("#E69F00","#56B4E9")) +
  ylim(c(0,.4))



plot <- ggdraw() +
  draw_plot(p1, x = 0, y = 0.05, width = .33, height = .95) +
  draw_plot(p2, x = 0.33, y = 0, width = .33, height = 1) +
  draw_plot(p3, x = 0.66, y = 0.05, width = .33, height = .95) +
  draw_plot_label(label = c("A", "B", "C"), size = 20, x = c(0, .33,.66), y = c(1, 1,1)) +
  draw_plot_label(label = c("Low spread", "Medium spread", "High spread"), size = 25, x = c(0.12, .42,.77), y = c(1, 1,1))

plot
#ggsave("Figure_SM5.pdf", plot = plot, dpi = 600, width=22)