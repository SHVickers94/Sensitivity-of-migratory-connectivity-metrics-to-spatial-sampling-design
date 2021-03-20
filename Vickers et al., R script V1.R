#### Scripts for Vickers et al. Sensitivity of migratory connectivity metrics to spatial sampling design ####

# This script file simulates breeding and non-breeding locations of a migratory popualiton with three different strengths of migratory connectiivty
# This is done for a single contiguous populaiton scenario and a patchy populaiton scenario
# Several spatial sampling designs are then applied to assess potential bias

# Certain elements of this script were conducted on a high permance cluster. 
# As such, these elements may take considerable time to run.

# Code for individual figures is not included but are available on request.

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
if (!require('ggspatial')) install.packages('ggspatial'); library(ggspatial)
if (!require('viridis')) install.packages('viridis'); library(viridis)
if (!require('ade4')) install.packages('ade4'); library(ade4)

if (!require('cruts')) install.packages('cruts'); library(cruts)
if (!require('rasterVis')) install.packages('rasterVis'); library(rasterVis)
if (!require('ncdf4')) install.packages('ncdf4'); library(ncdf4)
if (!require('tictoc')) install.packages('tictoc'); library(tictoc)
if (!require('sf')) install.packages('sf'); library(sf)
if (!require('concaveman')) install.packages('concaveman'); library(concaveman)
if (!require('plotrix')) install.packages('plotrix'); library(plotrix)
if (!require('remotes')) install.packages('remotes'); library(remotes)
if (!require('RFunctionsSN')) install.packages('RFunctionsSN'); library(RFunctionsSN)
if (!require('RColorBrewer')) install.packages('RColorBrewer'); library(RColorBrewer)
if (!require('scales')) install.packages('scales'); library(scales)
if (!require('gridExtra')) install.packages('gridExtra'); library(gridExtra)
if (!require('rgeos')) install.packages('rgeos'); library(rgeos)
if (!require('NLMR')) install.packages('NLMR'); library(NLMR)
if (!require('rnaturalearth')) install.packages('rnaturalearth'); library(rnaturalearth)
if (!require('rnaturalearthdata')) install.packages('rnaturalearthdata'); library(rnaturalearthdata)

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
  data_sheet[1,3] <- n2
  
  # Sampling the population
  
  # Delimit sampling areas (breeding range):
  
  for(t in 1:3){ # loop through spread scenarios
    
    # Create Big box - whole study area (no sampling outside this)
    bigbox <- bigboxsize[t] # sensible range 0.5 to 1.5. This would be added and small box fixed to allow big box to change but fix small box size.
    bigboxdims <- c(bigbox*diff(brd_lonrange),bigbox*diff(brd_latrange))
    midbrd <- c(mean(brd_lonrange),mean(brd_latrange))#midpoint of breeding range
    midwin <- c(mean(win_lonrange),mean(win_latrange))
    area1 <- c(midbrd[1]-(0.5*bigboxdims[1]),midbrd[1]+(0.5*bigboxdims[1]),
               midbrd[2]-(0.5*bigboxdims[2]),midbrd[2]+(0.5*bigboxdims[2]))
    
    # #### Mantel methodlogy for inds within big box extent
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
    
    # individuals for spread scenario
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
    data_sheet[q,4] <- Spread[t]
    data_sheet[q,5] <- Area[t]
    
    # individuals for area scenario
    bigsamp_n<- ifelse(fixN==T,min(sampleN,length(biginds)),floor(sampleprop*length(biginds))) #number of individuals to sample
    bigsamp_inds <- biginds[sample(1:length(biginds),bigsamp_n,replace=F)] #randomly sample the pop
    bigsamp_brd <- data.frame(brd_locs[bigsamp_inds,])
    bigsamp_win <- data.frame(win_locs[bigsamp_inds,]) 
    brddistssubset <- dist(bigsamp_brd[1:2]) # 
    windistssubset <- dist(bigsamp_win[1:2])
    mantelsubset <-  mantel.rtest(brddistssubset,windistssubset,nrepet=1)
    data_sheet[f,1] <- mantelsubset$obs #extract mantel statistic 
    data_sheet[f,3] <- length(bigsamp_inds)
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

data_sheet <- data.frame(matrix(nrow=6,ncol=3))
names(data_sheet) <- c("Mantel_MC1", "Mantel_MC2", "Mantel_MC3")

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
  data_sheet[1,3] <- n2
  
  # Sampling the population
  
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
    
    # individuals for spread scenario
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
    data_sheet[t+1,5] <- Area[t]
    q <- q+1
  }
  data_sheet_list[[g]] <- data_sheet  
}

# 100 replicates were completed and the data compiled into a single .csv file available in the file depository

# Section 6: Realistic population scenario ----
# Produce range map and patchy population density

# One species is done at a time, and is done for three species range maps. Relevant .shp files in repository. Process shown for just one species here and must be repeated and results produced are stored at the end of the process for inclusion in final figure.

blr2 <- st_read("Passerculus henslowii.shp") 

summerrangemap <- subset(blr2, SEASONA!=3)
winterrangemap <- subset(blr2, SEASONA!=2)

summerrangemap2 <- st_combine(summerrangemap)
winterrangemap2 <- st_combine(winterrangemap)

world <- ne_countries(scale = "medium", returnclass = "sf")

r <- nlm_gaussianfield(2000, 2000, resolution = 0.000001, autocorr_range = 10,  mag_var = 100, nug = 10, mean = .001,rescale=T)

# Create summer populaton
crs(r) <- crs(summerrangemap2)
r2 <- projectRaster(r, crs = crs(summerrangemap2))

extent(r2) <- c(st_bbox(summerrangemap2)[1],st_bbox(summerrangemap2)[3],st_bbox(summerrangemap2)[2], st_bbox(summerrangemap2)[4])

sum1 <- as(summerrangemap2,"Spatial")
sum2 <- mask(r2, sum1)

(f1 <- ggplot() +  
    geom_sf(data=world) +
    layer_spatial(sum2) +
    scale_fill_viridis_c(na.value = NA) +
    #geom_sf(data=summerrangemap2,aes(fill="#F8766D"), col="#F8766D", alpha=0.6) +
    #geom_sf(data=winterrangemap2,aes(fill="#F8766D"), col="blue", fill="blue", alpha=0.6) +
    coord_sf(xlim = c(-18,34), ylim = c(10,57), expand = FALSE) +
    theme_void() +
    theme(panel.border = element_rect(colour = "black", fill=NA, size=1), legend.position = "right", legend.text = element_text(size = 25)) +
    labs(fill = "", colour = ""))

p <- rasterToPoints(sum2)
p2 <- data.frame(p)

# buld in values below 0.2 being 0 (increases patchiness of population)
p2$layer <- ifelse(p2$layer < 0.5, 0, p2$layer)

samp_idx <- sample(seq_len(nrow(p2)), 50000, prob=p2$layer)
p3 <- p2[samp_idx,]

ggplot(p3, aes(x=x, y=y)) +
  geom_point(size=0.01)

(f2 <- ggplot() +  
    geom_sf(data=world) +
    geom_point(data=p3, aes(x=x, y=y), size=.5, colour="#E69F00") +
    #geom_sf(data=summerrangemap2,aes(fill="#F8766D"), col="#F8766D", alpha=0.6) +
    #geom_sf(data=winterrangemap2,aes(fill="#F8766D"), col="blue", fill="blue", alpha=0.6) +
    coord_sf(xlim = c(-18,34), ylim = c(10,57), expand = FALSE) +
    theme_void() +
    theme(panel.border = element_rect(colour = "black", fill=NA, size=1), legend.position = "right", legend.text = element_text(size = 25)) +
    labs(fill = "", colour = ""))

#### Create winter population

crs(r) <- crs(winterrangemap2)
r2 <- projectRaster(r, crs = crs(winterrangemap2))

extent(r2) <- c(st_bbox(winterrangemap2)[1],st_bbox(winterrangemap2)[3],st_bbox(winterrangemap2)[2], st_bbox(winterrangemap2)[4])
win1 <- as(winterrangemap2,"Spatial")
win2 <- mask(r2, win1)

q <- rasterToPoints(win2)
q2 <- data.frame(q)

q2$layer <- ifelse(q2$layer < 0.5, 0, q2$layer)

samp_idx <- sample(seq_len(nrow(q2)), 50000, prob=q2$layer)
q3 <- q2[samp_idx,]

ggplot(q3, aes(x=x, y=y)) +
  geom_point(size=0.01)

# Create varying levels of migrtory connectivity

p3$Season <- 'Breeding'
q3$Season <- 'Non-breeding'

p3 <- p3[order(p3$y),]
p3$yrank <- seq_len(nrow(p3))
p3 <- p3[order(p3$x),]
p3$xrank <- seq_len(nrow(p3))

q3 <- q3[order(q3$y),]
q3$yrank <- seq_len(nrow(q3))
q3 <- q3[order(q3$x),]
q3$xrank <- seq_len(nrow(q3))

q3.4 <- data.frame()
q3.0 <- q3
p3$Ind <- seq_len(nrow(p3))
p3rownames <- sample(rownames(data.frame(p3)), nrow(p3))
xspan <- 25000

# spans of 1000, 13,000, 25,000

for (i in 1:nrow(p3)){
  ind <- p3[rownames(data.frame(p3)) == p3rownames[i],]
  
  q3.0$xdiff <- abs(q3.0$xrank - ind$xrank)
  q3.0$ydiff <- abs(q3.0$yrank - ind$yrank)
  q3.0$diff <- q3.0$xdiff + q3.0$ydiff
  q3.0$avdiff <- (q3.0$xdiff + q3.0$ydiff)/2
  
  q3.1 <- q3.0[q3.0$xdiff < xspan,]
  ifelse(nrow(q3.1)==0, try(q3.1 <- q3.0[q3.0$diff == min(q3.0$diff),]),try(q3.1 <- q3.0[q3.0$xdiff < xspan,]))
  
  q3.2 <- q3.1[sample(nrow(q3.1), 1),] 
  
  q3.3 <- q3.2[sample(nrow(q3.2), 1),]
  q3.3$Ind <- ind$Ind
  q3.4 <- rbind(q3.4, q3.3)
  
  q3.0 <- q3.0[!rownames(q3.0) == rownames(q3.3), ]
}

q3.4 <- q3.4[order(q3.4$Ind),]

rows <- sample(seq(1:50000), 10000)
rownames(q3.4) <- seq(1:50000)
rownames(p3) <- seq(1:50000)

p3.0.1 <- p3[rows,]
q3.4.1 <- q3.4[rows,]

brddists <- dist(data.frame(p3.0.1[1:2]))  
windists <- dist(data.frame(q3.4.1[1:2]))

wholepop <- mantel.rtest(brddists,windists,nrepet=1)

test <- join(data.frame(p3),q3.4, by="Ind")

names(test)[1:6] <- c("long","lat","val","Season2","yrank2","xrank2")

(f3 <- ggplot() +
    geom_sf(data=world) +
    geom_segment(data=test, aes(x = x, y = y, xend = long, yend = lat, group=1), size=0.2) +
    geom_point(data=test, aes(x=x, y=y), colour="#56B4E9", size=2) +
    geom_point(data=test, aes(x=long, y=lat), colour="#E69F00", size=2) +
    coord_sf(xlim = c(-18,34), ylim = c(10,57), expand = FALSE) +
    theme_void())

# Create levels of mig connectivity: weakest connectivty

q4 <- q3[sample(row.names(q3), nrow(q3)),]

brddists <- dist(data.frame(p3[1:2]))  
windists2 <- dist(data.frame(q4[1:2]))
mantel.rtest(brddists,windists2,nrepet=1)

test2 <-  cbind(p3,q4[1:5])
test2$row <- seq_len(nrow(test2))
# Sample the population

samplingrast <- raster(nrow=50,ncol=50) # 50x50 ok for 50k inds? 30x30 for 10k

crs(samplingrast) <- crs(summerrangemap2)
samplingrast2 <- projectRaster(samplingrast, crs = crs(summerrangemap2))

extent(samplingrast2) <- c(st_bbox(summerrangemap2)[1],st_bbox(summerrangemap2)[3],st_bbox(summerrangemap2)[2], st_bbox(summerrangemap2)[4])

pointcount = function(r, pts){
  # make a raster of zeroes like the input
  r2 = r
  r2[] = 0
  # get the cell index for each point and make a table:
  counts = table(cellFromXY(r,pts))
  # fill in the raster with the counts from the cell index:
  r2[as.numeric(names(counts))] = counts
  return(r2)
}

r2 <- pointcount(samplingrast2, p3)
r2.5 <- r2
r2.5[r2.5 == 0] <- NA

(f4 <- ggplot() +  
    geom_sf(data=world) +
    layer_spatial(r2.5) +
    scale_fill_viridis_c(na.value = NA) +
    #geom_sf(data=summerrangemap2,aes(fill="#F8766D"), col="#F8766D", alpha=0.6) +
    #geom_sf(data=winterrangemap2,aes(fill="#F8766D"), col="blue", fill="blue", alpha=0.6) +
    coord_sf(xlim = c(-18,34), ylim = c(10,57), expand = FALSE) +
    theme_void() +
    theme(panel.border = element_rect(colour = "black", fill=NA, size=1), legend.position = "right", legend.text = element_text(size = 25)) +
    labs(fill = "", colour = ""))

r3 <- as.data.frame(r2)

minval <- min(r3[order(-r3$layer)[1:20],])

r3[-((order(-r3$layer)[1:20])),] <- NA

r4 <- r2

values(r4)[values(r4) < minval] = NA
values(r4) <- r3$layer

(f5 <- ggplot() +  
    geom_sf(data=world) +
    layer_spatial(r4) +
    scale_fill_viridis_c(na.value = NA) +
    geom_sf(data=summerrangemap2, col="#E69F00", alpha=0) +
    #geom_sf(data=winterrangemap2,aes(fill="#F8766D"), col="blue", fill="blue", alpha=0.6) +
    coord_sf(xlim = c(-18,34), ylim = c(10,57), expand = FALSE) +
    theme_void() +
    theme(panel.border = element_rect(colour = "black", fill=NA, size=1), legend.position = "right", legend.text = element_text(size = 25)) +
    labs(fill = "", colour = ""))


r5 <- rasterToPolygons(r4)
r6 <- fortify(r5, region='layer')

(f5 <- ggplot() +  
    geom_sf(data=world) +
    geom_polygon(data=r5, aes(x = long, y = lat, group = id), col="black", fill="#E69F00", alpha=1) +
    geom_sf(data=summerrangemap2, col="#E69F00", alpha=0) +
    #geom_sf(data=winterrangemap2,aes(fill="#56B4E9"), col="#56B4E9", fill="#56B4E9", alpha=0.6) +
    coord_sf(xlim = c(-99,-72), ylim = c(25,47), expand = FALSE) +
    theme_void() +
    theme(panel.border = element_rect(colour = "black", fill=NA, size=1), legend.position = "right", legend.text = element_text(size = 25)) +
    labs(fill = "", colour = ""))

data <- data.frame(matrix(ncol=36, nrow=100))
names(data) <- c('MC3', 'MC4', 'MC5', 'MC6', 'MC7', 'MC8', 'MC9', 'MC10', 'MC11', 'MC12', 'MC13', 'MC14', 'MC15', 'MC16', 'MC17', 'MC18', 'MC19', 'MC20', 'D3', 'D4', 'D5', 'D6', 'D7', 'D8', 'D9', 'D10', 'D11', 'D12', 'D13', 'D14', 'D15', 'D16', 'D17', 'D18', 'D19', 'D20')

for (t in 3:20){
  
  for (i in 1:100){
    
    plots <- sample(seq(1:20), t)
    
    r6 <- r5[plots,]
    try(coordinates(p3)<-c("x","y"), silent=T)
    crs(p3) <- crs(r5)
    
    inds <- as.numeric(which(over(p3,r6)>1)) 
    
    sampleN <- 200
    fixN <- T
    
    samp_n<- ifelse(fixN==T,min(sampleN,length(inds)),floor(sampleprop*length(inds))) #number of individuals to sample
    samp_inds <- inds[sample(1:length(inds),samp_n,replace=F)]
    
    samp_brd <- data.frame(p3[samp_inds,]) 
    samp_win <- data.frame(q3.4[samp_inds,])

    brddistssubset <- dist(samp_brd[1:2])
    windistssubset <- dist(samp_win[1:2])
    mantelsubset <-  mantel.rtest(brddistssubset,windistssubset,nrepet=1)
    data[i, t-2] <- mantelsubset$obs
    
    centroids <- getSpPPolygonsLabptSlots(r6)
    dst <- pointDistance(centroids, lonlat=TRUE)
    data[i, t+16] <- sum(dst, na.rm=T)/(((t-1)*t)/2)
    
  }
}

# median distance
data2 <- melt(data[1:18])
data3 <- melt(data[19:36])
names(data2) <- c('Sites', 'Mantel')
data2$dist <- data3$value

plot <- ggdraw() +
  draw_plot(f1, x = 0.015, y = .66, width = .5, height = .33) +
  draw_plot(f2, x = .5, y = .66, width = .5, height = .33) +
  draw_plot(f3, x = 0, y = .33, width = .5, height = .33) +
  draw_plot(f4, x = .5, y = .33, width = .5, height = .33) +
  draw_plot(f5, x = 0.25, y = 0, width = .5, height = .33) +
  draw_plot_label(label = c("A", "B", "C", "D","E"), size = 30, x = c(0, .49, 0, .49, 0.25), y = c(1,1,.66,.66,.33))
plot

ggsave("realisticmethods_fig.pdf", plot = plot, dpi = 600, width=25, height=25)

#### AQUATIC WARBLER

#highAQUWA <- data2
#highAQUWAbrddists <- brddists
#highAQUWAwindists <- windists
highAQUWAwholepop <- mantel.rtest(highAQUWAbrddists,highAQUWAwindists,nrepet=1)

#mediumAQUWA <- data2
#mediumAQUWAbrddists <- brddists
#mediumAQUWAwindists <- windists
mediumAQUWAwholepop <- mantel.rtest(mediumAQUWAbrddists,mediumAQUWAwindists,nrepet=1)

#lowAQUWA <- data2
#lowAQUWAbrddists <- brddists
#lowAQUWAwindists <- windists
lowAQUWAwholepop <- mantel.rtest(lowAQUWAbrddists,lowAQUWAwindists,nrepet=1)

(plot1.1 <- ggplot(highAQUWA, aes(y=Mantel, x=Sites, colour=dist)) +
    geom_jitter(size=2, width=.25) +
    geom_hline(yintercept = highAQUWAwholepop$obs, colour="black", size=1.1, linetype="dashed") +
    #geom_smooth(highAQUWA, mapping=aes(y=Mantel, x=as.numeric(Sites)), colour='black') +
    labs(x='Sampling sites', y='Mantel score') +
    scale_x_discrete(limits=levels(lowAQUWA$Sites), breaks=levels(lowAQUWA$Sites)[seq(1,18, by=2)], labels=seq(3,20,2)) +
    theme_pubr() +
    ylim(-0.1,1) +
    #scale_color_viridis() +
    scale_color_gradientn(colours = c("red", "grey", "blue"), name="Mean distance\nbetween sites")+
    theme(legend.position = 'none', axis.title = element_blank(), axis.text = element_blank())) 


(plot2.1 <- ggplot(mediumAQUWA, aes(y=Mantel, x=Sites, colour=dist)) +
    geom_jitter(size=2, width=.25) +
    geom_hline(yintercept = mediumAQUWAwholepop$obs, colour="black", size=1.1, linetype="dashed") +
    #geom_smooth(mediumAQUWA, mapping=aes(y=Mantel, x=as.numeric(Sites)), colour='black') +
    labs(x='Sampling sites', y='Mantel score') +
    scale_x_discrete(limits=levels(lowAQUWA$Sites), breaks=levels(lowAQUWA$Sites)[seq(1,18, by=2)], labels=seq(3,20,2)) +
    theme_pubr() +
    ylim(-0.1,1) +
    #scale_color_viridis() +
    scale_color_gradientn(colours = c("red", "grey", "blue"), name="Mean distance\nbetween sites")+
    theme(legend.position = 'none', axis.title = element_blank(), axis.text = element_blank())) 


(plot3.1 <- ggplot(lowAQUWA, aes(y=Mantel, x=Sites, colour=dist)) +
    geom_jitter(size=2, width=.25) +
    geom_hline(yintercept = lowAQUWAwholepop$obs, colour="black", size=1.1, linetype="dashed") +
    #geom_smooth(lowAQUWA, mapping=aes(y=Mantel, x=as.numeric(Sites)), colour='black') +
    labs(x='Sampling sites', y='Mantel score') +
    scale_x_discrete(limits=levels(lowAQUWA$Sites), breaks=levels(lowAQUWA$Sites)[seq(1,18, by=2)], labels=seq(3,20,2)) +
    theme_pubr() +
    ylim(-0.1,1) +
    #scale_color_viridis() +
    scale_color_gradientn(colours = c("red", "grey", "blue"), name="Mean distance\nbetween sites")+
    theme(legend.position = 'none', axis.title.x = element_text(size=30), axis.title.y = element_blank(), axis.text.x = element_text(size=15), axis.text.y = element_blank())) 

#AQUWApolys <- r5
#AQUWAsummap <- summerrangemap2
#AQUWAwinmap <- winterrangemap2

(sampmapAQUWA <- ggplot() +  
    geom_sf(data=world) +
    geom_sf(data=AQUWAsummap, col="#E69F00", col="#E69F00", fill="#E69F00", alpha=0.6) +
    geom_polygon(data=AQUWApolys, aes(x = long, y = lat, group = id), col="black", fill="black", alpha=1) +
    geom_sf(data=AQUWAwinmap,aes(fill="#56B4E9"), col="#56B4E9", fill="#56B4E9", alpha=0.6) +
    coord_sf(xlim = c(-18,34), ylim = c(10,57), expand = FALSE) +
    theme_void() +
    theme(panel.border = element_rect(colour = "black", fill=NA, size=1), legend.position = "right", legend.text = element_text(size = 25)) +
    labs(fill = "", colour = ""))

#### FALCATED DUCK

#highFALDU <- data2
#highFALDUbrddists <- brddists
#highFALDUwindists <- windists
highFALDUwholepop <- mantel.rtest(highFALDUbrddists,highFALDUwindists,nrepet=1)

#mediumFALDU <- data2
#mediumFALDUwindists <- windists
#mediumFALDUbrddists <- brddists
mediumFALDUwholepop <- mantel.rtest(mediumFALDUbrddists,mediumFALDUwindists,nrepet=1)

#lowFALDU <- data2
#lowFALDUbrddists <- brddists
#lowFALDUwindists <- windists
lowFALDUwholepop <- mantel.rtest(lowFALDUbrddists,lowFALDUwindists,nrepet=1)


(plot1.2 <- ggplot(highFALDU, aes(y=Mantel, x=Sites, colour=dist)) +
    geom_jitter(size=2, width=.25) +
    geom_hline(yintercept = highFALDUwholepop$obs, colour="black", size=1.1, linetype="dashed") +
    #geom_smooth(highFALDU, mapping=aes(y=Mantel, x=as.numeric(Sites)), colour='black') +
    labs(x='Sampling sites', y='Mantel score') +
    scale_x_discrete(limits=levels(highFALDU$Sites), breaks=levels(lowAQUWA$Sites)[seq(1,18, by=2)], labels=seq(3,20,2)) +
    theme_pubr() +
    ylim(-0.1,1) +
    #scale_color_viridis() +
    scale_color_gradientn(colours = c("red", "grey", "blue"), name="Mean distance\nbetween sites")+
    theme(legend.position = 'none', axis.title = element_blank(), axis.text = element_blank())) 


(plot2.2 <- ggplot(mediumFALDU, aes(y=Mantel, x=Sites, colour=dist)) +
    geom_jitter(size=2, width=.25) +
    geom_hline(yintercept = mediumFALDUwholepop$obs, colour="black", size=1.1, linetype="dashed") +
    #geom_smooth(mediumFALDU, mapping=aes(y=Mantel, x=as.numeric(Sites)), colour='black') +
    labs(x='Sampling sites', y='Mantel score') +
    scale_x_discrete(limits=levels(mediumFALDU$Sites), breaks=levels(lowAQUWA$Sites)[seq(1,18, by=2)], labels=seq(3,20,2)) +
    theme_pubr() +
    ylim(-0.1,1) +
    #scale_color_viridis() +
    scale_color_gradientn(colours = c("red", "grey", "blue"), name="Mean distance\nbetween sites")+
    theme(legend.position = 'none', axis.title = element_blank(), axis.text = element_blank())) 


(plot3.2 <- ggplot(lowFALDU, aes(y=Mantel, x=Sites, colour=dist)) +
    geom_jitter(size=2, width=.25) +
    geom_hline(yintercept = lowFALDUwholepop$obs, colour="black", size=1.1, linetype="dashed") +
    #geom_smooth(lowFALDU, mapping=aes(y=Mantel, x=as.numeric(Sites)), colour='black') +
    labs(x='Sampling sites', y='Mantel score') +
    scale_x_discrete(limits=levels(lowFALDU$Sites), breaks=levels(lowAQUWA$Sites)[seq(1,18, by=2)], labels=seq(3,20,2)) +
    theme_pubr() +
    ylim(-0.1,1) +
    labs(x="   ") +
    #scale_color_viridis() +
    scale_color_gradientn(colours = c("red", "grey", "blue"), name="Mean distance\nbetween sites")+
    theme(legend.position = 'none', axis.title.y = element_blank(), axis.text.y = element_blank(), axis.title.x = element_text(size=30), axis.text.x = element_text(size=15))) 

#FALDUpolys <- r5
#FALDUsummap <- summerrangemap2
#FALDUwinmap <- winterrangemap2

(sampmapFALDU <- ggplot() +  
    geom_sf(data=world) +
    geom_sf(data=FALDUsummap, col="#E69F00", col="#E69F00", fill="#E69F00", alpha=0.6) +
    geom_polygon(data=FALDUpolys, aes(x = long, y = lat, group = id), col="black", fill="black", alpha=1) +
    geom_sf(data=FALDUwinmap,aes(fill="#56B4E9"), col="#56B4E9", fill="#56B4E9", alpha=0.6) +
    coord_sf(xlim = c(75,165), ylim = c(6,78), expand = FALSE) +
    theme_void() +
    theme(panel.border = element_rect(colour = "black", fill=NA, size=1), legend.position = "right", legend.text = element_text(size = 25)) +
    labs(fill = "", colour = ""))
g
#### HENSLOWS SPARROW

#highHENSP <- data2
#highHENSPbrddists <- brddists
#highHENSPwindists <- windists
highHENSPwholepop <- mantel.rtest(highHENSPbrddists,highHENSPwindists,nrepet=1)

#mediumHENSP <- data2
#mediumHENSPwindists <- windists
#mediumHENSPbrddists <- brddists
mediumHENSPwholepop <- mantel.rtest(mediumHENSPbrddists,mediumHENSPwindists,nrepet=1)

#lowHENSP <- data2
#lowHENSPbrddists <- brddists
#lowHENSPwindists <- windists
lowHENSPwholepop <- mantel.rtest(lowHENSPbrddists,lowHENSPwindists,nrepet=1)

(plot1.0 <- ggplot(highHENSP, aes(y=Mantel, x=Sites, colour=dist)) +
    geom_jitter(size=2, width=.25) +
    geom_hline(yintercept = highHENSPwholepop$obs, colour="black", size=1.1, linetype="dashed") +
    #geom_smooth(highHENSP, mapping=aes(y=Mantel, x=as.numeric(Sites)), colour='black') +
    labs(x='Sampling sites', y='Mantel score') +
    scale_x_discrete(limits=levels(highHENSP$Sites), breaks=levels(lowAQUWA$Sites)[seq(1,18, by=2)], labels=seq(3,20,2)) +
    theme_pubr() +
    ylim(-0.1,1) +
    labs(y=" ") +
    #scale_color_viridis() +
    scale_color_gradientn(colours = c("red", "grey", "blue"), name="Mean distance\nbetween sites")+
    theme(legend.position = 'none', axis.title.x = element_blank(), axis.text.x = element_blank(), axis.title.y = element_text(size=30), axis.text.y = element_text(size=15))) 

(plot2.0 <- ggplot(mediumHENSP, aes(y=Mantel, x=Sites, colour=dist)) +
    geom_jitter(size=2, width=.25) +
    geom_hline(yintercept = mediumHENSPwholepop$obs, colour="black", size=1.1, linetype="dashed") +
    #geom_smooth(mediumHENSP, mapping=aes(y=Mantel, x=as.numeric(Sites)), colour='black') +
    labs(x='Sampling sites', y='Mantel score') +
    scale_x_discrete(limits=levels(mediumHENSP$Sites), breaks=levels(lowAQUWA$Sites)[seq(1,18, by=2)], labels=seq(3,20,2)) +
    theme_pubr() +
    ylim(-0.1,1) +
    #scale_color_viridis() +
    scale_color_gradientn(colours = c("red", "grey", "blue"), name="Mean distance\nbetween sites")+
    theme(legend.position = 'none', axis.title.x = element_blank(), axis.text.x = element_blank(), axis.title.y = element_text(size=30), axis.text.y = element_text(size=15))) 

(plot3.0 <- ggplot(lowHENSP, aes(y=Mantel, x=Sites, colour=dist)) +
    geom_jitter(size=2, width=.25) +
    geom_hline(yintercept = lowHENSPwholepop$obs, colour="black", size=1.1, linetype="dashed") +
    #geom_smooth(lowHENSP, mapping=aes(y=Mantel, x=as.numeric(Sites)), colour='black') +
    labs(x='Sampling sites', y='Mantel score') +
    scale_x_discrete(limits=levels(lowHENSP$Sites), breaks=levels(lowAQUWA$Sites)[seq(1,18, by=2)], labels=seq(3,20,2)) +
    theme_pubr() +
    ylim(-0.1,1) +
    labs(x="   ", y= "  ") +
    #scale_color_viridis() +
    scale_color_gradientn(colours = c("red", "grey", "blue"), name="Mean distance\nbetween sites", breaks=c(min(lowHENSP$dist)+50000,max(lowHENSP$dist)-50000),labels=c("Low","High"))+
    theme(legend.position = 'none', legend.text = element_text(size=15), legend.title=element_text(size=30), axis.title.y = element_text(size=30), axis.text.y = element_text(size=15), axis.title.x = element_text(size=30), axis.text.x = element_text(size=15))) 

legend <- cowplot::get_legend(plot3.0.0)

#HENSPpolys <- r5
#HENSPsummap <- summerrangemap2
#HENSPwinmap <- winterrangemap2

(sampmapHENSP <- ggplot() +  
    geom_sf(data=world) +
    geom_sf(data=HENSPsummap, col="#E69F00", col="#E69F00", fill="#E69F00", alpha=0.6) +
    geom_polygon(data=HENSPpolys, aes(x = long, y = lat, group = id), col="black", fill="black", alpha=1) +
    geom_sf(data=HENSPwinmap,aes(fill="#56B4E9"), col="#56B4E9", fill="#56B4E9", alpha=0.6) +
    coord_sf(xlim = c(-99,-72), ylim = c(24.6,47.4), expand = FALSE) +
    theme_void() +
    theme(panel.border = element_rect(colour = "black", fill=NA, size=1), legend.position = "right", legend.text = element_text(size = 25)) +
    labs(fill = "", colour = ""))

#### COMPILE PLOT

ggarrange(sampmapHENSP, sampmapAQUWA, sampmapFALDU, plot1.0, plot1.1, plot1.2, plot2.0, plot2.1, plot2.2, plot3, plot3.1, plot3.2,ncol = 3, nrow = 4)

plot <- ggdraw() +
  draw_plot(sampmapHENSP, x = 0.07, y = .75, width = .3, height = .25) +
  draw_plot(sampmapAQUWA, x = .39, y = .75, width = .29, height = .25) +
  draw_plot(sampmapFALDU, x = .70, y = .75, width = .29, height = .25) +
  draw_plot(plot1.0, x = 0, y = .53, width = .38, height = .22) +
  draw_plot(plot1.1, x = 0.38, y = .53, width = .31, height = .22) +
  draw_plot(plot1.2, x = 0.69, y = .53, width = .31, height = .22) +
  draw_plot(plot2.0, x = 0, y = .31, width = .38, height = .22) +
  draw_plot(plot2.1, x = 0.38, y = .31, width = .31, height = .22) +
  draw_plot(plot2.2, x = .69, y = .31, width = .31, height = .22) +
  draw_plot(plot3.0, x = 0, y = 0.04, width = .38, height = .27) +
  draw_plot(plot3.1, x = 0.38, y = 0.04, width = .31, height = .27) +
  draw_plot(plot3.2, x = .69, y = 0.04, width = .31, height = .27) +
  draw_plot(legend, x = 0.1, y =0, width=0.3, height = 0.05) +
  draw_plot_label(label = c("A", "B", "C"), size = 30, x = c(0.067, 0.385, 0.695), y = c(.995,.995,.995))

# Section 7: Post analysis data handling ----

# Sampling processes with sections 3–5 are random and as such specific results will vary slightly.
# Specific output files (after compliling) are provided in the file repository for this section.

rm(list = ls())
gc()

## Run code chunks 1. and 2., after this each figure chunk can be ran ##

# 7.1. Spread and Area dataset production ----

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

# 7.2. Sample size dataset production ----

Data <- read.csv(file.choose()) # read in samplingdata.csv file

Meltdatabias <- reshape2::melt(Data[1:9], id='Sampled')

MantelMDB <- subset(Meltdatabias, variable %in% c("Mantel_MC1","Mantel_MC2","Mantel_MC3"))
CohenMDB <- subset(Meltdatabias, variable %in% c("Cohen_MC1","Cohen_MC2","Cohen_MC3"))

MDBa <- summaryBy(value ~ variable + Sampled, MantelMDB, FUN = c(mean, sd))
MDBb <- summaryBy(value ~ variable + Sampled, CohenMDB, FUN = c(mean, sd))