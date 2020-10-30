Hummingbird Movement
================
Gannon et al.
October 2020

## Overview

The final dataset will include the frequency of hummingbird movements
among different RFID readers across the HJA meadows in each of four
years. We combine information on distances between readers, cover-type
in which each reader was placed, and the proportion of intervening
forest.

#### Raw data

The raw data include six columns of data and one row for each record of
a PIT tag as passed through the reader. The fields include:

  - PIT tag ID (character string)

  - Date (yyyy-mm-dd)

  - Time (hh:mm:ss)

  - Seconds past midnight (integer)

  - Meadow complex ![\\in](https://latex.codecogs.com/png.latex?%5Cin
    "\\in") {CM,LOM,M1,M2}

  - Reader position within the complex
    ![\\in](https://latex.codecogs.com/png.latex?%5Cin "\\in")
    {N,S,E,W,C} (sometimes SE or SW, etc.)

![\~](https://latex.codecogs.com/png.latex?~ "~")

**Read in RFID-reader data and reader location data:**

``` r
#read in data
   vis <- read.table(here("Data","all_visits2.txt"), 
                   header = F, sep = "\t", as.is = T)
   names(vis) <- c("bird", "date", "time", "seconds", "complex", "reader")

# remove headquarters test visits
   vis <- vis[-which(vis$complex == "HQ"), ]

# combine C-upper and C-lower with C
  # These readers are in the same meadow, 
  # so we can consider visits as evidence the bird is in that meadow
  vis$reader[which(vis$reader == "CU")] <- "C"
  vis$reader[which(vis$reader == "CL")] <- "C"
  vis$reader[which(vis$reader == "CW")] <- "C"
   
# combine complex and reader to make unique identifier for each reader
  vis$reader <- paste(vis$complex, vis$reader, sep = "_")

# load data on reader locations 
  rlocs <- read.table(here("Data", "hja_reader_locations_centercomb.txt"), 
                      header = T, sep = "\t", as.is = T)
```

![\~](https://latex.codecogs.com/png.latex?~ "~")

Remove any duplicated observations. The data were sometimes downloaded
multiple times and resulted in duplicates.

``` r
  vis <- vis[which(!duplicated(vis[,c(1:3)])), ]
```

**Convert records to movements**

To, create a movement dataset, sort by `bird` then `date` then `time`.
Loop through to ask whether a bird was recorded at reader
![i](https://latex.codecogs.com/png.latex?i "i") at time
![t](https://latex.codecogs.com/png.latex?t "t"), then again at reader
![j](https://latex.codecogs.com/png.latex?j "j") at time
![t'\>t](https://latex.codecogs.com/png.latex?t%27%3Et "t'\>t") within
the same day.

``` r
# order the dataframe
  vis_sort <- vis[order(vis$bird, vis$date, vis$seconds), ]

# add a column for the year
  vis_sort$year <- str_extract(vis_sort$date, "[[:digit:]]{4}")
  
# turn date into day of the year
  vis_sort$date2 <- ymd(vis_sort$date)
  vis_sort$day <- yday(vis_sort$date)

# create variables for use later
  rdrs <- rlocs$reader
  nrdrs <- length(rdrs)
  
# create movement dataset
  birds_all <- unique(vis_sort$bird)
  yrs <- sort(unique(vis_sort$year))
  mv_columns <- c("reader1", "reader2", "year", "count")
  mv <- matrix(nrow = 1, ncol = length(mv_columns))
  colnames(mv) <- mv_columns
  mv <- as.data.frame(mv)

# Fill dataframe with movements
  for(y in yrs){
    df_temp <- subset(vis_sort, year==y)
    df_temp <- df_temp[order(df_temp$bird, df_temp$day, df_temp$seconds), ]
    ds <- 0
  # create a count for each possible movement
  # for each year
    cnts <- matrix(data=0, nrow = nrdrs, ncol = nrdrs)
     if(nrow(df_temp) > 1){
      for(i in 2:nrow(df_temp)){
        if(df_temp$date[i] == df_temp$date[i-1] &
           (df_temp$bird[i] == df_temp$bird[i-1]) &
            (df_temp$reader[i] != df_temp$reader[i-1])){
          #if all conditions hold, find reader numbers
            r1 <- which(rdrs == df_temp$reader[i-1])
            r2 <- which(rdrs == df_temp$reader[i])
          #add a count to that cell  
            cnts[r1,r2] <- cnts[r1,r2] + 1
        }
      }
     }
      # combine with movement dataframe
      if(y=="2014"){
        rems <- which(rlocs$meadow_complex == "CM" | rlocs$meadow_complex == "LOM")
        cnts_y <- cnts[-rems,-rems]
        rdrs_y <- rdrs[-rems]
      } else if(y=="2015"){
        rems <- which(rlocs$meadow_complex == "CM")
        rdrs_y <- rdrs[-rems]
        cnts_y <- cnts[-rems,-rems]
      } else{
        rdrs_y <- rdrs
        cnts_y <- cnts
      }
        df_temp2 <- data.frame(reader1=rep(rdrs_y, length(rdrs_y)),
                               reader2=rep(rdrs_y, each=length(rdrs_y)),
                               year=rep(y, length(rdrs_y)^2),
                               count=as.vector(cnts_y))
        mv <- rbind(mv, df_temp2)
    
  }
  
  mv <- mv[-1,]
```

**Remove rows for ![i \\to
i](https://latex.codecogs.com/png.latex?i%20%5Cto%20i "i \\to i")
movements.**

``` r
  rows_ii <- which(mv$reader1==mv$reader2)

  mv <- mv[-rows_ii, ]
```

### Adding in landscape variables

Now add in spatial data pertaining to each “edge”: distance between
readers, proportion of intervening forest, and cover types connected by
the edge (i.e. meadow-meadow, meadow-mix, meadow-forest, meadow-scrub).

**Edge length**

``` r
 # coordinate matrix for reader locations
  coords_r <- as.matrix(rlocs[,c("x","y")])
  rownames(coords_r) <- rlocs$reader
  ones <- matrix(1, nrow = nrow(coords_r), ncol = 1)

# compute squared euclidean distance matrix
  D_rsq <- diag(coords_r%*%t(coords_r))%*%t(ones) + 
    ones%*%t(diag(coords_r%*%t(coords_r))) - 2*(coords_r%*%t(coords_r))

# square root of values and label the rows and cols
  D_r <- sqrt(D_rsq)
  rownames(D_r) <- rlocs$reader
  colnames(D_r) <- rownames(D_r)
 
# vectorize the distance matrix 
  dist_list <- data.frame(reader1=rep(rdrs,nrdrs),
                          reader2=rep(rdrs, each=nrdrs),
                          dist.km=as.vector(D_r)/1000)
# add distance to movement data
  mv$edgename <- paste(mv$reader1, mv$reader2, sep = "-")
  dist_list$edgename <- paste(dist_list$reader1, dist_list$reader2, sep = "-")
  mv_wdist <- merge(mv, dist_list[,c(3,4)])
```

![\~](https://latex.codecogs.com/png.latex?~ "~")

**Cover types**

``` r
#habitat for reader 1
  mv_wdist <- merge(mv_wdist, rlocs[,c(1,6)], by.x="reader1", by.y="reader")
  names(mv_wdist)[which(names(mv_wdist)=="habitat")] <- "hab1"
  
#habitat for reader 2
  mv_wdist <- merge(mv_wdist, rlocs[,c(1,6)], by.x="reader2", by.y="reader")
  names(mv_wdist)[which(names(mv_wdist)=="habitat")] <- "hab2"
  
# connectivity reader 1
  mv_wdist <- merge(mv_wdist, rlocs[,c(1,10)], by.x="reader1", by.y="reader")
  names(mv_wdist)[which(names(mv_wdist) == "prop_for100")] <- "p_for1"
  
# connectivity reader 2
  mv_wdist <- merge(mv_wdist, rlocs[,c(1,10)], by.x="reader2", by.y="reader")
  names(mv_wdist)[which(names(mv_wdist) == "prop_for100")] <- "p_for2" 
```

![\~](https://latex.codecogs.com/png.latex?~ "~")

**Proportion of forest vs. non-forest cells intersected by each edge**

``` r
# Create points data
 utm_proj <- CRS("+proj=utm +zone=10 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")
 crs <- st_crs(utm_proj)
 
# create list to store start-end coords
  edges_coordlist <- vector(mode = "list")

  for(i in 1:nrdrs){
      for(j in 1:nrdrs){
        coords <- as.matrix(rbind(coords_r[i, ],
                                coords_r[j, ]))
        nm <- paste(rownames(coords_r)[i], 
                    rownames(coords_r)[j],
                    sep = "-")
        edges_coordlist[[nm]] <- coords
      }
    }
  
# convert to lines
  edges_aslines <- map(edges_coordlist, ~Line(.))
 
# Convert the many lines to a single lines object
  edges_ls <- map(1:length(edges_aslines), 
                          ~Lines(edges_aslines[[.]], ID=names(edges_aslines)[.]))
  
# create one spatial lines object
  edges_spl <- SpatialLines(edges_ls, proj4string = utm_proj)
  
# read in forest data
  forest <- raster(here("Data","forest.tif"))
  
# extract proportion of forest between the nodes
  beginCluster(n=2)
    prop_for <- raster::extract(forest, edges_spl, fun=mean, na.rm=T,
                      df=T, buffer=25)
  endCluster()  
  prop_for$ID <- names(edges_aslines)

  
#merge two datasets
   mv_wpropfor <- merge(mv_wdist, prop_for,
                        by.x="edgename",
                        by.y="ID")
# reorder
   cols.order <- c("reader1", "reader2", "year", "hab1", "hab2", 
                   "p_for1", "p_for2", "dist.km", "forest", "count")
   mv_noedgename <- mv_wpropfor[,cols.order]
   
  names(mv_noedgename) <- c("reader_i", "reader_j", "year", "hab_i", "hab_j",
                  "pfor_100m_i", "pfor_100m_j", "dist_ij", "pfor_ij", "count")
  
# now create a dyad identifier
  mv_noedgename$dyad <- NA
  for(i in 1:nrow(mv_noedgename)){
    dyad_i <- sort(c(mv_noedgename$reader_i[i], mv_noedgename$reader_j[i]))
    mv_noedgename$dyad[i] <- paste(dyad_i, collapse = '-')
  }
  
  mv_noedgename <- mv_noedgename[,c(1,2,11,3:10)]
```

**Including an offset for the number of days the readers were maintained
on the landscape for a given year**

``` r
  active_days <- data.frame(year=c(2014:2017),
                            active_days = c(11,64,54,50))

  mv_final <- merge(mv_noedgename, active_days)
```

### Write out the cleaned movement data

``` r
  mv_final <- mv_final[order(mv_final$year, mv_final$reader_i, mv_final$reader_j),]
 
  write.table(mv_final, 
              file = here("Data", "total_movements_by_year.tsv"), 
              row.names = F, 
              sep = "\t",
              quote = F)
```
