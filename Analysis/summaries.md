Pairing automated mark-recapture and social network models to explore
the effects of landscape configuration on hummingbird foraging patterns
- Data summaries
================
D. G. Gannon, A. S. Hadley, S. J. K. Frey
June 2021

### Load banding data

``` r
# load banding data
  caps <- read_csv(
    here("Data", "hja_hummingbird_captures_rfid.csv")
  )

# remove records of returning birds
  caps <- caps[-which(duplicated(caps$rfid)), ]
```

### Read in RFID-reader data

``` r
#read in data
   rds <- read.table(here("Data","all_reads.txt"), 
                   header = T, sep = "\t", as.is = T)

# decompose reader id into a complex id and location
  split <- as.data.frame(str_split(rds$site, "_", simplify = T))
  rds <- cbind(rds[,c(1:4)], split)
  names(rds) <- c("bird", "date", "time", "seconds", "complex", "reader")

# remove split for the sake of memory
  rm(split)

# combine C-upper and C-lower with C
  # These readers are in the same meadow, 
  # so we can consider visits as evidence the bird is in that meadow
  rds$reader[which(rds$reader == "CU")] <- "C"
  rds$reader[which(rds$reader == "CL")] <- "C"
  rds$reader[which(rds$reader == "CW")] <- "C"
  rds$reader[which(rds$reader == "CE")] <- "C"
  
# combine complex and reader to make unique identifier for each reader
  rds$reader <- paste(rds$complex, rds$reader, sep = "_")
  
### Sanity check ###
  #unique(rds$reader)

# order the dataframe
  rds_sort <- rds[order(rds$bird, rds$date, rds$seconds), ]

# add a column for the year
  rds_sort$year <- str_extract(rds_sort$date, "[[:digit:]]{4}")
  
### Sanity check ###
  #unique(rds_sort$year)
  
# turn date into day of the year
  rds_sort$date2 <- ymd(rds_sort$date)
  rds_sort$day <- yday(rds_sort$date)
```

**Load separate visits data**

``` r
# load separate visits data
  vis <-  read.table(file = here(
    "Data",
    "all_visits.txt"
  ), header = T, sep = "\t", as.is = T)

# create complex column
  split2 <- as.data.frame(str_split(vis$site, "_", simplify = T))
  vis <- cbind(vis[,c(1:5)], split2)
  names(vis) <- c("bird", "date", "time", "seconds", "year", "complex", "reader")
  
# combine center-upper and center-lower
  vis$reader[which(vis$reader == "CU")] <- "C"
  vis$reader[which(vis$reader == "CL")] <- "C"
  vis$reader[which(vis$reader == "CW")] <- "C"
  vis$reader[which(vis$reader == "CE")] <- "C" 
  
# relabel reader
  vis$reader <- paste(vis$complex, vis$reader, sep = "_")
  
# remove split for the sake of memory
  rm(split2)
```

![\~](https://latex.codecogs.com/png.latex?~ "~")

### Some summary statistics

``` r
# How many birds were returning birds over multiple days?
  returning_birds_d <- group_by(rds_sort, bird, date) %>%
    summarise(n=n()) %>%
    group_by(., bird) %>%
    summarise(days=n())

# merge with demographic data
  returning_birds_d$rfid <- 
    str_sub(
      returning_birds_d$bird, 1, 5
    )
  returning_birds_d <- merge(
    returning_birds_d,
    caps[,c("age", "sex", "rfid")],
    sort=F
  )

 sum((returning_birds_d$days > 1) & 
       (returning_birds_d$sex == "M"))
```

    ## [1] 12

``` r
 sum((returning_birds_d$days > 1) & 
       (returning_birds_d$sex == "F"))
```

    ## [1] 33

``` r
 sum((returning_birds_d$days > 1) & 
       (returning_birds_d$sex == "U"))
```

    ## [1] 6

``` r
# returning in multiple years
 returning_birds_y <- group_by(rds_sort, bird, year) %>%
    summarise(n=n()) %>%
    group_by(., bird) %>%
    summarise(years=n())
 
# merge with demographic data
  returning_birds_y$rfid <- 
    str_sub(
      returning_birds_y$bird, 1, 5
    )
  returning_birds_y <- merge(
    returning_birds_y,
    caps[,c("age", "sex", "rfid")],
    sort=F
  )

 sum((returning_birds_y$years > 1) & 
       (returning_birds_y$sex == "M"))
```

    ## [1] 1

``` r
 sum((returning_birds_y$years > 1) & 
       (returning_birds_y$sex == "F"))
```

    ## [1] 7

#### Age and sex ratios in total data

``` r
# get first five letters of RFID to match
  rds_sort$rfid <- str_sub(
    rds_sort$bird, 1, 5
  )
  
# merge with records
  rds_dem <- merge(
    rds_sort,
    caps,
    by="rfid",
    sort = F
  )
  
# male to female ratio in all visits
  birds_in_data <- unique(
    rds_dem[,c("bird", "sex", "age")]
  )
  
  sex_table_all <- data.frame(
    source=c("captures", "rfid_data", "weighted_rfid_data"),
    male= c(
      sum(caps$sex == "M"),
      sum(birds_in_data$sex == "M"),
      sum(rds_dem$sex == "M")
    ),
    female= c(
      sum(caps$sex == "F"),
      sum(birds_in_data$sex == "F"),
      sum(rds_dem$sex == "F")
    ),
    unknown= c(
      sum(caps$sex == "U"),
      sum(birds_in_data$sex == "U"),
      sum(rds_dem$sex == "U")
    )
  )
  
# Age distribution in the data
  age_table_all <- data.frame(
    source = c("captures", "rfid_data", "weighted_rfid_data"),
    HY = c(
      sum(caps$age == "HY"),
      sum(birds_in_data$age == "HY"),
      sum(rds_dem$age == "HY")
    ),
    AHY = c(
      sum(caps$age == "AHY"),
      sum(birds_in_data$age == "AHY"),
      sum(rds_dem$age == "AHY")
    )
  )
```

<table>

<caption>

Summary of sex ratios in each dataset

</caption>

<thead>

<tr>

<th style="text-align:left;">

source

</th>

<th style="text-align:right;">

male

</th>

<th style="text-align:right;">

female

</th>

<th style="text-align:right;">

unknown

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

captures

</td>

<td style="text-align:right;">

24

</td>

<td style="text-align:right;">

95

</td>

<td style="text-align:right;">

40

</td>

</tr>

<tr>

<td style="text-align:left;">

rfid\_data

</td>

<td style="text-align:right;">

13

</td>

<td style="text-align:right;">

40

</td>

<td style="text-align:right;">

10

</td>

</tr>

<tr>

<td style="text-align:left;">

weighted\_rfid\_data

</td>

<td style="text-align:right;">

47658

</td>

<td style="text-align:right;">

653306

</td>

<td style="text-align:right;">

3220

</td>

</tr>

</tbody>

</table>

<table>

<caption>

Summary of age ratios in each dataset

</caption>

<thead>

<tr>

<th style="text-align:left;">

source

</th>

<th style="text-align:right;">

HY

</th>

<th style="text-align:right;">

AHY

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

captures

</td>

<td style="text-align:right;">

43

</td>

<td style="text-align:right;">

114

</td>

</tr>

<tr>

<td style="text-align:left;">

rfid\_data

</td>

<td style="text-align:right;">

9

</td>

<td style="text-align:right;">

53

</td>

</tr>

<tr>

<td style="text-align:left;">

weighted\_rfid\_data

</td>

<td style="text-align:right;">

3123

</td>

<td style="text-align:right;">

701055

</td>

</tr>

</tbody>

</table>

#### How many birds made inter-complex movements?

``` r
# make sure data is in right order
  vis_sort <- vis[order(vis$bird, vis$date, vis$seconds), ]

# create empty dataframe
 mvs_by_bird <- as.data.frame(
   matrix(nrow = nrow(vis_sort), ncol = 3)
 )
 names(mvs_by_bird) <- c(
   "bird",
   "complex1",
   "complex2"
 )
 
 for(i in 2:nrow(vis_sort)){
   if(
     (vis_sort$bird[i] == vis_sort$bird[(i-1)]) &
     (vis_sort$date[i] == vis_sort$date[(i-1)]) &
     (vis_sort$reader[i] != vis_sort$reader[(i-1)])
   ){
     mvs_by_bird[i,] <- c(
       vis_sort$bird[i],
       vis_sort$complex[(i-1)],
       vis_sort$complex[i]
     )
   }
 }
 
 mvs_by_bird <- na.omit(mvs_by_bird)
 
# add demographic data
 mvs_by_bird$rfid <- 
   str_sub(
     mvs_by_bird$bird, 1, 5
   )
 mvs_by_bird <- merge(
   mvs_by_bird,
   caps[,c("age", "sex", "rfid")],
   sort=F
 )
 
 sum(mvs_by_bird$sex == "M")
```

    ## [1] 5

``` r
 sum(mvs_by_bird$sex == "F")
```

    ## [1] 2216

``` r
 sum(mvs_by_bird$age == "AHY")
```

    ## [1] 2221

``` r
 inter_comps <- 
   mvs_by_bird[which(mvs_by_bird$complex1 != mvs_by_bird$complex2), ]
 
 inter_mvs_by_bird <- group_by(inter_comps, bird) %>%
   summarise(
     n=n(),
     age=unique(age),
     sex=unique(sex)
   )
 
  kable(inter_mvs_by_bird,
        caption = "Summary of inter-complex (i.e., long distance) movements by bird.")
```

<table>

<caption>

Summary of inter-complex (i.e., long distance) movements by bird.

</caption>

<thead>

<tr>

<th style="text-align:left;">

bird

</th>

<th style="text-align:right;">

n

</th>

<th style="text-align:left;">

age

</th>

<th style="text-align:left;">

sex

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

00F27598596F0001

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

AHY

</td>

<td style="text-align:left;">

F

</td>

</tr>

<tr>

<td style="text-align:left;">

0CF27598596F0001

</td>

<td style="text-align:right;">

11

</td>

<td style="text-align:left;">

AHY

</td>

<td style="text-align:left;">

F

</td>

</tr>

<tr>

<td style="text-align:left;">

11F27598596F0001

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

AHY

</td>

<td style="text-align:left;">

M

</td>

</tr>

<tr>

<td style="text-align:left;">

24B27598596F0001

</td>

<td style="text-align:right;">

7

</td>

<td style="text-align:left;">

AHY

</td>

<td style="text-align:left;">

F

</td>

</tr>

<tr>

<td style="text-align:left;">

39EBC199F0870001

</td>

<td style="text-align:right;">

9

</td>

<td style="text-align:left;">

AHY

</td>

<td style="text-align:left;">

F

</td>

</tr>

<tr>

<td style="text-align:left;">

3BB27598596F0001

</td>

<td style="text-align:right;">

35

</td>

<td style="text-align:left;">

AHY

</td>

<td style="text-align:left;">

F

</td>

</tr>

<tr>

<td style="text-align:left;">

57327598596F0001

</td>

<td style="text-align:right;">

21

</td>

<td style="text-align:left;">

AHY

</td>

<td style="text-align:left;">

F

</td>

</tr>

<tr>

<td style="text-align:left;">

D9EBC199F0870001

</td>

<td style="text-align:right;">

275

</td>

<td style="text-align:left;">

AHY

</td>

<td style="text-align:left;">

F

</td>

</tr>

</tbody>

</table>
