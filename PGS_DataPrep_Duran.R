#######################################################################
#######################################################################
##     This script was created by Dr. Jen Cruz and modified by 
##    Zoe Duran as part of the Applied Population Ecology Class     ###
##                                                                  ###
##  Our study species is the Piute ground squirrel, Urocitellus      ## 
## mollis. Ground squirrels are widely distributed in sagebrush      ##
## habitats of the Great Basin and Columbia Plateau. We focus        ##
##  on their populations inside the birds of prey NCA.               ##
## In this script we prepare the data that we will use              ###
## to estimate site abundance of Piute ground squirrels.             ## 
##  Site locations were randomly selected within the study area     ###
## ensuring that they were at least 1.5km m apart for independence.  ##
## Sampling occurred during 2013-2019, 2021 and 2022. Mark-recapture 
## trapping       ##
## occurred at eight sites representing four general habitat types.  ##
## A robust study design was used with primary trapping periods      ##
## representing either year or season within year (depending on the  ##
## question) with at least 3 repeat secondary trapping periods within## 
## each primary. Closed population is assumed during secondary       ##
## periods.                                                          ##
##                                                                   ##
##   For our abundance data we create two data summaries:            ##
##   1) for year 1 only so that we can use closed population models.  #
##   2) for all years as a classic open ROBUST DESIGN.               ##
##                                                                   ## 
## Predictors:                                                       ##
# We use landcover data from the National Geospatial Data Asset (NGDA).#
# https://www.mrlc.gov. We downloaded sagebrush and annual           ##
# grasses (which we assumed are mainly cheatgrass) from 2007-2018    ##
# Data are available as a raster of 30 x 30 m cell resolution.       ##
# We summarized annual % cover surrounding our sites with a 200 m buffer #
# Habitat data are at the site X year level resolution.               #
# Climate data were downloaded from the Twins Falls Station           #
# including daily measures of min, mean and max Temperatures          # 
# and total precipitation. We extracted minimum Temperatures for      #
# Feb, when squirrels come out of hibernation and max temperatures    #
# for Apr-May, when consistently hot days may prevent them from       #
# foraging.                                                           #
##ZD- fill out with actual predictors                                 #
##Detection predictors: weather (wind, temp min, temp max, precip);   #
##secondary day (e.g., first, second or last day of trapping);        #
##individual variation (e.g., boldness)?, julian date, effort?        #
##Abundance predictor: site (e.g., habitat); year; season             #
##Survivorship predictor: boldness BLUP (or random effect?); site,    #
##year, season, age (juv vs ad), sex, weight?                         #
#######################################################################
##### Set up your workspace and load relevant packages -----------
# Clean your workspace to reset your R environment. #
rm( list = ls() )
# Check that you are in the right project folder
getwd()

# Install new packages from "CRAN" repository. # 
install.packages( "tidyverse" ) #actually a collection of packages 

# load packages:
library( tidyverse ) 
## end of package load ###############
###################################################################
#### Load or create data -----------------------------------------
# set directory where your data are. 
# I have it in a Data folder in my Rstudio project:
datadir <- paste( getwd(), "/GroundSquirrelCapt-Recapt/", sep = "" )
datadir

# load observed occurrences:
obs_df <- read.csv( file = paste( datadir, "CleanTrapData.csv", sep = ""),
                    header = TRUE )
#view
head( obs_df ); dim(obs_df)
str(obs_df)


#load predictor data; using trap_preds_ijk.csv for now
preddf <- read.csv( file = paste( datadir, "trap_preds_ijk.csv", sep = ""),
                    header = TRUE )
#view
head( preddf ); dim( preddf )
str(preddf)

########## end of data load ################################
#######################################################################

######## explore data #############

# Start viewing our observations dataframe:
head( obs_df ); dim( obs_df )
# Rename column (name to right of '=' is the old name)
obs_df <- obs_df %>%
  rename(Recap = Recapture..Y.N.)
obs_df <- obs_df %>%
  rename(year = Year)
# Our observations dataframe metadata:
# WebName = site id (I)
# OpenTime = time traps were opened
# CloseTime = time traps were shut
# ProcessTime = time the individual trapped was processed
# Trap ID = trap number; also associated with trap location which is static
# Sex = M/F
# BagWeight.g. = handling bag weight in g
# GrossWeight.g. = weight of handling bag plus squirrel in g
# NetWeight.g. = gross minus bag (i.e., squirrel weight)
# AgeClass = adult or juvenile (i.e., young of year)
# ReproductiveStatus = females- presence of nipples; males- testes distended
#    note: phenology of trapping later than breeding (no males with testes),
#    juv are not reproductive in their first year so juv repro status = NA
# Parasites = presence of parasites (Y/N/NR) --ZD check that FLEAS = Y?
# Recap = recapture where N = new animal, Y = previously tagged
# PitTagTD. = unique identifier for each animal
# BLUP = best linear unbiased predictor for boldness (>0 = more bold, 
#    <0 = less bold)
# Year = year of trapping efforts
# Poop = whether a fecal sample was collected
# Wiggle.time..s. = time spent in motion during handling bag test (use BLUP for
#    personality, not wiggle)
# date = date of trapping day
# ProcessHr = hour of day animal was processed
# jday = julian date
# month = month of year

# NOTE- throughout- NR = no record

# Let's define some parameters
# How many sites were sampled?
I <- length( unique( obs_df$WebName ) )
I

# What years were sampled (primary seasons)?
yrrange <- sort( unique( obs_df$year) )
yrrange

#How many years (i.e., primary seasons )?
T <- length( yrrange )
T


# How many replicate surveys
J <- 3
#view
I; yrrange; T; J

# For our 1st set of data we create a closed population dataframe #
# that only includes the last year regular data in 2019, and columns that are #
# relevant.

# Select obs_df:
closeddf <- obs_df %>% 
  #filter only rows for 2019:
  dplyr::filter( year == '2019') %>%
  #select desired columns to keep:
  dplyr::select( WebName, TrapID, Sex, NetWeight.g., AgeClass, Recap, 
                 PitTagID., BLUP, jday, month, year, Survey) 
#view resulting dataframe
head( closeddf ); dim( closeddf )

# How many recaptures across our 903 captures (903 determined by dim output)?
table(closeddf$Recap)
#result = 564 recap, 339 new

# We also check for missing values in recap and PIT tag #
colSums( is.na( closeddf[, c("Recap", "PitTagID.")]) )
#none are present

#### Checking our predictors ------------
# What about our predictors?
# view
head( preddf ); dim( preddf )

###########ZD ENDED HERE######################################################
# Before thinking of whether we can include predictors in our#
# models we should check their distribution and correlation #
# Why?
# We start by checking for outliers, skewed distribution etc #

# turn trapping effort into an integar rather than character
preddf$effort_mins = lubridate::time_length( preddf$effort, unit = "minutes" )

# create a vector with predictor names
prednames <- c("tempC_st", "wind_kmph_st", "effort_mins")
# loop over each to create histograms for each predictor:
for( p in 1:length(prednames) ){
  # create an object with the ggplot so that you can display it 
  # in a loop 
  a <- ggplot( preddf ) + #choose your data
    theme_bw( base_size = 15 ) + #choose a preset theme
    labs( x = prednames[p] ) + #label x axis using our predictor names
    geom_histogram( aes( get(prednames[p]) ), bins= 10 ) #plot histogram
  # display your plot object
  print( a )
}
# What do you note? Any apparent issues with these data?
# Answer:
#
# Let's plot how predictors vary annually:
for( p in 1:length(prednames) ){
  # We can also incorporate site variability for habitat:
  bp<- ggplot( preddf ) +
    theme_bw( base_size = 15 ) + #choose a preset theme
    theme( legend.position = "none" ) + #remove legend display
    labs( y = prednames[p] ) + #label x axis using our predictor names
    # plot each site individually
    geom_line( aes( x = year, y = get(prednames[p]),
                    color = as.factor(WebName) ), size = 1.5 )
  
  # Here we rely on a smoothing spline to get mean annual trends #
  # across all sites:
  cp<- ggplot( preddf, aes( x = year, y = get(prednames[p]) ) ) +
    theme_bw( base_size = 15 ) + #choose a preset theme
    labs( y = prednames[p] ) + #label x axis using our predictor names
    geom_smooth( size = 2 ) #smooth mean across all sites
  #display plots
  print( bp )
  print( cp )
}
# Any notable sites? 
# Answer:
#
# Now that we are happy with no outliers, we can check for #
# correlation among predictors. Why is this important?
cor( preddf[ , prednames] )
# Are there any predictors we need to worry about?
# What correlations would be worrisome?
### end predictor check ----------------

# Now that we are satisfied with our predictor data we can #
# append it to our new closeddf:
closeddf <- #select the columns we want to keep in preddf 
  preddf %>% dplyr::select( jday, WebName, year, Survey, all_of(prednames) ) %>%
  right_join( closeddf, by = c("WebName", "year","jday") )
# why did we use right_join instead of left_join?
# check output #note that I always check dimensions when joining
# dataframes. Sometimes we add or substract rows unintentionally 
# so it is good to check that your code had the desired result
head( closeddf); dim( closeddf )

##############CLOSED ENDS HERE#########################################
########################################################################
# We repeat the process for our robust design. we want to join #
# the two dataframes keeping relevant columns 
head( obs_df )
opendf <-  obs_df %>%
  #select desired columns to keep:
  dplyr::select( WebName, TrapID, Sex, NetWeight.g., AgeClass, Recap, 
                 PitTagID., BLUP, jday, month, year ) 
#check
head( opendf ); dim( opendf )
# We append predictors:
opendf <- preddf %>% 
  dplyr::select( WebName, year, all_of(prednames) ) %>%
  left_join( opendf, by = c("WebName", "year") )
#check
head( opendf ); dim( opendf )
# We also check for missing values in the response #
# as those are often not allowed in frequentist analyses
colSums( is.na( opendf[, c("Recap", "PitTagID.")]) )
###QUESTION- what is response in mark-recap dataset? "Recap"? "PITTagID"?
#What happens when PIT or recap is NA? wHow can I quick view these records
#and how do I remove them?

################################################################
##########    save relevant data and workspaces     ###########
# Check that you are in the right project folder
getwd()
#save closed dataframe in our data folder:
write.csv( closeddf, paste( getwd(),
                            "/GroundSquirrelCapt-Recapt/PGS_2019_closedf.csv", 
                            sep = "" ),  
           row.names = FALSE )

#save open dataframe in our data folder:
write.csv( opendf, paste( getwd(),"/GroundSquirrelCapt-Recapt/PGS_opendf.csv", 
                          sep = "" ),  
           row.names = FALSE )

# if you want to save any of the plots we produced #
# for your presentation or ms, you do it here too:
# Save the most recently viewed plot with ggsave() to define file type, resolution, 
# and plot dimensions:
ggsave("GroundSquirrelCapt-Recapt/AprMayTXYear.png", dpi=500, 
       height = 10, width = 15, units= "cm" )
# or if you saved it as an object:
#start by calling the file where you will save it
tiff( 'GroundSquirrelCapt-Recapt/FebTXYear.tiff',
      height = 10, width = 15, units = 'cm', 
      compression = "lzw", res = 400 )
#call the plot
bp
#turn off
dev.off()

# if you want to save your workspace, because you are still #
# working through it use the following command:
save.image( "DataPrepWorkspace.RData" )
########## End of saving section ##################################
################## Save your data and workspace ###################

############# END OF SCRIPT ########################################