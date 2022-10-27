#######################################################################
#######################################################################
####   This script was created by Dr. Jen Cruz as part of      ########
####          the Applied Population Ecology Class              #######
####                                                           ########
# Our study species is the Piute ground squirrel, Spermophilus        # 
# mollis at the Morley Nelson Birds of Prey Nature Conservation Area. #
#                                                                     #            
# 20 sites were randomly selected for trapping over three days, after #
# three days of pre-baiting. This approach is meant to increase       #
# trappability, but may not avoid trap-happiness.                     #
# Trapping occurred over multiple years but tags used for marking     #
# individuals were temporary, so they lasted during the primary season#
# but not between seasons.                                            # 
#######################################################################
##### Set up your workspace and load relevant packages -----------
# Clean your workspace to reset your R environment. #
rm( list = ls() )
# Check that you are in the right project folder
getwd()

# load packages:
library( tidyverse ) 
## end of package load ###############
###################################################################
#### Load or create data -----------------------------------------
# set directory where your data are. 
# I have it in a Data folder in my Rstudio project:
datadir <- paste( getwd(), "/simdata/", sep = "" )

# load observed occurrences:
ind_df <- read.csv( file = paste( datadir, "Indf.csv", sep = ""),
                    header = TRUE )
#view
head( ind_df ); dim(ind_df)

# This time we don't join our siteXyear predictors with the individual-
#level dataframe but we keep them separate so we don't load the 
# predictor dataframe

########## end of data load ################################
#######################################################################
######## explore individual capture data #############
# Start viewing our observations dataframe:
head( ind_df ); dim( ind_df )

#how many sites?
(M <- length( unique(ind_df$o.sites) ) )
# N=20

#which years (primary seasons) were sampled:
(yrrange <- sort( unique( ind_df$year ) ) )
# 2007-2018

#how many years were sampled:
(T <- length( yrrange ) )
# 12 years

#define number of repeat surveys
J <- 3 

#calculate total number of captures for each sex:
N_totcap <- ind_df %>% group_by( o.sites, year, sex ) %>% 
  summarise( captures = n()  ) %>% ungroup()
#check
head( N_totcap ); dim(N_totcap )
#plot totcap estimates:
ggplot( N_totcap, aes( x = year, y = captures, color = as.factor(o.sites) ) ) +
  geom_line( size = 2 ) +
  theme_bw( base_size = 15) +
  facet_wrap( ~sex, nrow = 2 )
#plot totcap estimates:
ggplot( N_totcap, aes( x = year, y = captures, color = as.factor(sex) ) ) +
  geom_line( size = 2 ) +
  theme_bw( base_size = 15) +
  facet_wrap( ~o.sites, nrow = 2 )

# Can we see any sex differences?
# Answer:
# Note- this is not naive captures, it is total captures (includes recaps), 
# but appears that sex caps are similar across site a year

# Plot totcap annual estimates of abundance
N_totcap %>% group_by( year ) %>% 
  summarise( total = sum(captures) ) %>% 
  ggplot(., aes( x = year, y = total )) +
  geom_line( size = 2 ) +
  theme_bw( base_size = 15 )

# What does this tell us?
# Answer:
# PGS total captures have decreased over time

# we can also summarize capture histories
table( ind_df$ch )
# what do you notice? 
# Answer:
# Most common capture history is "111" which indicates high trappability and 
# possibly 'trap-happy' effect (p and c may need behavioral model)

# Let's fix our capture histories
ind_df$ch <- apply( ind_df[,c('trap.j1', 'trap.j2', 'trap.j3') ], 
                    1, paste, collapse="" )

#check again
table( ind_df$ch )

# We still do have to select a year for the single season analysis.
# I randomly choose 2013 this time
ind_2013 <- ind_df %>% 
  filter( year == 2013 )
#check
head( ind_2013 ); dim( ind_2013 )

##########    save relevant data and workspaces     ###########
#save multi-season dataframe:
write.csv( ind_df, paste( getwd(),"/simdata/ind_multi.csv", sep = "" ),  
           row.names = FALSE )

#save single season dataframe:
write.csv( ind_2013, paste( getwd(),"/simdata/ind_2013.csv", sep = "" ),  
           row.names = FALSE )

# if you want to save your workspace, because you are still #
# working through it use the following command:
save.image( "MarkPrepWorkspace.RData" )
########## End of saving section ##################################
################## Save your data and workspace ###################

############# END OF SCRIPT ########################################