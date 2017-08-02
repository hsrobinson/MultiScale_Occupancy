#Simulation of multi scale occupancy

library(dplyr) #used to subsample and manipulate data
library(R2jags) #to communicate with jags
library(mcmcplots) 
library(tibble) #simple data frames
load.module("glm")

setwd("D:/Jaguar/Costa_Rica/Combined_Cam_Occ/r-code")

study_areas <- 5 #the number of study areas
cells <- 500 #the number of broad scale cells in the study area
subunits_per <- 5 #subunits per cell
subcells <- cells*subunits_per #the number of small scale cells in the study area
#set.seed(24) #if you want consistancy

##################################################################################################
#Occupancy is a hierarchical "top down" process and it is impossible for a subcell to be occupied#
#if a cell is not.  However to test the model's ablity to estimate beta coefs this truth is ignored
##################################################################################################

#create  spatial variables - normal distribution for continous variables with mean of 0 and standard deviation of 1

#broad scale first
cellndvi<-rnorm(n=cells, 0, 1)
cellprotect<-rnorm(n=cells, 0, 1)
cellroad<-rnorm(n=cells, 0,1) #format is (number of obs, mean, SD)

#of uniform if you prefer
#cellroad<-runif(n=cells, -1,1) #format is (number of obs, min, max)
#cellndvi<-runif(n=cells, -1, 1)
#cellprotect<-runif(n=cells, -1, 1)

#create beta values for psi at cell scale
alpha_cell.occ <- 0
beta_cellndvi.occ <- 0.2
beta_cellprotect.occ <- 0.4
beta_cellroad.occ <- 0.6

#create continous probability of occupancy based on above betas and random covariate values
cellocc.prob <- plogis(alpha_cell.occ + beta_cellndvi.occ*cellndvi + beta_cellprotect.occ*cellprotect + beta_cellroad.occ*cellroad) #plogis is slick single step probability

#occupancy is binomial so convert continous prob into 1s and 0s
true_cell.occ <- rbinom(n=cells, size=1, prob=cellocc.prob)

#to test how well a logistic regression can return created betas
#glm(true_cell.occ ~ cellndvi + cellprotect + cellroad, family = binomial)

#combine needed data into a single dataframe with cellID and studyID columns 
studyID <-rep(c(1:study_areas), each = cells/study_areas)
celldata <- as.tibble(cbind (studyID, true_cell.occ, cellndvi, cellprotect, cellroad))
celldata <- mutate(celldata, ID = row_number()) %>%
  arrange(ID)

#small scale or subunit next
#aditional covariate of cameras per small cell as detection covariate
cameras<- sample(1:3, subcells, replace=T) #random number between 1 and 3
subroad<-rnorm(n=subcells, 0,1) 
subndvi<-rnorm(n=subcells, 0, 1)
subprotect<-rnorm(n=subcells, 0, 1)
alpha_sub.occ <- 0
beta_subndvi.occ <- 0.2
beta_subprotect.occ <- 0.4
beta_subroad.occ <- 0.6
alpha_sub.p <- 0
betacams.p <- .25


#create continous probability of occupancy based on above betas and random covariate values
subocc.prob <-plogis(alpha_sub.occ + beta_subroad.occ*subroad + beta_subndvi.occ*subndvi + beta_subprotect.occ*subprotect)

#convert continous prob into 1s and 0s
true_sub.occ <- rbinom(n=subcells, size=1, prob=subocc.prob)

#to test how well a logistic regression can return created subbetas
#glm(true_sub.occ ~ subndvi + subprotect + subroad, family = binomial)

#create continous probability of detection in each subcell
det.prob <- plogis(alpha_sub.p + betacams.p*cameras)


#create an ID column; repeating series based on number of broad cells and subunits per 
#combine two dataframes and reference cells properly
ids <- expand.grid(ID = 1:cells, subID = 1:5) 

#  cbind is not the greatest here b/c it returns a matrix and matrices can only
#   contain one class of data
simdata <- cbind(ids, true_sub.occ, det.prob, cameras, subndvi, subprotect, subroad) %>%
  tibble::as_tibble() %>%
  full_join(., celldata, by = "ID") %>%
  rename(cellID = ID) %>%
  arrange(cellID, subID)

#correct true sub unit occupancy and detection probability based on true cell and subunit occupancy
#i.e. force 0s for unoccupied cells  
# turn on or off depending on desire to estimate subcell betas 
simdata$true_sub.occ<- simdata$true_cell.occ*simdata$true_sub.occ
simdata$det.prob <- simdata$true_sub.occ * simdata$det.prob


#SIMULATED SURVEY
prop <- 0.5 #percent of total cells surveyed
trials <- 7 # number of survey periods or trials

#create a blank array (y) with same structure as survey data
y <- array(NA, dim = c(nrow(simdata), trials))

#for loop filling blank array with detections and non-detections based on det.prob
for(i in 1:trials){
  y[,i] <- rbinom(n = subcells, size = 1, prob = simdata$det.prob)
}

#combine trials/detection history with simulated covariates and arrange columns for easier reading 
fulldata <- bind_cols(simdata, as.tibble(y)) %>%#using as.tibble to maintain column names of y
  select("studyID", "cellID", "subID", "true_cell.occ", "true_sub.occ", everything())

#select random sample of full data based on proportion of area or cells we sampled and sort on subID
#sampdata <- fulldata %>% group_by(cellID, subID) %>%
#  sample_frac( size= prop)


sampdata <- fulldata %>%
  tidyr::gather(trial, obs, V1:V7) %>%
  mutate(
    trial = as.integer(sub("[A-z]", "", trial)),
    obs = replace(obs, sample(1:n(), 0.5*n(), replace = F), NA)
        ) %>%
  arrange(cellID, subID, trial)



#clean-up
rm(list= ls()[!(ls() %in% c("fulldata", "sampdata" ))]); gc()

#################################################
#  Covariate creation
#  Cell level
cellcov <- fulldata %>% 
  select(cellID, cellndvi, cellprotect, cellroad) %>% 
  distinct()%>%
  select(-cellID)

#  Sub level
subndvi <- fulldata %>% 
  select(cellID, subID, subndvi) %>% 
  tidyr::spread(subID, subndvi) %>%
  select(-cellID)

subprotect <- fulldata %>% 
  select(cellID, subID, subprotect) %>% 
  tidyr::spread(subID, subprotect) %>%
  select(-cellID)

subroad <- fulldata %>% 
  select(cellID, subID, subroad) %>% 
  tidyr::spread(subID, subroad) %>%
  select(-cellID)

#  Trial level covariates
cameras <- array(NA, dim = c(cells, subunits_per, trials))  

for(i in 1:nrow(sampdata)){
  cameras[sampdata$cellID[i], sampdata$subID[i], sampdata$trial[i]] <- sampdata$cameras[i] 
}

#################################################################################################################
#NOW PASS SIMULATED/SAMPLED DATA ON TO JAGS FOR ANALYSIS
# Converting everything to Kery's language
#################################################################################################################

# Bundle and summarize data set
y <- data$y
str( win.data <- list(y = y, n.pond = dim(y)[1], n.samples = dim(y)[2], n.pcr = dim(y)[3],
                      covA = data$covA, covB = data$covB, covC = data$covC))
