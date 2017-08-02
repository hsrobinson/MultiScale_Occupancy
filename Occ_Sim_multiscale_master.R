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
cellroad<-rnorm(n=cells, 0,1) #format is (number of obs, mean, SD)
cellndvi<-rnorm(n=cells, 0, 1)
cellprotect<-rnorm(n=cells, 0, 1)

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
cellocc.prob <- plogis(alpha_cell.occ + beta_cellroad.occ*cellroad + beta_cellndvi.occ*cellndvi + beta_cellprotect.occ*cellprotect) #plogis is slick single step probability

#occupancy is binomial so convert continous prob into 1s and 0s
true_cell.occ <- rbinom(n=cells, size=1, prob=cellocc.prob)

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

#select random sample of y based on proportion of area or cells we sampled and sort on subID
sampdata <- fulldata %>%
  tidyr::gather(trial, obs, V1:V7) %>%
  mutate(
    trial = as.integer(sub("[A-z]", "", trial)),
    obs = replace(obs, sample(1:n(), 0.5*n(), replace = F), NA)
  ) %>%
  arrange(cellID, subID, trial)

#clean-up
#rm(list= ls()[!(ls() %in% c("fulldata", "sampdata" ))]); gc()

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
#################################################################################################################

sp.params = c("alpha_cell.occ","alpha_sub.occ", "alpha_sub.p", "beta_cellndvi.occ","beta_cellprotect.occ",
              "beta_cellroad.occ", "beta_subndvi.occ", "beta_subprotect.occ", "beta_subroad.occ",
              "beta_cams.p")

#  Create intial values for det_hist
det_init <- sampdata %>% 
    group_by(cellID, subID) %>% 
    mutate(
      hold = is.na(obs),
      init = as.numeric(any(obs == 1, na.rm = T)),
      init = replace(init, !hold, NA)
    ) %>% 
    .$init

cell_init <- sampdata %>% 
  group_by(cellID) %>% 
  summarise(
    init = max(obs, na.rm = T)
  ) %>% 
  .$init

sub_init <- sampdata %>% 
  group_by(cellID, subID) %>% 
  summarise(
    init = max(obs, na.rm = T)
  ) %>% 
  mutate(
    init = replace(init, init == -Inf, 0)
  ) %>%
  ungroup %>%
  tidyr::spread(subID, init) %>%
  select(-cellID)


#Specify the initial values
sp.inits = function() {
  list(
        alpha_cell.occ = runif(1, -1, 1),#n, min, max
        alpha_sub.occ = runif(1, -1, 1),
        alpha_sub.p = runif(1, -1, 1),
        beta_cellndvi.occ = runif(1, -1, 1),
        beta_cellprotect.occ = runif(1, -1, 1),
        beta_cellroad.occ = runif(1, -1, 1),
        beta_subndvi.occ = runif(1, -1, 1),
        beta_subprotect.occ = runif(1, -1, 1),
        beta_subroad.occ = runif(1, -1, 1),
        beta_cams.p = runif(1, -1, 1),
        det_hist = det_init,
        cellocc = cell_init,
        subocc = as.matrix(sub_init)
  )
}

#  Organize data for jags
jdat <- list(
  det_hist = sampdata$obs,
  cellID = sampdata$cellID,
  subID = sampdata$subID,
  trial = sampdata$trial,
  ncell = length(unique(sampdata$cellID)),
  nsub = length(unique(sampdata$subID)),
  ntrial = length(unique(sampdata$trial)),
  cellcov = cellcov,
  subndvi = subndvi,
  subprotect = subprotect,
  subroad = subroad,
  cameras = cameras,
  nobs = nrow(sampdata),
  det_init = det_init,
  cell_init = cell_init,
  sub_init = sub_init
  
)

jag.res <- jags.parallel(
  jdat, 
  sp.inits, 
  sp.params, 
  "modelworkinglong.txt", 
  n.chains=3, 
  n.iter=1000, 
  n.burnin=100, 
  n.thin=10
)

jag.res$BUGSoutput$summary
mcmcplot(jag.res)

#################################################################################################################
#NOW PASS SIMULATED/SAMPLED DATA ON TO JAGS FOR ANALYSIS
#################################################################################################################

#seperate detection history from covariate data and reshape for jags
det_hist = select(sampdata, studyID, cellID, subID, starts_with("V")) #takes advantage of the original naming in the for loop
det_hist = simplify2array(by(det_hist, det_hist$cellID))

#jags wants a list of variables that will be used in the model

sp.data <- list(det_hist=det_hist, cellndvi= sampdata$cellndvi, cellprotect=sampdata$cellprotect,
                cellroad=sampdata$cellroad, cameras=sampdata$cameras, subndvi=sampdata$subndvi, subprotect=sampdata$subprotect, 
                subroad=sampdata$subroad, cells= length(unique(sampdata$cellID)), subcells= length(unique(sampdata$subID)), trials=ncol(det_hist))

#Specify the parameters to be monitored
sp.params = c("alpha_cell.occ","alpha_sub.occ", "alpha_sub.p", "beta_cellndvi.occ","beta_cellprotect.occ",
              "beta_cellroad.occ", "beta_subndvi.occ", "beta_subprotect.occ", "beta_subroad.occ",
              "beta_cams.p")

#creat initial value for occupancy
occ_init <- apply(det_hist, 1, function(x) as.numeric(any(x > 0)) )

#Specify the initial values
sp.inits = function() {
  list(alpha_cell.occ = runif(1, -1, 1),#n, min, max
       alpha_sub.occ = runif(1, -1, 1),
       alpha_sub.p = runif(1, -1, 1),
       beta_cellndvi.occ = runif(1, -1, 1),
       beta_cellprotect.occ = runif(1, -1, 1),
       beta_cellroad.occ = runif(1, -1, 1),
       beta_subndvi.occ = runif(1, -1, 1),
       beta_subprotect.occ = runif(1, -1, 1),
       beta_subroad.occ = runif(1, -1, 1),
       beta_cams.p = runif(1, -1, 1),
       occ = occ_init)
      }



#Run the model and call the results
#jags text file must be in working directory
jag.res = jags(sp.data, sp.inits, sp.params, "modelworking.txt", 
               n.chains=3, n.iter=10000, n.burnin=1000, n.thin=10)



