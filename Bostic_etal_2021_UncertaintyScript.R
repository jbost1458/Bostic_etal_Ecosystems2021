#########################
#R script associated with Bostic et al., 2021 Ecosystems paper
#using bootstrapping to estimate mean and uncertainty of:
#D17O, NO3-Atm concentrations, and Processing Efficiency 
#########################


#Data needed for this analysis is at:
#insert link to github page


#For full details please see the Supplementary Information with Bostic et al., 2021 Ecosystems
#Note that values calculated using this script will vary slightly from what is reported in 
#the SI of Bostic et al., 2021 due to the randomness inherent in sampling with replacement (bootstrapping)
#These variations are minor and do not change the significance of regression relationships reported
#in the manuscript

#packages needed:
#tidyverse


library(tidyverse)



#Psuedo-code
#1 Import relevant data
#2 Estimate mean D17O and uncertainty for each watershed
#3 Estimate mean flow-weighted NO3-Atm concentration using D17O and uncertainty
#with flow-weighted total NO3 concentration and uncertainty
#4 Estimate processing efficiency and uncertainty using D17O and uncertainty
#with total No3 yields and wet NO3 deposition


#Uncertainty in the following calculations is derived from:
#1 Analytical uncertainty associated with isotope measurements
#2 Uncertainty in the exact D17O value of terrestrial NO3 (accounted for by varying beta
#in the D17O equation)
#3 Uncertainty in the D17O value of our deposition end-member
#4 Uncertainty in estimates of NO3 Total flow-weighted concentrations and yields 



#1 Import relevant data
#Import .csv with water quality data - named "isotope_NO3_Bostic2021.csv"
dat <- read.csv("path and name here")
dat <- read.csv("/Users/joelbostic 1/Desktop/PhD/Dissertation/ch1/Writing/Ecosystems_Round2/github_dataStuff/isotope_NO3_Bostic2021.csv")




#2 Estimate mean D17O and uncertainty for each watershed
#split data into a list - each site is a list
dat.list <- split(dat, dat$Site)

#set number of iterations for bootstrapping
iterations <- 10000


#Define the normal-uniform mixed probability density function
support <- seq(-10, 10, length = 1000)
pd.f <- function(x, mu, sigma, w){
  h <- 1/(1+w/(sqrt(2*3.14*sigma)))
  y <- vector()
  for(i in 1:length(support))  
    if(x[i] <= (mu-w/2)){
      y[i] <- (h/sqrt(2*pi*sigma)*exp(-(1/(2*sigma^2))*(x[i]-mu+w/2)^2))
    } else if(x[i] >= (mu+w/2)){
      y[i] <- (h/sqrt(2*pi*sigma)*exp(-(1/(2*sigma^2))*(x[i]-mu-w/2)^2))
    } else{
      
      y[i] <- h/sqrt(2*pi*sigma)
    }
  y
}

#sample with replacement function
samp_fun <- function(n) sample(
  support, 
  n, 
  TRUE,  
  probs 
)


#Estimate mean+uncertainty D17O over entire sampling period
#First generate bootstraps of individual D17O measurement uncertainty
#Second - sample from this uncertainty distribution to estimate mean+uncertainty D17O from each watershed

#create list to store D17O bootstrap distributions
D17O.dist <- list()

#set beta = 0.52
beta <- 0.52

#for loop to estimate D17O by site
for(i in 1:length(dat.list)){
  #temporary data frame of d17O and d18O data for each site
  temp.dat <- dat.list[[i]]
  #temporary matrix - each column represents collected streamwater samples
  #and rows represent bootstraps of each collected streamwater sample
  temp.mat <- data.frame(matrix(ncol = length(temp.dat$d17O), nrow = iterations))
  
  #for loop to 
  for(j in 1:length(temp.dat$d17O)){
    
    #calculate D17O for each sample given measured d17O and d18O
    temp.dat$D17O[j] <- (((1+temp.dat$d17O[j]/1000)/(1+temp.dat$d18O[j]/1000)^beta)-1)*1000
    
    #calculate "window" associated with each measurement
    min <- temp.dat$D17O[j]-1
    max <- temp.dat$D17O[j]+1
    
    #define mean value (measured D17O)
    mu <- temp.dat$D17O[j]
    
    #generate normal-uniform mixed distribution for each sample
    probs <- pd.f(x = support, mu = mu, sigma = 0.5, w = (max-min))
    probs <- probs/sum(probs)
    
    
    temp.mat[,j] <- samp_fun(iterations)
    
    
    
  }
  D17O.dist[[i]] <- temp.mat
  names(D17O.dist[[i]]) <- unique(dat.list[[i]]$Site)
}


###Note the next following step takes a few minutes to run
#Now use D17O.dist to estimate mean and standard deviation of D17O for each site
#create list for mean D17O estimates - 1 per watershed
D17Omean.list <- list()
for(i in 1:length(D17O.dist)){
  #data frame for watershed
  #each column = 1 sample, rows are 10,000 bootstrap samples from probability distribution
  tempdf <- D17O.dist[[i]]
  
  #vector to store watershed means
  site.mean <- vector()
  #vector to store 
  temp.site.mean <- vector()
  for(j in 1:iterations){
    for(k in 1:length(1:ncol(tempdf))){
      #for every column (sample), randomly sample 1 row, then store in temp.site.mean vector
      #repeat this 10,000 times (iterations)
      temp.site.mean[k] <- sample(tempdf[,k], 1, replace = T)
      
      
    }
    #calculate mean D17O, store in site.mean vector that houses 10,000 bootstraped means
    site.mean[j] <- mean(temp.site.mean)
  }
  #store the 10,000 estimtes of mean D17O in D17Omean.list
  D17Omean.list[[i]] <- site.mean
  #add in site name
  names(D17Omean.list[[i]]) <- names(D17O.dist[[i]][1])
  
}


#From D17O.mean.list the mean and standard deviation of D17O can be readily calculated for each watershed
#for example
D17O.dat <- data.frame(D17O.mean = sapply(D17Omean.list, mean), 
                                    D17O.sd = sapply(D17Omean.list, sd),
                                    Site = sapply(D17Omean.list, names)[1,1:14])
D17O.dat




#3 Estimate mean flow-weighted NO3-Atm concentration using D17O and uncertainty
#with flow-weighted total NO3 concentration and uncertainty

#first we calculate the mean fraction of atm No3 using the D17O end-members (terr and atm)
#terrestrial end member has mean = 0, sd = 0.3 with normal distribution
mean.terr <- 0
sd.terr <- 0.3
#deposition end member is weibull distribution with following parameters
dep.shape <- 11.59
dep.scale <- 26.76

#create list to store estimates of fraction of atmospheric no3
f.atmlist <- list()

for(i in 1:length(D17Omean.list)){
  #vector to store estimates of mean f.atm
  temp <- vector()
  
  #randomly sample with replacement from watershed mean D17O, use this with D17O end-members
  #to estimate fraction atm No3 distribution
  for(j in 1:iterations){
    temp[j] <- (sample(as.numeric(D17Omean.list[[i]]), 1, replace = T)-rnorm(1, mean = mean.terr, sd = sd.terr))/
      (rweibull(1, shape = dep.shape, scale = dep.scale)-rnorm(1, mean = mean.terr, sd = sd.terr))
    
  }
  
  #store 10,000 estimtes of f.atm in f.atmlist
  f.atmlist[[i]] <- temp
  
  #add site names
  names(f.atmlist[[i]]) <- names(D17O.dist[[i]][1])
}

#Next, use the estimates of the mean annual fraction of atmospheric nitrate
#multipled by estimates of mean annual NO3-total flow-weighted concentrations 
#to calculate NO3-Atm flow-weighted concentrations (mean and standard deviation)

#This step requires estimates of mean water year 2016-2017 flow-weighted No3-total 
#estimates. These were made using WRTDS-K and uncertainty was estimated using a block bootstrapping approach

#import NO3-Total concentration data - separate csv's exist for each watershed
#that include 800 bootstrap replicates of water year, flow-weighted mean nitrate total concentration

#save files from github (all files ending in "DWMC.csv") into single folder
#import files with next 2 lines of code
file_names <- list.files("path to folder here", full.names = T)

########DELETE###############
file_names <- list.files("/Users/joelbostic 1/Desktop/PhD/Dissertation/ch1/R/WRTDS_K_uncert_output/conc/", full.names = T)
########DELETE###############
no3.tot.conc <- do.call(rbind,lapply(file_names,read.csv))

#add site names
site.nm <- data.frame(id = unique(no3.tot.conc$id), site.new = c("ANT2", "ANT", "BIGR", "BLAC", "CAC", "CON", "DPRN", 
                                                                 "GEO", "GUN", "GWN", "MON", "TERR", "TOW", "WIL"))
no3.tot.conc <- left_join(no3.tot.conc, site.nm, by = "id")

#only keep WY2016 and 2017
no3.tot.conc <- no3.tot.conc %>% filter(WY == 2016 | WY == 2017)

#create new matrix for no3.total bootstrap replicates
#rows = average No3-total flow-weighted concentration for WY2016-2017
no3.total <- matrix(ncol = 801, nrow = 14)
for(j in 1:14){
  for(i in 2:801){
    no3.total[j, i] <- mean(no3.tot.conc[which(no3.tot.conc$site.new == unique(no3.tot.conc$site.new)[j]), i])
  }
}
no3.total <- as.data.frame(no3.total)
no3.total$site.new <- unique(no3.tot.conc$site.new)

#no3.total contains 800 bootstrap replicates of flow-weighted No3-total averages over water years 2016-2017

#now multiply f.atm by no3.total to estimate average flow-weighted concentrations of No3-atm during water
#years 2016-2017

#note that this will take a 5-10 minutes to run
no3.atm.conc <- list()
for(i in 1:length(f.atmlist)){
  temp <- vector()
  
  for(j in 1:iterations){
    temp[j] <- sample(as.numeric(no3.total[no3.total$site.new == names(f.atmlist[[i]][1]), 2:801]), 1, replace = T)*
      sample(f.atmlist[[i]], 1, replace = T)
    
  }
  no3.atm.conc[[i]] <- temp
  names(no3.atm.conc[[i]]) <- names(f.atmlist[[i]][1])
}

#mean and sd of NO3-Atm flow-weighted concentrations can be displayed by
data.frame(atm.NO3.mean = sapply(no3.atm.conc, mean), 
                           atm.NO3.sd = sapply(no3.atm.conc, sd),
                           Site = sapply(no3.atm.conc, names)[1,1:14])



#4 Estimate processing efficiency and uncertainty using D17O and uncertainty
#with total No3 yields and wet NO3 deposition

#To do this, we first must estimate yields of NO3-Atm in kg N/ha/yr, then use
#this distribution along with rates of No3 deposition in the same units
#to estimate processing efficiency

#import NO3 total load/yield data
#now for flux
#same as above for NO3 total concentration
file_names <- list.files("path to folder containing No3 total load files here", full.names = T)
file_names <- list.files("/Users/joelbostic 1/Desktop/PhD/Dissertation/ch1/R/WRTDS_K_uncert_output/flux/", full.names = T)

#import files
no3.tot.flux <- do.call(rbind,lapply(file_names,read.csv))

#assign names
site.nm <- data.frame(id = unique(no3.tot.flux$id), site.new = c("ANT2", "ANT", "BIGR", "BLAC", "CAC", "CON", "DPRN", 
                                                                 "GEO", "GUN", "GWN", "MON", "TERR", "TOW", "WIL"))
no3.tot.flux <- left_join(no3.tot.flux, site.nm, by = "id")

#create matrix that contains water years 2016-2017 average flux
#note that these are loads in units of kg/day
#we'll re-calculate to yields in units of kg/ha/yr below
no3.total.flux <- matrix(ncol = 801, nrow = 14)
for(j in 1:14){
  for(i in 2:801){
    no3.total.flux[j, i] <- mean(no3.tot.flux[which(no3.tot.flux$site.new == unique(no3.tot.flux$site.new)[j]), i])
  }
}
no3.total.flux <- as.data.frame(no3.total.flux)
no3.total.flux$site.new <- unique(no3.tot.flux$site.new)

#dataframe of site ID, watershed area in hectares used to convert loads to yields,
#and average wet NO3 deposition between WY 2016-2017 needed to calculate processing efficiency
site.area <- data.frame(Site = site.nm$site.new, 
                        ws.area = c(24200, 72800, 160, 560, 17300, 127900, 1620, 18800,
                                    41400, 8400, 44800, 570, 38300, 64000),
                        no3.dep = c(1.66, 1.54, 1.77, 1.67, 1.55, 1.65,
                                    1.47, 1.66, 1.63, 1.55, 1.63, 1.60, 1.53, 1.65))

#create list to hold processing efficiency data
PE <- list()

#note that this will take 10ish minutes to run
#calculate NO3-Atm yields by sampling from distribution of total No3 loads and fraction of atmospheric NO3
for(i in 1:length(no3.atm.conc)){
  temp <- vector()
  
  for(j in 1:iterations){
    
    #temp stores loads - converted from kg/d to kg/year
    temp[j] <- sample(as.numeric(no3.total.flux[no3.total.flux$site.new == names(f.atmlist[[i]][1]), 2:801]), 1, replace = T)*
      sample(f.atmlist[[i]], 1, replace = T)*365
    
  }
  #convert loads to yields (kg N/ha/yr) then calculate PE
  PE[[i]] <- (1-((temp/site.area$ws.area[i])/site.area$no3.dep))*100
  names(PE[[i]]) <- names(f.atmlist[[i]][1])
} 


#mean and sd of PE can be displayed by
data.frame(PE.mean = sapply(PE, mean), 
           PE.sd = sapply(PE, sd),
           Site = sapply(PE, names)[1,1:14])




