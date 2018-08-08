########################################################################
# Emanuel A. Fronhofer et al.
# Bottom-up and top-down control of dispersal across major organismal groups
#
# food chain model
# 
# August 2018
########################################################################
# clear workspace
rm(list=ls())

########################################################################
# load required packages
library(deSolve)

########################################################################
# define 2-patch food chain model
spatial_food_chain_model <- function(time,y,parms){
  with(as.list(c(y,parms)), {
    
    # if dispersal is CDD
    if(cdd == T){
      # set dispersal rate to zero
      m_N1 <- 0
      m_N2 <- 0
      # reset dispersal rate according to CDD rule
      if(R1 <= cRd | P1 >= cPd){
        m_N1 <- 1
      }
      if(R2 <= cRd | P2 >= cPd){
        m_N2 <- 1
      }
    }else{
      # else dispersal is RD
      m_N1 <- m_N
      m_N2 <- m_N
    }
    
    # PATCH 1
    # focal consumer dynamics: LV model with type I FR
    dN1 <- e_N * a_N * R1 * N1 - d_N * N1 - a_P * N1 * P1 - m_N1 * N1 + (1-dm_N)*m_N2 * N2
    
    # bottom up --- resource dynamics following a chemostat model
    dR1 <- omega * R0 - omega * R1 - a_N * R1 * N1 - m_R * R1 + (1-dm_R)*m_R * R2
    
    # top down --- top predator dynamics identical to focal consumer for simplicity
    dP1 <- e_P * a_P * N1 * P1 - d_P * P1 - m_P * P1 + (1-dm_P)*m_P * P2
    
    # PATCH 2
    # focal consumer dynamics: LV model with type I FR
    dN2 <- e_N * a_N * R2 * N2 - d_N * N2 - a_P * N2 * P2 - m_N2 * N2 + (1-dm_N)*m_N1 * N1
    
    # bottom up --- resource dynamics following a chemostat model
    dR2 <- omega * R0 - omega * R2 - a_N * R2 * N2 - m_R * R2 + (1-dm_R)*m_R * R1
    
    # top down --- top predator dynamics identical to focal consumer for simplicity
    dP2 <- e_P * a_P * N2 * P2 - d_P * P2 - m_P * P2 + (1-dm_P)*m_P * P1
    
    # retun function output
    list(c(dR1, dN1, dP1, dR2, dN2, dP2))
  })
}

########################################################################
# define parameters
parms <- c(
  omega = 0.5, # flow rate of abiotic resource
  R0 = 1000, # external resource density for inflow
  m_R = 0, #d ispersal rate for resource
  dm_R = 0, # dispersal cost for resource
  e_N = 0.1, # assimilation efficiency of focal species
  a_N = 0.01, # feeding rate of focal species
  d_N = 0.1, # mortality of focal species
  cdd = F, # swith that determines whether dispersal is CDD or RD
  m_N = 0, # dispersal rates of focal species if RD
  cRd = NA, # CDD resource thershold 
  cPd = NA, # CDD predation thershold
  dm_N = 0, # dispersal cost for focal species
  e_P = 0.005, # assimiliation efficiency for top predator
  a_P = 4, # feeding rate for top predator
  d_P = 0.1, # mortality for top predator
  m_P = 0, # dispersal rates for top predator
  dm_P = 0 # disperal mortality for top predator
)

# define initial conditions: patches start at different starting points
cinit=c(R1 = parms[["R0"]], N1 = 0.1, P1 = 0.1, R2 = parms[["R0"]], N2 = 1, P2 = 1)

# define length of modelled dynamics
t <- seq(0,100,len=1000)

########################################################################
# run model to get equilibirum densities
t1 <- c(0:1000)
out1 <- ode(y=cinit, times=t1, func=spatial_food_chain_model, parms=parms)

########################################################################
# function to obtain population density CV for different dispersal rates for RD
CV_RD <- function(disprate){
  
  # set disprate to correct value
  parms[["cdd"]] = FALSE
  parms[["m_N"]] = disprate
  
  out <- ode(y=cinit, times=t, func=spatial_food_chain_model, parms=parms)
  
  m <- mean(c(out[,"N1"],out[,"N2"]))
  cv <- sd(c(out[,"N1"],out[,"N2"]))/m
  return(c(m=m,cv=cv))
}

# run
all_dispRates <- seq(0,1,len=50)
res_all_RD <- sapply(all_dispRates,CV_RD)

########################################################################
# function to obtain population density CV for different dispersal rates for CDD
CV_CDD <- function(threshold_R,thershold_P){
  
  # set disprate to correct value
  parms[["cdd"]] = TRUE
  parms[["cRd"]] = threshold_R * mean(c(out1[length(t1),"R1"],out1[length(t1),"R2"]))
  parms[["cPd"]] = thershold_P * mean(c(out1[length(t1),"P1"],out1[length(t1),"P2"]))
  
  out <- ode(y=cinit, times=t, func=spatial_food_chain_model, parms=parms)
  
  m <- try(mean(c(out[,"N1"],out[,"N2"])),silent=T)
  
  cv <- try(sd(c(out[,"N1"],out[,"N2"]))/m,silent=T)
  
  if(!is.numeric(m)){m <- NA}
  if(!is.numeric(cv)){cv <- NA}
  return(c(m=m,cv=cv))
}

# run
all_thresholds_R <- seq(0,2,len=20)
all_thresholds_P <- seq(0,5,len=20)

all_combs <- expand.grid(all_thresholds_R,all_thresholds_P)

res_all_CDD <- do.call(mapply,c(CV_CDD,unname(all_combs)))

########################################################################
# redefine parameters for specific model run with CDD (relative to equil dens of pred and resource)
# this is the parameter combination that minimizes the CV
parms[["cRd"]] = all_combs[which.min(res_all_CDD["cv",]),1] * mean(c(out1[length(t1),"R1"],out1[length(t1),"R2"]))
parms[["cPd"]] = all_combs[which.min(res_all_CDD["cv",]),2] * mean(c(out1[length(t1),"P1"],out1[length(t1),"P2"]))
parms[["cdd"]] = T

########################################################################
# run model with CDD
out_minCV_CDD <- ode(y=cinit, times=t, func=spatial_food_chain_model, parms=parms)

########################################################################
# redefine parameters for specific model run with RD
# this is the parameter combination that minimizes the CV
parms[["cdd"]] = F
parms[["m_N"]] = all_dispRates[which.min(res_all_RD["cv",])]

########################################################################
# run model with CDD
out_minCV_RD <- ode(y=cinit, times=t, func=spatial_food_chain_model, parms=parms)

########################################################################
# plotting all results
x11(width=5,height=8)
par(mfrow=c(3,1),bty="l",mar=c(0.5,0.5,0.5,0.5),oma=c(4,4,0.5,2),cex=1)

# top predator dynamics
plot(out_minCV_CDD[,"P1"]~out_minCV_CDD[,"time"],type="n",ylim=range(c(out_minCV_RD[,"P1"],out_minCV_RD[,"P2"],out_minCV_CDD[,"P1"],out_minCV_CDD[,"P2"])),xaxt="n")
lines(out_minCV_RD[,"P1"]~out_minCV_RD[,"time"],lwd=2,col="lightpink")
lines(out_minCV_RD[,"P2"]~out_minCV_RD[,"time"],lwd=2,col="lightpink",lty=2)
lines(out_minCV_CDD[,"P1"]~out_minCV_CDD[,"time"],lwd=2,col="red")
lines(out_minCV_CDD[,"P2"]~out_minCV_CDD[,"time"],lwd=2,col="red",lty=2)
legend("topright",fill=c("red","lightpink"),legend=c("CDD","RD"),bty="o",bg="white",box.col="white")

# focal species dynamics
plot(out_minCV_CDD[,"N1"]~out_minCV_CDD[,"time"],type="n",ylim=range(c(out_minCV_RD[,"N1"],out_minCV_RD[,"N2"],out_minCV_CDD[,"N1"],out_minCV_CDD[,"N2"])),xaxt="n")
lines(out_minCV_RD[,"N1"]~out_minCV_RD[,"time"],lwd=2,col="grey")
lines(out_minCV_RD[,"N2"]~out_minCV_RD[,"time"],lwd=2,col="grey",lty=2)
lines(out_minCV_CDD[,"N1"]~out_minCV_CDD[,"time"],lwd=2,col="black")
lines(out_minCV_CDD[,"N2"]~out_minCV_CDD[,"time"],lwd=2,col="black",lty=2)
legend("topright",fill=c("black","grey"),legend=c("CDD","RD"),bty="o",bg="white",box.col="white")

# resource dynamics
plot(out_minCV_CDD[,"R1"]~out_minCV_CDD[,"time"],type="n",ylim=range(c(out_minCV_RD[,"R1"],out_minCV_RD[,"R2"],out_minCV_CDD[,"R1"],out_minCV_CDD[,"R2"])))
lines(out_minCV_RD[,"R1"]~out_minCV_RD[,"time"],lwd=2,col="lightblue")
lines(out_minCV_RD[,"R2"]~out_minCV_RD[,"time"],lwd=2,col="lightblue",lty=2)
lines(out_minCV_CDD[,"R1"]~out_minCV_CDD[,"time"],lwd=2,col="darkblue")
lines(out_minCV_CDD[,"R2"]~out_minCV_CDD[,"time"],lwd=2,col="darkblue",lty=2)
legend("topright",fill=c("darkblue","lightblue"),legend=c("CDD","RD"),bty="o",bg="white",box.col="white")

mtext(side=1,at=0.5,outer=T,"Time",line=2.5,cex=1.2)
mtext(side=2,at=0.5,outer=T,"Population size",line=2.5,cex=1.2)
mtext(side=4,at=1/6,outer=T,"Resource",line=0,cex=1.2)
mtext(side=4,at=3/6,outer=T,"Focal species",line=0,cex=1.2)
mtext(side=4,at=5/6,outer=T,"Predator",line=0,cex=1.2)

#dev.copy2pdf(file="dynamics_all.pdf")

########################################################################
# calculate putput of CV and COV reduction
# get CVs for focal species
cv_RD_N <- min(res_all_RD["cv",])
cv_CDD_N <- min(res_all_CDD["cv",],na.rm=T)

# get CVs for resource
cv_RD_R <- sd(c(out_minCV_RD[,"R1"],out_minCV_RD[,"R2"]))/mean(c(out_minCV_RD[,"R1"],out_minCV_RD[,"R2"]))
cv_CDD_R <- sd(c(out_minCV_CDD[,"R1"],out_minCV_CDD[,"R2"]))/mean(c(out_minCV_CDD[,"R1"],out_minCV_CDD[,"R2"]))

# get CVs for top predator
cv_RD_P <- sd(c(out_minCV_RD[,"P1"],out_minCV_RD[,"P2"]))/mean(c(out_minCV_RD[,"P1"],out_minCV_RD[,"P2"]))
cv_CDD_P <- sd(c(out_minCV_CDD[,"P1"],out_minCV_CDD[,"P2"]))/mean(c(out_minCV_CDD[,"P1"],out_minCV_CDD[,"P2"]))

# get covariances between the patches for the focal species
cov_RD_N <- cov(out_minCV_RD[,"N1"],out_minCV_RD[,"N2"])
cov_CDD_N <- cov(out_minCV_CDD[,"N1"],out_minCV_CDD[,"N2"])

# get covariances between the patches for the resource
cov_RD_R <- cov(out_minCV_RD[,"R1"],out_minCV_RD[,"R2"])
cov_CDD_R <- cov(out_minCV_CDD[,"R1"],out_minCV_CDD[,"R2"])

# get covariances between the patches for the top predator
cov_RD_P <- cov(out_minCV_RD[,"P1"],out_minCV_RD[,"P2"])
cov_CDD_P <- cov(out_minCV_CDD[,"P1"],out_minCV_CDD[,"P2"])

########################################################################
# output of relative reduction of CV and COV
print(paste("CV reduction N: ",round((cv_RD_N - cv_CDD_N)/cv_RD_N*100,digits=0)," %"))
print(paste("CV reduction R: ",round((cv_RD_R - cv_CDD_R)/cv_RD_R*100,digits=0)," %"))
print(paste("CV reduction P: ",round((cv_RD_P - cv_CDD_P)/cv_RD_P*100,digits=0)," %"))

print(paste("COV reduction N: ",round((cov_RD_N - cov_CDD_N)/cov_RD_N*100,digits=0)," %"))
print(paste("COV reduction R: ",round((cov_RD_R - cov_CDD_R)/cov_RD_R*100,digits=0)," %"))
print(paste("COV reduction P: ",round((cov_RD_P - cov_CDD_P)/cov_RD_P*100,digits=0)," %"))

########################################################################