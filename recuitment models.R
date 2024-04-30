#Beverton-Holt recruitment model
#https://bookdown.org/Joseph_Hightower/bayesianfish/Recruit.html
#A Bayesian Introduction to Fish Population Analysis Joseph E. Hightowe
#The process of recruitment starts with eggs, but egg production is difficult to estimate directly so it is common to use spawning stock biomass (Bs) as a proxy. 
#One of the most commonly used curves relating spawning stock size and recruitment is the Beverton-Holt equation


SpBio <- seq(from=10, to=200, by=5) #stock biomass
a <- 0.001 #recruitmen parameters
b <- 0.1
R <- 1/(a+b/SpBio) #Beverton-Holt equation
plot(SpBio, R)

#This curve approaches zero as Bs approaches zero and approaches 1/?? as Bs approaches infinity. 
#It is often a good descriptor of observed stock-recruitment data in that recruitment is low when spawning stock size is close to zero but appears to be without trend at higher spawning stock sizes.


# fitting a Beverton-Holt model to simulated stock-recruitment data -----

rm(list=ls()) # Clear Environment

A <- 10 # max.age - Arbitrary. If changed, redo maturity and meanwt vectors
Y <- 12 # year - Arbitrary. If changed, redo AnnualF vector

AnnualF <- c(0.59, 0.41, 0.74, 0.91, 0.86, 0.74, 1.07, 0.9, 0.87, 1.1, 0.93)
# Arbitrary, increasing trend to create contrast in SpBio. Redo if Y changes
AnnualM <- 0.2

# F by age, using separable model:
#Obtaining year- and age- specific fishing mortality rates as a product of a year-specific measure of fishing intensity and an age-specific vulnerability to fishing
Fishery_k <- 1 # Slope for logistic function
Fishery_a_50 <- 2 # Age at 0.5 selectivity
SelF <- array(data=NA, dim=A) # Fishery selectivity pattern
F_a <- array(data=NA, dim=c(A,(Y-1))) # Matrix for F by age and year
Z_a <- array(data=NA, dim=c(A,(Y-1))) # Matrix for Z by age and year
S_a <- array(data=NA, dim=c(A,(Y-1))) # Matrix for survival rate by age and year

for (a in 1:A){
  SelF[a] <- 1/(1+exp(-Fishery_k*(a-Fishery_a_50)))
} #a

plot(SelF, main="fishing mortality rates by age")

for (y in 1:(Y-1)) {
  for (a in 1:A){
    F_a[a,y] <- AnnualF[y] * SelF[a]
    Z_a[a,y] <- F_a[a,y]+AnnualM
    S_a[a,y] <- exp(-Z_a[a,y])
  } #a
} #y

N <- array(data=NA, dim=c(A, Y))

BH_alpha <- 1/1E5 # Parameter defining Beverton-Holt asymptotic recruitment
BH_beta <- 0.1 # Arbitrary
SD_rec <- 0.25 # Low level of lognormal error in recruitment

# Set up year-1 vector based on asymptotic recruitment and year-1 survival rates
Mat <- c(0, 0, 0.2, 0.4, 0.8, 1, 1, 1, 1, 1) # Arbitrary maturity schedule. Redo if max age (A) changes
Wgt <- c(0.01, 0.05, 0.10, 0.20, 0.40, 0.62, 0.80, 1.01, 1.30, 1.56) # Mean wt. Redo if A changes
SpBio <- array(data=NA, dim=(Y-1)) # Spawning biomass vector
# Year-1 N (arbitrarily) from asymptotic recruitment, year-1 survival.
N[1,1] <- 1/BH_alpha*rlnorm(n=1, 0, SD_rec)
for (a in 2:A){
  N[a,1] <- rbinom(n=1, trunc(N[(a-1),1]), S_a[(a-1),1])
}#a
for (y in 2:Y){
  SpBio[y-1] <- sum(N[,y-1]*Wgt*Mat) # Spawning stock biomass
  N[1,y] <- 1/(BH_alpha+BH_beta/SpBio[y-1])*rlnorm(n=1, 0, SD_rec) #Beverton-Holt ecuation
  for (a in 2:A){
    N[a,y] <- rbinom(1, trunc(N[(a-1), (y-1)]), S_a[(a-1),(y-1)])
  } #a
} #y
plot(SpBio,N[1,2:Y], ylab="Recruits", xlab="Spawning stock biomass")

# Generate survey data on recruitment and spawning stock
Survey_q <- 1E-3 # Arbitrary catchability coefficient for survey
SD_survey <- 0.2 # Arbitrary ln-scale SD
Exp_ln_B <- log(Survey_q*SpBio)
Survey_B <- rlnorm(n=Y, meanlog=Exp_ln_B, sdlog=SD_survey)

plot(SpBio, Survey_B[1:Y-1])

Exp_ln_R <- log(Survey_q*N[1,2:Y])
Survey_R <- c(NA,rlnorm(n=(Y-1), meanlog=Exp_ln_R, sdlog=SD_survey))

plot(N[1,], Survey_R)



#Size-based egg-per-recruit model ()-----

#
#formula three-parameter Gaussian curve (model parameters are A, u, and sigma2) assuming a lognormal error structure:
#E(L)=Ae^[-(L-u)^2]/(s2)
  
  
L =c(0.00000, 11.33991, 15.19475, 16.50514, 16.95059, 17.10202, 17.15349, 17.17099)
A= 1000
u=5 #size start producing eggs
s2= 0.08

n.egg=A*exp((-(L-u)^2)/s2)

plot(n.egg)

#----
# Modelo de madurez sexual:

#probability of sexual maturity at length Lt

s=0.5#slope parameter
L50 = 10 # size 50% fish are mature
L =c(0.00000, 11.33991, 15.19475, 16.50514, 16.95059, 17.10202, 17.15349, 17.17099) #length at age class form 0 to 7

pm=1/(1+(exp(-s*(L-L50))))
plot(pm, xlab = "age class", ylab = "p(sexual maturity)")

#recruitment model n of eggs at length-----

#stares et al 2007 (study on atlantic cod)
#The number of eggs produced per recruit was calculated as ameasure of individual egg production rate

b=0.1
a=100 #eggs per length production paramter 

Fa=a*L^b
plot(Fa)

Na=c(100, 60, 45, 38, 30, 25, 20, 18) #where Na is numbers-at-age (Na=N(a-1)*e^-Z) starting with 1recruit at age 0, jthe maximum age in the age-structured popu-lation analyses used in the TEP estimates above 

#The number of eggs produced per recruit calculated as a measure of individual egg production rate
Ne=Na*pm*Fa

plot(Ne)


#A*e^bL
#A and b parameters affecting fecundity rate 

A=0.8
b=0.7

L=seq(0, 20, 1)

n.egg = A*exp(b*L)
plot(n.egg, )

