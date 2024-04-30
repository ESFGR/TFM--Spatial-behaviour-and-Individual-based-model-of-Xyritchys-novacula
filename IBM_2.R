#Individual based model Raor #2

#Parameters definition taking results from  Dimitriadis & Fournari-Konstantinidou (2018)
#Length-Weight relationship --> W = a*L^b (weight in g)
a <- 0.00955
b <- 3.00

# Parameters for the Von Bertalanffy model for length growth
k <- 1.079     # Growth rate 
Linf= 175 #  : 
t0= -0.73
agemax=7
t=seq(0,agemax, 1)

lt=Linf*(1 - exp(-k *(t-t0)))

plot(t, lt, type = "l", xlab = "Age (years)", ylab = "Lenght (mm)")
#Mortality----
f=-0.05
z=1.3
mortn=exp(-z*t)
mortf=exp(-f*t)
plot(t, mort, type = "l", xlab = "Age (years)", ylab = "Mortality (prob)")
plot(t, mortf, type = "l", xlab = "Age (years)", ylab = "Fishing Mortality (prob)")

pmn=0.008 #natural mortality parameter
pmf=0.002 #fishing mortality parameter


#Mortality functions:
# Beverton and Holt natural mortality model function
beverton_holt_Mn <- function(L, M) {
  # L: Length in cm
  # M: natural mortality parameter
  return(M * L / (1 + M * L))
}

# Beverton and Holt total mortality (z) model function
beverton_holt_Z <- function(L_inf, k, L, Lmin) {
  
  #L_inf and K = Von Bertalanffy parameters
  #L = mean lenght of fish
  #Lmin=min length under exploitation
 
   Z=k*((L_inf-Lmin)/(L-Lmin))
  
  return(Z)
} #this is not currently applied in the simulation


#Build matrix----
M<-matrix(NA, 20, 7)
M[,1]=0.1 #population size matrix, initial size=0.1cm

remove(i)

for (i in 2:7) {
  
  t0=i-1 #previous age class
  
  Lt<-Linf*(1 - exp(-k * i)) #lenght at age i; von Berttalannfy function
  gt1=20-sum(is.na(M[,t0]==FALSE)) #n of individuals at previos age class
  nextgen<-rep(Lt, gt1) #vector of n individual lenghts at age class i before appliying mortality
  
  # Apply natural mortality using Beverton and Holt model
  Mn <- beverton_holt_Mn(Lt, pmn)
  
  # Apply fishing mortality starting from the third year
  if (i > 2) { #it can also be Lt>20cm to apply fihing pressure from given length
    Mf <- pmf*gt1 #fishing mortality parameter is density dependent
  } else {
    Mf <- 0
  }
  
  #nº survivals at age class i
  surv=rbinom(1, gt1, (1-(Mn+Mf)))
  
  nextgen[sum(surv):gt1]<-NA #eliminate dead indiviuduals with NAs
  
  M[,i]<-c(nextgen, rep(NA, 20-length(nextgen))) #fill with NAs to length 20 and join the matrix
  
}

M

#Biomass
B<-rep(NA, 7)

for (i in 1:7) {
  
  L=sum(na.omit(M[,i])) #total length of age class i
  B[i]<-a*L^b #LWR formula

}

data.frame(B) #biomass (g) per class age

