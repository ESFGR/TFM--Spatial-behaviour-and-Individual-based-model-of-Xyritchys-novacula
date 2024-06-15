#Individual based model Raor #7
#integrating natural - fishing mortality models (by Miquel Palmer)
#adding simulated behaviour (home range HR and Chronotype) induced mortality

#Author: Pep Alós & Enrique Sánchez-Fabrés
#Project Metaraor
#Supervisor: Josep Alós

#Last edit: 16/05/2024

rm(list = ls())
#Formulas - functions ----

##Von Bertalanffy model for length growth: L_inf: asintotic length, K: growth rate (y^-1), age, age0=t0=null growth theoric age
von_ber_fun<-function(L_inf, K, AGE, AGE0){ 
  L_inf*(1 - exp(-K *(AGE-AGE0)))
}

plot(x=seq(0, 7, 1), von_ber_fun(17.87, 1.079, seq(0, 7, 1), 0), ylab = "Length (cm)", xlab="Age (years)")

#LWR
lwr_fun<- function(a, b, L){
  a*L^b
}

####Mortality ###
#######

####NAtural mortality ######

# Gislason's model 
M_fun <- function(M, A, s) { #the natural mortality at age A of the individual in years. 
  M*exp(-s * (A))  # Inverse logistic function
}

#M is a baseline mortality rate adjusted for specific species and environmental conditions. s: parameter controlling the rate of decrease in natural mortality, which in Gislason et al., (2010) is given a value of s=0.06. 

#plot(M_fun(0.1, seq(0, 7, 1), 0.06), ylab = "Natural mortality prob", xlab = "Age")

#Lorenzen (1996): natural mortality based on length, with von bertalanffy paramters
M_funL <- function(M, L_inf, K){ 
  M*(L_inf^(1/3))*K
}

### Behaviour induced natural mortality (home range and chronotype)###
# Shr and Schron are parameters controlling the rate of change in mortality due to HR and chronotype individual scores, respectively.
M_hr_fun<-function(M, HR, Shr){
  M * (1 - exp(-Shr * HR))
}
#plot(x=seq(0,100, 1), M_hr_fun(0.1,HR=0.1*(seq(0,100, 1)), 0.2), xlab = "HR score", ylab = "HR natural mortality rate")

M_chron_fun<-function(M, Chron, Schron){
  M * (1 - exp(-Schron * Chron))
}

##length dependent fishing mortality: f(length)

####f(length)
f_length_fun<-function(s, L, Lx){ #L=lengths (from growth model), s & Lx= slope and length of inflexion for fishing mortality function
  1-(1/(1 + exp (s * (L - Lx))))
}

### Behaviour induced fishing mortality (activity and chronotype)
####f(beha)
f_beha_fun<-function(s, Behavoiur_score, be_in){ #L=lengths (from growth model), s & be_in= slope and behaviour of inflexion for behaviour induiced fishing mortality function
  1-(1/(1 + exp (s * (Behavoiur_score - be_in))))
}
#####f(chron)
f_chron_fun<-function(s, crhonotype_score, be_in){ #L=lengths (from growth model), s & be_in= slope and behaviour of inflexion for behaviour induiced fishing mortality function
  1-(1/(1 + exp (s * (crhonotype_score - be_in))))
}

#lt1=von_ber_fun(linf1, k1, age, age01); lt1[1]<-1
#testM<-M_fun(age, lt1)
#plot(prob_M_fun(testM,age,lt1), main = "Natural surv. prob")

#Fishing mortality
prob_F_fun<-function(F_max, F_length, F_beha, F_chron, AGE){ ## F_max=maximum fishing mortality rate; M=natural mortality
    1-(exp(-(F_max*(F_length+F_beha+F_chron))*(AGE)))         
}

#test1=(1-exp(-(0.1*(0.5+0.5+0.5))*seq(0,7, 1))) 
#plot(x=seq(0,7, 1), test1, xlab = "age", main = "Overall fishing mortality rate")
#natural mortality surv prob

prob_M_fun<-function(M_max, M, M_hr, M_chron){ ## M_max=maximum natural mortality rate; M=natural mortality
  1-(exp(-(M_max*(M+M_hr+M_chron))))         
}

#duda: incluir la edad? como en la funcion de mortalidad por pesca? de forma inversa?

#test1=(1-exp(-(0.1*(0.5+0.1+0.2))))
#plot(test1, xlab = "age")


####Recruitment model: Beverholt
recruit_fun<- function(alfa, beta, SpBio){ #alpha is the maximum recruitment rate achievable (denisty independent paramters) 
  alfa*SpBio/1+(beta*SpBio) #beta  reflects the strength of density-dependent effects
}


#Population parameters----
ind=100 #number of initial indviduals
agemax=7
age= seq(0,agemax, 1) #age classes
#LWR paramters: males a=0.002; b=3.383 | Females a=0.003; b= 3.32
a1=0.0025
b1=3.35

#Von Bertalanffy parameters (Mekni et al., 2018)
linf1= 17.18 #cm assympotic length
k1= 1.079 #growth rate
age01= 0#null growth theoric age

#Natural mortality paramters
Mslope=0.3 # parameter for disminution rate of natural mortality
#M <- 0.002  # natural mortality max value

#f(length) parameters
fslope1=0.5 #slope and inflection of a logistic function (or any other function)
#lt1=#obtained from Von Bert. model
lx=10 #length at which fishing mortality starts

#f(behaviour) paramters
hr_inf=0.60 #behaviour score at which fishing mortality starts
ht_inm=0.60  #behaviour score at which natural mortality starts
hr_slopef=0.1 #slope for behaviour induced fishing morality
hr_slopem=0.1 #slope for behaviour induced natural morality

#f(chron) paramters
chron_inf=0.80 #chron score at which fishing mortality starts
chron_inm=0.70  #chron score at which natural mortality starts
chron_slopef=0.1 #slope for chron induced fishing morality
chron_slopem=0.1 #slope for chron induced natural morality

#simulating behaviour scores
hr=runif(ind, 0, 1)
plot(hr)
summary(hr)
boxplot(hr)

#simulating behaviour scores
chrtype=runif(ind, 0, 1)
plot(chrtype)
summary(chrtype)
boxplot(chrtype)

#surv. prob. paramters 
Mmax=0.3#baseline natural mortality rate 
fmax1=0.1 # maximum fishing mortality rate 

#Recruitment model
ar <- 0.1 #recruitmen density independent parameter proportional to fecundity
#The units of ar are "recruitment per spawner" 
br <- 2 #b is a density-dependent parameter that is proportional to both fecundity and density-dependent mortality
#The Beverton-Holt model is based on the assumptions that juvenile competition results in a mortality rate that is linearly dependent upon the number of fish alive in the cohort at any time and that predators are always present. The Beverton-Holt model is appropriate "if there is a maximum abundance imposed by food availability or space, or if the predator can adjust its predatory activity immediately to changes in prey abundance" (Wootton 1990, 
biomassr<- 5 #biomass per recruit
srvr<-0.1 #surv prob of recruits



#R <- a*SpBio/((1+(b*SpBio))) #Beverton-Holt equation
#plot(SpBio, R)

####### code pep IBM

#numbe rof individuals
ind=100
#years to simualte
IniYear=2024
Years=10
#ind adultos
lt=von_ber_fun(linf1, k1, age, age01) #vector of lengths at age according VB equation and paramters
wt=lwr_fun(a1, b1, lt)

F<-data.frame(L=rep(0.1,ind),W=rep(biomassr,ind),A=rep(0,ind),ID=round(seq(1:ind)))
assign(paste("F",IniYear,sep=""),F)
str(F2024)
str(F)
range(F$ID)
#ind adultos parametros
par<-data.frame(k=rep(k1,ind),hr=hr,chrtype=chrtype,ID=round(seq(1:ind)))
assign(paste("par",IniYear,sep=""),par)
str(par2024)


#Bioyear
BioyearAll<-sum(F$W)

for (years in c(IniYear+1):c(IniYear+Years)){
  #years=2026
  
  Ftemp<-get(paste("F",years-1,sep=""))
  partemp<-get(paste("par",years-1,sep=""))
  
  #remove vectors from previous year
  rm(list=paste("F",years-1,sep=""))
  rm(list=paste("par",years-1,sep=""))
  #for new ids
  #range(Ftemp$ID)
  IDRini<-max(Ftemp$ID)+1
  
  #survival to natural mortality
  M<-M_fun(Mmax, s = Mslope, A = Ftemp$A)
  str(Ftemp$L)
  Mhr<-M_hr_fun(Mmax, partemp$hr, hr_slopem)
  Mchron<-M_chron_fun(Mmax, partemp$chrtype, chron_slopem)
  
  #survival to fihisng mortality
  fl<-f_length_fun(s = fslope1 ,L = Ftemp$L, Lx = lx)
  
  fhr<-f_beha_fun(s = 0.01,Behavoiur_score = partemp$hr,be_in = hr_inf)
  
  fchron<-f_chron_fun(s = chron_slopef,crhonotype_score = partemp$chrtype,be_in = chron_inf)
  
  
  #assgining survival probabilities
  #surv.prob.M<-runif(nrow(Ftemp), 0.9, 1)
  surv.prob.M<-1-prob_M_fun(Mmax, M = M, Mhr, Mchron) # como aplicar mortalidad a edad 0?
  #plot(surv.prob.M)

  surv.prob.M[Ftemp$A==0] <- runif(length (Ftemp$A[Ftemp$A==0]), 0.7, 0.99)##asignacion aleatoria de prob supevivencia a reclutas de edad 0
  #range(surv.prob.M)
  
  surv.prob.f<-1-prob_F_fun(fmax1, fl, fhr,fchron,Ftemp$A)
  #range(surv.prob.f)
  #surv.prob.f<-runif(nrow(Ftemp), 0.9,1)
  surv.prob.f[Ftemp$L <= lx] <- 1 #fishing mortality only applies from min fishing length
  #range(surv.prob.f)
  
  #apploying survival probabilities
  survM_status <- rbinom(n = nrow(partemp), size = 1, prob = surv.prob.M) #binomial varaible for each individual (0=die, 1=survive)
  survM_index <- which(survM_status == 1) #extract index of survivals
  partemp <- partemp[survM_index, ] #substracts death individuals
  Ftemp <- Ftemp[survM_index, ]

  survF_status <- rbinom(n = nrow(partemp), size = 1, prob = surv.prob.f) #binomial varaible for each individual (0=die, 1=survive)
  survF_index <- which(survF_status == 1) #extract index of survivals
 partemp <- partemp[survF_index, ] #substracts death individuals
  Ftemp <- Ftemp[survF_index, ]
  
  partemp <- partemp[which(Ftemp$A<=agemax),]
  Ftemp <- subset(Ftemp,A<=agemax) #ind older than age max die
  
  #calculating n of survivers for each age class
  #for (i in 1:length(age)){ # time steps = years from 0 to 8
    #prob1=
    #s2=1-prob2[i+1]
    #s1=1-prob2[i]
    #p=(s1-s2)/s1 #-- # p for not surviving at t+1 given you did survive at t
    #temp=sum(runif(n1[i],0,1)>p)
    #n1=c(n1,temp) # p for not surviving at t+1 given you did survive at t
  #}
  
  #age
  Ftemp$A=Ftemp$A+1
  #VB
  Ftemp$L=von_ber_fun(linf1,k1,AGE=Ftemp$A, age01)
  #weight - LWR function
  Ftemp$W=lwr_fun(a1,b1,Ftemp$L)
  #Biomas of the year
  #table(Ftemp$W)
  SpBio<-subset(Ftemp, Ftemp$A>1)$W
  Bioyear<-sum(SpBio)
  BioyearAll<-c(BioyearAll,Bioyear)
  #R of the next year
  rbiomass<-round(recruit_fun(ar,br,Bioyear)) #biomass of recruits
  R<-round(rbiomass/biomassr) #number of recruits

  
  #recxruits of the next year
  if(R>0){  
  FtempR<-data.frame(L=rep(0.1,R),
                     W=rep(br,R), 
                     A=rep(0,R), #age= -1?)
                     ID=seq(IDRini,(IDRini+R-1)))
  
  #str(FtempR)
  #Recruits parameters
  #simulating behaviour scores
  hrR=runif(R, 0, 1)
  chrtypeR=runif(R, 0, 1)
  
  parR<-data.frame(k=rep(k1,R),
                   hr=hrR,
                   chrtype=chrtypeR, 
                   ID=seq(IDRini,IDRini+R-1))
  #join DBs
  Ftemp<-rbind(Ftemp,FtempR)
  
  partemp<-rbind(partemp, parR)
  }
  #str(partemp)
  #str(parR)
  #save for the year
  assign(paste("F",years,sep=""),Ftemp)
  assign(paste("par",years,sep=""),partemp)
  
}

plot(seq(IniYear,IniYear+length(BioyearAll)-1),BioyearAll,type="l", xlab = "Year", ylab = "Biomass")



plot(BioyearAll)

boxplot(surv.prob.f)
boxplot(surv.prob.M)

table(F2040$A)
table(F2059$A)
table(F2059$W)
w2050<-sum(F2050$W)
table(F2070$W)
boxplot(F2028$A)
w2070<-sum(F2070$W)

recruit
recruit_fun(10, -0.9, 1000)/5
recruit_fun(0.5, 0.8, 100)
#in recl año +1
#L2024=rep(1,10000)
#W2024=rep(10,10000)
#A2024=rep(-1,1)

table(F2030$A)
table(F2031$A)
table(F2032$A)
table(F2035$A)
table(F2034$A)

fmax2=0.5
plusyears=IniYear+Years
for (years in c(plusyears+1):c(plusyears+Years)){
  #years=2032
  
  Ftemp<-get(paste("F",years-1,sep=""))
  partemp<-get(paste("par",years-1,sep=""))
  
  #for new ids
  #range(Ftemp$ID)
  IDRini<-max(Ftemp$ID)+1
  
  #survival to natural mortality
  M<-M_fun(Mmax, s = Mslope, L = Ftemp$L)
  str(Ftemp$L)
  Mhr<-M_hr_fun(Mmax, partemp$hr, hr_slopem)
  Mchron<-M_chron_fun(Mmax, partemp$chrtype, chron_slopem)
  
  #survival to fihisng mortality
  fl<-f_length_fun(s = fslope1 ,L = Ftemp$L, Lx = lx)
  
  fhr<-f_beha_fun(s = hr_slopef,Behavoiur_score = partemp$hr,be_in = hr_inf)
  
  fchron<-f_chron_fun(s = chron_slopef,crhonotype_score = partemp$chrtype,be_in = chron_inf)
  
  
  #assgining survival probabilities
  #surv.prob.M<-runif(nrow(Ftemp), 0.9, 1)
  surv.prob.M<-1-prob_M_fun(Mmax, M = M, Mhr, Mchron) # como aplicar mortalidad a edad 0?
  #plot(surv.prob.M)
  
  surv.prob.M[Ftemp$A==0] <- runif(length (Ftemp$A[Ftemp$A==0]), 0.7, 0.99)##asignacion aleatoria de prob supevivencia a reclutas de edad 0
  #range(surv.prob.M)
  
  surv.prob.f<-1-prob_F_fun(fmax2, fl, fhr,fchron,Ftemp$A)
  #range(surv.prob.f)
  #surv.prob.f<-runif(nrow(Ftemp), 0.9,1)
  surv.prob.f[Ftemp$L <= lx] <- 1 #fishing mortality only applies from min fishing length
  #range(surv.prob.f)
  
  #apploying survival probabilities
  survM_status <- rbinom(n = nrow(partemp), size = 1, prob = surv.prob.M) #binomial varaible for each individual (0=die, 1=survive)
  survM_index <- which(survM_status == 1) #extract index of survivals
  partemp <- partemp[survM_index, ] #substracts death individuals
  Ftemp <- Ftemp[survM_index, ]
  
  
  survF_status <- rbinom(n = nrow(partemp), size = 1, prob = surv.prob.f) #binomial varaible for each individual (0=die, 1=survive)
  survF_index <- which(survF_status == 1) #extract index of survivals
  partemp <- partemp[survF_index, ] #substracts death individuals
  Ftemp <- Ftemp[survF_index, ]
  
  partemp <- partemp[which(Ftemp$A<=agemax),]
  Ftemp <- subset(Ftemp,A<=agemax) #ind older than age max die
  
  
  
  #calculating n of survivers for each age class
  #for (i in 1:length(age)){ # time steps = years from 0 to 8
  #prob1=
  #s2=1-prob2[i+1]
  #s1=1-prob2[i]
  #p=(s1-s2)/s1 #-- # p for not surviving at t+1 given you did survive at t
  #temp=sum(runif(n1[i],0,1)>p)
  #n1=c(n1,temp) # p for not surviving at t+1 given you did survive at t
  #}
  
  #age
  Ftemp$A=Ftemp$A+1
  #VB
  Ftemp$L=von_ber_fun(linf1,k1,AGE=Ftemp$A, age01)
  #weight - LWR function
  Ftemp$W=lwr_fun(a1,b1,Ftemp$L)
  #Biomas of the year
  #table(Ftemp$W)
  SpBio<-subset(Ftemp, Ftemp$A>1)$W
  Bioyear<-sum(SpBio)
  BioyearAll<-c(BioyearAll,Bioyear)
  #R of the next year
  rbiomass<-round(recruit_fun(ar,br,Bioyear)) #biomass of recruits
  R<-round(rbiomass/biomassr) #number of recruits
  
  
  #recxruits of the next year
  
  FtempR<-data.frame(L=rep(0.1,R),
                     W=rep(br,R), 
                     A=rep(0,R), #age= -1?)
                     ID=seq(IDRini,(IDRini+R-1)))
  
  #str(FtempR)
  #Recruits parameters
  #simulating behaviour scores
  hrR=runif(R, 0, 1)
  chrtypeR=runif(R, 0, 1)
  
  parR<-data.frame(k=rep(k1,R),
                   hr=hrR,
                   chrtype=chrtypeR, 
                   ID=seq(IDRini,IDRini+R-1))
  #join DBs
  Ftemp<-rbind(Ftemp,FtempR)
  
  partemp<-rbind(partemp, parR)
  
  #str(partemp)
  #str(parR)
  #save for the year
  assign(paste("F",years,sep=""),Ftemp)
  assign(paste("par",years,sep=""),partemp)
  
  
}

plot(seq(plusyears,plusyears+length(BioyearAll)-1),BioyearAll,type="l", xlab = "Year", ylab = "Biomass")

plot()