#Individual based model Raor #8
#integrating natural - fishing mortality models (by Miquel Palmer)
#adding simulated behaviour (home range HR and Chronotype) induced mortality

#Author: Pep Alós & Enrique Sánchez-Fabrés
#Colaborator: Miquel Palmer
#Project Metaraor
#Supervisor: Josep Alós

#Last edit: 02/06/2024

rm(list = ls())

#Population dynamics model: B(t+1) = B(t) + G + R + Z
#Where: B(t+1) and B(t) are the stock biomass at times t and t-1, respectively
#G is the is the biomass growth
#R is the biomass of new recruits
#Z is the total mortality or biomass loss, given by the sum of natural mortality (M) and fishing mortality (F): Z = M + F

####Formulas - functions ##########

####Growth####
#############

###Von Bertalanffy model for length growth###
#L_inf: asintotic length, K: growth rate (y^-1), age, age0=t0=null growth theoric age
von_ber_fun<-function(L_inf, K, AGE, AGE0){ 
  L_inf*(1 - exp(-K *(AGE-AGE0)))
}


#plot(x=seq(0, 7, 1), von_ber_fun(17.87, 1.079, seq(0, 7, 1), 0), ylab = "Length (cm)", xlab="Age (years)", main = "Von Bertalanffy growth model for pearly razorfish", sub = "k=1.079 | Linf=17.87")

#Length Weight relationship funcion (LWR)
lwr_fun<- function(a, b, L){ #a and b are LWD parameters and
  a*L^b                       #L is length
}

####Mortality######
###################

#The ammount of biomass lost will be determined by the individuals that dont survive each year. 

####Natural mortality (M)######
#In this approach, the total survival probability (Stotal) of each individual will be determined by a base survival probability (Sbase), result of the base natural mortality rate, and modified by the survival probability associated to individual's home range (HR) and chronotype (chron):

#Stotal = Sbase + SHR + Schron

#hence adding spatial mortality induced by behaviour 

#####Age dependant natural mortality (Ma): Lorenzen 2022 
#The Lorenzen allometric and length-inverse natural mortality model is designed to estimate the natural mortality rate as a function of age (a). This model assumes that mortality decreases with increasing age in a manner that can be described using the von Bertalanffy growth parameters:

Ma_fun<-function(M_inf, K, AGE, AGE0){
  
  M_inf*(1-exp(-K*(AGE-AGE0)^-1))
}

#M_inf is the natural mortality rate at L_inf. It can be estimated from the mortality (Mref) at a given reference age: recruitment age/length (AGEr)

Minf_fun<-function(M_ref, K, AGEr, AGE0){
  
  M_ref/(1-exp(-K*(AGEr-AGE0)^-1))
}

#Based on the Ergon et al. (2018), for each fish, the base probability of survive each year will be defined by natural mortality rate at previous age class (Ma1) and at next age class (Ma2) as:

Sbase_fun<-function(Ma1, Ma2){
  
  exp(-(Ma1+Ma2)/2)
}

#The survival probability induced by individual's HR and chron is defined a alogistic function as follows:

SHR_fun<-function(s_HR, score_HR, inf_HR){ 
  
  1/(1+exp(s_HR*(score_HR-inf_HR)))
}#Where s_HR and inf_HR are the slope and the inflexion point of the function, respectively


Schron_fun<-function(s_chron, score_chron, inf_chron){
  
  1/(1+exp(s_chron*(score_chron-inf_chron)))
}

#plot(1/(1+exp(0.1*(seq(0,1,0.1)-0.7))), x=seq(0,1,0.1), xlab="score" ,ylab="surv. prob")

#The hypothesis is that each individual can die independiently by any of the three sources. 

#SHR and Scrhon are modifiers (values ranging from 0 to 1) of Sbase (ranging between Minf for big individuals (~Linf) and virtualy 0 for individuals of age class 0)

####Fishing mortality (F) ####

#Fishing mortality is the rate at which fish are removed from a population due to fishing activities, including commercial, recreational, and subsistence fishing. Fishing efforts influence this type of mortality and is directly controllable through management practices such as gear restrictions, size limits, and quotas (Beverton and Holt, 1957). 

#For the present model, the fishing mortality rate F is influenced by the size of the fish. 

####Length dependent fihsing mortality (Flength)####
f_length_fun<-function(s, L, Lx){ #L=lengths (from growth model), s & Lx= slope and length of inflexion for fishing mortality function
  1-(1/(1 + exp (s * (L - Lx))))
}

#Similarly to natural mortality, the individual survival probability to fishing will be influenced by their spatial behaviour traits (HR and chronotype). Same functions than defined for natural mortality. 

#therefore, the total survival probability with respect to fishing will also be obtained by 
#Stotal_f = S_Flength * SHR * Schron

#where S_Flength is, similarly to the age-dependant natural mortality,= exp-((Fl1+Fl2)/2). (same function created above can be applied)


####Recruitment (R) model: Beverholt and Holt

#Two possible approach: based on n of individuals, or based on spawning biomass (SpBio):

recruit_fun<- function(alfa, beta, SpBio){ #alpha is the maximum recruitment rate achievable (denisty independent paramters) 
  alfa*SpBio/(1+(beta*SpBio)) #beta  reflects the strength of density-dependent effects
}

#based on n of individuals:
#recruit_fun<- function(alfa, beta, N){ #alpha is the maximum recruitment rate achievable (denisty independent paramters) 
 # alfa*N/1+(beta*N) #beta  reflects the strength of density-dependent effects
#}


#####Population parameters######

agemax=7
age= seq(1,agemax, 1) #age classes
#LWR paramters: males a=0.002; b=3.383 | Females a=0.003; b= 3.32 (source: fishbase)
a1=0.0025
b1=3.35
ind=1000
####Growth######
#Von Bertalanffy parameters (Mekni et al., 2018)
linf1= 17.18 #cm assympotic length
k1= 1.079 #growth rate
age01= 0#null growth theoric age

# Ajustar los márgenes y luego graficar
par(mgp = c(2, 0.5, 0))  # Reducir espacio entre las etiquetas de los ejes y el gráfico

lt<-von_ber_fun(linf1, k1, age, age01) #length at age 
plot(lt, ylab = "Length (cm)", xlab = "Age class (years)", main = "Growth function of pearly razorfish (Von Bertalanffy)", 
     #sub = paste("k=", k1, " | ", "L_inf=", linf1, " | ", "source: Mekni et al. (2018)")
     ); lines(lt, col="red")

#LWR
wt=lwr_fun(a1, b1, lt)

#####Natural mortality#####
a_rec = 1 #reference age of recruitment
m_ref = 0.5 #ref mortality rate at a_ref

#Minf: natural mortality rate at age of reference (recruitment size)
m_inf<- Minf_fun(M_ref = m_ref, K = k1, AGEr = a_rec, AGE0 = age01)
#exp(-m_inf)
#Natural mortality rate at age given by age dependent function from Lorenzen 2022
Ma<-Ma_fun(M_inf = m_inf, K = k1, AGE = age, AGE0 = age01)

plot(Ma, xlab = "Age class (year)", main = "Age-dependant natural mortality rate (Ma)", ylab = "Ma"); lines(Ma, col="blue")


#f(length) parameters
fslope1=0.1 #slope logistic function of length dependant fishing mortality
#lt1=#obtained from Von Bert. model
lx=15 #inflection length for fishing mortality

#length dependant fishing mortality rate
F_length<-f_length_fun(s=fslope1, L=lt, Lx=lx)

plot(F_length, main = "Length-dependant F",  ylab="F", xlab = "Age class (years)", sub = paste("s=", fslope1, " | ", "Linflection=", lx, "(cm)") ); lines(F_length)

#check function shape
#plot(f_length_fun(fslope1, seq(0, 100, 1), 40), type = "l", col="darkblue", lwd=2, ylab = "F(l)", xlab = "Size")



####survival prob paramters####
#s(behaviour) paramters
hr_inf=0.02 #behaviour score inflection for fishing mortality
hr_inm=0.02  #behaviour score inflection for natural mortality
hr_slopef=0.1 #slope for behaviour induced fishing morality
hr_slopem=0.1 #slope for behaviour induced natural morality

#s(chron) paramters
chron_inf=0.02 #chron score inflection for fishing mortality
chron_inm=0.02  #chron score inflection for natural mortality
chron_slopef=0.1 #slope for chron induced fishing morality
chron_slopem=0.1 #slope for chron induced natural morality

#simulating behaviour scores

# function to normalize data
normalize <- function(x) {
  return ((x - min(x)) / (max(x) - min(x)))
}


####home range##################
#mean KUD95 (m2) from Alos et al 2016 ####
kuds<-c(0.2, 0.22, 0.16, 0.22, 0.17, 0.26,0.17,0.16,0.39, 0.17, 0.32,0.16, 0.18, 0.17,0.16,0.17,0.16,0.17,0.19, 0.18)

#sd KUD96 from Alos et al 2016
sdkuds<-c(0.02,0.03,0,0.03,0.01,0.08,0.01,0.01,0.11,0.02,0.04,0.01,0.01,0.02,0.02,0.01,0.02,0,0.01,0.03,0.01,0.01,0)


plot(density(kuds), xlab="Home range KUD (m2)", main="Home range of wild pearly razorfish")

hr<-rnorm(ind, mean = mean(kuds), sd = mean(sdkuds))
#?rnorm
plot(density(hr), main="Home range of simulated pearly razorfish")

hr<-normalize(hr)

plot(density(hr), main="HR scores of simulated pearly razorfish")

plot(hr)
summary(hr)
boxplot(hr)

#######Chronotype#########


# Alós et al 2017: The mean and s.d. of awakening time (min relative to sunrise) and rest onset (min relative to sunset) were 117.9±85.1 min and 5±8.1 min, respectively
mean_awMreal=117.9
sd_awMreal=85.1
#Awakening time (min relative to sunrise)
awMreal<-c(100.2,81.4,211.1,114.9,29.3,18.2,118.3,130.4,74.7,211.8,271.0,192.0,144.6,194.1)#data from Alos et al. 2017
awSDreal<-c(49.4,59.1,105.2,140.0,13.3,18.1,33.1,32.4,27.0,16.2,39.9,56.0,33.2,14.3)#data from Alos et al. 2017

plot(density(awMreal), xlab="min relative to sunrise", main="Awakening time of wild pearly razorfish")
hist(awMreal, xlab="min relative to sunrise", main="", breaks=15)

awMsim<-rnorm(100, mean = mean_awMreal, sd=sd_awMreal)

plot(density(awMsim), main="Sim. awaking time")
#the lowest values of awMsim correspond to the individuals awaking earlier
#this individuals are expected to be more vulnerable to fishing mortality
#to express the score of the chronotype the awaking values are normalized to range from 0 to 1
#the resulting values are rested to one to express that a higher score will correspond to higher risk chronotype
chrtype<-1-normalize(awMsim) 

plot(density(chrtype), main=NA, xlab="Sim. chronotype scores")
hist(chrtype, breaks = 50)

###Recruitment model: Beverton and Holt####

ar <- 20 #recruitmen density independent parameter proportional to fecundity
#The units of ar are "recruitment per spawner" 
br <- 0.05 #b is a density-dependent parameter that is proportional to both fecundity and density-dependent mortality
#The Beverton-Holt model is based on the assumptions that juvenile competition results in a mortality rate that is linearly dependent upon the number of fish alive in the cohort at any time and that predators are always present. The Beverton-Holt model is appropriate "if there is a maximum abundance imposed by food availability or space, or if the predator can adjust its predatory activity immediately to changes in prey abundance" (Wootton 1990)

#biomass per recruit
biomassr<- wt[1] ## biomass of individuals age=1
plot(recruit_fun(ar, br, seq(0,1000, 1))/biomassr, ylab = "Nº of recruits", xlab = "Spawning biomass (g)", main = "Beverton and Holt recruitment function", sub = paste("alpha=", ar, " | ", "beta=", br)) #Beverton-Holt equation
#plot(SpBio, R)
# Graficar sin los valores de los ejes
# Ajustar los márgenes y luego graficar
par(mgp = c(1, 0.0, 0))  # Reducir espacio entre las etiquetas de los ejes y el gráfico


plot(recruit_fun(ar, br, seq(0, 1000, 1)) / biomassr,
     ylab = "Nº of recruits", 
     xlab = "Spawning biomass (g)",
     main = "Beverton and Holt recruitment function",
     xaxt = "n",  # Desactivar los valores del eje x
     yaxt = "n")  # Desactivar los valores del eje y


####### code pep IBM

#numbe rof individuals
#ind=100
#years to simualte
IniYear=2024
Years=100
#ind adultos

#the starting population is formed by recruits. Based the beverton and holt model, the recruits are considered as those individuals entering the fishery. In other words, individuals bigger than the minimum length of fishing (lx = 10cm). 
#Given the Von Bertalanffy model, the individuals of age class 2 (1 year) are 11cm long:
F<-data.frame(L=rep(lt[1],ind),W=rep(wt[1],ind),A=rep(age[1],ind),ID=round(seq(1:ind)))
assign(paste("F",IniYear,sep=""),F)
#str(F2024)
#str(F)
#range(F$ID)


#ind adultos parametros
par<-data.frame(k=rep(k1,ind),hr=hr,chrtype=chrtype,ID=round(seq(1:ind)))
assign(paste("par",IniYear,sep=""),par)
#str(par2024)


#rnorm(100, mean = k1, sd=0.2) if individual k were to be assignated...

#Bioyear
BioyearAll<-sum(F$W)
rm(years)


for (years in (IniYear+1):(IniYear+Years)){
  #years=2025
  
  Ftemp<-get(paste("F",years-1,sep=""))
  partemp<-get(paste("par",years-1,sep=""))
  
  #remove vectors from previous year
  #rm(list=paste("F",years-1,sep=""))
  #rm(list=paste("par",years-1,sep=""))
  rm(Ftemp2)
  rm(partemp2)
  rm(FtempR)
  rm(parR)
  #for new ids
  #range(Ftemp$ID)
  IDRini<-max(Ftemp$ID)+1
  
  #Module for obtaining survivals 
  #survival probability to natural and fishing mortality are calculated and appplied separately
  Ftemp2=Ftemp
  partemp2=partemp
  rm(i)
  for (i in 1:nrow(Ftemp)) {
    #i=1
    ###natural mortality###
    itemp<-Ftemp[i,]
    id<-itemp$ID
    a<-itemp$A
    #k<-partemp[i]$K
    ma1<-Ma_fun(m_inf, k1, a, age01)
    ma2<-Ma_fun(m_inf, k1, (a+1), age01)
    
    #Base natural survival probability
    Sbase<-Sbase_fun(Ma1 = ma1, Ma2 = ma2)
    
    ### Chron and HR survival probabilities
    itemp2<-partemp[i,]
    chron<-itemp2$chrtype 
    hr1<-itemp2$hr
    
    #applying Schorn and SHR functions
    Schron<-Schron_fun(chron_slopem, chron, chron_inm)
    SHR<-SHR_fun(hr_slopem, hr1, hr_inm)
    
    #Sbehaviour=(Schron*SHR)*2
    
    Stotal=Sbase#*Schron*SHR #always results in very low surv prob
    
    #Stotal=Sbase*Sbehaviour
    
    #condition to determine if fish survives to natural mortality
    if (runif(1,0,1)>=Stotal){#runif generates a random value between 0 and 1 
      Ftemp2<-Ftemp2[!Ftemp2$ID==id,]#if it is > than survival prob,
      partemp2<-partemp2[!partemp2$ID==id,]# the individual dies and is removed from population
     } else {
     ###Fishing mortality###
      l<-itemp$L
      l2<-von_ber_fun(linf1, k1, (a+1), age01)
      f1<-f_length_fun(fslope1, l, lx)#similarly as for natural mortality
      f2<-f_length_fun(fslope1, l2, lx)
    
      Sbase_f<-Sbase_fun(f1, f2)
    
      #applying Schron and SHR functions (for fishing)
      Schron_f<-Schron_fun(chron_slopef, chron, chron_inf)
      SHR_f<-SHR_fun(hr_slopef, hr1, hr_inf)
    
      Stotal_f=Sbase_f*Schron_f*SHR_f
    
      #condition if it is > than survival prob,
      if (runif(1,0,1)>=Stotal_f){
        Ftemp2<-Ftemp2[!Ftemp2$ID==id,]
        partemp2<-partemp2[!partemp2$ID==id,]
       #the individual dies and is removed from population
      }
    }
  }
  
  Ftemp=Ftemp2
  partemp=partemp2
  #partemp <- partemp[which(Ftemp$A<agemax),]
  #Ftemp <- subset(Ftemp,A<agemax) #ind older than age max die
  
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
  
  #Recruits of the next year
  rbiomass<-recruit_fun(ar,br,Bioyear) #biomass of recruits
  R<-floor(rbiomass/biomassr) #number of recruits
  
  #recruits of the next year
  if(R>0){  
  FtempR<-data.frame(L=rep(lt[1],R),
                     W=rep(biomassr,R), 
                     A=rep(a_rec,R),
                     ID=seq(IDRini,(IDRini+R-1)))
  
  #str(FtempR)
  #Recruits parameters
  #simulating recruits behaviour scores, sampled from spawners
  hrR=sample(partemp$hr, size = R, replace = T)
  chrtypeR=sample(partemp$chrtype, size = R, replace = T)
 
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

length(BioyearAll)
par(mfrow=c(1,1))
# Reducir espacio entre las etiquetas de los ejes y el gráfico
par(mgp = c(2, 1, 0))


plot(seq(IniYear,IniYear+length(BioyearAll)-1),BioyearAll,type="l", xlab = "Year", ylab = "Biomass (g)", sub = paste("Mref=", m_ref, " | ar=", ar, "; br=",br, " | Lx=",lx ), main = "Ma+Fl()" )


plot(seq(IniYear+1,IniYear+length(BioyearAll)-1),BioyearAll[2:length(BioyearAll)],type="l", xlab = "Year", ylab = "Biomass", sub = paste("Mref=", m_ref, " | ar=", ar, "; br=",br, " | Lx=",lx ), main = "Ma+Fl()")




table(F2034$A)

table(F2124$A)

#####further projections ########
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