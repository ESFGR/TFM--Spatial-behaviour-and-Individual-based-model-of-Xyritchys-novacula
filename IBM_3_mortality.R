#Individual based model Raor #3
#integrating natural - fishing mortality models (by Miquel Palmer)

#Author: Enrique Sánchez-Fabrés
#Project Metaraor
#Supervisor: Josep Alós

#Last edit: 16/04/2024

#Formulas - functions ----

##Von Bertalanffy model for length growth
von_ber_fun<-function(L_inf, K, AGE, AGE0){ 
  L_inf*(1 - exp(-K *(AGE-AGE0)))
}

##Survival probability with length dependent fishing mortality: f(length)

####f(length)
f_length_fun<-function(s, L, Lx){ #L=lengths (from growth model), s & Lx= slope and length of inflexion for fishing mortality function
  1-(1/(1 + exp (s * (L - Lx))))
}

###Survival probability
prob_fun<-function(F_max, F_length, M, AGE){ ## F_max=maximum fishing mortality rate; M=natural mortality
(1 - exp(-(F_max*F_length+M)*AGE)) 
}

#Population simulations----

agemax=7
age= seq(0,agemax, 1)

##scenario 1 #parameters from Mekni et al. (2018) -----

#LWR paramters: males a=0.002; b=3.383 | Females a=0.003; b= 3.32
a1=0.0025
b1=3.35

#Von Bertalanffy parameters

linf1= 17.18#cm assympotic length
k1= 1.079 #growth rate
age01= 0#null growth theoric age

lt1=von_ber_fun(linf1, k1, age, age01)

#growth graph
plot(lt1, xlab= "Age class (year)", ylab = "Length (cm)", main = "Von Bertalanffy growth" , sub = paste("linf=17.18", "k=1.079", "t0=0", "|"," Mekni et al. (2018)")) 
lines(lt1, col="red")

#f(length) parameters
fslope1=0.5 #slope and inflection of a logistic function (or any other function)
#lt1=#obtained from Von Bert. model
lx=10 #length at which fishing mortality starts


lexample=seq(0, 20, 1)

f_length=f_length_fun(fslope1, lexample, lx)
plot(f_length, ylab = "Fishing mortality rate", xlab = "Length (cm)",  main = "Length-dependent fishing mortality ", sub=paste("s=",fslope1, "|", "Lx=", lx, " | Flength=1-(1/(1 + e^(s*(L-Lx))"))


#surv. prob. paramters #in mekni et al. (2018): m=0.949 and f =0.35...
m1=0.1#natural mortality
fmax1=0.35 # maximum fishing mortality rate (f(length) will move between 0 (small length) and f_max)

#Build matrix--

n1=100 #initial number of fish

M<-matrix(NA, n1, length(age)) #8 age clases: 0, 1, 2, 3, 4, 5, 6, 7

remove(i)

for (i in 0:length(age)) { 
  nextgen<-rep(lt1[i], n1) #vector of n individual lenghts at age class "i" before appliying mortality
  M[, i]<-nextgen #adding length to class age i at matrix 
  
}

M[,1]=0.1 #population size matrix, initial size=0.1cm, at first age class (0)
M#matrix of population sizes per age class before applying mortality
  

#applying mortality
# DISCRETIZATION model#3 Length-dependent fishing mortality: f(length)

f_length1=f_length_fun(fslope1, lt1, lx) #length-dependent fishing mortality rate at age class "i"
#plot(f_length1)

prob1=prob_fun(fmax1, f_length1,m1, age)#surv. probability

plot(prob1, xlab = "Age class", ylab = "Surv. Probability", main="p(surv)1 - e^(-(F_max*F_length+m)*t)", sub = paste("fmax=",fmax1, "", "m=", m1))


for (i in 1:length(age)){ # time steps = years from 0 to 8
  s2=1-prob1[i+1]
  s1=1-prob1[i]
  p=(s1-s2)/s1 #-- # p for not surviving at t+1 given you did survive at t
  temp=sum(runif(n1[i],0,1)>p)
  n1=c(n1,temp) # p for not surviving at t+1 given you did survive at t
}

plot(n1/100, sub = "Survival curve", ylab = "prop. of survivals", xlab = "Age class")

remove(i)

for (i in 2:length(age)) {
  
  M[(n1[i]+1):n1[1], i] <- NA
  
}

M #final matrix of population  

#Biomass
B1<-rep(NA, 8)

for (i in 1:8) {
  
  L=sum(na.omit(M[,i])) #total length of age class i
  B1[i]<-a1*L^b1 #LWR formula
  
}

data.frame(B1) #biomass (g) per class age

plot(B1/1000, ylab = "Biomass (kg)", xlab = "Age class")
lines(B1/1000)  

#Scenario 2. Parameters from Bataglia et al., 2010 -----

#LWR paramters

a2=0.0139
b2=2.9326

#Von Bertalanffy parameters
linf2= 17.50#cm assympotic length
k2= 0.8 #growth rate
age02= -0.75#null growth theoric age

lt2=von_ber_fun(linf2, k2, age, age02)

#growth graph
plot(lt2, xlab= "Age class (year)", ylab = "Length (cm)", main = "Von Bertalanffy growth" , sub = paste("linf=",linf2, "k=",k2,"t0=",age02," Bataglia et al. (2010)"))
lines(lt2, col="red")

#f(length) parameters
fslope2=0.5 #slope and inflection of a logistic function (or any other function)
lx=10 #length at which fishing mortality starts


lexample=seq(0, 20, 1)

f_length=f_length_fun(fslope2, lexample, lx)
plot(f_length, ylab = "Fishing mortality rate", xlab = "Length (cm)",  main = "Length-dependent fishing mortality | Flength=1-(1/(1 + e^(s*(L-Lx))", sub=paste("s=",fslope2, "|", "Lx=", lx))


#
m2=0.1#natural mortality
fmax2=0.35 # maximum fishing mortality rate (f(length) will move between 0 (small length) and f_max)

#Build matrix--

n2=100 #initial number of fish

M2<-matrix(NA, n2, length(age)) #8 age clases: 0, 1, 2, 3, 4, 5, 6, 7

remove(i)

for (i in 0:length(age)) { 
  nextgen<-rep(lt2[i], n2) #vector of n individual lenghts at age class "i" before appliying mortality
  M2[, i]<-nextgen #adding length to class age i at matrix 
  
}

M2[,1]=0.1 #population size matrix, initial size=0.1cm, at first age class (0)
M2#matrix of population sizes per age class before applying mortality


#applying mortality
# DISCRETIZATION model#3 Length-dependent fishing mortality: f(length)

f_length2=f_length_fun(fslope2, lt2, lx) #length-dependent fishing mortality rate at age class "i"
#plot(f_length1)

prob2=prob_fun(fmax2, f_length2,m2, age)#surv. probability
plot(prob2, xlab = "Age class", ylab = "Surv. Probability", sub = paste("fmax=",fmax2, "", "m=", m2), main="p(surv)1 - e^(-(F_max*F_length+m)*t)",)

for (i in 1:length(age)){ # time steps = years from 0 to 8
  s2=1-prob2[i+1]
  s1=1-prob2[i]
  p=(s1-s2)/s1 #-- # p for not surviving at t+1 given you did survive at t
  temp=sum(runif(n[i],0,1)>p)
  n2=c(n2,temp) # p for not surviving at t+1 given you did survive at t
}

plot(n2/100, sub = paste("fmax=",fmax2, "", "m=", m2), ylab = "prop. of survivals", xlab = "Age class")

remove(i)

for (i in 2:length(age)) {
  
  M2[(n2[i]+1):n2[1], i] <- NA
  
}

M2 #final matrix of population  

#Biomass
B2<-rep(NA, 8)

for (i in 1:8) {
  
  L=sum(na.omit(M2[,i])) #total length of age class i
  B2[i]<-a2*L^b2 #LWR formula
  
}

data.frame(B2) #biomass (g) per class age

plot(B2/1000, ylab = "Biomass (kg)", xlab = "Age class", main=paste("LWR paramters a=", a2, "", "b=", b2))
lines(B2/1000)

#Scenario 3.  -----

#LWR paramters from Dimitriadis and Fournari-Konstantinidou 2018 

a3=0.003
b3=3.075

#Von Bertalanffy parameters
linf3= 17.50#cm assympotic length
k3= 0.8 #growth rate from fishbase
age03= -0.75#null growth theoric age

lt3=von_ber_fun(linf3, k3, age, age03)

#growth graph
plot(lt3, xlab= "Age class (year)", ylab = "Length (cm)", main = "Von Bertalanffy growth" , sub = paste("linf=",linf3, "k=",k3,"t0=",age03," Dimitriadis and Fournari-Konstantinidou (2018) "))
lines(lt3, col="red")

#f(length) parameters
fslope3=0.5 #slope and inflection of a logistic function (or any other function)
lx=10 #length at which fishing mortality starts


lexample=seq(0, 20, 1)

f_length=f_length_fun(fslope3, lexample, lx)
plot(f_length, ylab = "Fishing mortality rate", xlab = "Length (cm)",  main = "Length-dependent fishing mortality | Flength=1-(1/(1 + e^(s*(L-Lx))", sub=paste("s=",fslope1, "|", "Lx=", lx))


#
m3=0.1#natural mortality
fmax3=0.35 # maximum fishing mortality rate (f(length) will move between 0 (small length) and f_max)

#Build matrix--

n3=100 #initial number of fish

M3<-matrix(NA, n3, length(age)) #8 age clases: 0, 1, 2, 3, 4, 5, 6, 7

remove(i)

for (i in 0:length(age)) { 
  nextgen<-rep(lt3[i], n3) #vector of n individual lenghts at age class "i" before appliying mortality
  M3[, i]<-nextgen #adding length to class age i at matrix 
  
}

M3[,1]=0.1 #population size matrix, initial size=0.1cm, at first age class (0)
M3#matrix of population sizes per age class before applying mortality


#applying mortality
# DISCRETIZATION model#3 Length-dependent fishing mortality: f(length)

f_length3=f_length_fun(fslope3, lt3, lx) #length-dependent fishing mortality rate at age class "i"
#plot(f_length1)

prob3=prob_fun(fmax3, f_length3,m3, age)#surv. probability
plot(prob3, xlab = "Age class", ylab = "Surv. Probability", sub = paste("fmax=",fmax3, "", "m=", m3), main="p(surv)1 - e^(-(F_max*F_length+m)*t)",)


for (i in 1:length(age)){ # time steps = years from 0 to 8
  s2=1-prob3[i+1]
  s1=1-prob3[i]
  p=(s1-s2)/s1 #-- # p for not surviving at t+1 given you did survive at t
  temp=sum(runif(n[i],0,1)>p)
  n3=c(n3,temp) # p for not surviving at t+1 given you did survive at t
}

plot(n3/100, sub = paste("fmax=",fmax3, "", "m=", m3), ylab = "prop. of survivals")
lines(n3/100)
remove(i)

for (i in 2:length(age)) {
  
  M3[(n3[i]+1):n3[1], i] <- NA
  
}

M3 #final matrix of population  

#Biomass
B3<-rep(NA, 8)

for (i in 1:8) {
  
  L=sum(na.omit(M3[,i])) #total length of age class i
  B3[i]<-a*L^b #LWR formula
  
}

data.frame(B3) #biomass (g) per class age

plot(B3/1000, ylab = "Biomass (kg)", xlab = "Age class", main=paste("LWR paramters a=", a3, "", "b=", b3))
lines(B3/1000)

#Scenario 4.  Cardinale, M., Colloca, F., & Ardizzone, G. D. (1998). ----

#LWR paramters from Fishbase 

a4=0.003
b4=3.075

#Von Bertalanffy parameters -- Cardinale, M., Colloca, F., & Ardizzone, G. D. (1998).
		
linf4= 23.441#cm assympotic length
k4= 0.289 #growth rate
age04= -0.5915#null growth theoric age

lt4=von_ber_fun(linf4, k4, age, age04)

#growth graph
plot(lt4, xlab= "Age class (year)", ylab = "Length (cm)", main = "Von Bertalanffy growth" , sub = paste("linf=",linf4, "k=",k4,"t0=",age04," Cardinale, Colloca, & Ardizzone (1998) "))
lines(lt4, col="red")

#f(length) parameters
fslope4=0.5 #slope and inflection of a logistic function (or any other function)
lx=10 #length at which fishing mortality starts


lexample=seq(0, 20, 1)

f_length=f_length_fun(fslope4, lexample, lx)
plot(f_length, ylab = "Fishing mortality rate", xlab = "Length (cm)",  main = "Length-dependent fishing mortality | Flength=1-(1/(1 + e^(s*(L-Lx))", sub=paste("s=",fslope1, "|", "Lx=", lx))


#
m4=0.1#natural mortality
fmax4=0.35 # maximum fishing mortality rate (f(length) will move between 0 (small length) and f_max)

#Build matrix--

n4=100 #initial number of fish

M4<-matrix(NA, n4, length(age)) #8 age clases: 0, 1, 2, 3, 4, 5, 6, 7

remove(i)

for (i in 0:length(age)) { 
  nextgen<-rep(lt4[i], n4) #vector of n individual lenghts at age class "i" before appliying mortality
  M4[, i]<-nextgen #adding length to class age i at matrix 
  
}

M4[,1]=0.1 #population size matrix, initial size=0.1cm, at first age class (0)
M4#matrix of population sizes per age class before applying mortality


#applying mortality
# DISCRETIZATION model#3 Length-dependent fishing mortality: f(length)

f_length4=f_length_fun(fslope4, lt4, lx) #length-dependent fishing mortality rate at age class "i"
#plot(f_length1)

prob4=prob_fun(fmax4, f_length4,m4, age)#surv. probability
plot(prob4, xlab = "Age class", ylab = "Surv. Probability", sub = paste("fmax=",fmax4, "", "m=", m4), main="p(surv)1 - e^(-(F_max*F_length+m)*t)",)


for (i in 1:length(age)){ # time steps = years from 0 to 8
  s2=1-prob4[i+1]
  s1=1-prob4[i]
  p=(s1-s2)/s1 #-- # p for not surviving at t+1 given you did survive at t
  temp=sum(runif(n4[i],0,1)>p)
  n4=c(n4,temp) # p for not surviving at t+1 given you did survive at t
}

plot(n4/100, sub = paste("fmax=",fmax4, "", "m=", m4), ylab = "prop. of survivals")
lines(n4/100)

remove(i)

for (i in 2:length(age)) {
  
  M4[(n4[i]+1):n4[1], i] <- NA
  
}

M4 #final matrix of population  

#Biomass
B4<-rep(NA, 8)

for (i in 1:8) {
  
  L=sum(na.omit(M4[,i])) #total length of age class i
  B4[i]<-a*L^b #LWR formula
  
}

data.frame(B4) #biomass (g) per class age

plot(B4/1000, ylab = "Biomass (kg)", xlab = "Age class", main=paste("LWR paramters a=", a4, "", "b=", b4))
line(B4/1000)


#escenario 1 plots ----

par(mfrow=c(2,2))

plot(lt1, xlab= "Age class (year)", ylab = "Length (cm)", main = "Von Bertalanffy growth (1)" , sub = paste("linf=",linf1, "k=", k1, "t0=", age01,"|"," Mekni et al. (2018)")) 
lines(lt1, col="red")

plot(f_length, ylab = "Fishing mortality rate", xlab = "Length (cm)",  main = "Length-dependent fishing mortality ", sub=paste("s=",fslope1, "|", "Lx=", lx, " | Flength=1-(1/(1 + e^(s*(L-Lx))"))

plot(prob1, xlab = "Age class", ylab = "Surv. Probability (1)", main="p(surv)=1-e^(-(F_max*F_length+m)*t)", sub = paste("fmax=",fmax1, "", "m=", m1))
lines(prob1)

plot(n1/100, main = "Survival curve (1)", ylab = "prop. of survivals", xlab = "Age class")
lines(n1/100)

par(mfrow=c(1,1))
plot(B1/1000, ylab = "Biomass (kg)", xlab = "Age class", sub  =paste("LWR paramters a=", a1, "", "b=", b1), main="Total biomass (1)")
lines(B1/1000)

#escenario 2 plots

par(mfrow=c(2,2))
#growth graph
plot(lt2, xlab= "Age class (year)", ylab = "Length (cm)", main = "Von Bertalanffy growth (2)" , sub = paste("linf=",linf2, "k=",k2,"t0=",age02," |Bataglia et al. (2010)"))
lines(lt2, col="red")

plot(f_length, ylab = "Fishing mortality rate", xlab = "Length (cm)",  main = "Length-dependent fishing mortality ", sub=paste("s=",fslope2, "|", "Lx=", lx, " | Flength=1-(1/(1 + e^(s*(L-Lx))"))

plot(prob2, xlab = "Age class", ylab = "Surv. Probability (2)", main="p(surv)=1-e^(-(F_max*F_length+m)*t)", sub = paste("fmax=",fmax2, "", "m=", m2))
lines(prob2)

plot(n2/100, main = "Survival curve (2)", ylab = "prop. of survivals", xlab = "Age class")
lines(n2/100)

par(mfrow=c(1,1))
plot(B2/1000, ylab = "Biomass (kg)", xlab = "Age class", main="Total biomass (2)",sub=paste("LWR paramters a=", a2, "", "b=", b2))
lines(B2/1000)

#Escenario 3 plots

par(mfrow=c(2,2))
#growth graph
plot(lt3, xlab= "Age class (year)", ylab = "Length (cm)", main = "Von Bertalanffy growth (3)" , sub = paste("linf=",linf3, "k=",k3,"t0=",age03,"|Dimitriadis 2018"))
lines(lt3, col="red")


plot(f_length, ylab = "Fishing mortality rate", xlab = "Length (cm)",  main = "Length-dependent fishing mortality ", sub=paste("s=",fslope3, "|", "Lx=", lx, " | Flength=1-(1/(1 + e^(s*(L-Lx))"))

plot(prob3, xlab = "Age class", ylab = "Surv. Probability (3)", main="p(surv)=1-e^(-(F_max*F_length+m)*t)", sub = paste("fmax=",fmax3, "", "m=", m3))
lines(prob3)

plot(n3/100, main = "Survival curve (3)", ylab = "prop. of survivals", xlab = "Age class")
lines(n3/100)

par(mfrow=c(1,1))

plot(B3/1000, ylab = "Biomass (kg)", xlab = "Age class", main="Total biomass (3)",sub=paste("LWR paramters a=", a3, "", "b=", b3))
lines(B3/1000)

#Escenario 4 - plots

par(mfrow=c(2,2))
#growth graph
plot(lt4, xlab= "Age class (year)", ylab = "Length (cm)", main = "Von Bertalanffy growth (4)" , sub = paste("linf=",linf4, "k=",k4,"t0=",age04,"|Cardinale (1998) "))
lines(lt4, col="red")


plot(f_length, ylab = "Fishing mortality rate", xlab = "Length (cm)",  main = "Length-dependent fishing mortality ", sub=paste("s=",fslope4, "|", "Lx=", lx, " | Flength=1-(1/(1 + e^(s*(L-Lx))"))

plot(prob4, xlab = "Age class", ylab = "Surv. Probability (4)", main="p(surv)=1-e^(-(F_max*F_length+m)*t)", sub = paste("fmax=",fmax4, "", "m=", m4))
lines(prob4)

plot(n4/100, main = "Survival curve (4)", ylab = "prop. of survivals", xlab = "Age class")
lines(n4/100)

par(mfrow=c(1,1))

plot(B4/1000, ylab = "Biomass (kg)", xlab = "Age class", main="Total biomass (4)",sub=paste("LWR paramters a=", a4, "", "b=", b4))
lines(B4/1000)

