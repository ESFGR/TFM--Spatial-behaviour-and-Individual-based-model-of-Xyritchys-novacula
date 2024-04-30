#Individual based model Raor #5: no fihing mortality projections

#Author: Enrique Sánchez-Fabrés
#Project Metaraor
#Supervisor: Josep Alós

#Last edit: 29/04/2024

#Formulas - functions ----

##Von Bertalanffy model for length growth: L_inf: asintotic length, K: growth rate (y^-1), age, age0=t0=null growth theoric age
von_ber_fun<-function(L_inf, K, AGE, AGE0){ 
  L_inf*(1 - exp(-K *(AGE-AGE0)))
}

##Survival probability with length dependent fishing mortality: f(length)

####f(length)
f_length_fun<-function(s, L, Lx){ #L=lengths (from growth model), s & Lx= slope and length of inflexion for fishing mortality function
  1-(1/(1 + exp (s * (L - Lx))))
}

### Behaviour induced mortality

####f(length)
f_beha_fun<-function(s, Behavoiur, be_in){ #L=lengths (from growth model), s & be_in= slope and behaviour of inflexion for behaviour induiced fishing mortality function
  1-(1/(1 + exp (s * (Behavoiur - be_in))))
}

###Survival probability
prob_fun<-function(F_max, F_length, F_beha, M, AGE){ ## F_max=maximum fishing mortality rate; M=natural mortality
  (1 - exp(-(F_max*(F_length+F_beha)+M)*AGE)) 
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

#f(behaviour) paramters
be_in1=60 #behaviour score at which fishing mortality starts
be_slope1=0.07

#simulating behaviour scores

be_scores1=runif(100, 0, 100)
plot(be_scores1)
summary(be_scores1)
boxplot(be_scores1)

#plot function

be_scores_example=seq(0,100,1)

f_beha=f_beha_fun(be_slope1, be_scores_example, be_in1 )
plot(f_beha)


#surv. prob. paramters #in mekni et al. (2018): m=0.949 and f =0.35...
m1=0.1#natural mortality
fmax1=0.1 # maximum fishing mortality rate (f(length) will move between 0 (small length) and f_max)

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


f_beha1=f_beha_fun(fslope1, be_scores1, be_in1 )
plot(f_beha1)

prob1=prob_fun(fmax1, f_length1,f_beha1, m1, age)#surv. probability

#plot(prob1, xlab = "Age class", ylab = "Surv. Probability", main="p(surv)1 - e^(-(F_max*F_length+m)*t)", sub = paste("fmax=",fmax1, "", "m=", m1))


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

Mbiomass= a1*M^b1 #matriz de masas (g): #LWR formula


B1<-rep(NA, 8)

for (i in 1:8) {
  
  B1[i]=sum(na.omit(Mbiomass[,i]))
  
}


data.frame(B1) #biomass (g) per class age

plot(B1, ylab = "Biomass (g)", xlab = "Age class")
lines(B1)  


#Recruitment model
SpBio <- B1 #age class stock biomass
ar <- 10 #recruitmen parameters
br <- 0.9
R <- ar*SpBio/((1+(br*SpBio))) #Beverton-Holt equation
plot(SpBio, R)



#Proyecting population-----

#scearios with no fishing mortality
fmax_0 = 0
f_length_0 = 0
f_beha_0 = 0

prob2=prob_fun(F_max = fmax_0, F_length = f_length_0, F_beha = f_beha_0, AGE = age, M=m1)

# yearly total stock biomass
biomasa_total <- numeric(11) #change according to the years to proyect in the loop
biomasa_total[1]<- sum(B1)

# yearly total stock biomass per age class
age_class_biomass<- matrix(NA, nrow = length(age), ncol = 11)
age_class_biomass[, 1]<- B1

#number of individuals per year
n_population<-numeric(11)
n_population[1]=sum(!is.na(M))


M2=Mbiomass

remove(i)
# Simulación de la evolución de la población durante 10 años
for (year in 2:11) {
  
  M2 <- cbind(NA, M2[, 1:7]) #new population matrix
  
  # Cálculo del reclutamiento utilizando la ecuación de Beverton-Holt
  #recruitmen parameters, with varaibility among years 
  ar <- rnorm(1, mean = 15, sd = 1) #(every iteration gets a different value)
  br <- rnorm(10, mean = 0.8, sd = 0.2)
  
  #Beverton-Holt equation: number or recruits
  R <- ar*SpBio/((1+(br*SpBio)))
  
  # adjust matrix dimensions to add recruits
  num_reclutas <- round(sum(R))
  num_filas <- nrow(M2)
  if (num_reclutas > num_filas) {
    M2 <- rbind(Mbiomass, matrix(NA, nrow = num_reclutas - num_filas, ncol = length(age)))
  }
  
  # Adding first age clas recruits
  M2[1:num_reclutas, 1] <- 0.1
  
  n1=nrow(M2)
  
  #calculating n of survivers for each age class
  for (i in 1:length(age)){ # time steps = years from 0 to 8
    s2=1-prob2[i+1]
    s1=1-prob2[i]
    p=(s1-s2)/s1 #-- # p for not surviving at t+1 given you did survive at t
    temp=sum(runif(n1[i],0,1)>p)
    n1=c(n1,temp) # p for not surviving at t+1 given you did survive at t
  }
  
  
  # Fill M2 with survivals at each age class and their respective length
  
  for (i in 2:length(age)) {
    M2[,i]<-lt1[i]
    M2[(n1[i]+1):n1[1], i] <- NA
  }
  
  # yearly stock biomass
  
  Mbiomass= a1*M2^b1 #matriz de masas (g): #LWR formula
  
  biomasa_total[year] <- sum(a1 * colSums(na.omit(M2)^b1)) # yearly total biomass evolution
  
  age_class_biomass[, year] = a1 * colSums(na.omit(M2)^b1)
  
  #number of infividuals: 
  n_population[year] = sum(!is.na(M2))
  
}

# Visualizar la biomasa total del stock para cada año
data.frame(biomasa_total)

plot(biomasa_total, ylab="Total biomass (g)", xlab = "Year", ylim = c(0,30000), main = paste("Population proyection with Fm=0"))
lines(biomasa_total)

data.frame(n_population)
plot(n_population, ylab="Number of individuals", xlab = "Year", ylim = c(0,1500), main = paste("Population proyection: 10 years"))
lines(n_population)
