rm()
setwd("c:/Users/Usuario/Desktop/cosas beca jae/codigos/")

load("IBM_5_pop.RData")
#Biomass
Mbiomass= a1*M^b1 #matriz de masas (g): #LWR formula
w1<-Mbiomass[1,]

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


#scenarios with fishing mortality--------

y=100

# yearly total stock biomass
biomasa_total <- numeric(y+1) #change according to the years to proyect in the loop
biomasa_total[1]<- sum(B1)



#number of individuals per year
n_population<-numeric(y+1)
n_population[1]=sum(!is.na(M))

# initial population matrix
Mb0=Mbiomass #of biomass at age
Ml0=M #length at age

age_num<-numeric(length(age))
for (year in 1:length(age)) {
  age_num[year]=sum(!is.na(Mb0[,year])) 
}


#Matrix surv. prob. ----

#matrix of home range (hr) behaviour and chronotype (chr) dependant survival probabilities
be_scores0 <- runif(length(Mb0), 0, 100) #individual scores for hr 
f_beha0 = f_beha_fun(s1, be_scores0, be_in1)#individual fm prob dependant on hr behaviour

chron_scores0 <- runif(length(Mb0), 0, 100)
f_chron0 = f_chron_fun(s1, chron_scores0, be_in1)


Mhr0 <- matrix(f_beha0, nrow = nrow(Mb0), ncol = ncol(Mb0), byrow = TRUE)
Mchr0 <- matrix(f_beha0, nrow = nrow(Mb0), ncol = ncol(Mb0), byrow = TRUE)

for (i in 2:length(age)) {
  Mhr0[(age_num[i]+1):nrow(Mhr0), i] <- NA
  Mchr0[(age_num[i]+1):nrow(Mchr0), i] <- NA
}


Msurv0 = matrix(NA, nrow = nrow(Mb0), ncol = ncol(Mb0)) #empty matrix of probabilities

#for (a in 1:length(age)) {
#  for (i in 1:nrow(Mb0)) {
#    Msurv0[i,a]<-prob_fun(fmax1,f_length[a],Mhr0[i,a], Mchr0[i,a], m1, AGE = a)
#  }
#}

identical(na.omit(Mhr0),na.omit(Mchr0), na.omit(Mb0))
dim(Msurv0)
#?identical

Mhr=Mhr0
Mchr=Mchr0
Msurv=Msurv0
Mb=Mb0


#Msurv #surv. prob. matrix
#year=3
#remove(y)
# Simulación de la evolución de la población durante 10 años
for (year in 2:(y+1)) {
  
  # Cálculo del reclutamiento utilizando la ecuación de Beverton-Holt
  #recruitmen parameters, with varaibility among years 
  ar <- rnorm(8, mean = 15, sd = 4) #(every iteration gets a different value)
  br <- rnorm(8, mean = 0.2, sd = 0.0)
  
  SpBio <- sum(Mb, na.rm = T)
  #Beverton-Holt equation: number or recruits
  R <- ar*SpBio/((1+(br*SpBio)))
  
  #elder than 7 day and new recruits column is added
  Mb <- cbind(NA, Mb[, 1:7]) #new population matrix 
  Mhr <- cbind(NA, Mhr[, 1:7]) #new hr scores matrix
  Mchr <- cbind(NA, Mchr[, 1:7]) #new hr scores matrix
  
  # adjust matrix dimensions to add recruits
  num_reclutas <- sum(floor(R)) #floor rounds to lower integer
  num_filas <- nrow(Mb)
  if (num_reclutas > num_filas) {
    Mb <- rbind(Mb, matrix(NA, nrow = num_reclutas - num_filas, ncol = length(age)))
    Mhr <- rbind(Mhr, matrix(NA, nrow = num_reclutas - num_filas, ncol = length(age)))
    Mchr <- rbind(Mchr, matrix(NA, nrow = num_reclutas - num_filas, ncol = length(age)))
    Msurv <- rbind(Msurv, matrix(NA, nrow = num_reclutas - num_filas, ncol = length(age)))
  } else if (num_reclutas < num_filas) {
    Mb <- Mb[1:num_reclutas, ]
    Mhr <- Mhr[1:num_reclutas, ]
    Mchr <- Mchr[1:num_reclutas, ]
    Msurv <- Msurv[1:num_reclutas, ]
  }
  
  
  
  identical(na.omit(Mhr), na.omit(Mchr))
  dim(Mb)
  #updating sizes at age
  # Adding first age class recruits
  Mb[1:num_reclutas, 1] <- w1[1] #biomass at class age 1
  
  for (b in 2:length(age)) {
    for (i in 1:nrow(Mb)) {
      if (!is.na(Mb[i, b])) {
        Mb[i, b] <- w1[b]
      }
    }
  }
  
  # obtaining home range and chronotype behaviour scores
  be_scoresi <- runif(num_reclutas, 0, 100)
  chron_scoresi <- runif(num_reclutas, 0, 100)
  
  # Calcular la probabilidad de supervivencia para los nuevos reclutas
  f_beha <- f_beha_fun(fslope1, be_scoresi, be_in1)
  f_chr <- f_chron_fun(fslope1, chron_scoresi, chron_in1)
  prob_r <- prob_fun(fmax1, f_length1[1], f_beha, f_chr, m1, AGE = 1) 
  
  plot(prob_r)
  
 
  Mhr[,1]<-f_beha
  Mchr[,1]<-f_chr
  
  #fill surv prob matrix
  for (a in 2:length(age)){ # time steps = years from 0 to 8
    #i=1
    for (n in 1:nrow(Mb)) {
      Msurv[n,a]<-prob_fun(fmax1,f_length[a],Mhr[n,a], Mchr[n,a], m1, AGE = a)
    }
  }
  
  Msurv[,1]<-prob_r #adding recruits surv. prob.
  
  identical(na.omit(Mhr), na.omit(Mchr))
  sum(is.na(Mhr))
  sum(is.na(Mchr))
  sum(is.na(Msurv))
  #applying mortality
    
  #Palmer approach...
    #s2=1-prob1[i+1]
    #s1=1-prob1[i]
    #p=(s1-s2)/s1 #-- # p for not surviving at t+1 given you did survive at t
    #temp=sum(runif(n1[i],0,1)>p)
    #n1=c(n1,temp) # p for not surviving at t+1 given you did survive at t
  vsurv<-rep( NA, length(age))
  #binomial probability distribution
  for (a in 2:length(age)) {
    n1<-sum(!is.na(Mb[,a]))
    d=0
    for (r in 1:n1) {
      if(rbinom(1, 1, Msurv[r,a])==1){ #1 equals success in the binom prob
        d=d+1 #count of survivals
      }
    }
    vsurv[a]<-d  #number of survivals per age class
  }
  
  vsurv[1]<- sum(!is.na(Mb[,1]))
  

  # Fill Matrix with survivals at each age class and their respective length
  
  for (i in 2:length(age)) {
    Mb[,i]<-w1[i]
    Mb[(vsurv[i]+1):vsurv[1], i] <- NA
  }
  
  # yearly stock biomass
  
  biomasa_total[year] <- sum(na.omit(Mb)) # yearly total biomass evolution
  
  age_class_biomass[, year] = colSums(na.omit(Mb))
  
  #number of infividuals: 
  n_population[year] = sum(!is.na(Mb))
  
}

# Visualizar la biomasa total del stock para cada año
data.frame(biomasa_total)

plot(biomasa_total, ylab="Total biomass (g)", xlab = "Year", ylim = c(0,6000), main = paste("Population proyection with Fm=0"))
lines(biomasa_total)

data.frame(n_population)
plot(n_population, ylab="Number of individuals", xlab = "Year", ylim = c(0,1600), main = paste("Population proyection: 10 years"), type = "l")
#lines(n_population)
