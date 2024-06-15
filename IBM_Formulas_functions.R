#Formulas - functions ----

##Von Bertalanffy model for length growth: L_inf: asintotic length, K: growth rate (y^-1), age, age0=t0=null growth theoric age
von_ber_fun<-function(L_inf, K, AGE, AGE0){ 
  L_inf*(1 - exp(-K *(AGE-AGE0)))
}

#LWR
lwr_fun<- function(a, b, L){
  a*L^b
}

####Mortality ###
#######

#NAtural mortality ###
M_fun <- function(M, L, s) {
  M*exp(-s * (L))  # Inverse logistic function
}
#plot(M_fun(Mmax,lt, Mslope), ylab = "Natural mortality prob", xlab = "Length")

### Behaviour induced natural mortality (activity and chronotype)
M_hr_fun<-function(M, HR, Shr){
  M * (1 - exp(-Shr * HR))
}

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

#Fishing mortality
prob_F_fun<-function(F_max, F_length, F_beha, F_chron, AGE){ ## F_max=maximum fishing mortality rate; M=natural mortality
  1-(exp(-(F_max*(F_length+F_beha+F_chron))*(AGE)))         
}


#natural mortality surv prob
prob_M_fun<-function(M_max, M, M_hr, M_chron){ ## M_max=maximum natural mortality rate; M=natural mortality
  1-(exp(-(M_max*(M+M_hr+M_chron))))         
}

####Recruitment model: Beverholt
recruit_fun<- function(alfa, beta, SpBio){
  alfa*SpBio/1+(beta*SpBio) 
}

#####Plotting functions: #####

