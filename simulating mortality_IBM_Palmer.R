#------------------------------------------------------------------------------
# Simulating mortality in an IBM context
# M. Palmer
#------------------------------------------------------------------------------------------------------------------------------------

# last update: 18-4-2024


# 1: simplest case: only fishing mortality at constant rate
age=seq(0,100)
f=0.03
prob=(1 - exp(-f*age)) # survival prob = 1-e^(-f*age)
plot(age,prob)
#  for adding stochasticity: (1 - exp(-f)) should be compared with rnorm(0,1)
# This is equivalent to rnbinom
# example:
n=10000 # initial number of fish
for (i in 1:100){ # time steps
  temp=sum(runif(n[i],0,1)>(1 - exp(-f)))
  n=c(n,temp) # number of surviving fish
}
# theoretical versus simulation  
plot(age,n[1]-n, ylab="Number of dead fish") #cómo cambia el número de peces sobrevivientes a lo largo del tiempo según la simulación.
lines(age,n[1]*prob,col="red") #número teórico de peces que sobreviven en función de la edad, utilizando la probabilidad de supervivencia teórica calculada previamente y el número inicial de peces

# 2: adding natural mortality
m=0.01
prob=(1 - exp(-(f+m)*age)) # survival prob
plot(age,prob)

# 3: Length-dependent fishing mortality: f(length)
m=0.01
f_max=0.05 # maximum fishing mortality rate (f(length) will move between 0 (small length) and f_max)
length=seq(0,100) # here you should add your preferred growth model (that is, length=f(age))
# slope and inflection of a logistic function (or any other function)
f_slope=0.1 
f_inflection=50

f_length=1-(1/(1 + exp (f_slope * (length - f_inflection))))
plot(length,f_length, main = "Length-dependent fishing mortality")

prob=(1 - exp(-(f_max*f_length+m)*age)) # survival prob
plot(age,prob, main="Survival")

# 4: It is also possible to model a non constant natural mortality rate (Lorenzen 2022).

# 5: Differences in vulnerability
# You can deal with fish displaying different fishing mortality rate by setting 
#  different f_max. The simplest case: two values for f_max: more vulnerable fish 
# (=larger overlap between fish daily activity period and fisher daily
#  activity period) will display larger f_max
#
#...however, fish with larger daily activity period will have larger
# foraging time, thus they will display larger energy assimilation and they will growth 
# faster, thus it will modifiy also the dependency between (length-dependent) fishing mortality and age...

# 6: a final note of caution: If you are planing to use a large time step (say one year), 
# you should take into account that mortality cannot be considered constant
# within a given time step (length is changing continuously; it is not reliable that
# fish growth suddenly only at the last instant of the year).

# DISCRETIZATION model#1
age=seq(0,100)
f=0.03
prob=(1 - exp(-f*age)) # survival prob
n=10000 # initial number of fish
for (i in 1:100){ # time steps
  s2=1-prob[i+1]
  s1=1-prob[i]
  p=(s1-s2)/s1 # p for not surviving at t+1 given you did survive at t
  temp=sum(runif(n[i],0,1)>p) #n individuos sobreviven en t
  n=c(n,temp) # number of surviving fish
}
plot(age,n[1]-n, ylab="Number of dead fish")
lines(age,n[1]*prob,col="red",lwd=3)

# DISCRETIZATION model#3
m=0.01
f_max=0.05 # maximum fishing mortality rate (f(length) will move between 0 (small length) and f_max)
length=seq(0,100) # here you should add your preferred growth model (that is, length=f(age))
f_slope=0.5 # slope and inflection of a logistic function (or any other function)
f_inflection=50
f_length=1-(1/(1 + exp (f_slope * (length - f_inflection))))
plot(length,f_length)
prob=(1 - exp(-(f_max*f_length+m)*age)) # survival prob

n=10000 # initial number of fish
for (i in 1:100){ # time steps
  s2=1-prob[i+1]
  s1=1-prob[i]
  p=(s1-s2)/s1 #-- # p for not surviving at t+1 given you did survive at t
  temp=sum(runif(n[i],0,1)>p)
  n=c(n,temp) # p for not surviving at t+1 given you did survive at t
}
plot(age,n[1]-n)
lines(age,n[1]*prob,col="red",lwd=3)

# DISCRETIZATION model#3 at different dt (emulating, e.g., dt=1 year)
n2=10000 # initial number of fish
age2=seq(0,100,10)
length2=age2
f_length2=1-(1/(1 + exp (f_slope * (length2 - f_inflection))))
prob2=(1 - exp(-(f_max*f_length2+m)*age2)) # survival prob
for (i in 2:length(age2)){ # time steps
  s2=1-prob2[i]
  s1=1-prob2[i-1]
  p=(s1-s2)/s1 # probabilidad de no supervivencia 
  temp=sum(runif(n2[i-1],0,1)>p) #Se simula la mortalidad generando números aleatorios uniformemente distribuidos entre 0 y 1 (runif). Si estos números aleatorios son mayores que la probabilidad de no supervivencia (p), entonces los individuos correspondientes sobreviven y se suman a la cuenta temporal de mortalidad (temp).
  n2=c(n2,temp) # p for not surviving at t+1 given you did survive at t
}
plot(age2,n[1]-n2,xlab="Age",ylab="Number of dead fish",cex=1.2,pch=19,col="blue")
points(age,n[1]-n)
lines(age,n[1]*prob,col="red",lwd=3)
legend("bottomright",c("theoretic","simulation dt=1","simulation dt=10"),
       lwd=c(2,NA,NA),col=c("red","black","blue"),pch=c(NA,1,19))
