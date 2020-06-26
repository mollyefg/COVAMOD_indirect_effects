### SEIR model of vaccination scenarios

setwd("C:/Users/mgalla9/Dropbox/Emory/Coronavirus/code")
setwd("C:/Users/molly/Dropbox/Emory/Coronavirus/code")

#install.packages('deSolve','viridis')
require(deSolve)
require(viridis)

#########################		

### scenarios: no vaccine; vaccine that reduces mortality; vaccine that reduces transmission; vaccine that reduces both

SEIR<-function(t,y,p){
  
  # specify the initial conditions, for each age class
  
  S = y[1:nAgeClasses]
  E = y[(nAgeClasses+1):(2*nAgeClasses)]
  Ic = y[(2*nAgeClasses+1):(3*nAgeClasses)]
  Isc = y[(3*nAgeClasses+1):(4*nAgeClasses)]
  R = y[(4*nAgeClasses+1):(5*nAgeClasses)]
  M = y[(5*nAgeClasses+1):(6*nAgeClasses)] 
    
  Sv = y[(6*nAgeClasses+1):(7*nAgeClasses)]
  Ev = y[(7*nAgeClasses+1):(8*nAgeClasses)]
  Icv = y[(8*nAgeClasses+1):(9*nAgeClasses)]
  Iscv = y[(9*nAgeClasses+1):(10*nAgeClasses)]
  
  dX = c()	# vector to store the outputs
  
  with(as.list(p),{
    
    for(i in 1:nAgeClasses){
      
      vaccination_rate = c()
      if(t < tau){                    # optional time delay until vaccination starts
        vaccination_rate = c(0, 0)
      }else{
        vaccination_rate = lambda
      }
      
      # equation for dS.dt:
      dX[i] = - beta[i]*(S[i]/N)*(Ic[i] + alpha*Isc[i]) - betav[i]*(S[i]/N)*(Icv[i] + alpha*Iscv[i]) - vaccination_rate[i]*S[i]/N
      
      #equation for dE.dt: 	
      dX[i+nAgeClasses] = beta[i]*(S[i]/N)*(Ic[i] + alpha*Isc[i]) +  betav[i]*(S[i]/N)*(Icv[i] + alpha*Iscv[i]) - sigma*E[i]
      
      #equation for dIc.dt:
      dX[i + 2*nAgeClasses] = nu[i]*sigma*E[i] - gamma*Ic[i]
      
      #equation for dIsc.dt:
      dX[i + 3*nAgeClasses] = (1 - nu[i])*sigma*E[i] - gamma*Isc[i]
      
      #equation for dR.dt:
      dX[i + 4*nAgeClasses] = gamma*(Ic[i] + Isc[i] + Icv[i] + Iscv[i])
            
      #equation for dM.dt:
      dX[i + 5*nAgeClasses] = gamma*(rhoc[i]*(Ic[i]+Icv[i]) + rhosc[i]*(Isc[i]+Iscv[i])) 
      
      ##### now the equations for the vaccinated individuals:
      # equation for dSv.dt:
      dX[i + 6*nAgeClasses] = vaccination_rate[i]*S[i]/N - beta[i]*(Sv[i]/N)*(Ic[i] + alpha*Isc[i]) - betav[i]*(Sv[i]/N)*(Icv[i] + alpha*Iscv[i])
      
      # equation for dEv.dt:
      dX[i + 7*nAgeClasses] = beta[i]*(Sv[i]/N)*(Ic[i] + alpha*Isc[i]) + betav[i]*(Sv[i]/N)*(Icv[i] + alpha*Iscv[i]) - sigma*Ev[i]
      
      #equation for dIcv.dt:
      dX[i + 8*nAgeClasses] = nuv[i]*sigma*Ev[i] - gamma*Icv[i]
      
      #equation for dIscv.dt:
      dX[i + 9*nAgeClasses] = (1 - nuv[i])*sigma*Ev[i] - gamma*Iscv[i]
      
    } # closes the loop over age classes	
    
    return(list(dX)) 	
  }) # stops the parameter list		
}	# closes the function



###########################################################################		Define parameters & conditions

vaccineRate = 5000 

# define the vaccination scenarios:
scenarios = 4   

deaths = array(dim = c(scenarios, 4))   # the 4 columns are: perfect vaccination, all effects; perfect vaccination, direct effects only; 
                                        # all-or-nothing vaccination, all effects; and all-or-nothing vaccination, direct effects only

##### Universal values:

nAgeClasses = 1
alpha = 1; sigma = 1/4; gamma = 1/5; beta = 0.4; nu = 0.35
rhoc = 0.05; rhosc = 0.001

minT = 0; maxT = 365; dt = 1
t = seq(from=minT,to=maxT,by=dt)  # creating the vector of time for output
tau = 0 # delay before vaccination starts; prior to this, lambda = 0

E0 = 200 #100
Ic0 = 140
Isc0 = 260
R0 = 0 
M0 = 0
S0 = 1000000 - E0 - Ic0 - Isc0 - R0 - M0


initialConditions = c(S0, 
                      E0,
                      Ic0, Isc0, 
                      R0, M0,
                      0, 0, 0, 0)	# ordered by S[1]...S[nAgeClasses], then E[1]...E[nAgeClasses],etc

maxN = sum(initialConditions)
fixedN = maxN


###### Values for different scenarios:
reduction = 0.9	# reduce nu, beta by a this percent
factor = 1 - reduction

### 1: lambda = 0
lambda = 0; betav = beta; nuv = nu
p = list(N = fixedN, beta, betav, alpha, sigma, nu, nuv, gamma, rhoc, rhosc, lambda, tau)		 # creating the LIST of parameter values
out1 = ode(y=initialConditions, times=t, func=SEIR, parms=p, method = 'ode45') # run the ode solver

### 2: vaccination reduces clinical infections by... 90%?
lambda = vaccineRate; betav = beta; nuv = factor*nu
p = list(N = fixedN, beta, betav, alpha, sigma, nu, nuv, gamma, rhoc, rhosc, lambda, tau)		 # creating the LIST of parameter values
out2 = ode(y=initialConditions, times=t, func=SEIR, parms=p, method = 'ode45') # run the ode solver

### 3: vaccination reduces transmission by... 90%?
lambda = vaccineRate; betav = factor*beta; nuv = nu
p = list(N = fixedN, beta, betav, alpha, sigma, nu, nuv, gamma, rhoc, rhosc, lambda, tau)		 # creating the LIST of parameter values
out3 = ode(y=initialConditions, times=t, func=SEIR, parms=p, method = 'ode45') # run the ode solver

### 3: vaccination reduces both clinical infections AND transmission by... 90%?
lambda = vaccineRate; betav = factor*beta; nuv = factor*nu
p = list(N = fixedN, beta, betav, alpha, sigma, nu, nuv, gamma, rhoc, rhosc, lambda, tau)		 # creating the LIST of parameter values
out4 = ode(y=initialConditions, times=t, func=SEIR, parms=p, method = 'ode45') # run the ode solver


##### set up for plots:
par(mfrow = c(1, 2))
par(mai = c(0.9, 1, 0.5, 0.3))

end = length(t)

cols = viridis(n=4, begin = 0, end = .85)
legendnames = c("no vaccination", "direct protection", "indirect protection", "total effects")

##### plot the death dynamics:

dailyDeath = array(dim = c(4, end-1))
dailyDeath[1,] = c(out1[2:end,7]-out1[1:(end-1),7])
dailyDeath[2,] = c(out2[2:end,7]-out2[1:(end-1),7])
dailyDeath[3,] = c(out3[2:end,7]-out3[1:(end-1),7])
dailyDeath[4,] = c(out4[2:end,7]-out4[1:(end-1),7])

plot(dailyDeath[1,], xlab = "day", ylab = "deaths per day", xlim = c(0,250), cex.axis = 1.5, cex.lab = 1.5, type = "l", col = cols[1], lwd = 3)
lines(dailyDeath[2,], type = "l", col = cols[2], lwd = 3)
lines(dailyDeath[3,], type = "l", col = cols[3], lwd = 3)
lines(dailyDeath[4,], type = "l", col = cols[4], lwd = 3)
legend("topright", legend = legendnames, col = cols, lwd = 2, cex = 1.5, bty = "n")

##### plot the number of deaths averted for each scenario:

deaths = c(out1[end,7], out2[end,7], out3[end,7], out4[end,7])

averted = array(dim = c(1, length(deaths)))
averted[1,] = (deaths[1] - deaths)
averted; sum(averted[,2] + averted[,3])	# the total deaths averted is slightly less than the sum of these two, due to overlap (check this)

barplot(averted, beside = T, space = c(0, .1), names.arg = legendnames, main = "", xlab = "vaccination scenario", ylab = "number of deaths averted", col=cols,
            cex.axis = 1.5, cex.lab = 1.5, cex.names = 1.15)


##### end
