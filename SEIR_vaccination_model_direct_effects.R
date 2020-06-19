### Age-structured SEIR of vaccination scenarios

setwd("C:/Users/mgalla9/Dropbox/Emory/Coronavirus/code")
setwd("C:/Users/molly/Dropbox/Emory/Coronavirus/code")

#install.packages('deSolve','viridis')
require(deSolve)
require(viridis)

#########################			First, create and name the functions: base model - "SEIR" 
# perfect vaccination - "Perfect"
# all-or-nothing vaccination - "AllOrNothing"

######### scenarios to consider: 

####  Scenario 1:

# perfect vaccination: all vaccinated individuals are fully and permanently immune; 
# still assuming no waning immunity from natural infection as well -
# so you could run this 'Perfect' function with the vaccination rate lambda set to zero
# to get the same result as the SEIR function

# to determine the magnitude of direct effects of vaccination - the number of vaccinated individuals who would have died -
# versus indirect effects - the number of unvaccinated individuals who are protected by herd immunity -
# run simulations where vaccinated individuals behave like unvaccinated ones, but can't die of infection; so they contribute to transmission
# but are protected from death. The difference between the number of deaths in a model where vaccinated indivuals do not transmit
# and a model where they do, then, is the indirect effect of vaccination on mortality

# parameter 'f' here is set to zero if vaccination is providing group-level benefits, and to 1 if it is not

Perfect_direct<-function(t,y,p){
  
  # specify the initial conditions, for each age class
  
  S = y[1:nAgeClasses]
  E = y[(nAgeClasses+1):(2*nAgeClasses)]
  Ic = y[(2*nAgeClasses+1):(3*nAgeClasses)]
  Isc = y[(3*nAgeClasses+1):(4*nAgeClasses)]
  Rc = y[(4*nAgeClasses+1):(5*nAgeClasses)]
  Rsc = y[(5*nAgeClasses+1):(6*nAgeClasses)]
  M = y[(6*nAgeClasses+1):(7*nAgeClasses)] 
    
  Sv = y[(7*nAgeClasses+1):(8*nAgeClasses)]
  Ev = y[(8*nAgeClasses+1):(9*nAgeClasses)]
  Icv = y[(9*nAgeClasses+1):(10*nAgeClasses)]
  Iscv = y[(10*nAgeClasses+1):(11*nAgeClasses)]
  Rcv = y[(11*nAgeClasses+1):(12*nAgeClasses)]
  Rscv = y[(12*nAgeClasses+1):(13*nAgeClasses)]
  
  dX = c()	# vector to store the outputs
  
  with(as.list(p),{
    
    for(i in 1:nAgeClasses){
      
      vaccination_rate = c()
      if(t < tau){                    # optional time delay until vaccination starts
        vaccination_rate = c(0, 0)
      }else{
        vaccination_rate = lambda
      }
      
      
      # calculate the sum of contacts for the two infectious classes
      contactSumIc = 0; contactSumIsc = 0
      for(j in 1:nAgeClasses){
        contactSumIc = contactSumIc + C[i,j]*(Ic[j] + f*Icv[j])     # if f = 0, vaccinated individuals do not contribute to transmission
        contactSumIsc = contactSumIsc + C[i,j]*(Isc[j] + f*Iscv[j])
      }
      
      # equation for dS.dt:
      dX[i] = omegac[i]*Rc[i] + omegasc[i]*Rsc[i] - beta[i]*(S[i]/N)*(contactSumIc + alpha*contactSumIsc) - vaccination_rate[i]*S[i]/N
      
      #equation for dE.dt: 	
      dX[i+nAgeClasses] = beta[i]*(S[i]/N)*(contactSumIc + alpha*contactSumIsc) - sigma*E[i]
      
      #equation for dIc.dt:
      dX[i + 2*nAgeClasses] = nu[i]*sigma*E[i] - gamma*Ic[i]
      
      #equation for dIsc.dt:
      dX[i + 3*nAgeClasses] = (1 - nu[i])*sigma*E[i] - gamma*Isc[i]
      
      #equation for dRc.dt:
      dX[i + 4*nAgeClasses] = gamma*Ic[i] + nu[i]*vaccination_rate[i]*S[i]/N - omegac[i]*Rc[i]  
      
      #equation for dRsc.dt:
      dX[i + 5*nAgeClasses] = gamma*Isc[i] + (1 - nu[i])*vaccination_rate[i]*S[i]/N - omegasc[i]*Rsc[i]  
      
      #equation for dM.dt:
      dX[i + 6*nAgeClasses] = gamma*(rhoc[i]*Ic[i] + rhosc[i]*Isc[i]) 
      
      ##### now the equations for the vaccinated individuals:
      # equation for dSv.dt:
      dX[i + 7*nAgeClasses] = vaccination_rate[i]*S[i]/N + omegac[i]*Rcv[i] + omegasc[i]*Rscv[i] - beta[i]*(Sv[i]/N)*(contactSumIc + alpha*contactSumIsc) 
      
      # equation for dEv.dt:
      dX[i + 8*nAgeClasses] = beta[i]*(Sv[i]/N)*(contactSumIc + alpha*contactSumIsc) - sigma*Ev[i]
      
      #equation for dIcv.dt:
      dX[i + 9*nAgeClasses] = nu[i]*sigma*Ev[i] - gamma*Icv[i]
      
      #equation for dIscv.dt:
      dX[i + 10*nAgeClasses] = (1 - nu[i])*sigma*Ev[i] - gamma*Iscv[i]
      
      #equation for dRcv.dt:
      dX[i + 11*nAgeClasses] = gamma*Icv[i] - omegac[i]*Rcv[i]  
      
      #equation for dRscv.dt:
      dX[i + 12*nAgeClasses] = gamma*Iscv[i] - omegasc[i]*Rscv[i]  
      
    } # closes the loop over age classes	
    
    return(list(dX)) 	
  }) # stops the parameter list		
}	# closes the function


### similar to the perfect vaccination scenario, but with primary vaccine failure incorporated

AllOrNothing_direct<-function(t,y,p){
  
  # specify the initial conditions, for each age class
  
  S = y[1:nAgeClasses]
  EF = y[(nAgeClasses+1):(2*nAgeClasses)]
  E = y[(2*nAgeClasses+1):(3*nAgeClasses)]
  Ic = y[(3*nAgeClasses+1):(4*nAgeClasses)]
  Isc = y[(4*nAgeClasses+1):(5*nAgeClasses)]
  Rc = y[(5*nAgeClasses+1):(6*nAgeClasses)]
  Rsc = y[(6*nAgeClasses+1):(7*nAgeClasses)]
  M = y[(7*nAgeClasses+1):(8*nAgeClasses)] 
  
  Sv = y[(8*nAgeClasses+1):(9*nAgeClasses)]
  Ev = y[(9*nAgeClasses+1):(10*nAgeClasses)]
  Icv = y[(10*nAgeClasses+1):(11*nAgeClasses)]
  Iscv = y[(11*nAgeClasses+1):(12*nAgeClasses)]
  Rcv = y[(12*nAgeClasses+1):(13*nAgeClasses)]
  Rscv = y[(13*nAgeClasses+1):(14*nAgeClasses)]
  
  dX = c()	# vector to store the outputs
  
  with(as.list(p),{
    
    for(i in 1:nAgeClasses){
      
      vaccination_rate = c()
      if(t < tau){
        vaccination_rate = c(0, 0)
      }else{
        vaccination_rate = lambda
      }
      
      
      # calculate the sum of contacts for the two infectious classes
      contactSumIc = 0; contactSumIsc = 0
      for(j in 1:nAgeClasses){
        contactSumIc = contactSumIc + C[i,j]*(Ic[j] + f*Icv[j])
        contactSumIsc = contactSumIsc + C[i,j]*(Isc[j] + f*Iscv[j])
      }
      
      # equation for dS.dt:
      dX[i] = omegac[i]*Rc[i] + omegasc[i]*Rsc[i] - beta[i]*(S[i]/N)*(contactSumIc + alpha*contactSumIsc) - vaccination_rate[i]*S[i]/N
      
      #e quation for dEF.dt
      dX[i+nAgeClasses] = phi[i]*lambda[i]*(S[i]/N) - beta[i]*(EF[i]/N)*(contactSumIc + alpha*contactSumIsc)
      
      # equation for dE.dt
      dX[i+2*nAgeClasses] = beta[i]*((S[i]+EF[i])/N)*(contactSumIc + alpha*contactSumIsc) - sigma*E[i]
      
      #equation for dIc.dt:
      dX[i + 3*nAgeClasses] = nu[i]*sigma*E[i] - gamma*Ic[i]
      
      #equation for dIsc.dt:
      dX[i + 4*nAgeClasses] = (1 - nu[i])*sigma*E[i] - gamma*Isc[i]
      
      #equation for dRc.dt:
      dX[i + 5*nAgeClasses] = gamma*Ic[i] - omegac[i]*Rc[i]  
      
      #equation for dRsc.dt:
      dX[i + 6*nAgeClasses] = gamma*Isc[i] - omegasc[i]*Rsc[i]  
      
      #equation for dM.dt:
      dX[i + 7*nAgeClasses] = gamma*(rhoc[i]*Ic[i] + rhosc[i]*Isc[i]) 
      
      ###
      # equation for dSv.dt:
      dX[i + 8*nAgeClasses] = (1 - phi[i])*vaccination_rate[i]*S[i]/N + omegac[i]*Rcv[i] + omegasc[i]*Rscv[i] - beta[i]*(Sv[i]/N)*(contactSumIc + alpha*contactSumIsc) 
      
      # equation for dEv.dt:
      dX[i + 9*nAgeClasses] = beta[i]*(Sv[i]/N)*(contactSumIc + alpha*contactSumIsc) - sigma*Ev[i]
      
      #equation for dIcv.dt:
      dX[i + 10*nAgeClasses] = nu[i]*sigma*Ev[i] - gamma*Icv[i]
      
      #equation for dIscv.dt:
      dX[i + 11*nAgeClasses] = (1 - nu[i])*sigma*Ev[i] - gamma*Iscv[i]
      
      #equation for dRcv.dt:
      dX[i + 12*nAgeClasses] = gamma*Icv[i] - omegac[i]*Rcv[i]  
      
      #equation for dRscv.dt:
      dX[i + 13*nAgeClasses] = gamma*Iscv[i] - omegasc[i]*Rscv[i]  
      
    } # closes the loop over age classes	
    
    return(list(dX)) 	
  }) # stops the parameter list		
}	# closes the function


###########################################################################		End of functions

###########################################################################		Define parameters & conditions

# proportion of the population in each class:
classProportions = c(0.80, 0.20)  # 80% adults, 20% children

popStructure = T    # if true, class 1 and class 2 differ in contact rate, risk of death, risk of clinical infection, risk of vaccine failure

vaccineRate = 5000 

# define the vaccination scenarios:
scenarios = 4   
vacc_vals = array(dim = c(scenarios, 2))  # four different situations, two different types of vaccines

vacc_vals[1,] = c(0, 0)	# no vaccination
vacc_vals[2,] = c(vaccineRate*classProportions[1], vaccineRate*classProportions[2])	# vaccinate evenly
vacc_vals[3,] = c(vaccineRate, 0)	# vaccinate 'adults' - 
vacc_vals[4,] = c(0, vaccineRate)	# vaccinate 'children' - 

deaths = array(dim = c(scenarios, 4))   # the 4 columns are: perfect vaccination, all effects; perfect vaccination, direct effects only; 
                                        # all-or-nothing vaccination, all effects; and all-or-nothing vaccination, direct effects only

for(v in 1:4){

lambda = vacc_vals[v,]

##### Universal values:

nAgeClasses = 2
alpha = 0.80; sigma = 0.3; gamma = 0.18; beta = c(.44, .44) 

# with population structure;
rhoc = c(0.02, 0.001)
rhosc = c(0, 0)
nu = c(0.5, 0.1)
contacts = c(0.8, .9, .9, 1)

phi = c(0.12, 0.04)   # for primary vaccine failure (all-or-nothing)
#theta = c(0.75, 0.90) # degree of risk reduction as a result of vaccination (leaky)

if(popStructure == F){  # if we want to simulate an instructured population, where the classes still exist but behave identically
  # no population structure:
  rho = rep(rho[1]*classProportions[1] + rho[2]*classProportions[2], nAgeClasses)
  nu = rep(nu[1]*classProportions[1] + nu[2]*classProportions[2], nAgeClasses)
  contacts = rep(contacts[1]*classProportions[1]^2 + contacts[4]*classProportions[2]^2 + contacts[2]*classProportions[1]*classProportions[2]+ contacts[3]*classProportions[1]*classProportions[2], nAgeClasses^2)
  phi = rep(phi[1]*classProportions[1] + phi[2]*classProportions[2], nAgeClasses)   # for primary vaccine failure
  theta = rep(theta[1]*classProportions[1] + theta[2]*classProportions[2], nAgeClasses)
}

C = matrix(contacts, nrow = nAgeClasses, ncol = nAgeClasses)

minT = 0; maxT = 365*2; dt = 1

t = seq(from=minT,to=maxT,by=dt)  # creating the vector of time for output

E0 = 500 #100
Ic0 = 0
Isc0 = 0
R0 = 50000
S0 = 1000000 - E0 - R0

###### Values for different scenarios:

# for vaccination models:
tau = 0 # delay before vaccination starts; prior to this, lambda = 0

# finally, for waning immunity:
waningOff = c(0, 0)	# waning immunity rate
#waningOn = c(0.005, 0.005)	# waning immunity rate - don't worry about this for nwo


##### Perfect vaccination model:

initialConditions = c(S0*classProportions[1], S0*classProportions[2], 
                      E0*classProportions[1], E0*classProportions[2],
                      0, 0, 0, 0, 
                      R0*classProportions[1], R0*classProportions[2],
                      0, 0, 0, 0)	# ordered by S[1]...S[nAgeClasses], then E[1]...E[nAgeClasses],etc

maxN = sum(initialConditions)*classProportions[1]

fixedN = sum(initialConditions)

initialConditions_direct = c(initialConditions, rep(0, 6*nAgeClasses))

# for all vaccine effects, f = 0; for direct effects only, f = 1
p_all = list(N = fixedN, beta, C, alpha, sigma, nu, gamma, rhoc, rhosc, lambda, tau, phi, f = 0, omegac=waningOff, omegasc = waningOff)		 # creating the LIST of parameter values
p_direct = list(N = fixedN, beta, C, alpha, sigma, nu, gamma, rhoc, rhosc, lambda, tau, phi, f = 1, omegac=waningOff, omegasc = waningOff)

# run the simulations (perfect vaccination):

Perfect_allout = ode(y=initialConditions_direct, times=t, func=Perfect_direct, parms=p_all, method = 'ode45') # run the ode solver
Perfect_directout = ode(y=initialConditions_direct, times=t, func=Perfect_direct, parms=p_direct, method = 'ode45') # run the ode solver

##### All-Or-Nothing vaccination model: 

newinitialConditions = c(S0*classProportions[1], S0*classProportions[2], 
                      0, 0,    
                      E0*classProportions[1], E0*classProportions[2],
                      0, 0, 0, 0, 
                      R0*classProportions[1], R0*classProportions[2],
                      0, 0, 0, 0)	# ordered by S[1]...S[nAgeClasses], then EF (failed vaccine), then E[1]...E[nAgeClasses],etc

initialConditions_AON = c(newinitialConditions, rep(0, 6*nAgeClasses))

AON_allout = ode(y=initialConditions_AON, times=t, func=AllOrNothing_direct, parms=p_all, method = 'ode45') # run the ode solver

AON_directout = ode(y=initialConditions_AON, times=t, func=AllOrNothing_direct, parms=p_direct, method = 'ode45') # run the ode solver

deaths[v,1] = sum(Perfect_allout[731,14:15])
deaths[v,2] = sum(Perfect_directout[731,14:15])
deaths[v,3] = sum(AON_allout[731,16:17])
deaths[v,4] = sum(AON_directout[731,16:17])

}

rowN = c("no vaccination", "vaccinate evenly", "vaccinate 'adults' ", "vaccinate 'children'")
colN = c("Perfect", "Perfect, direct", "All-Or-Nothing", "All-Or-Nothing, direct")

rownames(deaths) = rowN
colnames(deaths) = colN
deaths

totalAvert_perfect = c(deaths[1,1] - deaths[,1])    # total deaths without vaccination, minus number of deaths under each scenario
directAvert_perfect = c(deaths[1,1] - deaths[,2])   # total deaths without vaccination, minus deaths with direct effects *only*
indirectAvert_perfect = c(deaths[,2] - deaths[,1])  # deaths with direct effects only, minus deaths with all effects

totalAvert_aon = c(deaths[1,3] - deaths[,3])    # total deaths without vaccination, minus number of deaths under each scenario
directAvert_aon = c(deaths[1,3] - deaths[,4])   # total deaths without vaccination, minus deaths with direct effects *only*
indirectAvert_aon = c(deaths[,4] - deaths[,1])  # deaths with direct effects only, minus deaths with all effects

### make the first plot, for perfect vaccination:
par(mai = c(1, 1, 1, .02))
cols = viridis(n=3, begin = 0, end = 0.5)
legendnames = c("total", "direct effect", "indirect effect")
aversion_perfect = rbind(totalAvert_perfect, directAvert_perfect, indirectAvert_perfect)
barplot(aversion_perfect, ylim = c(0, max(aversion_perfect, aversion_aon)), beside = T,
          legend = legendnames, args.legend = list(x=15.5, y = max(aversion), cex = 1.5, bty = "n"), 
            main = "", xlab = "Perfect Vaccination", ylab = "Number of deaths averted", col=cols,
            cex.axis = 1.5, cex.lab = 1.5, cex.names = 1.3)

### make the second plot, for all-or-nothing vaccination:
cols = viridis(n=3, begin = 0, end = 0.5)
legendnames = c("total", "direct effect", "indirect effect")
aversion_aon = rbind(totalAvert_aon, directAvert_aon, indirectAvert_aon)
barplot(aversion_aon, ylim = c(0, max(aversion_perfect, aversion_aon)), beside = T,
        legend = legendnames, args.legend = list(x=15.5, y = max(aversion), cex = 1.5, bty = "n"), 
        main = "", xlab = "All-or-Nothing Vaccination", ylab = "Number of deaths averted", col=cols,
        cex.axis = 1.5, cex.lab = 1.5, cex.names = 1.3)
###
