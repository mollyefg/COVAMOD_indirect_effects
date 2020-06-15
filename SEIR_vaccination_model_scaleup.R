### Age-structured SEIR of vaccination scenarios

setwd("C:/Users/mgalla9/Dropbox/Emory/Coronavirus/code")
setwd("C:/Users/molly/Dropbox/Emory/Coronavirus/code")
#install.packages('deSolve', 'viridis')
require(deSolve); require(viridis)

#########################			First, create and name the functions: base model - "SEIR" 
													# perfect vaccination - "Perfect"
													# all-or-nothing vaccination - "AllOrNothing"
													# leaky vaccination - "Leaky"
SEIR<-function(t,y,p){

	# specify the initial conditions, for each age class

	S = y[1:nAgeClasses]
	E = y[(nAgeClasses+1):(2*nAgeClasses)]
	Ic = y[(2*nAgeClasses+1):(3*nAgeClasses)]
	Isc = y[(3*nAgeClasses+1):(4*nAgeClasses)]
	R = y[(4*nAgeClasses+1):(5*nAgeClasses)]
	
	dX = c()	# vector to store the outputs

	with(as.list(p),{

	for(i in 1:nAgeClasses){

		# calculate the sum of contacts for the two infectious classes
		contactSumIc = 0; contactSumIsc = 0
		for(j in 1:nAgeClasses){
			contactSumIc = contactSumIc + C[i,j]*Ic[j]
			contactSumIsc = contactSumIsc + C[i,j]*Isc[j]
			}

		# equation for dS.dt:
		dX[i] = omega[i]*R[i] -beta[i]*S[i]*contactSumIc - alpha*beta[i]*S[i]*contactSumIsc - mu*S[i] 	

		#equation for dE.dt: 	
		dX[i+nAgeClasses] = beta[i]*S[i]*contactSumIc + alpha*beta[i]*S[i]*contactSumIsc - (sigma + mu)*E[i]

		#equation for dIc.dt:
		dX[i + 2*nAgeClasses] = nu[i]*sigma*E[i] - (gamma + rho[i] + mu)*Ic[i]

		#equation for dIsc.dt:
		dX[i + 3*nAgeClasses] = (1 - nu[i])*sigma*E[i] - (gamma + mu)*Isc[i]

		#equation for dR.dt:
		dX[i + 4*nAgeClasses] = gamma*Ic[i] + gamma*Isc[i] - omega[i]*R[i] - mu*R[i]


		} # closes the loop over age classes	
		
	return(list(dX)) 	
		}) # stops the parameter list		
	}	# closes the function


######### scenarios to consider: 

####  Scenario 1:

# perfect vaccination: all vaccinated individuals are fully and permanently immune; 
# still assuming no waning immunity from natural infection as well -
# so you could run this 'Perfect' function with the vaccination rate lambda set to zero
# to get the same result as the SEIR function

Perfect<-function(t,y,p){

	# specify the initial conditions, for each age class

	S = y[1:nAgeClasses]
	E = y[(nAgeClasses+1):(2*nAgeClasses)]
	Ic = y[(2*nAgeClasses+1):(3*nAgeClasses)]
	Isc = y[(3*nAgeClasses+1):(4*nAgeClasses)]
	R = y[(4*nAgeClasses+1):(5*nAgeClasses)]
	
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
			contactSumIc = contactSumIc + C[i,j]*Ic[j]
			contactSumIsc = contactSumIsc + C[i,j]*Isc[j]
			}

		# equation for dS.dt:
		dX[i] = omega[i]*R[i] - beta[i]*S[i]*contactSumIc - alpha*beta[i]*S[i]*contactSumIsc - vaccination_rate[i]*S[i] - mu*S[i]

		#equation for dE.dt: 	
		dX[i+nAgeClasses] = beta[i]*S[i]*contactSumIc + alpha*beta[i]*S[i]*contactSumIsc - (sigma + mu)*E[i]

		#equation for dIc.dt:
		dX[i + 2*nAgeClasses] = nu[i]*sigma*E[i] - (gamma + rho[i] + mu)*Ic[i]

		#equation for dIsc.dt:
		dX[i + 3*nAgeClasses] = (1 - nu[i])*sigma*E[i] - (gamma + mu)*Isc[i]

		#equation for dR.dt:
		dX[i + 4*nAgeClasses] = gamma*Ic[i] + gamma*Isc[i] + vaccination_rate[i]*S[i] - omega[i]*R[i] - mu*R[i] 


		} # closes the loop over age classes	
		
	return(list(dX)) 	
		}) # stops the parameter list		
	}	# closes the function


####  Scenario 2:
# Primary vaccine failure, also referred to as ``all-or-nothing'' vaccines or vaccines with some failure to take
# represent a scenario in which some proportion of vaccinated hosts receive perfect immunity, 
# but a subset of vaccinated hosts receive no protection whatsoever (Vf) and are just as easily infected as unvaccinated, susceptible hosts. 
# These vaccinated, yet still susceptible, hosts receive their own compartment to represent an inability to get revaccinated, 
# yet otherwise behave identically as $S_1$ hosts. Measles and rubella are an example of ``all-or-nothing'' vaccines \citep{Farrington2003}. 

AllOrNothing<-function(t,y,p){

	# specify the initial conditions, for each age class

	S = y[1:nAgeClasses]
	Vf = y[(nAgeClasses+1):(2*nAgeClasses)]
	E = y[(2*nAgeClasses+1):(3*nAgeClasses)]
	Ic = y[(3*nAgeClasses+1):(4*nAgeClasses)]
	Isc = y[(4*nAgeClasses+1):(5*nAgeClasses)]
	R = y[(5*nAgeClasses+1):(6*nAgeClasses)]
	
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
			contactSumIc = contactSumIc + C[i,j]*Ic[j]
			contactSumIsc = contactSumIsc + C[i,j]*Isc[j]
			}

		# equation for dS.dt:
		dX[i] = omega[i]*R[i] - beta[i]*S[i]*contactSumIc - alpha*beta[i]*S[i]*contactSumIsc - mu*S[i] - vaccination_rate[i]*S[i]

		# equation for dVf.dt:
		dX[i+nAgeClasses] = phi[i]*vaccination_rate[i]*S[i] - beta[i]*Vf[i]*contactSumIc - alpha*beta[i]*Vf[i]*contactSumIsc - mu*Vf[i]

		#equation for dE.dt: 	
		dX[i+2*nAgeClasses] = beta[i]*(S[i]+Vf[i])*contactSumIc + alpha*beta[i]*S[i]*contactSumIsc - (sigma + mu)*E[i]

		#equation for dIc.dt:
		dX[i + 3*nAgeClasses] = nu[i]*sigma*E[i] - (gamma + rho[i] + mu)*Ic[i]

		#equation for dIsc.dt:
		dX[i + 4*nAgeClasses] = (1 - nu[i])*sigma*E[i] - (gamma + mu)*Isc[i]

		#equation for dR.dt:
		dX[i + 5*nAgeClasses] = gamma*Ic[i] + gamma*Isc[i] + (1 - phi[i])*vaccination_rate[i]*S[i] - omega[i]*R[i] - mu*R[i]


		} # closes the loop over age classes	
		
	return(list(dX)) 	
		}) # stops the parameter list		
	}	# closes the function
	

####  Scenario 3:
# Leaky vaccines, or vaccines that fail in degree, are when vaccinated hosts are less likely to be infected upon contact with infected hosts 
# relative to unvaccinated hosts. 
# This could be modeled by adding some proportionality constant that represents vaccinated hosts' degree of protection from transmission. 
# Upon infection, however, they behave like unvaccinated hosts which are infected. 
# Pertussis and malaria vaccines are examples of leaky vaccines

Leaky<-function(t,y,p){

	# specify the initial conditions, for each age class

	S = y[1:nAgeClasses]
	V = y[(nAgeClasses+1):(2*nAgeClasses)]
	E = y[(2*nAgeClasses+1):(3*nAgeClasses)]
	Ic = y[(3*nAgeClasses+1):(4*nAgeClasses)]
	Isc = y[(4*nAgeClasses+1):(5*nAgeClasses)]
	R = y[(5*nAgeClasses+1):(6*nAgeClasses)]
	
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
			contactSumIc = contactSumIc + C[i,j]*Ic[j]
			contactSumIsc = contactSumIsc + C[i,j]*Isc[j]
			}

		# equation for dS.dt:
		dX[i] = omega[i]*(V[i] + R[i]) -beta[i]*S[i]*contactSumIc - alpha*beta[i]*S[i]*contactSumIsc - mu*S[i] - vaccination_rate[i]*S[i]

		# equation for dV.dt:
		dX[i+nAgeClasses] = vaccination_rate[i]*S[i] - (1 - theta[i])*beta[i]*V[i]*contactSumIc - alpha*(1 - theta[i])*beta[i]*V[i]*contactSumIsc - omega[i]*V[i] -mu*V[i]

		#equation for dE.dt: 	
		dX[i+2*nAgeClasses] = beta[i]*(S[i]+(1 - theta[i])*V[i])*contactSumIc + alpha*beta[i]*(S[i]+(1 - theta[i])*V[i])*contactSumIsc - (sigma + mu)*E[i]

		#equation for dIc.dt:
		dX[i + 3*nAgeClasses] = nu[i]*sigma*E[i] - (gamma + rho[i] + mu)*Ic[i]

		#equation for dIsc.dt:
		dX[i + 4*nAgeClasses] = (1 - nu[i])*sigma*E[i] - (gamma + mu)*Isc[i]

		#equation for dR.dt:
		dX[i + 5*nAgeClasses] = gamma*Ic[i] + gamma*Isc[i] - omega[i]*R[i] - mu*R[i] 

		} # closes the loop over age classes	
		
	return(list(dX)) 	
		}) # stops the parameter list		
	}	# closes the function


###########################################################################		End of functions
###########################################################################		

vacc_vals = array(dim = c(3, 2))
vacc_vals[1,] = c(0.01, 0.01)	# vaccinate evenly
vacc_vals[2,] = c(0.02, 0.0)	# vaccinate 'adults' - 
vacc_vals[3,] = c(0.0, 0.02)	# vaccinate 'children' - 

totalDeath = array(dim = c(3, 8))

ptm <- proc.time()
for(v in 1:3){

lambda = vacc_vals[v,]

##### Universal values:

nAgeClasses = 2
beta = c(0.2, 0.5); alpha = 0.2; mu = 0; sigma = 0.2; gamma = 0.18; rho = c(0.005, 0.00015)

nu = c(0.1, 0.4)

contacts = rep(5e-6, 4)
C = matrix(contacts, nrow = nAgeClasses, ncol = nAgeClasses)

minT = 0; maxT = 365; dt = 1

t = seq(from=minT,to=maxT,by=dt)  # creating the vector of time for output

S0 = 1000000 #330000000
E0 = 500 #100
Ic0 = 0
Isc0 = 0
R0 = 100
V0 = 0

# proportion of the population in each class:
classProportions = c(0.8, 0.2)  # 80% adults, 20% children

###### Values for different scenarios:

# for vaccination models:
tau = 20 # delay before vaccination starts; prior to this, lambda = 0

# for primary vaccine failure
phi = c(0.15, 0.02) # vaccine failure rate, phi

# for 'leaky' vaccines
theta = c(0.60, 0.90) # degree of risk reduction as a result of vaccination

# finally, for waning immunity:
waningOff = c(0, 0)	# waning immunity rate
waningOn = c(0.005, 0.005)	# waning immunity rate


##### Base model:
	
initialConditions = c(S0*classProportions[1], S0*classProportions[2], 
			    E0*classProportions[1], E0*classProportions[2],
			    0, 0, 0, 0, 0, 0)	# ordered by S[1]...S[nAgeClasses], then E[1]...E[nAgeClasses],etc

maxN = sum(initialConditions)*classProportions[1]


p = list(beta, C, alpha, mu, sigma, nu, gamma, rho, omega=waningOff)		 # creating the LIST of parameter values
p_wane = list(beta, C, alpha, mu, sigma, nu, gamma, rho, omega = waningOn)


# run the simulations:
SEIR_out = ode(y=initialConditions, times=t, func=SEIR, parms=p, method = 'ode45') # run the ode solver
sum(SEIR_out[max(t),2:11]) + sum(rho[1]*SEIR_out[,6]) + sum(rho[2]*SEIR_out[,7]) - sum(initialConditions)

SEIR_wane = ode(y=initialConditions, times=t, func=SEIR, parms=p_wane, method = 'ode45') # run the ode solver


##### Perfect vaccination model:

p = list(beta, C, alpha, mu, sigma, nu, gamma, rho, lambda, tau, omega=waningOff)		 # creating the LIST of parameter values
p_wane = list(beta, C, alpha, mu, sigma, nu, gamma, rho, lambda, tau, omega=waningOn)		 # creating the LIST of parameter values

# run the simulation:

Perfect_out = ode(y=initialConditions, times=t, func=Perfect, parms=p, method = 'ode45') # run the ode solver
Perfect_wane = ode(y=initialConditions, times=t, func=Perfect, parms=p_wane, method = 'ode45') # run the ode solver


##### All-Or-Nothing vaccination model:

initialConditions_AllOrNothing = c(S0*classProportions[1], S0*classProportions[2], 
			    0, 0,
			    E0*classProportions[1], E0*classProportions[2],
			    0, 0, 0, 0, 0, 0)	# ordered by S[1]...S[nAgeClasses], Vf[1]...Vf[nAgeClasses], then E[1]...E[nAgeClasses],etc

p = list(beta, C, alpha, mu, sigma, nu, gamma, rho, lambda, tau, phi, omega=waningOff)		 # creating the LIST of parameter values
p_wane = list(beta, C, alpha, mu, sigma, nu, gamma, rho, lambda, tau, phi, omega=waningOn)		 # creating the LIST of parameter values

# run the simulation:

AllOrNothing_out = ode(y=initialConditions_AllOrNothing, times=t, func=AllOrNothing, parms=p, method = 'ode45') # run the ode solver
AllOrNothing_wane = ode(y=initialConditions_AllOrNothing, times=t, func=AllOrNothing, parms=p_wane, method = 'ode45') # run the ode solver


##### Leaky vaccination model:

initialConditions_Leaky = c(S0*classProportions[1], S0*classProportions[2], 
			    0, 0, 
			    E0*classProportions[1], E0*classProportions[2],
			    0, 0, 0, 0, 0, 0)	# ordered by S[1]...S[nAgeClasses], V[1]...V[nAgeClasses], then E[1]...E[nAgeClasses],etc

p = list(beta, C, alpha, mu, sigma, nu, gamma, rho, lambda, tau, theta, omega=waningOff)		 # creating the LIST of parameter values
p_wane = list(beta, C, alpha, mu, sigma, nu, gamma, rho, lambda, tau, theta, omega=waningOn)		 # creating the LIST of parameter values

# run the simulation:	

Leaky_out = ode(y=initialConditions_Leaky, times=t, func=Leaky, parms=p, method = 'ode45') # run the ode solver
Leaky_wane = ode(y=initialConditions_Leaky, times=t, func=Leaky, parms=p_wane, method = 'ode45') # run the ode solver



###########################################################################		Plot the results

cols = c("black", "gray", "red", "hotpink", "blue", "green")
outnames = c("Base SEIR Model, No Vaccination or Waning Immunity", 
		"Base SEIR Model, Waning Immunity, No Vaccination", 
		"Perfect Vaccination, No Waning Immunity", 
		"Perfect Vaccination, Waning Immunity", 
		"All or Nothing Vaccination, No Waning Immunity", 
		"All or Nothing Vaccination, Waning Immunity", 
		"Leaky Vaccination, No Waning Immunity", 
		"Leaky Vaccination, Waning Immunity") 

outarray = array(dim = c(length(outnames), length(t), ncol(Leaky_wane)))
outarray[1,,1:ncol(SEIR_out)] = SEIR_out
outarray[2,,1:ncol(SEIR_wane)] = SEIR_wane
outarray[3,,1:ncol(Perfect_out)] = Perfect_out
outarray[4,,1:ncol(Perfect_wane)] = Perfect_wane
outarray[5,,1:ncol(AllOrNothing_out)] = AllOrNothing_out
outarray[6,,1:ncol(AllOrNothing_wane)] = AllOrNothing_wane
outarray[7,,1:ncol(Leaky_out)] = Leaky_out
outarray[8,,1:ncol(Leaky_wane)] = Leaky_wane


### for i = 1, no waning immunity; for i = 2, immunity wanes

par(mfrow = c(2, 2))
old_mar = par("mar")
par(mar = c(5.1, 5, 4.1, 2.1))

totalIc_1 = c(); totalIc_2 = c()
deaths_1 = c(); deaths_2 = c()

for(i in 1:2){
	j = c()
	names = c()
	if(i==1){
		j = c(1, 3, 5, 7)
		}else{
		j = c(2, 4, 6, 8)
		}

	names = c(outnames[j[1]], outnames[j[2]], outnames[j[3]], outnames[j[4]])

	# plot the base model:
	output = outarray[j[1],,]

	totalIc_1 = c(totalIc_1, dt*sum(nu[1]*sigma*output[,4]))
	totalIc_2 = c(totalIc_2, dt*sum(nu[2]*sigma*output[,5]))
	deaths_1 = c(deaths_1, dt*sum(rho[1]*output[,6]))
	deaths_2 = c(deaths_2, dt*sum(rho[2]*output[,7]))

	plot(t, output[,2], xlab = "time (days)", ylab = "number of individuals", lty = "solid", col = cols[1],
		lwd = 2,  type = "l", cex.axis = 1.5, cex.lab = 1.5, ylim = c(0,maxN+1), main = names[1])
	lines(t, output[,3], lty = "dashed", col = cols[1], lwd = 2,  type = "l")
	lines(t, output[,4], lty = 'solid', col = cols[2], lwd = 2, type = "l")
	lines(t, output[,5], lty = 'dashed', col = cols[2], lwd = 2, type = "l")
	lines(t, output[,6], lty = 'solid', col = cols[3], lwd = 2, type = "l")
	lines(t, output[,7], lty = 'dashed', col = cols[3], lwd = 2, type = "l")
	lines(t, output[,8], lty = 'solid', col = cols[4], lwd = 2, type = "l")
	lines(t, output[,9], lty = 'dashed', col = cols[4], lwd = 2, type = "l")
	lines(t, output[,10], lty = 'solid', col = cols[5], lwd = 2, type = "l")
	lines(t, output[,11], lty = 'dashed', col = cols[5], lwd = 2, type = "l")

	legend("right", legend = c("S", "E", "Ic", "Isc", "R"), col = cols, lwd = 2,  cex = 1.2, bty = "n")
	legend(35, 300, legend = c("age class 1", "age class 2"), col = 'black', lty = c(1, 2), lwd = 2, cex = 1.2, bty = "n")

	### end base model figure


	##### Scenario 1: perfect vaccination
	output = outarray[j[2],,]

	totalIc_1 = c(totalIc_1, dt*sum(nu[1]*sigma*output[,4]))
	totalIc_2 = c(totalIc_2, dt*sum(nu[2]*sigma*output[,5]))
	deaths_1 = c(deaths_1, dt*sum(rho[1]*output[,6]))
	deaths_2 = c(deaths_2, dt*sum(rho[2]*output[,7]))

	plot(t, output[,2], xlab = "time (days)", ylab = "number of individuals", lty = "solid", col = cols[1],
		lwd = 2,  type = "l", cex.axis = 1.5, cex.lab = 1.5, ylim = c(0,maxN+1), main = names[2])
	lines(t, output[,3], lty = "dashed", col = cols[1], lwd = 2,  type = "l")
	lines(t, output[,4], lty = 'solid', col = cols[2], lwd = 2, type = "l")
	lines(t, output[,5], lty = 'dashed', col = cols[2], lwd = 2, type = "l")
	lines(t, output[,6], lty = 'solid', col = cols[3], lwd = 2, type = "l")
	lines(t, output[,7], lty = 'dashed', col = cols[3], lwd = 2, type = "l")
	lines(t, output[,8], lty = 'solid', col = cols[4], lwd = 2, type = "l")
	lines(t, output[,9], lty = 'dashed', col = cols[4], lwd = 2, type = "l")
	lines(t, output[,10], lty = 'solid', col = cols[5], lwd = 2, type = "l")
	lines(t, output[,11], lty = 'dashed', col = cols[5], lwd = 2, type = "l")

	legend("right", legend = c("S", "E", "Ic", "Isc", "R"), col = cols, lwd = 2,  cex = 1.2, bty = "n")
	legend(35, 300, legend = c("age class 1", "age class 2"), col = 'black', lty = c(1, 2), lwd = 2, cex = 1.2, bty = "n")


	##### Scenario 2: all-or-nothing vaccination
	output = outarray[j[3],,]

	totalIc_1 = c(totalIc_1, dt*sum(nu[1]*sigma*output[,6]))
	totalIc_2 = c(totalIc_2, dt*sum(nu[2]*sigma*output[,7]))
	deaths_1 = c(deaths_1, dt*sum(rho[1]*output[,8]))
	deaths_2 = c(deaths_2, dt*sum(rho[2]*output[,9]))

	plot(t, output[,2]+output[,4], xlab = "time (days)", ylab = "number of individuals", lty = "solid", col = cols[1],
		lwd = 2,  type = "l", cex.axis = 1.5, cex.lab = 1.5, ylim = c(0,maxN+1), main = names[3])
	lines(t, output[,3]+output[,5], lty = "dashed", col = cols[1], lwd = 2,  type = "l")
	lines(t, output[,6], lty = 'solid', col = cols[2], lwd = 2, type = "l")
	lines(t, output[,7], lty = 'dashed', col = cols[2], lwd = 2, type = "l")
	lines(t, output[,8], lty = 'solid', col = cols[3], lwd = 2, type = "l")
	lines(t, output[,9], lty = 'dashed', col = cols[3], lwd = 2, type = "l")
	lines(t, output[,10], lty = 'solid', col = cols[4], lwd = 2, type = "l")
	lines(t, output[,11], lty = 'dashed', col = cols[4], lwd = 2, type = "l")
	lines(t, output[,12], lty = 'solid', col = cols[5], lwd = 2, type = "l")
	lines(t, output[,13], lty = 'dashed', col = cols[5], lwd = 2, type = "l")

	legend("right", legend = c("S+Vf", "E", "Ic", "Isc", "R"), col = cols, lwd = 2,  cex = 1.2, bty = "n")
	legend(35, 300, legend = c("age class 1", "age class 2"), col = 'black', lty = c(1, 2), lwd = 2, cex = 1.2, bty = "n")


	##### Scenario 3: leaky vaccination
	output = outarray[j[4],,]

	totalIc_1 = c(totalIc_1, dt*sum(nu[1]*sigma*output[,6]))
	totalIc_2 = c(totalIc_2, dt*sum(nu[2]*sigma*output[,7]))
	deaths_1 = c(deaths_1, dt*sum(rho[1]*output[,8]))
	deaths_2 = c(deaths_2, dt*sum(rho[1]*output[,9]))

	plot(t, output[,2], xlab = "time (days)", ylab = "number of individuals", lty = "solid", col = cols[1],
		lwd = 2,  type = "l", cex.axis = 1.5, cex.lab = 1.5, ylim = c(0,maxN+1), main = names[4])
	lines(t, output[,3], lty = "dashed", col = cols[1], lwd = 2,  type = "l")

	lines(t, output[,4], lty = 'solid', col = cols[6], lwd = 2, type = "l")	# V[1]
	lines(t, output[,5], lty = 'dashed', col = cols[6], lwd = 2, type = "l")# V[2]

	lines(t, output[,6], lty = 'solid', col = cols[2], lwd = 2, type = "l")
	lines(t, output[,7], lty = 'dashed', col = cols[2], lwd = 2, type = "l")
	lines(t, output[,8], lty = 'solid', col = cols[3], lwd = 2, type = "l")
	lines(t, output[,9], lty = 'dashed', col = cols[3], lwd = 2, type = "l")
	lines(t, output[,10], lty = 'solid', col = cols[4], lwd = 2, type = "l")
	lines(t, output[,11], lty = 'dashed', col = cols[4], lwd = 2, type = "l")
	lines(t, output[,12], lty = 'solid', col = cols[5], lwd = 2, type = "l")
	lines(t, output[,13], lty = 'dashed', col = cols[5], lwd = 2, type = "l")

	legend("right", legend = c("S", "E", "Ic", "Isc", "R", "V"), col = cols, lwd = 2,  cex = 1.2, bty = "n")
	legend(35, 300, legend = c("age class 1", "age class 2"), col = 'black', lty = c(1, 2), lwd = 2, cex = 1.2, bty = "n")

	}	# end loop over i (immune waning off/on)

############################################ End plots


### store death terms for outer plot:

totalDeath[v,] = 100*(deaths_1 + deaths_2)/sum(initialConditions)

} # close loop over vaccination options

proc.time() - ptm

###

cols = viridis(n=3, begin = 0, end = 0.5)
cols = rep(c(cols[1], cols[2], cols[3]), 4)

legendnames = c("Vaccinate evenly", "Vaccinate class 1 ('adults')", "Vaccinate class 2 ('children')")
deaths = totalDeath
colnames(deaths) = rep(c("None", "Perfect", "All-or-Nothing", "Leaky"), nAgeClasses)

pdf("example_results.pdf", width=11, height = 5)
par(mfrow = c(1, 2))
barplot(deaths[,1:4], ylim = c(0, max(deaths)), beside = T, legend = legendnames, args.legend = list(x=17, y = 0.55, cex = 0.85, bty = "n"), main = "Permanent Immunity", xlab = "Vaccination Scenario", ylab = "Population percent mortality", col=cols)
barplot(deaths[,5:8], ylim = c(0, max(deaths)), beside = T, main = "Waning Immunity", xlab = "Vaccination Scenario", ylab = "Population percent mortality", col=cols)
dev.off()

### end
