### SEIR model of vaccination scenarios

#setwd("C:/Users/mgalla9/Dropbox/Emory/Coronavirus/code")
setwd("C:/Users/molly/Dropbox/Emory/Coronavirus/code")

source("put_fig_letter.r")

#install.packages('deSolve','viridis','jpeg')
require(deSolve)
require(viridis)
require(jpeg)


#########################

fixedN = 1e7	# total population size

# loop over vaccination rates (per capita per day)
v = seq(0, .03, .001)
vplot = v[11] # for the snapshot plots - B and C are plotted for a given value of v

###### set up for the plots: 4 panels

cols = viridis(n=3, begin = 0, end = .7)	# 3 colors, one for each vaccine strategy (none, vaccine 1, vaccine 2)

pdf("4panelresults.pdf", width=11.36, height = 7.5)	# uncomment this to generate the .pdf file

par(mfrow = c(2, 2))
par(mai = c(0, 0, 0, 0))
plot(0, 0, xlim = c(0, 20), ylim = c(0, 20), cex = 0, xlab = "", ylab = "", xaxt = "n", yaxt = "n", bty = "n")	# empty plot to hold the top left quadrant for a schematic

schematic = readJPEG("schematic.jpg")

#now open a plot window with coordinates
#specify the position of the image through bottom-left and top-right coords
rasterImage(schematic,2,2,18, 20)
put.fig.letter(label="A", location="topleft", cex = 2.5, font=2)
par(mai = c(0.9, 1, 0.5, 0.5))


############################		onward to the actual modeling:

### scenarios: no vaccine; vaccine 1 that reduces nu 80%, beta 20%; vaccine 2 reduces nu 20%, beta 80%

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
      dX[i] = - beta[i]*(S[i]/N)*(Ic[i] + Isc[i]) - betav[i]*(S[i]/N)*(Icv[i] + Iscv[i]) - vaccination_rate[i]*S[i]
     
      #equation for dE.dt: 	
      dX[i+nAgeClasses] = beta[i]*(S[i]/N)*(Ic[i] + Isc[i]) +  betav[i]*(S[i]/N)*(Icv[i] + Iscv[i]) - sigma*E[i]
      
      #equation for dIc.dt:
      dX[i + 2*nAgeClasses] = nu[i]*sigma*E[i] - gamma*Ic[i]
      
      #equation for dIsc.dt:
      dX[i + 3*nAgeClasses] = (1 - nu[i])*sigma*E[i] - gamma*Isc[i]
      
      #equation for dR.dt:
      #dX[i + 4*nAgeClasses] = gamma*(Ic[i] + Isc[i] + Icv[i] + Iscv[i])
      dX[i + 4*nAgeClasses] = gamma*((1 - rhoc[i])*(Ic[i] +  Icv[i]) + Isc[i] + Iscv[i])

      #equation for dM.dt:
      #dX[i + 5*nAgeClasses] = gamma*(rhoc[i]*(Ic[i]+Icv[i]) + rhosc[i]*(Isc[i]+Iscv[i])) 
      dX[i + 5*nAgeClasses] = gamma*rhoc[i]*(Ic[i]+Icv[i]) 

      ##### now the equations for the vaccinated individuals:
      # equation for dSv.dt:
      dX[i + 6*nAgeClasses] = vaccination_rate[i]*S[i] - beta[i]*(Sv[i]/N)*(Ic[i] + Isc[i]) - betav[i]*(Sv[i]/N)*(Icv[i] + Iscv[i])

      # equation for dEv.dt:
      dX[i + 7*nAgeClasses] = beta[i]*(Sv[i]/N)*(Ic[i] + Isc[i]) + betav[i]*(Sv[i]/N)*(Icv[i] + Iscv[i]) - sigma*Ev[i]
      
      #equation for dIcv.dt:
      dX[i + 8*nAgeClasses] = nuv[i]*sigma*Ev[i] - gamma*Icv[i]
      
      #equation for dIscv.dt:
      dX[i + 9*nAgeClasses] = (1 - nuv[i])*sigma*Ev[i] - gamma*Iscv[i]
      
    } # closes the loop over age classes	
    
    return(list(dX)) 	
  }) # stops the parameter list		
}	# closes the function


##### Universal values:

nAgeClasses = 1
sigma = 1/4; gamma = 1/5; beta = 0.6; nu = 0.35
rhoc = 0.05; rhosc = 0 
minT = 0; maxT = 365; dt = 1

t = seq(from=minT,to=maxT,by=dt)  # creating the vector of time for output

tau = 0 # delay before vaccination starts; prior to this, lambda = 0

percentImmune = 0.2	# how many are already immune?

# the relative value of E0, Ic0, and Isc0 matters; don't change them without checking the equations in the supplement
E0 = 200 
Ic0 = 140
Isc0 = 260
R0 =  percentImmune*fixedN
M0 = 0
S0 = fixedN - E0 - Ic0 - Isc0 - R0 - M0

initialConditions = c(S0, 
                      E0,
                      Ic0, Isc0, 
                      R0, M0,
                      0, 0, 0, 0)	# ordered by S[1]...S[nAgeClasses], then E[1]...E[nAgeClasses],etc

maxN = sum(initialConditions)
maxN == fixedN	# check to be sure this is true


###### Values for different scenarios:

reduction = c(0.2, 0.8)	# reduce nu, beta by this percent
factor = 1 - reduction	


# storage:

opts = c("total", "direct", "indirect")
storeDeaths = array(dim = c(length(v), length(opts), 3))	# array of vaccination rates, total/direct/indirect deaths, and 
										# scenarios: no vaccine, vaccine 1, and vaccine 2

percentAverted = array(dim = c(length(v), 4)) # array of % deaths averted by V1 direct, V1 total, V2 direct, and V2 total

###########################################################################		open the loop over vaccination rates:

for(c in 1:length(v)){

vaccineRate = v[c] 

########### simulate the different scenarios:


### 0: no vaccination; lambda = 0
lambda = 0; betav = beta; nuv = nu
p = list(N = fixedN, beta, betav, sigma, nu, nuv, gamma, rhoc, rhosc, lambda, tau)		 # creating the LIST of parameter values
total0 = ode(y=initialConditions, times=t, func=SEIR, parms=p, method = 'ode45') # run the ode solver

# measure direct effects of vaccination only: 
betav = beta; 
p = list(N = fixedN, beta, betav, sigma, nu, nuv, gamma, rhoc, rhosc, lambda, tau)	
direct0 = ode(y=initialConditions, times=t, func=SEIR, parms=p, method = 'ode45') # run the ode solver

# measure indirect effects of vaccination only:
nuv = nu; 
p = list(N = fixedN, beta, betav, sigma, nu, nuv, gamma, rhoc, rhosc, lambda, tau)
indirect0 = ode(y=initialConditions, times=t, func=SEIR, parms=p, method = 'ode45') # run the ode solver


### 1: vaccination reduces clinical infections by a lot, transmission by a little
lambda = vaccineRate; betav = factor[1]*beta; nuv = factor[2]*nu
p = list(N = fixedN, beta, betav, sigma, nu, nuv, gamma, rhoc, rhosc, lambda, tau)		 # creating the LIST of parameter values
total1 = ode(y=initialConditions, times=t, func=SEIR, parms=p, method = 'ode45') # run the ode solver

# measure direct effects of vaccination only: 
betav = beta; 
p = list(N = fixedN, beta, betav, sigma, nu, nuv, gamma, rhoc, rhosc, lambda, tau)	
direct1 = ode(y=initialConditions, times=t, func=SEIR, parms=p, method = 'ode45') # run the ode solver

# measure indirect effects of vaccination only:
betav = factor[1]*beta; nuv = nu
p = list(N = fixedN, beta, betav, sigma, nu, nuv, gamma, rhoc, rhosc, lambda, tau)	
indirect1 = ode(y=initialConditions, times=t, func=SEIR, parms=p, method = 'ode45') # run the ode solver


### 2: vaccination reduces transmission by a lot, clinical infections by a little
lambda = vaccineRate; betav = factor[2]*beta; nuv = factor[1]*nu
p = list(N = fixedN, beta, betav, sigma, nu, nuv, gamma, rhoc, rhosc, lambda, tau)		 # creating the LIST of parameter values
total2 = ode(y=initialConditions, times=t, func=SEIR, parms=p, method = 'ode45') # run the ode solver

# measure direct effects of vaccination only: 
betav = beta; 
p = list(N = fixedN, beta, betav, sigma, nu, nuv, gamma, rhoc, rhosc, lambda, tau)	
direct2 = ode(y=initialConditions, times=t, func=SEIR, parms=p, method = 'ode45') # run the ode solver

# measure indirect effects of vaccination only:
betav = factor[2]*beta; nuv = nu
p = list(N = fixedN, beta, betav, sigma, nu, nuv, gamma, rhoc, rhosc, lambda, tau)	
indirect2 = ode(y=initialConditions, times=t, func=SEIR, parms=p, method = 'ode45') # run the ode solver


end = length(t)
normalize = fixedN/100000	# normalize values to "per 100k individuals"


##### for the big array of deaths under different scenarios

deaths = c(direct0[end,7], indirect0[end,7], total0[end,7], direct1[end,7], indirect1[end,7], total1[end,7], direct2[end,7], indirect2[end,7], total2[end,7])
deathArray = matrix(deaths, nrow = 3, ncol = 3, byrow=T)

storeDeaths[c,,] = deathArray	

#deathArray: rows = no vaccine, v1, v2; columns = 'direct'(nuv reduced), 'indirect'(betav reduced), & total
# percent deaths averted:
percentAverted[c,1] = 100*round((deathArray[1,1] - deathArray[2,1])/deathArray[1,1],4)
percentAverted[c,2] = 100*round((deathArray[1,1] - deathArray[2,3])/deathArray[1,1],4)
percentAverted[c,3] = 100*round((deathArray[1,1] - deathArray[3,1])/deathArray[1,1],4)
percentAverted[c,4] = 100*round((deathArray[1,1] - deathArray[3,3])/deathArray[1,1],4)

if(vaccineRate == vplot){		# here we plot C and D for a given vaccination rate:

	##### for the dynamic death plots:

	dailyDeath = array(dim = c(3, end-1))
	dailyDeath[1,] = c(total0[2:end,7]-total0[1:(end-1),7])
	dailyDeath[2,] = c(total1[2:end,7]-total1[1:(end-1),7])
	dailyDeath[3,] = c(total2[2:end,7]-total2[1:(end-1),7])
	dailyDeath = dailyDeath/normalize

	cumDeath = array(dim = c(3, end))
	cumDeath[1,] = total0[,7]
	cumDeath[2,] = total1[,7]
	cumDeath[3,] = total2[,7]
	cumDeath = cumDeath/normalize

	##### set up for plots:

	legendnames = c("no vaccine", "vaccine 1", "vaccine 2")

	## B: 
	par(mai = c(0.9, 1, 0.25, 1))
	types = c("dotted", "solid")

	plot(dailyDeath[1,], main ="", cex.main = 1.4, xlab = "day", ylab = "daily deaths per 100k", xlim = c(0,300), cex.axis = 1.5, cex.lab = 1.5, type = "l", lty = types[1], col = cols[1], lwd = 4)
	lines(dailyDeath[2,], type = "l", lty = types[1], col = cols[2], lwd = 4)
	lines(dailyDeath[3,], type = "l", lty = types[1], col = cols[3], lwd = 4)
	legendCTitles = c("daily deaths", "total deaths")
	legend(100, 32, legend = legendnames, title = legendCTitles[1], col = cols, lwd = 2, lty = types[1], cex = 1, y.intersp = .75, x.intersp = .5, text.width = 50, bty = "n")
	legend(200, 32, legend = legendnames, title = legendCTitles[2], col = cols, lwd = 2, lty = types[2], cex = 1, y.intersp = .75, x.intersp = .5, text.width = 50, bty = "n")

	par(new=T)	# add a second y-axis to the plot:
	plot(cumDeath[1,], xlim = c(0,300), xlab = "", ylab = "", xaxt = "null", yaxt = "null", type = "l", lty = types[2], lwd = 3, col = cols[1])
	lines(cumDeath[2,],type = "l", lty = types[2], lwd = 3, col = cols[2])
	lines(cumDeath[3,],type = "l", lty = types[2], lwd = 3, col = cols[3])

	axis(4, cex.axis = 1.5)
	mtext(side = 4, line = 3, "total deaths per 100k", cex = 1.3)
	#text(250, 300, , cex = 1.)	 
	put.fig.letter(label="B", location="topleft", cex = 2.5, font=2)

	## C:
	
	par(mai = c(0.9, 1, .25, 0.5))
	colors = c(cols[2], cols[2], cols[3], cols[3])
	namesC = c("direct effect", "total effect", "direct effect", "total effect")
	barplot(percentAverted[c,], ylim = c(0, 1.05*(max(percentAverted[c,]))), col = colors, density = c(50, 1000, 50, 1000),
		ylab = "percent deaths averted", space = c(0.05, 0.05, 0.1, 0.05), cex.main = 1.4,
			cex.axis = 1.5, cex.lab = 1.5,cex.names = 1.2, names.arg = namesC)
	#mtext(side = 1, line = -13, at  = 1.1, cex = 1.5, "vaccine 1") 
	#mtext(side = 1, line = -13, at  = 3.2, cex = 1.5, "vaccine 2") 
	mtext(side = 1, line = 2.5, at  = 1.05, cex = 1.5, "vaccine 1") 
	mtext(side = 1, line = 2.5, at  = 3.2, cex = 1.5, "vaccine 2") 

	put.fig.letter(label="C", location="topleft", cex = 2.5, font=2)

	
	par(mai = c(0.9, 1, 0.25, 0.5))

	}	# close the 'if' statement

}# close the loop over vaccination rates

## D:
plot(v, percentAverted[,2], col = cols[2], ylim = c(0, 100), xlab = "daily per capita vaccination rate", 
	ylab = "percent deaths averted", cex.axis = 1.5, cex.lab = 1.5, type = "l", lwd = 3)
lines(v, percentAverted[,4], col = cols[3], type = "l", lwd = 3)
legend("bottomright", legend = legendnames[2:3], col = cols[2:3], lwd = 3, lty = types[2], cex = 1.5, bty = "n")

put.fig.letter(label="D", location="topleft", cex = 2.5, font=2)

dev.off()

##### end
