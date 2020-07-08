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

#pdf("4panelresults.pdf", width=11.36, height = 7.5)	# uncomment this to generate the .pdf file

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
      dX[i + 4*nAgeClasses] = gamma*(Ic[i] + Isc[i] + Icv[i] + Iscv[i])
            
      #equation for dM.dt:
      dX[i + 5*nAgeClasses] = gamma*(rhoc[i]*(Ic[i]+Icv[i]) + rhosc[i]*(Isc[i]+Iscv[i])) 
      
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
rhoc = 0.05; rhosc = 0.001

minT = 0; maxT = 365; dt = 1

t = seq(from=minT,to=maxT,by=dt)  # creating the vector of time for output

tau = 0 # delay before vaccination starts; prior to this, lambda = 0

percentImmune = 0.2	# how many are already immune?

# the relative values of E0, Ic0, and Isc0 matters; don't change them without checking the equations in the supplement
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

reduction = c(0.2, 0.8)	# reduce nu, beta by a this percent
factor = 1 - reduction

# storage:

opts = c("total", "direct", "indirect")
storeDeaths = array(dim = c(length(v), length(opts), 3))	# array of vaccination rates, total/direct/indirect deaths, and 
										# scenarios: no vaccine, vaccine 1, and vaccine 2


###########################################################################		open the loop over vaccination rates:

for(c in 1:length(v)){

vaccineRate = v[c] 

##### should clean this up before publication, just for readability

### 0: lambda = 0
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

##### for the dynamic death plot:

dailyDeath = array(dim = c(3, end-1))
dailyDeath[1,] = c(total0[2:end,7]-total0[1:(end-1),7])
dailyDeath[2,] = c(total1[2:end,7]-total1[1:(end-1),7])
dailyDeath[3,] = c(total2[2:end,7]-total2[1:(end-1),7])

dailyDeath = dailyDeath/normalize

##### for the big array of deaths under different scenarios

deaths = c(direct0[end,7], indirect0[end,7], total0[end,7], direct1[end,7], indirect1[end,7], total1[end,7], direct2[end,7], indirect2[end,7], total2[end,7])
deathArray = matrix(deaths, nrow = 3, ncol = 3, byrow=T)

storeDeaths[c,,] = deathArray	

if(vaccineRate == vplot){		# here we plot C and D for a given vaccination rate:
	##### set up for plots:

	legendnames = c("no vaccination", "vaccine 1", "vaccine 2")

	## B: 
	plot(dailyDeath[1,], xlab = "day", ylab = "deaths per day per 100k", xlim = c(0,250), cex.axis = 1.5, cex.lab = 1.5, type = "l", col = cols[1], lwd = 3)
	lines(dailyDeath[2,], type = "l", col = cols[2], lwd = 3)
	lines(dailyDeath[3,], type = "l", col = cols[3], lwd = 3)
	legend("topright", legend = legendnames, col = cols, lwd = 2, cex = 1.5, bty = "n")

	text(175, 15, paste("daily per capita\n vaccination rate of", lambda), cex = 1.2)	 
	put.fig.letter(label="B", location="topleft", cex = 2.5, font=2)

	## C:
	# plot total deaths per 100k: direct effects only and total effects:
	colors = c(cols[1], cols[1], cols[2], cols[2], cols[3], cols[3])

	plotArray = rbind(deathArray[,1], deathArray[,3])/normalize	# an array with the total deaths (row 1)
										# and then deaths if you have direct effects only (row 2)
	newRow = plotArray[1,] - plotArray[2,]
	plotArray = rbind(plotArray,newRow)				# now add a row for the difference between these
	plotArray = plotArray[2:3,]

	# I'm sorry about this but making a stacked barplot in R with different colors for each bar is way more complicated than it should be
	newArray = array(dim = c(6, 3))
	zeros = rbind(0, 0)
	newArray[1:2,] = cbind(plotArray[1:2,1], zeros, zeros)
	newArray[3:4,] = cbind(zeros, plotArray[1:2,2], zeros)
	newArray[5:6,] = cbind(zeros, zeros, plotArray[1:2,3])

	# and finally, the plot:

	grays = c("black", "black")	# use this for the legend
	barplot(newArray, ylim = c(0, 1.14*max(newArray)), legend = F, beside = F, cex.names = 1.5,  names.arg = legendnames, main = "",
		ylab = "total deaths per 100k", col=colors, cex.axis = 1.5, cex.lab = 1.5, density = c(1000, 50))
	axis(side = 1, line = 2, at = c(1.9), paste("daily per capita vaccination rate of", lambda), cex.axis = 1.5, tick = F)
	legend(1.75, 1555, legend = c("total effect", "direct effect alone"), cex = 1.5, bty = "n", fill = grays, density =c(1000, 50))
	put.fig.letter(label="C", location="topleft", cex = 2.5, font=2)

	}	# close the 'if' statement

}# close the loop over vaccination rates


## D:
plot(v, storeDeaths[,1,3]/normalize, col = cols[1], ylim = c(0, max(storeDeaths/normalize)), xlab = "daily per capita vaccination rate", 
	ylab = "total deaths per 100k", cex.axis = 1.5, cex.lab = 1.5, type = "l", lwd = 3)
lines(v, storeDeaths[,2,3]/normalize, col = cols[2], type = "l", lwd = 3)
lines(v,storeDeaths[,3,3]/normalize, col = cols[3], type = "l", lwd = 3)
put.fig.letter(label="D", location="topleft", cex = 2.5, font=2)

#dev.off()

### note that the sum of indirect and direct effects is not equal to the total effect, and this seems to be more pronounced at higher vaccination rates:

averted = c(deathArray[1,3] - deathArray[2,3], deathArray[1,3] - deathArray[3,3])
directAvert = c(deathArray[1,3] - deathArray[2,1], deathArray[1,3] - deathArray[3,1])
indirectAvert = c(deathArray[1,3] - deathArray[2,2], deathArray[1,3] - deathArray[3,2])
totalAvert = c(directAvert[1] + indirectAvert[1], directAvert[2] + indirectAvert[2])
diff = totalAvert - averted	# raw numbers for the difference between direct+indirect deaths averted, minus total deaths actually averted

##### end
