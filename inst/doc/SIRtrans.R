## ----setup, include=FALSE, cache=FALSE--------------------------------------------------
library(knitr)
# set global chunk options
opts_chunk$set(fig.path='figure/minimal-', fig.align='center', fig.show='hold', cache=TRUE)
options(formatR.arrow=TRUE,width=90,tidy=FALSE)

## ----eval=TRUE, cache=TRUE--------------------------------------------------------------
#Choose initial S and I population here
S <- 140
I <- 10
beta = .5/(S+I)
gamma = .1
N <- 1000 #number of MC realizations: increase N in practice
t.end <- .5 #time interval length

#monte carlo estimate and standard error
tpm.MC <- getTrans.MC(N, t.end, S, I, beta, gamma)
sd.MC <- sqrt( (tpm.MC)*(1 - tpm.MC)/N )

## ---------------------------------------------------------------------------------------
# now, calculate probabilities using generating functions
gridLength = 256
system.time( tpm1 <- getTransProbsODE(t.end, gridLength, 
                                      beta, gamma, S, I)[1:(S+I),1:(S+I)] )
system.time( tpm2 <- getTransProbsClosed(t.end, gridLength, 
                                         beta, gamma, S, I)[1:(S+I),1:(S+I)]  )


## ---------------------------------------------------------------------------------------
# total errors:
sum(abs(tpm1-tpm2)) # compare the two methods
sum(abs(tpm1-tpm.MC)) #compare to true model Monte Carlo probabilities

## ----marginalizedSIRplots---------------------------------------------------------------
#marginalize over susceptibles: c
infectiveProbs.MC <- colSums(tpm.MC)
infectiveSD.MC <- sqrt( (infectiveProbs.MC)*(1 - infectiveProbs.MC)/(N) )
#check that this SD is correct...
lower <- infectiveProbs.MC - infectiveSD.MC*1.96
upper <- infectiveProbs.MC + infectiveSD.MC*1.96
infectiveProbs.FFT <-round(colSums(tpm2),4)

require(plotrix)
plot(seq(S+I), xlim = c(0,50), infectiveProbs.FFT, pch = 3, 
     col = 'purple', main = "Probabilities of ending with x infectives")
plotCI(seq(S+I), xlim = c(0,50), infectiveProbs.MC, pch = 16, 
       col = 4, ui = upper, li = lower, add=TRUE)

#marginalize other way:
suscepProbs.MC <- rowSums(tpm.MC)
suscepSD.MC <- sqrt( (suscepProbs.MC)*(1 - suscepProbs.MC)/(N) )
lower <- suscepProbs.MC - suscepSD.MC*1.96
upper <- suscepProbs.MC + suscepSD.MC*1.96
suscepProbs.FFT <-round(rowSums(tpm2),4)

plot(seq(S+I), xlim = c(110,160), suscepProbs.FFT, pch = 3, 
     col = 'purple', main = "Probabilities of ending with x susceptibles")
plotCI(seq(S+I), xlim = c(110,160), suscepProbs.MC, pch = 16, 
       col = 4, ui = upper, li = lower, add=TRUE)

## ----multiBDsetup-----------------------------------------------------------------------
library(MultiBD)
tList  <- c(.1, .2, .25, .3 ,.35, .4, .5, .6, .7, .8, .9, 1)
gridLength = 128
a0 = 110 # S_0
b0 = 15 # I_0
A = 0
B = gridLength - 1
alpha = 3.2 #3.2 #this is death rate
beta = .025 #.019 #this is transition or infection rates
nSim = 4000 #number of MC simulations

brates1=function(a,b){0}
drates1=function(a,b){0}
brates2=function(a,b){0}
drates2=function(a,b){alpha*b}
trans=function(a,b){beta*a*b}

## ----makeTPMs, cache=TRUE---------------------------------------------------------------
#indexed by time, type of computation, and dimensions of the tpm
tpmArray <- array(NA, dim= c(length(tList),3, 52, 25 )) #store a subset of transition probabilities 

for(i in 1:length(tList)){
  t.end <- tList[i]
  system.time( tpm.Closed <- getTransProbsClosed(t.end, gridLength, 
                                           beta, alpha, a0, b0) ) 
  tpm1 = tpm.Closed[1:(a0+1),] #using 2-type branching approximation
  
  #using continued fractions via MultiBD
  system.time( tpm2 <- dbd_prob(t.end, a0, b0, drates1, brates2, drates2, trans,
                                           a=A, B))#, computeMode=2))
  #MC simulation "ground truth"
  tpm.MC <- getTrans.MC(nSim, t.end, a0, b0, beta, alpha)
  tpm3 <- tpm.MC[1:(a0+1), ]

  #store subset of matrices containing about 99 percent of the mass:
  tpmArray[i,1,,] <- tpm1[60:(a0+1),1:25]
  tpmArray[i,2,,] <- tpm2[60:(a0+1),1:25]
  tpmArray[i,3,,] <- tpm3[60:(a0+1),1:25]
}

## ---------------------------------------------------------------------------------------
#for example, look at the ones with t.end = .5
small1 <- tpmArray[5,1,,]
small2 <- tpmArray[5,2,,]
small3 <- tpmArray[5,3,,]

#they comprise most of transition probability mass:
sum(small1); sum(small2); sum(small3)

# mean errors
mean(abs(small1- small3 ) ) #2-type vs MC
mean(abs(small2 - small3) ) #Continued Frac vs MC

# scaled heatmap images to compare tpm visually
par(mfrow=(c(3,1)))
image(small1, main = "Two-type branching approximation")
image(small2, main = "Continued Fraction expansion")
image(small3, main = "Monte Carlo estimates")

## ---------------------------------------------------------------------------------------
par(mfrow=(c(3,1)))
image(tpmArray[12,1,,], main = "Two-type branching approximation")
image(tpmArray[12,2,,], main = "Continued Fraction expansion")
image(tpmArray[12,3,,], main = "Monte Carlo estimates")

## ----transitionCompare------------------------------------------------------------------
library(plotrix)
inds <- t(which(tpmArray[7,2,,] >= sort(tpmArray[7,2,,], decreasing=T)[16], arr.ind=TRUE))
#ind1 <- sample(52,25, replace=T); ind2 <- sample(25,25,replace=T)
par(mfrow = c(4,4), oma = c(5,4,0,0) + 0.1, mar = c(0,0,1,1) + 0.1)
for(i in 1:16){
plot(tList, tpmArray[,2,inds[1,i], inds[2,i] ], pch = 17, col = 'red',
ylim = c(0,max(tpmArray[,,inds[1,i], inds[2,i]])),
yaxt = 'n', xlab = "dt")
MCp <- tpmArray[,3,inds[1,i], inds[2,i] ] #MC prob
plotCI(tList, MCp, pch = 4, col = 'green', ui=MCp+1.96*sqrt(MCp*(1-MCp)/nSim),
li=MCp-1.96*sqrt(MCp*(1-MCp)/nSim), add = TRUE)
points(tList, tpmArray[,1,inds[1,i], inds[2,i] ], col='purple', pch = 16)
}

## ----largerPopulation,cache=TRUE--------------------------------------------------------
tList  <- c( .5, 1)
gridLength = 256
a0 = 235 # S_0
b0 = 15 # I_0
A = 0
B = gridLength - 1
alpha = 3.2 #3.2 #this is death rate
beta = .025 #.019 #this is transition or infection rates
nSim <- 10000
tpmArray <- array(NA, dim= c(length(tList),2, (a0+1), 240 )) #store a subset of transition probabilities 

for(i in 1:length(tList)){
  t.end <- tList[i]
  system.time( tpm2 <- dbd_prob(t.end, a0, b0, drates1, brates2, drates2, trans,
                                           a=A, B))#, computeMode=2))
  #MC simulation "ground truth"
  tpm.MC <- getTrans.MC(nSim, t.end, a0, b0, beta, alpha)
  tpm3 <- tpm.MC[1:(a0+1), ]

  #store subset of matrices containing about 99 percent of the mass:
  tpmArray[i,1,,] <- tpm2[1:(a0+1),1:240]
  tpmArray[i,2,,] <- tpm3[1:(a0+1),1:240]
}

## ---------------------------------------------------------------------------------------
par(mfrow=(c(2,1)))
image(tpmArray[1,1,,1:60], main = "Continued Fraction approximation, t=.5")
image(tpmArray[1,2,,1:60], main = "Monte Carlo estimates")

par(mfrow=(c(2,1)))
image(tpmArray[2,1,,1:60], main = "Continued Fraction approximation, t=1")
image(tpmArray[2,2,,1:60], main = "Monte Carlo estimates")

