hypothesisTest = function(alpha = 0.05, testStatistic, observedOutcomes,
                          covariates, treatmentAllocation,
                          numIters = 500,  nboot = 500)
{
  
  referenceDistribution = generateReferenceDistribution(numIters = numIters,
                                                        nboot = nboot,
                                                        testStatistic = testStatistic,
                                                        observedOutcomes = observedOutcomes,
                                                        covariates = covariates,
                                                        treatmentAllocation = treatmentAllocation)
  
  # The Gaussian prepivoted test statistic
  observedTestStatistic = GaussianPrepivTestStat(testStatistic = testStatistic,
                                                 observedOutcomes = observedOutcomes,
                                                 covariates = covariates,
                                                 treatmentAllocation = treatmentAllocation,
                                                 nboot = nboot)
  
  # p-Value of the test & rejection of the null
  pValue = mean(referenceDistribution > observedTestStatistic)
  reject = (pValue <= alpha)
  
  return(list(rejectNull = reject, pValue = pValue))
}

GaussianPrepivTestStat = function(testStatistic, observedOutcomes, covariates, treatmentAllocation, nboot)
{
  N = length(treatmentAllocation) # Number of units
  nt = sum(treatmentAllocation) # Number of treated units
  nc = N - nt # Number of control units
  
  sortedTreat = c(rep(TRUE, nt), rep(FALSE, nc))
  
  Ytreated = observedOutcomes[treatmentAllocation, ,drop = FALSE] # Treated outcomes
  Xtreated = covariates[treatmentAllocation, ,drop = FALSE] # Treated covariates
  
  Ycontrol = observedOutcomes[!treatmentAllocation, ,drop = FALSE] # Control outcomes
  Xcontrol = covariates[!treatmentAllocation, ,drop = FALSE] # Control covariates
  
  # The test statistic without prepivoting
  unprepivStat = testStatistic(diffInMeans = colMeans(as.matrix(Ytreated)) - colMeans(as.matrix(Ycontrol)),
                               observedOutcomes = observedOutcomes, covariates = covariates,
                               treatmentAllocation = treatmentAllocation) 
  
  # Neyman's classical variance estimator for the difference in means
  varEst = (var(Ytreated)/nt) + (var(Ycontrol)/nc) 
  
  # Gaussian prepivot by using a parametric bootstrap on the difference in means
  bootstrapDistribution = numeric(nboot)
  for(b in 1:nboot)
  {
    gaussianDraw = mvtnorm::rmvnorm(n = 1, sigma = as.matrix(varEst))
    bootstrapDistribution[b] = testStatistic(diffInMeans = gaussianDraw,
                                             observedOutcomes = Yobs, covariates = covariates,
                                             treatmentAllocation = treatmentAllocation)
  }
  prepivResult = mean(bootstrapDistribution <= unprepivStat) # The prepivoted statistic
  
  return(prepivResult)
}

generateReferenceDistribution = function(numIters = 500, testStatistic, observedOutcomes,
                                         covariates, treatmentAllocation, nboot = 500)
{
  referenceDistribution = numeric(numIters)
  for(i in 1:numIters)
  {
    pseudoTreatmentAllocation = sample(treatmentAllocation, replace = FALSE)
    referenceDistribution[i] = GaussianPrepivTestStat(testStatistic = testStatistic, observedOutcomes = observedOutcomes,
                                                      covariates = covariates, treatmentAllocation = pseudoTreatmentAllocation,
                                                      nboot = nboot)
  }
  return(referenceDistribution)
}

DiM_GaussianPrepiv = function(diffInMeans, observedOutcomes, covariates, treatmentAllocation)
{
  return(sqrt(N)*abs(diffInMeans))
}

alpha = 0.05 # the level of the test


###########################################
#              Generate Data              #
###########################################
# Some algorithm hyperparameters
numIters = 100
nboot = 100

# First we make some synthetic data for our experiments
N = 1000 # the number of units
p = .3 # porportion of treated units
nt = floor(N*p) # number of treated
nc = N - nt # number of control

# Handy variables
treatind = c(rep(T, nt), rep(F, nc))

# Superpopulation level parameters
SdT = 2
SdC =1
SigmaX = cbind(c(1,.5, .2), c(.5,1,.3), c(.2,.3,1))
betaT = c(1,2,3)*.2
betaC = c(-1,-2, 4)*.8

# Taubar is the average additive treatment effect
taubar = 0

# Make covariates and potential outcomes (Gaussian model with noisey linear relations)
X = mvtnorm::rmvnorm(N, c(0,0,0), SigmaX)
Yttemp = X%*%betaT + rnorm(N, 0, SdT)
Yc = X%*%betaC + rnorm(N, 0, SdC) # The control POs
Yt = Yttemp - mean(Yttemp - Yc) + taubar # The treated POs (this centering forces Neyman's null to hold)

Z = sample(treatind, replace = FALSE) # Select a treatment allocation for a CRE design

# This is what the experimenter gets to see
Yobs = Z*Yt + (1 - Z)*Yc
Ytreat = Yt[Z]
Ycontrol = Yc[!Z]
Xtreat = X[Z,]
Xcontrol = X[!Z,]
Xtreat = data.frame(Xtreat)
Xcontrol = data.frame(Xcontrol)


###########################################
#     Start hypothesis testing.           #
###########################################
# Gaussian prepivot for finite population inference
DIM_Gaussian = hypothesisTest(alpha = alpha, 
                              testStatistic = DiM_GaussianPrepiv, 
                              observedOutcomes = Yobs, covariates = data.frame(X),
                              treatmentAllocation = Z, 
                              numIters = numIters, nboot = nboot)


if(DIM_Gaussian$rejectNull)
{
  cat('Rejected the null.\n')
}else
{
  cat('Failed to reject the null.\n')
}
