library(foreign) # This package is used to read in the data from the ALO paper (stored in a .dta file format)

set.seed(1) # For reproducible results
# Some hyperparamters
numPerms = 5000
nMC = 10000


################################################################################
# Pull in the ALO data set
################################################################################
# To download the data go to "https://www.aeaweb.org/articles?id=10.1257/app.1.1.136"
ALO_RawData <- read.dta("<SPECIFY YOUR FILEPATH HERE>/STARdatapost/STAR_public_use.dta")

ALO_males = ALO_RawData[!ALO_RawData$female,] # Subset to only males

experimentalData = ALO_males[, c("GPA_year1", "GPA_year2" ,"sfsp", "ssp", "gpa0")]
experimentalData = experimentalData[complete.cases(experimentalData), ]

# Examine only those who were offered services or services & support
experimentalData = experimentalData[experimentalData$sfsp + experimentalData$ssp == 1, ]

treatment = as.logical(experimentalData$sfsp) # Treatment is "financial incentive offered" given that services were also offered 
N = length(treatment)
n1 = sum(treatment)
n0 = N - n1

observedOutcomes = data.frame(experimentalData[, c('GPA_year1', 'GPA_year2')])
covariates = data.frame(experimentalData[,'gpa0', drop = FALSE])

################################################################################
# Define Test Statistics
################################################################################

# the Euclidean norm of the difference in means (scaled at sqrt(N))
norm_function = function(diffInMeans, observedOutcomes, covariates, treatmentAllocation)
{
  return(sqrt(N)*sqrt(sum(diffInMeans^2))) # Euclidean norm of the difference in means
}

# the multivariate studentized statistic
multivarStudentized_function = function(diffInMeans, observedOutcomes, covariates, treatmentAllocation)
{
  # diffInMeans must be a col vector
  diffInMeans = matrix(diffInMeans, ncol = 1)
  
  # Studentizing with the Neyman covariance estimator
  N = length(treatmentAllocation)
  n1 = sum(treatmentAllocation)
  n0 = N - n1
  V_neyman = N*((cov(observedOutcomes[treatmentAllocation, ]) / n1) + (cov(observedOutcomes[!treatmentAllocation, ])/ n0))
  
  return(as.numeric(N*t(diffInMeans)%*%solve(V_neyman)%*%diffInMeans))
}

# the multivariate POOLED studentized statistic
multivarPooled_function = function(diffInMeans, observedOutcomes, covariates, treatmentAllocation)
{
  # diffInMeans must be a col vector
  diffInMeans = matrix(diffInMeans, ncol = 1)
  
  # Studentizing with the Neyman covariance estimator
  N = length(treatmentAllocation)
  n1 = sum(treatmentAllocation)
  n0 = N - n1
  V_pool = ((N / n0) + (N / n1))*
    (((n1 - 1)*cov(observedOutcomes[treatmentAllocation, ]) + 
        (n0 - 1)*cov(observedOutcomes[!treatmentAllocation, ])) / (n1 + n0 - 2))
  
  return(as.numeric(N*t(diffInMeans)%*%solve(V_pool)%*%diffInMeans))
}

# The maximum absolute t-statistic
maxT_function = function(diffInMeans, observedOutcomes, covariates, treatmentAllocation)
{
  N = length(treatmentAllocation)
  n1 = sum(treatmentAllocation)
  n0 = N - n1
  
  V_neyman = N*((cov(observedOutcomes[treatmentAllocation, ]) / n1) + (cov(observedOutcomes[!treatmentAllocation, ])/ n0))
  
  tStat = sqrt(N)*abs(diffInMeans)/diag(V_neyman)
  
  return(max(tStat))
}

################################################################################
# Compute the observed test statistics
################################################################################

Ytreated = observedOutcomes[treatment, ,drop = FALSE] # treated outcomes
Xtreated = covariates[treatment, ,drop = FALSE] # treated covariates

Ycontrol = observedOutcomes[!treatment, ,drop = FALSE] # control outcomes
Xcontrol = covariates[!treatment, ,drop = FALSE] # control covariates


# Compute the non-prepivoted test statistics (without any regression adjustment)
MaxT = maxT_function(diffInMeans = colMeans(as.matrix(Ytreated)) - colMeans(as.matrix(Ycontrol)),
                     observedOutcomes = observedOutcomes, 
                     covariates = covariates,
                     treatmentAllocation = treatment) # the test statistic without prepivoting

Euclidean = norm_function(diffInMeans = colMeans(as.matrix(Ytreated)) - colMeans(as.matrix(Ycontrol)),
                          observedOutcomes = observedOutcomes, 
                          covariates = covariates,
                          treatmentAllocation = treatment) # the test statistic without prepivoting

Studentized = multivarStudentized_function(diffInMeans = colMeans(as.matrix(Ytreated)) - colMeans(as.matrix(Ycontrol)),
                                           observedOutcomes = observedOutcomes, 
                                           covariates = covariates,
                                           treatmentAllocation = treatment) # the test statistic without prepivoting

Pooled = multivarPooled_function(diffInMeans = colMeans(as.matrix(Ytreated)) - colMeans(as.matrix(Ycontrol)),
                                 observedOutcomes = observedOutcomes, 
                                 covariates = covariates,
                                 treatmentAllocation = treatment) # the test statistic without prepivoting

# Compute the non-prepivoted test statistics (with regression adjustment on the high-school gpa)
lmT = lm(cbind(GPA_year1, GPA_year2)~gpa0, data = data.frame(experimentalData), subset = treatment)
lmC = lm(cbind(GPA_year1, GPA_year2)~gpa0, data = data.frame(experimentalData), subset = !treatment)

predsT = predict(lmT, covariates)
predsC = predict(lmC, covariates)

regAdjEstimator = as.matrix(colMeans(predsT - predsC))

regAdj_MaxT = maxT_function(diffInMeans = regAdjEstimator,
                            observedOutcomes = observedOutcomes, 
                            covariates = covariates,
                            treatmentAllocation = treatment) # the test statistic without prepivoting

regAdj_Euclidean = norm_function(diffInMeans = regAdjEstimator,
                                 observedOutcomes = observedOutcomes, 
                                 covariates = covariates,
                                 treatmentAllocation = treatment) # the test statistic without prepivoting

regAdj_Studentized = multivarStudentized_function(diffInMeans = regAdjEstimator,
                                                  observedOutcomes = observedOutcomes, 
                                                  covariates = covariates,
                                                  treatmentAllocation = treatment) # the test statistic without prepivoting

regAdj_Pooled = multivarPooled_function(diffInMeans = regAdjEstimator,
                                        observedOutcomes = observedOutcomes, 
                                        covariates = covariates,
                                        treatmentAllocation = treatment) # the test statistic without prepivoting

# Compute the prepivoted test statistics (using the Neyman covariance estimator for prepivoting)
varEst = (cov(Ytreated)/n1) + (cov(Ycontrol)/n0) # Neyman's classical variance estimator for the difference in means

# Gaussian prepivot by using Monte-Carlo simulation on the difference in means
gaussianDraw = mvtnorm::rmvnorm(n = nMC, sigma = as.matrix(varEst))
MaxT_MC = apply(gaussianDraw, 1, maxT_function, observedOutcomes = observedOutcomes, covariates = covariates, treatmentAllocation = treatment)
Euclidean_MC = apply(gaussianDraw, 1, norm_function, observedOutcomes = observedOutcomes, covariates = covariates, treatmentAllocation = treatment)
Studentized_MC = apply(gaussianDraw, 1, multivarStudentized_function, observedOutcomes = observedOutcomes, covariates = covariates, treatmentAllocation = treatment)
Pooled_MC = apply(gaussianDraw, 1, multivarPooled_function, observedOutcomes = observedOutcomes, covariates = covariates, treatmentAllocation = treatment)

# Prepivoted statistics
MaxT_Prepiv_obs = mean(MaxT_MC <= MaxT) 
Euclidean_Prepiv_obs = mean(Euclidean_MC <= Euclidean) 
Studentized_Prepiv_obs = pchisq(Studentized, df = 2) #
# Studentized_Prepiv_obs = mean(Studentized_MC <= Studentized) 
# plot(ecdf(pchisq(Studentized_MC, df = 2)))
# abline(a=0, b =1)
Pooled_Prepiv_obs = mean(Pooled_MC <= Pooled)

regAdj_MaxT_Prepiv_obs = mean(MaxT_MC <= regAdj_MaxT) 
regAdj_Euclidean_Prepiv_obs = mean(Euclidean_MC <= regAdj_Euclidean) 
regAdj_Studentized_Prepiv_obs = pchisq(regAdj_Studentized, df = 2) #
# regAdj_Studentized_Prepiv_obs =  mean(Studentized_MC <= regAdj_Studentized) 
regAdj_Pooled_Prepiv_obs = mean(Pooled_MC <= regAdj_Pooled)

################################################################################
# Enumerate the reference distribution by permuting treatment labels
################################################################################
temp = c(rep(TRUE, n1), rep(FALSE, n0))

MaxT_withoutPrepiv = numeric(numPerms)
Euclidean_withoutPrepiv = numeric(numPerms)
Studentized_withoutPrepiv = numeric(numPerms)
Pooled_withoutPrepiv = numeric(numPerms)

regAdj_MaxT_withoutPrepiv = numeric(numPerms)
regAdj_Euclidean_withoutPrepiv = numeric(numPerms)
regAdj_Studentized_withoutPrepiv = numeric(numPerms)
regAdj_Pooled_withoutPrepiv = numeric(numPerms)

MaxT_Prepiv = numeric(numPerms)
Euclidean_Prepiv = numeric(numPerms)
Studentized_Prepiv = numeric(numPerms)
Pooled_Prepiv = numeric(numPerms)

regAdj_MaxT_Prepiv = numeric(numPerms)
regAdj_Euclidean_Prepiv = numeric(numPerms)
regAdj_Studentized_Prepiv = numeric(numPerms)
regAdj_Pooled_Prepiv = numeric(numPerms)

for (i in 1:numPerms)
{
  W = sample(temp)
  Ytreated = observedOutcomes[W, ,drop = FALSE] # treated outcomes
  Xtreated = covariates[W, ,drop = FALSE] # treated covariates
  
  Ycontrol = observedOutcomes[!W, ,drop = FALSE] # control outcomes
  Xcontrol = covariates[!W, ,drop = FALSE] # control covariates
  
  
  # Compute the non-prepivoted test statistics (without any regression adjustment)
  MaxT_withoutPrepiv[i] = maxT_function(diffInMeans = colMeans(as.matrix(Ytreated)) - colMeans(as.matrix(Ycontrol)),
                                        observedOutcomes = observedOutcomes, 
                                        covariates = covariates,
                                        treatmentAllocation = W) # the test statistic without prepivoting
  
  Euclidean_withoutPrepiv[i] = norm_function(diffInMeans = colMeans(as.matrix(Ytreated)) - colMeans(as.matrix(Ycontrol)),
                                             observedOutcomes = observedOutcomes, 
                                             covariates = covariates,
                                             treatmentAllocation = W) # the test statistic without prepivoting
  
  Studentized_withoutPrepiv[i] = multivarStudentized_function(diffInMeans = colMeans(as.matrix(Ytreated)) - colMeans(as.matrix(Ycontrol)),
                                                              observedOutcomes = observedOutcomes, 
                                                              covariates = covariates,
                                                              treatmentAllocation = W) # the test statistic without prepivoting
  
  Pooled_withoutPrepiv[i] = multivarPooled_function(diffInMeans = colMeans(as.matrix(Ytreated)) - colMeans(as.matrix(Ycontrol)),
                                                    observedOutcomes = observedOutcomes, 
                                                    covariates = covariates,
                                                    treatmentAllocation = W) # the test statistic without prepivoting
  
  # Compute the non-prepivoted test statistics (with regression adjustment on the high-school gpa)
  lmT = lm(cbind(GPA_year1, GPA_year2)~gpa0, data = data.frame(experimentalData), subset = W)
  lmC = lm(cbind(GPA_year1, GPA_year2)~gpa0, data = data.frame(experimentalData), subset = !W)
  
  predsT = predict(lmT, covariates)
  predsC = predict(lmC, covariates)
  
  regAdjEstimator = as.matrix(colMeans(predsT - predsC))
  
  regAdj_MaxT_withoutPrepiv[i] = maxT_function(diffInMeans = regAdjEstimator,
                                               observedOutcomes = observedOutcomes, 
                                               covariates = covariates,
                                               treatmentAllocation = W) # the test statistic without prepivoting
  
  regAdj_Euclidean_withoutPrepiv[i] = norm_function(diffInMeans = regAdjEstimator,
                                                    observedOutcomes = observedOutcomes, 
                                                    covariates = covariates,
                                                    treatmentAllocation = W) # the test statistic without prepivoting
  
  regAdj_Studentized_withoutPrepiv[i] = multivarStudentized_function(diffInMeans = regAdjEstimator,
                                                                     observedOutcomes = observedOutcomes, 
                                                                     covariates = covariates,
                                                                     treatmentAllocation = W) # the test statistic without prepivoting
  
  regAdj_Pooled_withoutPrepiv[i] = multivarPooled_function(diffInMeans = regAdjEstimator,
                                                           observedOutcomes = observedOutcomes, 
                                                           covariates = covariates,
                                                           treatmentAllocation = W) # the test statistic without prepivoting
  
  # Compute the prepivoted test statistics (using the Neyman covariance estimator for prepivoting)
  varEst = (cov(Ytreated)/n1) + (cov(Ycontrol)/n0) # Neyman's classical variance estimator for the difference in means
  
  # Gaussian prepivot by using Monte-Carlo simulation on the difference in means
  gaussianDraw = mvtnorm::rmvnorm(n = nMC, sigma = as.matrix(varEst))
  MaxT_MC = apply(gaussianDraw, 1, maxT_function, observedOutcomes = observedOutcomes, covariates = covariates, treatmentAllocation = W)
  Euclidean_MC = apply(gaussianDraw, 1, norm_function, observedOutcomes = observedOutcomes, covariates = covariates, treatmentAllocation = W)
  Studentized_MC = apply(gaussianDraw, 1, multivarStudentized_function, observedOutcomes = observedOutcomes, covariates = covariates, treatmentAllocation = W)
  Pooled_MC = apply(gaussianDraw, 1, multivarPooled_function, observedOutcomes = observedOutcomes, covariates = covariates, treatmentAllocation = W)
  # Prepivoted statistics
  MaxT_Prepiv[i] = mean(MaxT_MC <= MaxT_withoutPrepiv[i]) 
  Euclidean_Prepiv[i] = mean(Euclidean_MC <= Euclidean_withoutPrepiv[i]) 
  Studentized_Prepiv[i] = pchisq(Studentized_withoutPrepiv[i], df = 2) #
  # Studentized_Prepiv[i] = mean(Studentized_MC <= Studentized_withoutPrepiv[i]) 
  # plot(ecdf(pchisq(Studentized_MC, df = 2)))
  # abline(a=0, b =1)
  Pooled_Prepiv[i] = mean(Pooled_MC <= Pooled_withoutPrepiv[i])
  
  regAdj_MaxT_Prepiv[i] = mean(MaxT_MC <= regAdj_MaxT_withoutPrepiv[i]) 
  regAdj_Euclidean_Prepiv[i] = mean(Euclidean_MC <= regAdj_Euclidean_withoutPrepiv[i])
  regAdj_Studentized_Prepiv[i] = pchisq(regAdj_Studentized_withoutPrepiv[i], df = 2)#
  # regAdj_Studentized_Prepiv[i] = mean(Studentized_MC <= regAdj_Studentized_withoutPrepiv[i]) 
  regAdj_Pooled_Prepiv[i] = mean(Pooled_MC <= regAdj_Pooled_withoutPrepiv[i])
}

################################################################################
# Compute p-Values
################################################################################
# non-prepivoted test stat p-values
npMaxT_noRegAdj = mean(MaxT_withoutPrepiv >= MaxT)
npEuclidean_noRegAdj = mean(Euclidean_withoutPrepiv >= Euclidean)
npStudentized_noRegAdj = mean(Studentized_withoutPrepiv >= Studentized)
npPooled_noRegAdj = mean(Pooled_withoutPrepiv >= Pooled)


npMaxT_RegAdj = mean(regAdj_MaxT_withoutPrepiv >= regAdj_MaxT)
npEuclidean_RegAdj = mean(regAdj_Euclidean_withoutPrepiv >= regAdj_Euclidean)
npStudentized_RegAdj = mean(regAdj_Studentized_withoutPrepiv >= regAdj_Studentized)
npPooled_RegAdj = mean(regAdj_Pooled_withoutPrepiv >= regAdj_Pooled)


# prepivoted test stat p-values
pMaxT_noRegAdj = mean(MaxT_Prepiv >= MaxT_Prepiv_obs)
pEuclidean_noRegAdj = mean(Euclidean_Prepiv >= Euclidean_Prepiv_obs)
pStudentized_noRegAdj = mean(Studentized_Prepiv >= Studentized_Prepiv_obs)
pPooled_noRegAdj = mean(Pooled_Prepiv >= Pooled_Prepiv_obs)


pMaxT_RegAdj = mean(regAdj_MaxT_Prepiv >= regAdj_MaxT_Prepiv_obs)
pEuclidean_RegAdj = mean(regAdj_Euclidean_Prepiv >= regAdj_Euclidean_Prepiv_obs)
pStudentized_RegAdj = mean(regAdj_Studentized_Prepiv >= regAdj_Studentized_Prepiv_obs)
pPooled_RegAdj = mean(regAdj_Pooled_Prepiv >= regAdj_Pooled_Prepiv_obs)

################################################################################
#                              Summary of results
################################################################################

# Quickly output table for LaTex
cat('\\begin{tabular}{ccccc}
               & \\multicolumn{2}{c}{No Prepiv.}       & \\multicolumn{2}{c}{Prepiv.}          \\\\
        Base Statistic & Without Adjustment & With Adjustment & Without Adjustment & With Adjustment \\\\
       $T_{||\\cdot||_{2}}$  &', npEuclidean_noRegAdj, '&', npEuclidean_RegAdj,'&', pEuclidean_noRegAdj, '&', pEuclidean_RegAdj,'  \\\\
       $T_{\\chi^{2}}$  &', npStudentized_noRegAdj, '&', npStudentized_RegAdj,'&', pStudentized_noRegAdj, '&', pStudentized_RegAdj,'  \\\\
       $T_{Pool} $     &', npPooled_noRegAdj, '&', npPooled_RegAdj,'&', pPooled_noRegAdj, '&', pPooled_RegAdj,'  \\\\
       $T_{|max|}$     &', npMaxT_noRegAdj, '&', npMaxT_RegAdj,'&', pMaxT_noRegAdj, '&', pMaxT_RegAdj,'  \\\\
       \\end{tabular}')
