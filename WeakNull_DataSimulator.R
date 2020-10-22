# Import some helpful libraries
library(mvtnorm)
library(MASS)
library(emulator)
library(CompQuadForm)

################################################################################
#########                  Some setup                               ############
################################################################################

K= 25
rhoT = 0 # Correlation for treated outcomes
rhoC = .95 # Correlation for control outcomes

# Generate an equicorrelated covariance matrix
SigmaT = matrix(rhoT, K, K)
diag(SigmaT) = 1
#SigmaT = SigmaT

# Generate an equicorrelated covariance matrix
SigmaC = 1*matrix(rhoC, K, K)
diag(SigmaC) = 1
#SigmaC =1*SigmaC

N=300 # The number of individuals in the population
nrand = 5000 # The number of simulations to run
nperm = 500 # The number of pseudo-treatment-allocations used to generate the reference distribution
ngauss = 500 # The number of samples of Gaussian random variables used in Monte-Carlo computation of the prepivoted statistic

# The average treatment effect for each outcome
#   When set to 0 Neymans's weak null holds
#   (This can be modified to perform power simulations)
taubar = rep(.2, K)
if(all(taubar == 0))
{
  cat("Neyman's weak null holds.\n")
}else
{
  cat('The avergae treatment effect vector is: \n', taubar, "\n")
}

# Initialize lots of empty arrays to hold results
statprepool <- statprechi <- statpremax <- statpool <- statchi <- statmax <- rep(0, nrand)
pvalchiLS <- pvalmaxLS <- pvalpoolLS<- pvalprepool <- pvalprechi <- pvalpremax <- pvalpool <- pvalchi <- pvalmax <- rep(0, nrand)
reject = rep(0, nrand)
refprepool <- refprechi <- refpremax <- refpool<- refchi <- refmax <- matrix(0, nperm, nrand)


nt = round(.2*N) # The number of treated units
nc = N - nt # The number of control units
treatind = c(rep(T, nt), rep(F, nc)) # A handy variable for generating treatment allocations

################################################################################
#########                  Perform main simulations                 ############
################################################################################

# Loop over many simulations
for(i in 1:nrand)
{
  # Generate some potential outcomes
  Yttemp = rmvnorm(N, rep(0, K), SigmaT)
  Yc2 = rmvnorm(N, rep(0, K), SigmaC)
  Yc = Yc2
  Yt = Yttemp
  #Yc = Yttemp
  for(k in 1:K)
  {
    Yt[,k] = Yttemp[,k] - mean(Yttemp[,k] - Yc[,k]) + taubar[k]
  }

  # Select a treatment allocation to be experimentally observed
  treat = sample(treatind)

  # Given the treatment allocation, these are what the experimenter gets to see
  Ytreat = Yt[treat,]
  Ycontrol = Yc[!treat,]

  # Means and covariance estimators
  mvec = colMeans(Ytreat) - colMeans(Ycontrol) # The difference in means
  Cov_unpooled = cov(Ytreat)/nt + cov(Ycontrol)/nc # A Neyman-style covariance estimator
  Cov_pooled = ((cov(Ytreat)*(nt-1) + cov(Ycontrol)*(nc-1))/(nt+nc-2))*(1/nt + 1/nc) # A pooled covariance estimator
  se = sqrt(diag(Cov_unpooled)) # Standard errors computed from the diagonal of the Neyman covariance estimator
  Corr_unpooled = cov2cor(Cov_unpooled)# Correlation matrix for the Neyman covariance estimator

  # Compute the UNPREPIVOTED statistics
  statmax[i] = max(abs(mvec)/se) # Max coordinate-wise t-statistic
  statpool[i] = t(mvec)%*%solve(Cov_pooled)%*%mvec # A Wald-style statistic using the pooled covariance estimator
  statchi[i] = t(mvec)%*%solve(Cov_unpooled)%*%mvec # A Wald-style statistic using the un-pooled covariance estimator

  # Simulates Gaussian data for use in Monte-Carlo computation of prepivoted statistics
  ZZ = mvrnorm(ngauss, rep(0, K), Cov_unpooled) # Samples from a centered K-dimensional Gaussian with covariance Cov_unpooled
  pp = apply(ZZ, 1, function(y){quad.form.inv(Cov_pooled, y)} ) # Computes the Wald-statistic with pooled covariance estimator on the Gaussian samples
  ZZstu = abs(ZZ)*matrix(1/sqrt(diag(Cov_unpooled)), nrow = ngauss, ncol = K, byrow = T) # ``Studentizes'' each of the Gaussian samples
  pp2 = apply(ZZstu, 1, max) # Computes the max t-statistic on the Gaussian samples

  # Compute the PREPIVOTED statistics (via Monte-Carlo approximation)
  statprechi[i] = pchisq(statchi[i], K) # The prepivoted statistic corresponding to statchi (See Example 2 of https://arxiv.org/abs/2002.06654)
  statprepool[i] = mean(pp <= statpool[i]) # The prepivoted statistic corresponding to statpool (See Example 2 of https://arxiv.org/abs/2002.06654)
  statpremax[i] = mean(pp2 <= statmax[i]) # The prepivoted statistic corresponding to statmax (See Example 3 of https://arxiv.org/abs/2002.06654)

  # Now we generate the reference distributions for the UNPREPIVOTED and PREPIVOTED statistics
  Yperm = rbind(Ytreat, Ycontrol) # A handy way to record the observed outcomes
  # Initialize lots of empty arrays to hold results
  permpremax <- permprechi <- permprepool <- permpool <- permchi <- permmax <- rep(0, nperm)
  # Generate the reference distributions as if Fisher's sharp null holds
  for(j in 1:nperm)
  {
    # Generate a pseudo-treatment-allocation
    pseudoTreatmentAlloc = sample(treatind)

    # Assuming that Fisher's sharp null holds these are what the experimenter would have observed under the allocation pseudoTreatmentAlloc
    Ytreatperm = Yperm[pseudoTreatmentAlloc,]
    Ycontrolperm = Yperm[!pseudoTreatmentAlloc,]

    # Means and covariance estimators (just like above)
    mveCov_poolederm = colMeans(Ytreatperm) - colMeans(Ycontrolperm)
    Cov_unpooled_perm = cov(Ytreatperm)/nt + cov(Ycontrolperm)/nc
    Cov_pooled_perm = ((cov(Ytreatperm)*(nt-1) + cov(Ycontrolperm)*(nc-1))/(nt+nc-2))*(1/nt + 1/nc)
    seperm = sqrt(diag(Cov_unpooled_perm))
    Corr_unpooled_perm = cov2cor(Cov_unpooled_perm)

    # The UNPREPIVOTED statistics based upon Ytreatperm and Ycontrolperm
    permmax[j] = max(abs(mveCov_poolederm)/seperm)
    permpool[j] = t(mveCov_poolederm)%*%solve(Cov_pooled_perm)%*%mveCov_poolederm
    permchi[j] = t(mveCov_poolederm)%*%solve(Cov_unpooled_perm)%*%mveCov_poolederm

    # Simulates Gaussian data for use in Monte-Carlo computation of prepivoted statistics
    ZZ = mvrnorm(ngauss, rep(0, K), Cov_unpooled_perm)
    ZZstu = abs(ZZ)*matrix(1/sqrt(diag(Cov_unpooled_perm)), nrow = ngauss, ncol = K, byrow = T)
    pp = apply(ZZ, 1, function(y){quad.form.inv(Cov_pooled_perm, y)} )
    pp2 = apply(ZZstu, 1, max)

    # The PREPIVOTED statistics based upon Ytreatperm and Ycontrolperm
    permprechi[j] = pchisq(permchi[j], K)
    permprepool[j] = mean(pp <= permpool[j])
    permpremax[j] = mean(pp2 <= permmax[j])
  }

  # p-Values for UNPREPIVOTED statistics
  pvalpool[i] = (1 + sum(permpool >= statpool[i]))/(1+nperm)
  pvalchi[i] = (1 + sum(permchi >= statchi[i]))/(1+nperm)
  pvalmax[i] = (1 + sum(permmax >= statmax[i]))/(1+nperm)

  # p-Values for PREPIVOTED statistics
  pvalprepool[i] = (1 + sum(permprepool >= statprepool[i]))/(1+nperm)
  pvalprechi[i] = (1 + sum(permprechi >= statprechi[i]))/(1+nperm)
  pvalpremax[i] = (1 + sum(permpremax >= statpremax[i]))/(1+nperm)

  # Record the reference distribtions
  refpool[,i] = permpool
  refprepool[,i] = permprepool
  refchi[,i] = permchi
  refprechi[,i] = permprechi
  refmax[,i] = permmax
  refpremax[,i] = permpremax

  print(i) # Progress counter
}

cat("\n\nHotelling, UNPOOOLED\n")
print(mean(pvalchi[1:(i-1)] <= .25))
print(mean(pvalprechi[1:(i-1)] <= .25))
print(mean(1-statprechi <= .25)) # Used to generate LS column data in Table 2


cat("Hotelling, POOOLED\n")
print(mean(pvalpool[1:(i-1)] <= .25))
print(mean(pvalprepool[1:(i-1)] <= .25))
print(mean(1-statprepool <= .25)) # Used to generate LS column data in Table 2

cat("MAX t-STAT\n")
print(mean(pvalmax[1:(i-1)] <= .25))
print(mean(pvalpremax[1:(i-1)] <= .25))
print(mean(1-statpremax <= .25)) # Used to generate LS column data in Table 2



# To save results.
filename = paste('Power_multN300weak_tau', tau, sep ="")
save(file = filename, pvalmax, pvalpremax, pvalpool, pvalprepool, pvalchi, pvalprechi, statmax, statpremax, statpool, statprepool, statchi, statprechi, refmax, refpremax, refpool, refprepool, refchi, refprechi)
