# Import some helpful libraries
library(mvtnorm)
library(emulator)
library(MASS)

################################################################################
#########                  Some setup                               ############
################################################################################

nrand = 100#5000 # The number of simulations to run
nperm = 100 #1000 # The number of pseudo-treatment-allocations used to generate the reference distribution
ngauss = 1000 # The number of samples of Gaussian random variables used in Monte-Carlo computation of the prepivoted statistic

N= 1000 #50 # Number of individuals
nt = round(.2*N) # The number of treated units
nc = N - nt # The number of control units
treatind = c(rep(T, nt), rep(F, nc)) # A handy variable for generating treatment allocations

SdT = 10 # Standard deviation for  potential outcomes
# SdC = 1 # Standard deviation for control potential outcomes

SigmaX = cbind(c(1,.8, .2), c(.8,1,.3), c(.2,.3,1)) # Covariance matrix for covariates

# Vector for generating outcomes as noisy linear model of covariates
betaT = c(1,2,3)*.2


# Initialize lots of empty arrays to hold results
statabs = rep(0, nrand)
statstu = rep(0, nrand)
statpre = rep(0, nrand)
pvalpre =  pvalstu = pvalabs = rep(0, nrand)
reject = rep(0, nrand)
refdistabs = refdistpre = refdiststu = matrix(0, nperm, nrand)

tol = 1 # The ``tolerance" of the rerandomization balance criterion

################################################################################
#########                  Perform main simulations                 ############
################################################################################

# Loop over many simulations
i = 0 # Number of iterations so far
while(i < nrand)
{
  # Generate covariates
  X = rmvnorm(N, c(0,0,0), SigmaX)
  SighatX = cov(X)*(1/nt + 1/nc) # Covariance of the generated covariates

  # Generate potential outcomes as a noisy linear model of the covariates
  Yc = X%*%betaT - (rexp(N, 1/SdT) - SdT)
  Yt = Yc # Build the potential outomes to satisfy Fisher's sharp null

  # Select a treatment allocation to be experimentally observed
  treat = sample(treatind)

  # Given the treatment allocation, these are what the experimenter gets to see
  Ytreat = Yt[treat]
  Ycontrol = Yc[!treat]
  Xtreat = X[treat,]
  Xcontrol = X[!treat,]

  # Means and covariance estimators
  Vhat = cov(cbind(Ytreat, Xtreat))/nt + cov(cbind(Ycontrol, Xcontrol))/nc # Joint covariance estimator
  Vhatx = cov(cbind(Xtreat))/nt + cov(cbind(Xcontrol))/nc # Covariance estimator for just covariates
  mX = colMeans(Xtreat) - colMeans(Xcontrol) # Difference in means of covariates
  mah = (t(mX)%*%solve(SighatX)%*%t(t(mX))) # Mahalanobis distance between treated and control groups (based upon covariates)
  if(mah <= tol) # if the Mahalanobis distance doesn't exceed the tolerance then the treatment allocation is used; otherwise we re-draw a new treatment allocation
  {
    i = i + 1 # Increment the counter
    mvec =mean(Ytreat) - mean(Ycontrol) # Difference in means
    vpool = var(Ytreat)/nt + var(Ycontrol)/nc # Neyman variance estimator
    se = sqrt(vpool) # Standard error

    # Compute the UNPREPIVOTED statistics
    statabs[i] = abs(mvec) # The absolute difference in means statistic
    statstu[i] = abs(mvec)/se # The studentized absolute difference in means statistic

    # Simulates Gaussian data for use in Monte-Carlo computation of prepivoted statistics
    ZZ = mvrnorm(ngauss, rep(0, 4), Vhat) # Samples from a centered Gaussian with covariance matrix given by Vhat
    sSighat = solve(SighatX) # Invert the covariance of the generated covariates
    pp = apply(ZZ[,-1], 1, function(y){quad.form(sSighat, y)} ) # Compute the Mahalanobis distance based upon the Gaussian simulations

    # Compute the PREPIVOTED statistics (via Monte-Carlo approximation) while accounting for the rerandomization (by enforcing Mahalanobis distance satisfying the tolerance)
    statpre[i] = mean(abs(ZZ[,1]) <= statabs[i] & pp <= tol)/mean(pp <= tol) # The prepivoted absolute difference in means statistic

    # Now we generate the reference distributions for the UNPREPIVOTED and PREPIVOTED statistics
    Yperm = c(Ytreat, Ycontrol) # A handy way to record the observed outcomes
    Xperm = rbind(Xtreat, Xcontrol) # A handy way to record the observed covariates
    # Initialize lots of empty arrays to hold results
    permdistpre <- permdiststu <- permdistabs <- rep(0, nperm)
    j = 0 # Number of iterations so far
    while(j < nperm)
    {
      # Generate a pseudo-treatment-allocation
      pseudoTreatmentAlloc = sample(treatind)

      # Assuming that Fisher's sharp null holds these are what the experimenter would have observed under the allocation pseudoTreatmentAlloc
      Ytreatperm = Yperm[pseudoTreatmentAlloc]
      Ycontrolperm = Yperm[!pseudoTreatmentAlloc]
      Xtreatperm = Xperm[pseudoTreatmentAlloc,]
      Xcontrolperm = Xperm[!pseudoTreatmentAlloc,]

      # Means and Mahalanobis distance (just like above)
      mXperm = colMeans(Xtreatperm) - colMeans(Xcontrolperm)
      mahperm = (t(mXperm)%*%solve(SighatX)%*%t(t(mXperm)))
      if(mahperm <= tol) # if the Mahalanobis distance doesn't exceed the tolerance then the treatment allocation is used; otherwise we re-draw a new treatment allocation
      {
        j = j+1 # Increment the counter

        # Means and covariance estimators (just like above)
        mvecperm = mean(Ytreatperm) - mean(Ycontrolperm)
        vpoolperm = var(Ytreatperm)/nt + var(Ycontrolperm)/nc
        seperm = sqrt(vpoolperm)
        Vhatperm = cov(cbind(Ytreatperm, Xtreatperm))/nt + cov(cbind(Ycontrolperm, Xcontrolperm))/nc
        Vhatxperm = cov(cbind(Xtreatperm))/nt + cov(cbind(Xcontrolperm))/nc

        # The UNPREPIVOTED statistics based upon Ytreatperm and Ycontrolperm
        permdistabs[j] = abs(mvecperm)
        permdiststu[j] = abs(mvecperm)/seperm

        # Simulates Gaussian data for use in Monte-Carlo computation of prepivoted statistics
        ZZ = mvrnorm(ngauss, rep(0, 4), Vhatperm)
        sSighat = solve(SighatX)
        pp = apply(ZZ[,-1], 1, function(y){quad.form(sSighat, y)} )

        # The PREPIVOTED statistic based upon Ytreatperm and Ycontrolperm
        permdistpre[j] = mean(abs(ZZ[,1]) <= permdistabs[j] & pp <= tol)/mean(pp <= tol)
      }

    }

    # p-Values
    pvalabs[i] = (1 + sum(permdistabs >= statabs[i]))/(1+nperm)
    pvalstu[i] = (1 + sum(permdiststu >= statstu[i]))/(1+nperm)
    pvalpre[i] = (1 + sum(permdistpre >= statpre[i]))/(1+nperm)

    # Reference distributions
    refdistabs[,i] = permdistabs
    refdiststu[,i] = permdiststu
    refdistpre[,i] = permdistpre

    print(i) # Progress counter
  }

}

################################################################################
#########                  Output the results                       ############
################################################################################

print(mean(pvalabs[1:(i-1)] <= .1))
print(mean(pvalstu[1:(i-1)] <= .1))
print(mean(pvalpre[1:(i-1)] <= .1))
print(mean(1-statpre[1:(i-1)] <= .1))

################################################################################
#########                    Save the results                       ############
################################################################################
save(file = "rerandSharp", pvalabs, pvalstu, pvalpre, statabs, statstu, statpre, refdistabs, refdiststu, refdistpre)
