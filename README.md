## This is an R package for Bayesian rank analysis


### Installation
library(devtools)

install_github("li-xinran/BayesRankAnalysis")

### An example illustrating the use of the package

#### simulate inidividual ranking list
```{r }
M = 10  ## number of rankers
N = 50  ## number of ranked items
L = 3   ## number of covariates
rho=0.5   ## correlation for covariates
CovMat=diag(L) ## covariance matrix for the covariates
for(i in 1:(L-1)){
  for(j in (i+1):L){
    CovMat[i,j]=rho^(abs(i-j))
    CovMat[j,i]=rho^(abs(i-j))
  }
}
X.mat = rmvnorm(N, mean = rep(0, L), sigma = CovMat) ## covariate matrix 
beta.true = c(3,2,1)
mu.true = rowSums( X.mat^2 ) + as.vector(  X.mat %*% beta.true )  ## true evaluation score
sigma.true = 5  ## noise level
Z.real = t( rmvnorm(M, mean = mu.true, sigma = sigma.true^2 * diag(N) ) ) ## scores for all rankers
fullrank.real = apply(Z.real, 2, rank)  ## observed rank lists
```


library(devtools)
# install_github("li-xinran/BayesRankAnalysis", force = TRUE)
library(BayesRankAnalysis)

#### simulation setting ####
## number of rankers
M = 10
## number of ranked items
N = 50
## number of covariates
L = 4

## true beta and true alpha
alpha.true = rep(0,N) 
beta.true = c(3,2,-1,-0.5)

## noise level
sigma.true.vec = c(1, 5, 10, 20, 40)
sigma.true = 5

## correlation for covariates
rho=0.2
## covariates for the N items
CovMat=diag(L)
for(i in 1:(L-1)){
  for(j in (i+1):L){
    CovMat[i,j]=rho^(abs(i-j))
    CovMat[j,i]=rho^(abs(i-j))
  }
}
X.mat = rmvnorm(N, mean = rep(0, L), sigma = CovMat)

## true value for mu and rank
mu.true = alpha.true + as.vector( X.mat %*% beta.true )
rank.true = rank(mu.true)

## observed rank and toplist
Z.real = t( rmvnorm(M, mean = mu.true, sigma = sigma.true^2 * diag(N) ) )
fullrank.real = apply(Z.real, 2, rank)

## transform full ranking to partial ranking
pair.comp.ten = array(NA, dim = c(N, N, M))
for(j in 1:M){
  pair.comp.ten[,,j] = FullRankToPairComp( fullrank.real[,j] )
}

## standardized covariates
X.mat.sd = t( (t( X.mat ) - colMeans(X.mat)) / apply(X.mat, 2, sd) )

#### BAR without covariates ####
iter.max = 1000
iter.burn = 200
print.opt = 100

BAR.fit = BayesRankCovSimp(pair.comp.ten = pair.comp.ten, X.mat = matrix(NA, nrow =dim(pair.comp.ten)[1], ncol = 0), 
                              tau2.alpha = 1^2, nu.alpha = 3,
                              tau2.beta = 10^2, nu.beta = 3,
                              iter.max = iter.max, print.opt = print.opt)
BAR.fit$agg.rank = apply(BAR.fit$mu[, -c(1:iter.burn)], 1, mean)
RankDist(BAR.fit$agg.rank, rank.true)

BARC.fit = BayesRankCovSimp(pair.comp.ten = pair.comp.ten, X.mat = X.mat.sd, 
                                tau2.alpha = 1^2, nu.alpha = 3,
                                tau2.beta = 10^2, nu.beta = 3,
                                iter.max = iter.max, print.opt = print.opt)
BARC.fit$agg.rank = apply(BARC.fit$mu[, -c(1:iter.burn)], 1, mean)
RankDist(BARC.fit$agg.rank, rank.true)

BARCW.fit = BayesRankCovWeight(pair.comp.ten = pair.comp.ten, X.mat = X.mat.sd, 
                             tau2.alpha = 1^2, nu.alpha = 3,
                             tau2.beta = 10^2, nu.beta = 3,
                             iter.max = iter.max, print.opt = print.opt)
BARCW.fit$agg.rank = apply(BARCW.fit$mu[, -c(1:iter.burn)], 1, mean)
RankDist(BARCW.fit$agg.rank, rank.true)

rowMeans( BARCW.fit$weight.vec )


BARCM.fit = BayesRankCovMix(pair.comp.ten = pair.comp.ten, X.mat = X.mat.sd, 
                            tau2.alpha = 1^2, nu.alpha = 3,
                            tau2.beta = 10^2, nu.beta = 3,
                            gamma.a = 2, gamma.b = 4,
                            iter.max = iter.max, print.opt = print.opt)
BARCM.fit$agg.rank = apply(BARCM.fit$mu[, , -c(1:iter.burn)], 1, mean)
RankDist(BARCM.fit$agg.rank, rank.true)

for(j in 1:M){
  sum.j = table(BARCM.fit$cluster[j, -c(1:iter.burn)])
  BARCM.fit$cluster.map[j] = as.numeric( names( sum.j[which.max(sum.j)] ) )
}
BARCM.fit$cluster.map


BARCMW.fit = BayesRankCovMixWeight(pair.comp.ten = pair.comp.ten, X.mat = X.mat.sd, 
                                   tau2.alpha = 1^2, nu.alpha = 3,
                                   tau2.beta = 10^2, nu.beta = 3,
                                   gamma.a = 2, gamma.b = 4,
                                   iter.max = iter.max, print.opt = print.opt)
BARCMW.fit$agg.rank = apply(BARCMW.fit$mu[, , -c(1:iter.burn)], 1, mean)
RankDist(BARCMW.fit$agg.rank, rank.true)

rowMeans( BARCMW.fit$weight.vec )

for(j in 1:M){
  sum.j = table(BARCMW.fit$cluster[j, -c(1:iter.burn)])
  BARCMW.fit$cluster.map[j] = as.numeric( names( sum.j[which.max(sum.j)] ) )
}
BARCMW.fit$cluster.map


