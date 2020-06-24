## This is an R package for Bayesian rank analysis


### Installation
library(devtools)

install_github("li-xinran/BayesRankAnalysis")

### An example illustrating the use of the package

#### simulate inidividual ranking lists
```{r}
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
rank.true = rank(mu.true)  ## true ranking list
sigma.true = 5  ## noise level
Z.real = t( rmvnorm(M, mean = mu.true, sigma = sigma.true^2 * diag(N) ) ) ## scores for all rankers
fullrank.real = apply(Z.real, 2, rank)  ## observed ranking lists
```

#### Pre-steps for implementing Bayesian analysis 
```{r}
pair.comp.ten = array(NA, dim = c(N, N, M)) ## get pairwise comparison matrices from the ranking lists
for(j in 1:M){
  pair.comp.ten[,,j] = FullRankToPairComp( fullrank.real[,j] )
}
X.mat.sd = t( (t( X.mat ) - colMeans(X.mat)) / apply(X.mat, 2, sd) )  ## standardized covariates
iter.max = 1000   ## Gibbs sampler total iterations
iter.burn = 200   ## Gibbs sampler burn-in iterations
print.opt = 100  ## print a message every print.opt steps
```

#### Bayesian Analysis of Rank data without covariates 
```{r}
BAR.fit = BayesRankCovSimp(pair.comp.ten = pair.comp.ten, X.mat = matrix(NA, nrow =dim(pair.comp.ten)[1], ncol = 0), 
                              tau2.alpha = 1^2, nu.alpha = 3,
                              tau2.beta = 10^2, nu.beta = 3,
                              iter.max = iter.max, print.opt = print.opt)
                              
BAR.fit$agg.rank = apply(BAR.fit$mu[, -c(1:iter.burn)], 1, mean)  ## aggregated ranking list
RankDist(BAR.fit$agg.rank, rank.true)   ## Kendall tau distance between estimated and true ranking lists
```

#### Bayesian Analysis of Rank data with entities' Covariates
```{r}
BARC.fit = BayesRankCovSimp(pair.comp.ten = pair.comp.ten, X.mat = X.mat.sd, 
                                tau2.alpha = 1^2, nu.alpha = 3,
                                tau2.beta = 10^2, nu.beta = 3,
                                iter.max = iter.max, print.opt = print.opt)
                                
BARC.fit$agg.rank = apply(BARC.fit$mu[, -c(1:iter.burn)], 1, mean)  ## aggregated ranking list
RankDist(BARC.fit$agg.rank, rank.true)   ## Kendall tau distance between estimated and true ranking lists
```

#### Bayesian Analysis of Rank data with entities' Covariates and rankers' Weights
```{r}
BARCW.fit = BayesRankCovWeight(pair.comp.ten = pair.comp.ten, X.mat = X.mat.sd, 
                             tau2.alpha = 1^2, nu.alpha = 3,
                             tau2.beta = 10^2, nu.beta = 3,
                             iter.max = iter.max, print.opt = print.opt)
                             
BARCW.fit$agg.rank = apply(BARCW.fit$mu[, -c(1:iter.burn)], 1, mean)  ## aggregated ranking list
RankDist(BARCW.fit$agg.rank, rank.true)   ## Kendall tau distance between estimated and true ranking lists

rowMeans( BARCW.fit$weight.vec )  ## posterior means of weights for all rankers
```

#### Bayesian Analysis of Rank data with Covariates of entities and Mixture of rankers with different opinion
```{r}
BARCM.fit = BayesRankCovMix(pair.comp.ten = pair.comp.ten, X.mat = X.mat.sd, 
                            tau2.alpha = 1^2, nu.alpha = 3,
                            tau2.beta = 10^2, nu.beta = 3,
                            gamma.a = 2, gamma.b = 4,
                            iter.max = iter.max, print.opt = print.opt)
                            
BARCM.fit$agg.rank = apply(BARCM.fit$mu[, , -c(1:iter.burn)], 1, mean)  ## aggregated ranking list
RankDist(BARCM.fit$agg.rank, rank.true)   ## Kendall tau distance between estimated and true ranking lists

BARCM.fit$cluster.map = rep(NA, M)  ## get maximum a posterior estimates of cluster indicators for all rankers
for(j in 1:M){
  sum.j = table(BARCM.fit$cluster[j, -c(1:iter.burn)])
  BARCM.fit$cluster.map[j] = as.numeric( names( sum.j[which.max(sum.j)] ) )
}
BARCM.fit$cluster.map
```

#### Bayesian Analysis of Rank data with Covariates of entities and Mixture of rankers with different opinions and Weight
```{r}
BARCMW.fit = BayesRankCovMixWeight(pair.comp.ten = pair.comp.ten, X.mat = X.mat.sd, 
                                   tau2.alpha = 1^2, nu.alpha = 3,
                                   tau2.beta = 10^2, nu.beta = 3,
                                   gamma.a = 2, gamma.b = 4,
                                   iter.max = iter.max, print.opt = print.opt)
                                   
BARCMW.fit$agg.rank = apply(BARCMW.fit$mu[, , -c(1:iter.burn)], 1, mean)  ## aggregated ranking list
RankDist(BARCMW.fit$agg.rank, rank.true)   ## Kendall tau distance between estimated and true ranking lists

rowMeans( BARCMW.fit$weight.vec )  ## posterior means of weights for all rankers

BARCMW.fit$cluster.map = rep(NA, M)  ## get maximum a posterior estimates of cluster indicators for all rankers
for(j in 1:M){
  sum.j = table(BARCMW.fit$cluster[j, -c(1:iter.burn)])
  BARCMW.fit$cluster.map[j] = as.numeric( names( sum.j[which.max(sum.j)] ) )
}
BARCMW.fit$cluster.map
```


