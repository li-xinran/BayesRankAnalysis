#### Gibbs update for latent vector Z given pairwise comparison matrix of one individual ranker ####
### the prior for Z is normal with mean Z and variance 1/weight or equivalently standard deviance 1/sqrt(weight) ###
GibbsUpLatentGivenRankInd <- function(pair.comp, Z, mu, weight){
  up.order = sort( rowSums( pair.comp, na.rm = TRUE ), decreasing = FALSE, index.return = TRUE )$ix
  for(i in up.order){

    set1 = which( pair.comp[i, ] == 1)
    set0 = which( pair.comp[i, ] == 0)

    if(length(set1) > 0){
      upper = min(Z[set1])
    }else{
      upper = Inf
    }

    if(length(set0) > 0){
      lower = max(Z[set0])
    }else{
      lower = -Inf
    }

    Z[i] = rtruncnorm( 1, lower, upper, mean = mu[i], sd = 1/sqrt(weight) )
  }
  return(Z)
}

# rank.matrix = matrix(NA, nrow = 5, ncol = 2)
# rank.matrix[,1] = c(2, 1, NA, NA, NA)
# rank.matrix[,2] = c(NA, NA, 3, 1, 2)
# pair.comp = PartRankToPairComp(rank.matrix, overlap = FALSE)
#
# Z = -1 * c(2, 1, 3.5, 1.5, 2.5)
# mu = rep(0, 5)
# weight = 1
# Z.new = GibbsUpLatentGivenRank(pair.comp, Z, mu, weight)
# rank(Z.new)





### Gibbs update for a group of individual rankers ###
### pair.comp.ten[,,j]: pairwise comparison matrix for jth ranker ###
### Z.mat[,j]: latent variable vector for jth ranker ###
### mu: shared mean vector for this group of rankers ###
### weight.vec[j]: weight for jth ranker ###
GibbsUpLatentGivenRankGroup <- function(pair.comp.ten, Z.mat, mu, weight.vec = rep(1, ncol(Z.mat)), n.ranker = ncol(Z.mat) ){
  for(j in 1:n.ranker){
    Z.mat[,j] = GibbsUpLatentGivenRankInd(pair.comp.ten[,,j], Z.mat[,j], mu, weight = weight.vec[j])
  }
  return(Z.mat)
}

# mu = c(1, 2, 3, 4, 5)
# Z.mat = t( rmvnorm( 10, mean = mu, sigma = diag(5) ) )
# pair.comp.ten = array(NA, dim = c(5, 5, 10))
# for(j in 1:10){
#   pair.comp.ten[,,j] = PartRankToPairComp( as.matrix( rank( - Z.mat[,j] ), ncol = 1 ) )
# }
# GibbsUpLatentGivenRankGroup(pair.comp.ten, Z.mat, mu)






### Gibbs update for the shared mean mu ###
### Z.mat[,j]: latent variable vector for jth ranker ###
### weight.vec[j]: weight for jth ranker ###
### sigma2.alpha, sigma2.beta: prior parameters for mu = (alpha, beta) ###
### para.expan: whether use parameter expansion ###
GibbsUpMuGivenLatentGroup <- function(Z.mat, X.mat = matrix(NA, nrow = nrow(Z.mat), ncol = 0), weight.vec = rep(1, ncol(Z.mat)), sigma2.alpha = 2, sigma2.beta = 1, n.ranker = ncol(Z.mat), n.item = nrow(Z.mat), p.cov = ncol(X.mat), para.expan = TRUE){
  diagLambda = c( rep(sigma2.alpha, n.item), rep(sigma2.beta, p.cov) )
  V <- cbind( diag(n.item), X.mat )

  # Sigma.old = solve( diag(1/diagLambda, nrow = n.item + p.cov) + sum(weight.vec) * t(V) %*% V )

  Sigma.inv.eigen = eigen( diag(1/diagLambda, nrow = n.item + p.cov) + sum(weight.vec) * t(V) %*% V )
  Sigma = Sigma.inv.eigen$vectors %*% diag(1/Sigma.inv.eigen$values, nrow = n.item + p.cov, ncol = n.item + p.cov) %*% t(Sigma.inv.eigen$vectors)

  lambda = t(V) %*% rowSums( t( t(Z.mat) * weight.vec ) )
  eta = Sigma %*%  lambda

  if(para.expan){
    S = sum( colSums(Z.mat^2) * weight.vec ) - as.vector( t(lambda) %*% Sigma %*% lambda )
    theta = as.vector( sqrt( S/rchisq(1, df = n.item * n.ranker) ) )
  }else{
    theta = 1
  }

  # alpha.beta = as.vector( rmvnorm(1, mean = eta/theta, sigma = Sigma) )
  alpha.beta = as.vector( eta/theta + Sigma.inv.eigen$vectors %*% diag(1/sqrt(Sigma.inv.eigen$values), nrow = n.item + p.cov, ncol = n.item + p.cov) %*% rnorm(n.item + p.cov) )


  alpha = alpha.beta[c(1:n.item)]
  beta = alpha.beta[-c(1:n.item)]

  ### parameter move
  # alpha = alpha - mean(alpha) + mean( rnorm(n.item, mean = 0, sd = sqrt(sigma2.alpha)) )

  return(list(alpha = alpha, beta = beta, theta = theta))
}

# X.mat = matrix( rnorm(15), nrow = 5, ncol = 3 )
# # X.mat = matrix( NA, nrow = 5, ncol = 0 )
# mu = c(1, 2, 3, 4, 5)
# Z.mat = t( rmvnorm( 10, mean = mu, sigma = diag(5) ) )
# pair.comp.ten = array(NA, dim = c(5, 5, 10))
# for(j in 1:10){
#   pair.comp.ten[,,j] = PartRankToPairComp( as.matrix( rank( - Z.mat[,j] ), ncol = 1 ) )
# }
# GibbsUpMuGivenLatentGroup(Z.mat, X.mat, weight.vec = rep(1, ncol(Z.mat)), para.expan = TRUE)

### Gibbs update for individual weight ####
GibbsUpWeightInd <- function(Z, mu, weight.prior.value = c(0.5, 1, 2), weight.prior.prob = rep(1/length(weight.prior.value), length(weight.prior.value)), n.item = length(Z) ){
  n.value = length(weight.prior.value)

  log.post.prob = rep(NA, n.value)
  for(k in 1:n.value){
    log.post.prob[k] = log( weight.prior.prob[k] ) + n.item/2 * log( weight.prior.value[k] ) - weight.prior.value[k]/2 * sum( (Z - mu)^2 )
  }
  log.post.prob = log.post.prob - max(log.post.prob)
  post.prob = exp(log.post.prob)
  # post.prob/sum(post.prob)

  # post.prob = rep(NA, n.value)
  # for(k in 1:n.value){
  #   post.prob[k] = weight.prior.prob[k] * dmvnorm(Z, mean = mu, sigma = diag(n.item)/weight.prior.value[k] )
  # }
  # post.prob/sum(post.prob)

  weight <- weight.prior.value[ as.vector( rmultinom(1, 1, prob = post.prob) ) == 1 ]

  return(weight)
}


# Z = -1 * c(2, 1, 3.5, 1.5, 2.5)
# mu = -Z
# GibbsUpWeightInd(Z, mu, weight.prior.value = c(0.5, 1, 2) )

### Gibbs update for weights of a group of rankers with shared mean ####
GibbsUpWeightGroup <- function(Z.mat, mu, weight.prior.value = c(0.5, 1, 2), weight.prior.prob = rep(1/length(weight.prior.value), length(weight.prior.value)), n.item = nrow(Z.mat), n.ranker = ncol(Z.mat)){
  weight.vec = rep(NA, n.ranker)

  for(j in 1:n.ranker){
    weight.vec[j] = GibbsUpWeightInd(Z.mat[,j], mu, weight.prior.value = weight.prior.value, weight.prior.prob = weight.prior.prob, n.item = n.item )
  }

  return(weight.vec)
}

# mu = c(1, 2, 3, 4, 5)
# Z.mat = t( rmvnorm( 10, mean = mu, sigma = diag(5) ) )
# GibbsUpWeightGroup(Z.mat, mu, weight.prior.value = c(0.5, 1, 2) )


### marginal density of Z.mat for BARCM
LogMargDensity <- function(Z.mat = Z.mat, X.mat = X.mat, sigma2.alpha = sigma2.alpha, sigma2.beta = sigma2.beta, n.item = nrow(Z.mat), n.ranker = ncol(Z.mat), p.cov = ncol(X.mat), weight.vec = rep(1, ncol(Z.mat))){

  # print(n.item)
  # print(n.ranker)
  # print(p.cov)
  # print(weight.vec)
  #
  diagLambda = c( rep(sigma2.alpha, n.item), rep(sigma2.beta, p.cov) )
  V <- cbind( diag(n.item), X.mat )

  Sigma = solve( diag(1/diagLambda, nrow = n.item + p.cov) + sum(weight.vec) * t(V) %*% V )
  lambda = t(V) %*% rowSums( t( t(Z.mat) * weight.vec ) )
  # eta = Sigma %*%  lambda
  S = sum( colSums(Z.mat^2) * weight.vec ) - as.vector( t(lambda) %*% Sigma %*% lambda )

  logp = - n.item * n.ranker/2 * log(2*pi) - 1/2 * sum(log(diagLambda)) + 1/2 * n.item * sum( log(weight.vec) ) + 1/2 * log(det(Sigma)) - 1/2 * S
  # logp = - 1/2 * S
  return(logp)
}


### Gibbs update for gamma of the DP prior
GibbsUpDPgamma <- function(gamma, c.vec, a, b, n = length(c.vec)){
  k = length(unique(c.vec))
  # draw eta given gamma, k
  eta = rbeta(1, gamma+1, n)
  # draw gammg given eta, k
  pi.eta = (a+k-1)/(a+k-1 + n*(b - log(eta)))
  if(runif(1)<= pi.eta){
    gamma = rgamma(1, shape = a+k, rate = b-log(eta))
  }else{
    gamma = rgamma(1, shape = a+k-1, rate = b-log(eta))
  }
  return(gamma)
}


### Gibbs update for sigma2, given prior sigma2 ~ Scale-Inv-chi2(nu, tau2) and data iid ~ N(0, sigma2)
GibbsUpsigma2 <- function(x, nu, tau2){
  if(nu < Inf){
    n.x = length(x)

    nu.post = nu + n.x
    tau2.post = ( nu * tau2 + sum(x^2) )/(nu + n.x)

    sigma2 = tau2.post * nu.post/rchisq(1, df = nu.post)
  }else{
    sigma2 = tau2
  }

  return(sigma2)
}

# x = c(1, 2, 3, 4, 5)
# nu = 3
# tau2 = 1
# GibbsUpsigma2(x, nu = Inf, tau2 = 5)

