#' Bayesian Analysis of Rank data with entities' Covariates
#'
#' Implement the Bayesian model for rand data with ranked entities' covariates information.
#' @import truncnorm
#' @import mvtnorm
#' @param pair.com.ten An \eqn{N} by \eqn{N} by \eqn{M} pairwise comparison tensor for all \eqn{N} entities and \eqn{M} rankers, where the (\eqn{i},\eqn{j},\eqn{m}) element equals 1 if \eqn{i} is ranked higher than \eqn{j} by ranker \eqn{m}, 0 if \eqn{i} is ranker lower than \eqn{j}, and NA if the relation between \eqn{i} and \eqn{j} is missing. Note that the diagonal elements (\eqn{i},\eqn{i},\eqn{m})'s for all rankers should be set to NA as well.
#' @param X.mat An \eqn{N} by \eqn{L} covariate matrix for the \eqn{N} entities with \eqn{L} covariates.
#' @param tau2.alpha The scale parameter for the scaled inverse chi-squared prior on \eqn{\sigma^2_alpha}.
#' @param nu.alpha The d.f. for the scaled inverse chi-squared prior on \eqn{\sigma^2_alpha}.
#' @param tau2.beta The scale parameter for the scaled inverse chi-squared prior on \eqn{\sigma^2_beta}.
#' @param nu.beta The d.f. for the scaled inverse chi-squared prior on \eqn{\sigma^2_beta}.
#' @param iter.max Number of iterations for Gibbs sampler.
#' @param para.expan Logical variable for whether using parameter expansion in the Gibbs sampler.
#' @return A list containing posterior samples of all the missing evaluation scores for all rankers and all the model parameters.
#' @export
BayesRankCovSimp <- function(pair.comp.ten, X.mat = matrix(NA, nrow =dim(pair.comp.ten)[1], ncol = 0),
                             tau2.alpha = 5^2, nu.alpha = 3,
                             tau2.beta = 5^2, nu.beta = 3,
                             n.item = dim(pair.comp.ten)[1], n.ranker = dim(pair.comp.ten)[3], p.cov = ncol(X.mat),
                             iter.max = 5000, para.expan = TRUE, print.opt = 100,
                             initial.list = NULL){
  ## store MCMC draws
  draw = list(
    Z.mat = array(NA, dim = c(n.item, n.ranker, iter.max)),
    alpha = array(NA, dim = c(n.item, iter.max)),
    beta = array(NA, dim = c(p.cov, iter.max)),
    mu = array(NA, dim = c(n.item, iter.max)),
    sigma2.alpha = rep(NA, iter.max),
    sigma2.beta = rep(NA, iter.max)
  )

  if(is.null(initial.list)){
    ## initial values for Z
    Z.mat = matrix(NA, nrow = n.item, ncol = n.ranker)
    for(j in 1:n.ranker){
      Z.mat[sort( rowSums( pair.comp.ten[,,j], na.rm = TRUE ), decreasing = FALSE, index.return = TRUE )$ix, j] = (c(n.item : 1) - (1+n.item)/2)/sd(c(n.item : 1))
    }

    ## initial values for alpha, beta and thus mu
    alpha = rep(0, n.item)
    beta = rep(0, p.cov)
    mu = as.vector( alpha + X.mat %*% beta )

    ## initial values for sigma2.alpha and sigma2.beta
    sigma2.alpha = tau2.alpha
    sigma2.beta = tau2.beta

  }else{
    Z.mat = initial.list$Z.mat
    alpha = initial.list$alpha
    beta = initial.list$beta
    mu = as.vector( alpha + X.mat %*% beta )
    sigma2.alpha = initial.list$sigma2.alpha
    sigma2.beta = initial.list$sigma2.beta
  }



  ## store initial value
  draw$Z.mat[,,1] = Z.mat
  draw$alpha[,1] = alpha
  draw$beta[,1] = beta
  draw$mu[,1] = mu
  draw$sigma2.alpha[1] = sigma2.alpha
  draw$sigma2.beta[1] = sigma2.beta

  ## Gibbs iteration
  for(iter in 2:iter.max){

    # update Z.mat given (alpha, beta) or equivalently mu
    Z.mat = GibbsUpLatentGivenRankGroup(pair.comp.ten = pair.comp.ten, Z.mat = Z.mat, mu = mu, weight.vec = rep(1, n.ranker), n.ranker = n.ranker )

    # update (alpha, beta) or equivalently mu
    mean.para.update = GibbsUpMuGivenLatentGroup(Z.mat = Z.mat, X.mat = X.mat, weight.vec = rep(1, n.ranker), sigma2.alpha = sigma2.alpha, sigma2.beta = sigma2.beta, n.ranker = n.ranker, n.item = n.item, p.cov = p.cov, para.expan = para.expan)

    ### for check only
    Z.mat = Z.mat/mean.para.update$theta

    alpha = mean.para.update$alpha
    beta = mean.para.update$beta
    mu = as.vector( alpha + X.mat %*% beta )

    # update hyper para sigma2.alpha and sigma2.beta
    sigma2.alpha = GibbsUpsigma2(alpha, nu.alpha, tau2.alpha)
    if(p.cov > 0){
      sigma2.beta = GibbsUpsigma2(beta, nu.beta, tau2.beta)
    }


    # store value at this iteration
    draw$Z.mat[,,iter] = Z.mat
    draw$alpha[,iter] = alpha
    draw$beta[,iter] = beta
    draw$mu[,iter] = mu
    draw$sigma2.alpha[iter] = sigma2.alpha
    draw$sigma2.beta[iter] = sigma2.beta

    # print iteration number
    if(iter %% print.opt == 0){
      print(paste("Gibbs Iteration", iter))
      # print(c(sigma2.alpha, sigma2.beta))
    }

  }
  return(draw)
}



### BARCW model

#' Bayesian Analysis of Rank data with entities' Covariates and rankers' Weights.
#'
#' Implement the Bayesian model for rand data with ranked entities' covariates information and rankers' with varying qualities or weights.
#' @import truncnorm
#' @import mvtnorm
#' @param pair.com.ten An \eqn{N} by \eqn{N} by \eqn{M} pairwise comparison tensor for all \eqn{N} entities and \eqn{M} rankers, where the (\eqn{i},\eqn{j},\eqn{m}) element equals 1 if \eqn{i} is ranked higher than \eqn{j} by ranker \eqn{m}, 0 if \eqn{i} is ranker lower than \eqn{j}, and NA if the relation between \eqn{i} and \eqn{j} is missing. Note that the diagonal elements (\eqn{i},\eqn{i},\eqn{m})'s for all rankers should be set to NA as well.
#' @param X.mat An \eqn{N} by \eqn{L} covariate matrix for the \eqn{N} entities with \eqn{L} covariates.
#' @param tau2.alpha The scale parameter for the scaled inverse chi-squared prior on \eqn{\sigma^2_alpha}.
#' @param nu.alpha The d.f. for the scaled inverse chi-squared prior on \eqn{\sigma^2_alpha}.
#' @param tau2.beta The scale parameter for the scaled inverse chi-squared prior on \eqn{\sigma^2_beta}.
#' @param nu.beta The d.f. for the scaled inverse chi-squared prior on \eqn{\sigma^2_beta}.
#' @param weight.prior.value A vector for the support of the discrete prior on weight parameter.
#' @param weight.prior.prob A vector for the probability mass of the discrete prior on weight parameter.
#' @param iter.max Number of iterations for Gibbs sampler.
#' @param para.expan Logical variable for whether using parameter expansion in the Gibbs sampler.
#' @return A list containing posterior samples of all the missing evaluation scores for all rankers and all the model parameters.
#' @export
BayesRankCovWeight <- function(pair.comp.ten, X.mat = matrix(NA, nrow =dim(pair.comp.ten)[1], ncol = 0),
                               tau2.alpha = 5^2, nu.alpha = 3,
                               tau2.beta = 5^2, nu.beta = 3,
                               weight.prior.value = c(0.5, 1, 2), weight.prior.prob = rep(1/length(weight.prior.value), length(weight.prior.value)),
                               n.item = dim(pair.comp.ten)[1], n.ranker = dim(pair.comp.ten)[3], p.cov = ncol(X.mat),
                               iter.max = 5000, para.expan = TRUE, print.opt = 100,
                               initial.list = NULL){
  ## store MCMC draws
  draw = list(
    Z.mat = array(NA, dim = c(n.item, n.ranker, iter.max)),
    alpha = array(NA, dim = c(n.item, iter.max)),
    beta = array(NA, dim = c(p.cov, iter.max)),
    mu = array(NA, dim = c(n.item, iter.max)),
    weight.vec = array(NA, dim = c(n.ranker, iter.max) ),
    sigma2.alpha = rep(NA, iter.max),
    sigma2.beta = rep(NA, iter.max)
  )

  if(is.null(initial.list)){
    ## initial values for Z
    Z.mat = matrix(NA, nrow = n.item, ncol = n.ranker)
    for(j in 1:n.ranker){
      Z.mat[sort( rowSums( pair.comp.ten[,,j], na.rm = TRUE ), decreasing = FALSE, index.return = TRUE )$ix, j] = (c(n.item : 1) - (1+n.item)/2)/sd(c(n.item : 1))
    }

    ## initial values for alpha, beta and thus mu
    alpha = rep(0, n.item)
    beta = rep(0, p.cov)
    mu = as.vector( alpha + X.mat %*% beta )

    ## initial values for weights
    weight.vec = rep(1, n.ranker)

    ## initial values for sigma2
    sigma2.alpha = tau2.alpha
    sigma2.beta = tau2.beta
  }else{

    Z.mat = initial.list$Z.mat
    alpha = initial.list$alpha
    beta = initial.list$beta
    mu = as.vector( alpha + X.mat %*% beta )
    weight.vec = initial.list$weight.vec
    sigma2.alpha = initial.list$sigma2.alpha
    sigma2.beta = initial.list$sigma2.beta

  }

  ## store initial value
  draw$Z.mat[,,1] = Z.mat
  draw$alpha[,1] = alpha
  draw$beta[,1] = beta
  draw$mu[,1] = mu
  draw$weight.vec[,1] = weight.vec

  ## Gibbs iteration
  for(iter in 2:iter.max){

    # update Z.mat given (alpha, beta) or equivalently mu
    Z.mat = GibbsUpLatentGivenRankGroup(pair.comp.ten = pair.comp.ten, Z.mat = Z.mat, mu = mu, weight.vec = weight.vec, n.ranker = n.ranker )

    # update (alpha, beta) or equivalently mu
    mean.para.update = GibbsUpMuGivenLatentGroup(Z.mat = Z.mat, X.mat = X.mat, weight.vec = weight.vec, sigma2.alpha = sigma2.alpha, sigma2.beta = sigma2.beta, n.ranker = n.ranker, n.item = n.item, p.cov = p.cov, para.expan = para.expan)

    ### for check only
    Z.mat = Z.mat/mean.para.update$theta

    alpha = mean.para.update$alpha
    beta = mean.para.update$beta
    mu = as.vector( alpha + X.mat %*% beta )

    ### update weight
    weight.vec = GibbsUpWeightGroup(Z.mat = Z.mat, mu = mu, weight.prior.value = weight.prior.value, weight.prior.prob = weight.prior.prob, n.item = n.item, n.ranker = n.ranker)

    ### update sigma2
    sigma2.alpha = GibbsUpsigma2(alpha, nu.alpha, tau2.alpha)
    if(p.cov > 0){
      sigma2.beta = GibbsUpsigma2(beta, nu.beta, tau2.beta)
    }

    # store value at this iteration
    draw$Z.mat[,,iter] = Z.mat
    draw$alpha[,iter] = alpha
    draw$beta[,iter] = beta
    draw$mu[,iter] = mu
    draw$weight.vec[, iter] = weight.vec
    draw$sigma2.alpha[iter] = sigma2.alpha
    draw$sigma2.beta[iter] = sigma2.beta

    # print iteration number
    if(iter %% print.opt == 0){
      print(paste("Gibbs Iteration", iter))
      # print(table(weight.vec))
      # print(c(sigma2.alpha, sigma2.beta))
    }
  }
  return(draw)
}




### BARCM model
#' Bayesian Analysis of Rank data with Covariates of entities and Mixture of rankers with different opinion
#'
#' Implement the Bayesian model for rand data with ranked entities' covariates information and mixture of rankers with different ranking opinions.
#' @import truncnorm
#' @import mvtnorm
#' @param pair.com.ten An \eqn{N} by \eqn{N} by \eqn{M} pairwise comparison tensor for all \eqn{N} entities and \eqn{M} rankers, where the (\eqn{i},\eqn{j},\eqn{m}) element equals 1 if \eqn{i} is ranked higher than \eqn{j} by ranker \eqn{m}, 0 if \eqn{i} is ranker lower than \eqn{j}, and NA if the relation between \eqn{i} and \eqn{j} is missing. Note that the diagonal elements (\eqn{i},\eqn{i},\eqn{m})'s for all rankers should be set to NA as well.
#' @param X.mat An \eqn{N} by \eqn{L} covariate matrix for the \eqn{N} entities with \eqn{L} covariates.
#' @param tau2.alpha The scale parameter for the scaled inverse chi-squared prior on \eqn{\sigma^2_alpha}.
#' @param nu.alpha The d.f. for the scaled inverse chi-squared prior on \eqn{\sigma^2_alpha}.
#' @param tau2.beta The scale parameter for the scaled inverse chi-squared prior on \eqn{\sigma^2_beta}.
#' @param nu.beta The d.f. for the scaled inverse chi-squared prior on \eqn{\sigma^2_beta}.
#' @param gamma.a The shape parameter for the Gamma prior on the concentration parameter of Dirichlet process.
#' @param gamma.b The rate parameter for the Gamma prior on the concentration parameter of Dirichlet process.
#' @param iter.max Number of iterations for Gibbs sampler.
#' @param para.expan Logical variable for whether using parameter expansion in the Gibbs sampler.
#' @return A list containing posterior samples of all the missing evaluation scores for all rankers and all the model parameters.
#' @export
BayesRankCovMix <- function(pair.comp.ten, X.mat = matrix(NA, nrow =dim(pair.comp.ten)[1], ncol = 0),
                            tau2.alpha = 5^2, nu.alpha = 3,
                            tau2.beta = 5^2, nu.beta = 3,
                            gamma.a = 2, gamma.b = NULL,
                            n.item = dim(pair.comp.ten)[1], n.ranker = dim(pair.comp.ten)[3], p.cov = ncol(X.mat),
                            iter.max = 5000, para.expan = TRUE, print.opt = 100,
                            initial.list = NULL){
  ## store MCMC draws
  draw = list(
    Z.mat = array(NA, dim = c(n.item, n.ranker, iter.max)),
    alpha = array(NA, dim = c(n.item, n.ranker, iter.max)),
    beta = array(NA, dim = c(p.cov, n.ranker, iter.max)),
    mu = array(NA, dim = c(n.item, n.ranker, iter.max)),
    cluster = array(NA, dim = c(n.ranker, iter.max)),
    sigma2.alpha = rep(NA, iter.max),
    sigma2.beta = rep(NA, iter.max),
    gamma = rep(NA, iter.max)
  )

  if(is.null(initial.list)){
    ## initial values for Z
    Z.mat = matrix(NA, nrow = n.item, ncol = n.ranker)
    for(j in 1:n.ranker){
      Z.mat[sort( rowSums( pair.comp.ten[,,j], na.rm = TRUE ), decreasing = FALSE, index.return = TRUE )$ix, j] = (c(n.item : 1) - (1+n.item)/2)/sd(c(n.item : 1))
    }

    ## initial values for alpha, beta and thus mu
    alpha = array(0, dim = c(n.item, n.ranker))
    beta = array(0, dim = c(p.cov, n.ranker))
    mu = alpha + X.mat %*% beta

    ## initial value for cluster indicator
    # cluster = rep(1, n.ranker)
    cluster = sample(c(1:ceiling(log(n.ranker))), n.ranker, replace = TRUE )

    ## initial value for hyper para of DP
    sigma2.alpha = tau2.alpha
    sigma2.beta = tau2.beta
    if(is.null(gamma.b)){
      gamma = gamma.a
    }else{
      gamma = gamma.a/gamma.b
    }

  }else{
    Z.mat = initial.list$Z.mat
    alpha = initial.list$alpha
    beta = initial.list$beta
    mu = alpha + X.mat %*% beta
    cluster = initial.list$cluster
    sigma2.alpha = initial.list$sigma2.alpha
    sigma2.beta = initial.list$sigma2.beta
    gamma = initial.list$gamma
  }

  ## store initial value
  draw$Z.mat[,,1] = Z.mat
  draw$alpha[,,1] = alpha
  draw$beta[,,1] = beta
  draw$mu[,,1] = mu
  draw$cluster[, 1] = cluster
  draw$sigma2.alpha[1] = sigma2.alpha
  draw$sigma2.beta[1] = sigma2.beta
  draw$gamma[1] = gamma

  c.max = max(cluster)

  ## Gibbs iteration
  for(iter in 2:iter.max){
    ### Gibbs update for each clsuter
    alpha.all = c()
    beta.all = c()
    for(c in unique( cluster ) ){
      # update (alpha, beta) or equivalently mu
      mean.para.update = GibbsUpMuGivenLatentGroup(Z.mat = Z.mat[, cluster == c, drop=FALSE], X.mat = X.mat, weight.vec = rep(1, sum(cluster==c)), sigma2.alpha = sigma2.alpha, sigma2.beta = sigma2.beta, n.ranker = sum(cluster==c), n.item = n.item, p.cov = p.cov, para.expan = para.expan)

      alpha.c = mean.para.update$alpha
      beta.c = mean.para.update$beta
      mu.c = as.vector( alpha.c + X.mat %*% beta.c )

      alpha[, cluster==c] = alpha.c
      beta[, cluster==c] = beta.c
      mu[, cluster==c] = mu.c

      alpha.all = c(alpha.all, alpha.c)
      beta.all = c(beta.all, beta.c)

      ### for check only
      Z.mat[, cluster == c] = Z.mat[, cluster == c]/mean.para.update$theta

      # update Z.mat given (alpha, beta) or equivalently mu
      # mu.c = mu[, cluster==c, drop=FALSE][,1]
      Z.mat[, cluster == c] = GibbsUpLatentGivenRankGroup(pair.comp.ten = pair.comp.ten[,,cluster==c, drop=FALSE], Z.mat = Z.mat[, cluster==c, drop=FALSE], mu = mu.c, weight.vec = rep(1, sum(cluster==c)), n.ranker = sum(cluster==c) )
    }

    ### Gibbs update for hyper parameters of DP
    sigma2.alpha = GibbsUpsigma2(alpha.all, nu.alpha, tau2.alpha)
    if(p.cov > 0){
      sigma2.beta = GibbsUpsigma2(beta.all, nu.beta, tau2.beta)
    }
    if(!is.null(gamma.b)){
      gamma = GibbsUpDPgamma( gamma = gamma, c.vec = cluster, a = gamma.a, b = gamma.b, n = n.ranker )
    }

    ### Gibbs update for the cluster label
    # for(i in 1:n.ranker){
    for(i in sample(c(1:n.ranker)) ){
      candidate = c( unique( cluster[-i] ), c.max+1 )
      logp = rep(NA, length(candidate))
      for(c in candidate){
        if( sum(cluster[-i] == c) > 0 ){
          set = setdiff( which(cluster == c), i )
          logp[candidate==c] = log( sum(cluster[-i] == c) ) + LogMargDensity(Z.mat = Z.mat[, c(set, i), drop=FALSE], X.mat = X.mat, sigma2.alpha = sigma2.alpha, sigma2.beta = sigma2.beta) - LogMargDensity(Z.mat = Z.mat[, set, drop=FALSE], X.mat = X.mat, sigma2.alpha = sigma2.alpha, sigma2.beta = sigma2.beta)
        }else{
          logp[candidate==c] = log(gamma) + LogMargDensity(Z.mat = Z.mat[, i, drop=FALSE], X.mat = X.mat, sigma2.alpha = sigma2.alpha, sigma2.beta = sigma2.beta)
        }
      }
      logp = logp - max(logp)
      cluster[i] = candidate[which( as.vector( rmultinom(1, 1, prob = exp(logp)) ) == 1 )]
      if(cluster[i] == ( c.max+1)){
        c.max = c.max + 1
      }
    }

    ### store vlaues
    draw$Z.mat[,,iter] = Z.mat
    draw$alpha[,,iter] = alpha
    draw$beta[,,iter] = beta
    draw$mu[,,iter] = mu
    draw$cluster[, iter] = cluster
    draw$sigma2.alpha[iter] = sigma2.alpha
    draw$sigma2.beta[iter] = sigma2.beta
    draw$gamma[iter] = gamma

    # print iteration number
    if(iter %% print.opt == 0){
      print(paste("Gibbs Iteration", iter))
      # print(table(cluster))
      # print(c(sigma2.alpha, sigma2.beta))
      # print(gamma)
    }
  }
  return(draw)
}


### BARCMW model
#' Bayesian Analysis of Rank data with Covariates of entities and Mixture of rankers with different opinions and Weight
#'
#' Implement the Bayesian model for rand data with ranked entities' covariates information and mixture of rankers with different ranking opinions, where the rankers can also have different qualities or weights.
#' @import truncnorm
#' @import mvtnorm
#' @param pair.com.ten An \eqn{N} by \eqn{N} by \eqn{M} pairwise comparison tensor for all \eqn{N} entities and \eqn{M} rankers, where the (\eqn{i},\eqn{j},\eqn{m}) element equals 1 if \eqn{i} is ranked higher than \eqn{j} by ranker \eqn{m}, 0 if \eqn{i} is ranker lower than \eqn{j}, and NA if the relation between \eqn{i} and \eqn{j} is missing. Note that the diagonal elements (\eqn{i},\eqn{i},\eqn{m})'s for all rankers should be set to NA as well.
#' @param X.mat An \eqn{N} by \eqn{L} covariate matrix for the \eqn{N} entities with \eqn{L} covariates.
#' @param tau2.alpha The scale parameter for the scaled inverse chi-squared prior on \eqn{\sigma^2_alpha}.
#' @param nu.alpha The d.f. for the scaled inverse chi-squared prior on \eqn{\sigma^2_alpha}.
#' @param tau2.beta The scale parameter for the scaled inverse chi-squared prior on \eqn{\sigma^2_beta}.
#' @param nu.beta The d.f. for the scaled inverse chi-squared prior on \eqn{\sigma^2_beta}.
#' @param gamma.a The shape parameter for the Gamma prior on the concentration parameter of Dirichlet process.
#' @param gamma.b The rate parameter for the Gamma prior on the concentration parameter of Dirichlet process.
#' @param weight.prior.value A vector for the support of the discrete prior on weight parameter.
#' @param weight.prior.prob A vector for the probability mass of the discrete prior on weight parameter.
#' @param iter.max Number of iterations for Gibbs sampler.
#' @param para.expan Logical variable for whether using parameter expansion in the Gibbs sampler.
#' @return A list containing posterior samples of all the missing evaluation scores for all rankers and all the model parameters.
#' @export
BayesRankCovMixWeight <- function(pair.comp.ten, X.mat = matrix(NA, nrow =dim(pair.comp.ten)[1], ncol = 0),
                                  tau2.alpha = 5^2, nu.alpha = 3,
                                  tau2.beta = 5^2, nu.beta = 3,
                                  gamma.a = 2, gamma.b = NULL,
                                  weight.prior.value = c(0.5, 1, 2), weight.prior.prob = rep(1/length(weight.prior.value), length(weight.prior.value)),
                                  n.item = dim(pair.comp.ten)[1], n.ranker = dim(pair.comp.ten)[3], p.cov = ncol(X.mat),
                                  iter.max = 5000, para.expan = TRUE, print.opt = 100,
                                  initial.list = NULL){
  ## store MCMC draws
  draw = list(
    Z.mat = array(NA, dim = c(n.item, n.ranker, iter.max)),
    alpha = array(NA, dim = c(n.item, n.ranker, iter.max)),
    beta = array(NA, dim = c(p.cov, n.ranker, iter.max)),
    mu = array(NA, dim = c(n.item, n.ranker, iter.max)),
    cluster = array(NA, dim = c(n.ranker, iter.max)),
    weight.vec = array(NA, dim = c(n.ranker, iter.max) ),
    sigma2.alpha = rep(NA, iter.max),
    sigma2.beta = rep(NA, iter.max),
    gamma = rep(NA, iter.max)
  )

  if(is.null(initial.list)){
    ## initial values for Z
    Z.mat = matrix(NA, nrow = n.item, ncol = n.ranker)
    for(j in 1:n.ranker){
      Z.mat[sort( rowSums( pair.comp.ten[,,j], na.rm = TRUE ), decreasing = FALSE, index.return = TRUE )$ix, j] = (c(n.item : 1) - (1+n.item)/2)/sd(c(n.item : 1))
    }

    ## initial values for alpha, beta and thus mu
    alpha = array(0, dim = c(n.item, n.ranker))
    beta = array(0, dim = c(p.cov, n.ranker))
    mu = alpha + X.mat %*% beta

    ## initial value for cluster indicator
    # cluster = rep(1, n.ranker)
    cluster = sample(c(1:ceiling(log(n.ranker))), n.ranker, replace = TRUE )

    ## initial values for weights
    weight.vec = rep(1, n.ranker)

    ## initial values for hyper para of DP
    sigma2.alpha = tau2.alpha
    sigma2.beta = tau2.beta
    if(is.null(gamma.b)){
      gamma = gamma.a
    }else{
      gamma = gamma.a/gamma.b
    }
  }else{

    Z.mat = initial.list$Z.mat
    alpha = initial.list$alpha
    beta = initial.list$beta
    mu = alpha + X.mat %*% beta
    cluster = initial.list$cluster
    weight.vec = initial.list$weight.vec
    sigma2.alpha = initial.list$sigma2.alpha
    sigma2.beta = initial.list$sigma2.beta
    gamma = initial.list$gamma

  }



  ## store initial value
  draw$Z.mat[,,1] = Z.mat
  draw$alpha[,,1] = alpha
  draw$beta[,,1] = beta
  draw$mu[,,1] = mu
  draw$cluster[, 1] = cluster
  draw$weight.vec[, 1] = weight.vec
  draw$sigma2.alpha[1] = sigma2.alpha
  draw$sigma2.beta[1] = sigma2.beta
  draw$gamma[1] = gamma

  c.max = max(cluster)

  ## Gibbs iteration
  for(iter in 2:iter.max){
    ### Gibbs update for each clsuter
    alpha.all = c()
    beta.all = c()
    for(c in unique( cluster ) ){
      # update (alpha, beta) or equivalently mu
      mean.para.update = GibbsUpMuGivenLatentGroup(Z.mat = Z.mat[, cluster == c, drop=FALSE], X.mat = X.mat, weight.vec = weight.vec[cluster==c], sigma2.alpha = sigma2.alpha, sigma2.beta = sigma2.beta, n.ranker = sum(cluster==c), n.item = n.item, p.cov = p.cov, para.expan = para.expan)

      alpha.c = mean.para.update$alpha
      beta.c = mean.para.update$beta
      mu.c = as.vector( alpha.c + X.mat %*% beta.c )

      alpha[, cluster==c] = alpha.c
      beta[, cluster==c] = beta.c
      mu[, cluster==c] = mu.c

      alpha.all = c(alpha.all, alpha.c)
      beta.all = c(beta.all, beta.c)

      ### for check only
      Z.mat[, cluster == c] = Z.mat[, cluster == c]/mean.para.update$theta

      # update Z.mat given (alpha, beta) or equivalently mu
      # mu.c = mu[, cluster==c, drop=FALSE][,1]
      Z.mat[, cluster == c] = GibbsUpLatentGivenRankGroup(pair.comp.ten = pair.comp.ten[,,cluster==c, drop=FALSE], Z.mat = Z.mat[, cluster==c, drop=FALSE], mu = mu.c, weight.vec = weight.vec[cluster==c], n.ranker = sum(cluster==c) )
    }

    ### Gibbs update for weights
    for(i in 1:n.ranker){
      weight.vec[i] = GibbsUpWeightInd(Z = Z.mat[, i], mu[, i], weight.prior.value = weight.prior.value, weight.prior.prob = weight.prior.prob, n.item = n.item )
    }

    ### Gibbs update for hyper parameters of DP
    sigma2.alpha = GibbsUpsigma2(alpha.all, nu.alpha, tau2.alpha)
    if(p.cov > 0){
      sigma2.beta = GibbsUpsigma2(beta.all, nu.beta, tau2.beta)
    }
    if(!is.null(gamma.b)){
      gamma = GibbsUpDPgamma( gamma = gamma, c.vec = cluster, a = gamma.a, b = gamma.b, n = n.ranker )
    }

    ### Gibbs update for the cluster label
    # for(i in 1:n.ranker){
    for(i in sample(c(1:n.ranker)) ){
      candidate = c( unique( cluster[-i] ), c.max+1 )
      logp = rep(NA, length(candidate))
      for(c in candidate){
        if( sum(cluster[-i] == c) > 0 ){
          set = setdiff( which(cluster == c), i )
          logp[candidate==c] = log( sum(cluster[-i] == c) ) + LogMargDensity(Z.mat = Z.mat[, c(set, i), drop=FALSE], X.mat = X.mat, sigma2.alpha = sigma2.alpha, sigma2.beta = sigma2.beta, weight.vec = weight.vec[c(set, i)]) - LogMargDensity(Z.mat = Z.mat[, set, drop=FALSE], X.mat = X.mat, sigma2.alpha = sigma2.alpha, sigma2.beta = sigma2.beta, weight.vec = weight.vec[set])
        }else{
          logp[candidate==c] = log(gamma) + LogMargDensity(Z.mat = Z.mat[, i, drop=FALSE], X.mat = X.mat, sigma2.alpha = sigma2.alpha, sigma2.beta = sigma2.beta, weight.vec = weight.vec[i])
        }
      }
      logp = logp - max(logp)
      cluster[i] = candidate[which( as.vector( rmultinom(1, 1, prob = exp(logp)) ) == 1 )]
      if(cluster[i] == ( c.max+1)){
        c.max = c.max + 1
      }
    }

    ### store vlaues
    draw$Z.mat[,,iter] = Z.mat
    draw$alpha[,,iter] = alpha
    draw$beta[,,iter] = beta
    draw$mu[,,iter] = mu
    draw$cluster[, iter] = cluster
    draw$weight.vec[, iter] = weight.vec
    draw$sigma2.alpha[iter] = sigma2.alpha
    draw$sigma2.beta[iter] = sigma2.beta
    draw$gamma[iter] = gamma

    # print iteration number
    if(iter %% print.opt == 0){
      print(paste("Gibbs Iteration", iter))
      # print(table(cluster))
      # print(table(weight.vec))
      # print(c(sigma2.alpha, sigma2.beta))
      # print(gamma)
    }
  }
  return(draw)
}


