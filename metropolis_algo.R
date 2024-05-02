################# Metropolis algorithms ########################
# # # # # # # # # # # # # # # MODEL # # # # # # # # # # # # # # # # # # # #
# y : vector of predicted values (size n) #
# var : dataframe of explanatory variables (with intercept) (size : n*k) #
# X : matrix of n explanatory variables (size n*k) #
# theta : vector coefficients (size k) #
# sd : standard deviation of errors (additive errors) #
# param : c(theta, sd) (size : k+1) #
# # # # # # # # total model : y ~ normal ( X %*% theta , sd^2 ) # # # # #

likelihood <- function(y, X, param){
  theta = param[-length(param)]
  sd = param[length(param)]
  pred = X %*% theta
  singlelikelihoods = dnorm(y, mean = pred, sd = sd, log = T)
  sumll = sum(singlelikelihoods)
  return(sumll)
}
# PRIORS
# non informative priors on theta
# sd must be positive

prior <- function(param){
  theta = param[-length(param)]
  sd = param[length(param)]
  thetaprior = dnorm(theta, mean=0, sd=10, log = T)
  sdprior = dunif(sd, min=0, max=30, log = T)
  return(sum(thetaprior)+sdprior)
}

posterior <- function(y, X, param){
  return (likelihood(y, X, param) + prior(param))
}

# proposal = random walk : normal distribution
# if sdprop is too low : you won't be exploring all your parameters multidimentional space
# if sdprop is too high : propositions will be rejected too often

proposalfunction <- function(param, sdprop){
  return(rnorm(length(param) , mean = param, sd= sdprop))
}

run_metropolis_MCMC <- function(y, X, startvalue, iterations, sdprop){
  # chain : row i = set of parameters at iteration i
  chain = array(dim = c(iterations+1, length(startvalue)))
  colnames(chain) <- names(startvalue)
  lp = array(dim=iterations+1)
  chain[1,] = startvalue
  for (i in 1:iterations){
    proposal = proposalfunction(chain[i,], sdprop)
    probab = exp(posterior(y, X, proposal) - posterior(y, X, chain[i,]))
    # probab ranges between 0 and 1
    
    # probab > 1 <=> new proposal is more likely considering the data
    if (runif(1) < probab){
      chain[i+1,] = proposal
    }else{
      chain[i+1,] = chain[i,]
    }
    lp[i+1] = likelihood(y, X, chain[i+1,])
  }
  return(cbind(chain, lp))
}

############ Graphics ############
# Posterior distribution and trace plots
posterior_graph <- function(chain, ngraph){
  N <- dim(chain)[2]
  par(mfcol = c(2,ngraph))
  for (i in 1:N){
    hist(chain[-(1:burnIn),i],nclass=30,
         main=paste("Posterior of\n", colnames(chain)[i], "parameter"),
         xlab=NULL)
    plot(chain[-(1:burnIn),i], type = "l",
         main = paste("Posterior of\n", colnames(chain)[i], "parameter"),
         xlab=NULL, ylab=NULL)
  }
}

# Run and get param4
param4 <- function(y,x){
  y4 = log(y)
  var4 = data.frame(intercept=1,logagbt=log(x))
  # 
  # y4 = log(biomf$agb)
  # var4 = data.frame(intercept=1, logagbt=log(biomf$agbt))
  X4 = as.matrix(var4)
  
  start4 = c(0,0,1)
  names(start4) <- c(colnames(X4), "sd")
  chain4 = run_metropolis_MCMC(y=y4, X=X4, startvalue = start4, iterations = 110000,
                               sdprop = c(0.005,0.002,0.005))
  
  burnIn = 100000
  (acceptance = 1-mean(duplicated(chain4[-(1:burnIn),])))
  #posterior_graph(chain4[,1:3],3)
  
  #par(mfrow=c(1,3))
  #for (i in 1:3) acf(chain4[,i], lag.max = 500, main=colnames(chain4)[i])
  # avoid autocorrelation in the data : take 1 of every 1000 rows in chain
  # param_4 : 1001 sets of parameters for eq4 and their log likelihood (lp)
  param_4 <- data.frame(chain4[seq(burnIn+1, dim(chain4)[1], by=1000), ])
  
  return(param_4)
}
