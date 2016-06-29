# simulate data with measurement error neither from classical measurement error model nore from 
# our own model 2
# nX - number of cases
# nY - number of controls

# X.i.obs = X.i.true+e.X.i
# Y.j.obs = Y.j.true+e.Y.j
#
# log(X.i.true) ~ Normal(lambda, sigma_X^2), i=1, ..., nX
# log(Y.j.true) ~ Normal(lambda-\mu, sigma_Y^2), j=1, ..., nY
#
# log(e.X.i)~N(0, sigma.e.X^2)
# log(e.Y.i)~N(0, sigma.e.Y^2)
genSimDataModelIII=function(nX, nY, mu, lambda, sigma.X2, sigma.Y2, sigma.e.X, sigma.e.Y)
{
  x.true=exp(rnorm(nX, mean=lambda, sd=sqrt(sigma.X2)))
  y.true=exp(rnorm(nY, mean=lambda-mu, sd=sqrt(sigma.Y2)))

  e.X=exp(rnorm(nX, mean=0, sd=sigma.e.X))
  e.Y=exp(rnorm(nY, mean=0, sd=sigma.e.Y))

  x.obs=x.true+e.X
  y.obs=y.true+e.Y

  #####################
  # simulate replicated data
  #####################
  e.X.rep=exp(rnorm(nX, mean=0, sd=sigma.e.X))
  e.Y.rep=exp(rnorm(nY, mean=0, sd=sigma.e.Y))

  x.obs.rep=x.true+e.X.rep
  y.obs.rep=y.true+e.Y.rep

  #######################
  AUC.true = pnorm(mu/sqrt(sigma.X2+sigma.Y2))

  # form data frame with columns 
  #  'y' -- observations
  #  'subjID' -- subject ID
  #  'grp' -- group indicator
  #  'myrep' -- replication indicator
  nSubj=nX+nY
  datFrame=data.frame(y=c(x.obs, y.obs, x.obs.rep, y.obs.rep),
    subjID=c(seq(from=1, to=nSubj, by=1), seq(from=1, to=nSubj, by=1)),
    grp=c(rep(1, nX), rep(0, nY), rep(1, nX), rep(0, nY)),
    myrep=c(rep(1, nSubj), rep(2, nSubj)))
 

    res <- list(datFrame=datFrame, AUC.true = AUC.true)

  invisible(res)
}

