# created on Dec. 15, 2015
#  generate data from the simulation algorithm 
#  used by Reiser (2000). Each subject has exactly 2 observations

# X_{ik, obs}=X_{i,true}+\epsilon_{ik},
#   X_{i, true} ~ N(mu.X, sigma.X2)
#   \epsilon_{ik} ~ N(0, sigma.epsilon2)
#   i=1,\ldots, nX, k=1, 2
#
# Y_{jl, obs}=Y_{j,true}+\xi{jl}
#   Y_{j, true} ~ N(mu.Y, sigma.Y2)
#   \xi_{jl} ~ N(0, sigma.eta2)
#   j=1,\ldots, nY, l=1, 2

genSimDataReiser <-
function (nX = 100, nY = 100, sigma.X2 = 1, mu.X = 0.25, sigma.Y2 = 1, 
    mu.Y = 0, sigma.epsilon2 = 0.5, sigma.eta2 = 0.5) 
{
    # calculate AUC.true
    theta2 = (sigma.epsilon2 + sigma.eta2)/(sigma.X2 + sigma.Y2)
    mu.true = mu.X - mu.Y
    delta = mu.true/sqrt(sigma.X2 + sigma.Y2)
    AUC.true = pnorm(delta)

    # simulate 1st observation for a subject
    Xvec = rnorm(nX, mean = mu.X, sd = sqrt(sigma.X2))
    Yvec = rnorm(nY, mean = mu.Y, sd = sqrt(sigma.Y2))
    epsilon = rnorm(nX, mean = 0, sd = sqrt(sigma.epsilon2))
    eta = rnorm(nY, mean = 0, sd = sqrt(sigma.eta2))
    x.obs = Xvec + epsilon
    y.obs = Yvec + eta

    # simulate 2nd observation  for a subject
    epsilon = rnorm(nX, mean = 0, sd = sqrt(sigma.epsilon2))
    eta = rnorm(nY, mean = 0, sd = sqrt(sigma.eta2))
    x.obs.rep = Xvec + epsilon
    y.obs.rep = Yvec + eta

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
 
    res <- list(datFrame=datFrame, theta2 = theta2, mu.true = mu.true, 
        AUC.true = AUC.true)

    invisible(res)
}


