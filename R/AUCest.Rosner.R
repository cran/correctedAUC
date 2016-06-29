# created on Dec. 15, 2015
#  (1) obtain the estimate of corrected AUC and its confidence interval
#
# datFrame - data frame with at least 4 columns
#   'y' -- observations
#   'subjID' -- subject id
#   'grp' -- group id: 1 - case; 0 - control;
#   'myrep' -- replication id: 1, 2, 3...
AUCest.Rosner <- function (
  datFrame, 
  sidVar = "subjID", 
  obsVar = "y", 
  grpVar = "grp", 
  repVar = "myrep", 
  alpha = 0.05) 
{
    myrep=datFrame[, c(repVar)]
    grp = datFrame[, c(grpVar)]
    y = datFrame[, c(obsVar)]
    subjID = datFrame[, c(sidVar)]

    datFrame$myrep=myrep
    datFrame$grp=grp
    datFrame$y=y
    datFrame$subjID=subjID

    # use first observation to esstimate AUC
    dat1=datFrame[which(datFrame$myrep==1),]
    x.obs=dat1$y[which(dat1$grp==1)]
    y.obs=dat1$y[which(dat1$grp==0)]

    mu.mle = mu.mle.func(x = x.obs, y = y.obs, verbose = FALSE)$mu.mle
    res.MW = getMannWhitney(x = x.obs, y = y.obs)
    AUC.obs = res.MW$AUC
    var.AUC.obs = res.MW$varAUC.alt

    u.rep=sort(unique(datFrame$myrep))
    k0=max(u.rep, na.rm=TRUE) # ??

    # use the all observations to estimate ICC
    outcome.x=NULL
    outcome.y=NULL
    clusterid.x=NULL
    clusterid.y=NULL

    for(i in 1:k0)
    {
      dati=datFrame[which(datFrame$myrep==i),]
      x.obs.i=dati$y[which(dati$grp==1)]
      y.obs.i=dati$y[which(dati$grp==0)]

      z.obs.i = c(x.obs.i, y.obs.i)
      nX.i = length(x.obs.i)
      nY.i = length(y.obs.i)
      if(i==1)
      {
        nX=nX.i
        nY=nY.i
      }


      N.i = nX.i + nY.i
      H.obs.i = qnorm(rank(z.obs.i)/(N.i + 1))
      H.obs.x.i = H.obs.i[1:nX.i]
      H.obs.y.i = H.obs.i[-c(1:nX.i)]

      outcome.x=c(outcome.x, H.obs.x.i)
      clusterid.x=c(clusterid.x, 1:nX.i)

      outcome.y=c(outcome.y, H.obs.y.i)
      clusterid.y=c(clusterid.y, 1:nY.i)
    }

    ICC.x = ICCest(x = as.factor(clusterid.x), y = outcome.x)$ICC
    ICC.y = ICCest(x = as.factor(clusterid.y), y = outcome.y)$ICC

    tt1 = sqrt((1/ICC.x + 1/ICC.y)/2)
    tt2 = qnorm(AUC.obs) * tt1
    AUC.c = pnorm(tt2)
    za = qnorm(1 - alpha/2)
    a = qnorm(AUC.obs)
    b = tt1
    var.a = var.AUC.obs/(dnorm(a))^2
    a.low = a - za * sqrt(var.a)
    a.upp = a + za * sqrt(var.a)
    AUC.obs.low = pnorm(a.low)
    AUC.obs.upp = pnorm(a.upp)

    var.ICC.x = 2 * (1 - ICC.x)^2 * (1 + (k0 - 1) * ICC.x)^2
    var.ICC.x = var.ICC.x/(k0 * (k0 - 1) * (nX - 1))
    var.ICC.y = 2 * (1 - ICC.y)^2 * (1 + (k0 - 1) * ICC.y)^2
    var.ICC.y = var.ICC.y/(k0 * (k0 - 1) * (nX - 1))
    ph.pICC.x = -1/(2 * sqrt(2) * ICC.x^2 * sqrt(1/ICC.x + 1/ICC.y))
    ph.pICC.y = -1/(2 * sqrt(2) * ICC.y^2 * sqrt(1/ICC.x + 1/ICC.y))
    var.b = ph.pICC.x^2 * var.ICC.x + ph.pICC.y^2 * var.ICC.y
    ab = a * b
    se.ab = sqrt(b^2 * var.a + a^2 * var.b)
    c1 = ab - za * se.ab
    c2 = ab + za * se.ab
    AUC.c.low = pnorm(c1)
    AUC.c.upp = pnorm(c2)

    res = list(AUC.obs = AUC.obs, AUC.c = AUC.c, 
        ICC.x = ICC.x, ICC.y = ICC.y, mu.mle = mu.mle, AUC.obs.low = AUC.obs.low, 
        AUC.obs.upp = AUC.obs.upp, AUC.c.low = AUC.c.low, AUC.c.upp = AUC.c.upp)
    return(res)
}

##########################
# loglikelihood function of y sample
loglkh.func=function(mu, x, y)
{
  #x.s=sort(x[which(is.na(x)==FALSE)])
  x=x[which(is.na(x)==FALSE)]
  y=y[which(is.na(y)==FALSE)]
  x.s=sort(x)
  len=length(x)

  # value of ecdf of x sample at sorted x
  Fx.x=ecdf.func(x=x, z=x.s)
  H=qnorm(Fx.x)

  start=1
  end=len-1

  n1=sum(is.na(y)==FALSE & y<x.s[1], na.rm=TRUE)
  if(n1)
  {
    logLy=n1*log( pnorm( H[1] + mu ) )
  } else {
    logLy=0
  }

  #cat("logLy=", logLy, "\n")
  for(i in start:end)
  {
    n.i=sum(is.na(y)==FALSE & y>=x.s[i] & y <x.s[i+1])
    #cat("n.i=", n.i)
    if(is.na(n.i)==FALSE && n.i>0)
    {
      logLy=logLy+n.i*log( pnorm( H[i+1] + mu ) - pnorm( H[i] + mu ) )
    }
  }

  n2=sum(is.na(y)==FALSE & y>=x.s[len], na.rm=TRUE)
  if(n2)
  {
    logLy = logLy + n2*log( 1 - pnorm( H[len] + mu ) )
  }

  return(logLy)
}


##############
# first order derivative of loglikelihood function of y sample
dloglkh.func <-
function (mu, x, y) 
{
    x = x[which(is.na(x) == FALSE)]
    y = y[which(is.na(y) == FALSE)]
    x.s = sort(x)
    len = length(x)
    Fx.x = ecdf.func(x = x, z = x.s)
    H = qnorm(Fx.x)
    start = 1
    end = len - 1
    n1 = sum(is.na(y) == FALSE & y < x.s[1], na.rm = TRUE)
    dlogLy = n1 * dnorm(H[1] + mu)/pnorm(H[1] + mu)
    for (i in start:end) {
        n.i = sum(is.na(y) == FALSE & y >= x.s[i] & y < x.s[i + 
            1])
        dlogLy = dlogLy + n.i * (dnorm(H[i + 1] + mu) - dnorm(H[i] + 
            mu))/(pnorm(H[i + 1] + mu) - pnorm(H[i] + mu))
    }
    n2 = sum(is.na(y) == FALSE & y >= x.s[len], na.rm = TRUE)
    dlogLy = dlogLy + n2 * (-dnorm(H[len] + mu))/(1 - pnorm(H[len] + 
        mu))
    return(dlogLy)
}

# calculate ecdf for z based on x
#  modified the formula to avoid zero and one
ecdf.func=function(x, z, const=100)
{
  x.s=sort(x)
  len=length(x)
  start=1
  end=len-1

  cdf.z=rep(NA, length(z))
  pos1=which(is.na(z)==FALSE & z<x.s[1])
  if(length(pos1))
  {
    cdf.z[pos1]=1/(len*const)
  }

  pos2=which(is.na(z)==FALSE & z>=x.s[len])
  if(length(pos2))
  {
    cdf.z[pos2]=  1-1/(len*const)
  }

  for(i in start:end)
  {
    pos.i=which(is.na(z)==FALSE & z>=x.s[i] & z <x.s[i+1])
    cdf.z[pos.i]=i/len
  }

  return(cdf.z)

}

# mle of mu
mu.mle.func=function(x, y, alpha=0.05, verbose=TRUE)
{
  x=x[which(is.na(x)==FALSE)]
  y=y[which(is.na(y)==FALSE)]

  # initial estimate of mu
  mu.ini=est.mu.func(x, y)

  x.min=min(x, na.rm=TRUE)
  x.max=max(x, na.rm=TRUE)

  y.min=min(y, na.rm=TRUE)
  y.max=max(y, na.rm=TRUE)

  if(x.min>y.max || x.max < y.min)
  {
    return(mu.ini)
  }

  # to do maximization, set fnscale < 0
  res.optim=optim(par=mu.ini, fn=loglkh.func, gr=dloglkh.func, method="BFGS", 
    control=list(fnscale= -1), x=x, y=y, hessian =TRUE)

  # observed information  = - hessian matrix
  info= - res.optim$hessian
  v= 1/info
  s=sqrt(v)

  mu.mle=res.optim$par
  # 95% CI for mu
  za=qnorm(1-alpha/2)

  ###
  mu.L=mu.mle-za*s
  mu.U=mu.mle+za*s

  # 95% CI for AUC
  AUC.mle=pnorm(mu.mle/sqrt(2))
  AUC.L=pnorm(mu.L/sqrt(2))
  AUC.U=pnorm(mu.U/sqrt(2))

  if(verbose)
  {
   
    print(res.optim)

    cat("\nmu.ini=", mu.ini, ", mu.mle=", mu.mle, "\n")
    cat("\n95% CI for mu>>\n")
    print(c(mu.L, mu.U))

    cat("\nAUC.mle=", AUC.mle, "\n")
    cat("\n95% CI for AUC>>\n")
    print(c(AUC.L, AUC.U))
  }

  res=list(mu.ini=mu.ini, mu.mle=mu.mle, AUC.mle=AUC.mle, mu.CI=c(mu.L, mu.U), AUC.CI=c(AUC.L, AUC.U))
  return(res)
}

# estimate mu
est.mu.func=function(x, y)
{
  x=x[which(is.na(x)==FALSE)]
  y=y[which(is.na(y)==FALSE)]
  z=c(x,y)
  # ecdf based on x
  ecdf.x=ecdf.func(x=x, z=z)
  # ecdf based on y
  ecdf.y=ecdf.func(x=y, z=z)

  mu=mean(qnorm(ecdf.y) - qnorm(ecdf.x), na.rm=TRUE)
  return(mu)
}

# Pr(Y<X) assuming no tie
getMannWhitney=function(x, y)
{
  nX=length(x)
  nY=length(y)

  z=c(y, x)
  rz=rank(z)
  Ty=sum(rz[1:nY])

  U=nX*nY+nY*(nY+1)/2-Ty

  AUC=U/(nX*nY)

  ## under the null hypothesis
  EU=nX*nY/2
  varU=nX*nY*(nX+nY+1)/12
  z=(U-EU)/sqrt(varU)
  pval=1-pnorm(abs(z))+pnorm(-abs(z))

  # under alternative
  theta=AUC
  qn.theta=qnorm(AUC)

  part1=theta*(1-theta)+(nX+nY-2)
  part2=mnormt::pmnorm( x=c(qn.theta, qn.theta), mean=c(0,0), varcov=matrix(c(1,0.5, 0.5, 1), 2,2) )
  varAUC.alt=part1*(part2-theta^2)/(nX*nY)
  
  res=list(U=U, AUC=AUC, EU=EU, varU=varU, z=z, pval=pval, varAUC.alt=varAUC.alt)
  
  return(res)
}


