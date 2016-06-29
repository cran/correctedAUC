\name{genSimDataReiser}
\alias{genSimDataReiser}
\title{
Generate one simulated data set based on Reiser's (2000) model
}
\description{
Generate one simulated data set based on Reiser's (2000) model. The true
AUC will also be calculated.
}
\usage{
genSimDataReiser(
  nX = 100, 
  nY = 100, 
  sigma.X2 = 1, 
  mu.X = 0.25, 
  sigma.Y2 = 1, 
  mu.Y = 0, 
  sigma.epsilon2 = 0.5, 
  sigma.eta2 = 0.5) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nX}{
integer. number of cases.
}
  \item{nY}{
integer. number of controls.
}
  \item{sigma.X2}{
variance of the true value for cases.
}
  \item{mu.X}{
mean of the true value for cases.
}
  \item{sigma.Y2}{
variance of the true value for controls.
}
  \item{mu.Y}{
mean of the true value for controls.
}
  \item{sigma.epsilon2}{
variance of the random error term for cases.
}
  \item{sigma.eta2}{
variance of the random error term for controls.
}
}
\details{
Reiser's (2000) measurement error model is:
\deqn{
X_{ik, obs}=X_{i,true}+\epsilon_{ik},\\
   X_{i, true} \sim N\left(\mu_X, \sigma_X^2\right),\\
   \epsilon_{ik} \sim N\left(0, \sigma_{\epsilon}^2\right),\\
   i=1,\ldots, n_X, k=1, 2
}{%
 X_{ik, obs}=X_{i,true}+\epsilon_{ik},
   X_{i, true} ~ N(mu.X, sigma.X2),
   \epsilon_{ik} ~ N(0, sigma.epsilon2),
   i=1,\ldots, nX, k=1, 2
}
\deqn{
 Y_{jl, obs}=Y_{j,true}+\xi_{jl},\\
   Y_{j, true} \sim N\left(\mu_Y, \sigma_Y^2\right),\\
   \xi_{jl} \sim N(0, \sigma_{\eta}^2),\\
   j=1,\ldots, n_Y, l=1, 2
}{%
 Y_{jl, obs}=Y_{j,true}+\xi{jl},
   Y_{j, true} ~ N(mu.Y, sigma.Y2),
   \xi_{jl} ~ N(0, sigma.eta2),
   j=1,\ldots, nY, l=1, 2
}
}
\value{
A list of 4 elements: 
\item{datFrame}{A data frame with 4 elements: 
    \code{y}: observations;
    \code{subjID}: subject ID;
    \code{grp}: group indicator;
    \code{myrep}: replication indicator.
}
\item{theta2 }{\eqn{(sigma.epsilon2 + sigma.eta2)/(sigma.X2 + sigma.Y2)}}
\item{mu.true }{\eqn{mu.X-mu.Y}}
\item{AUC.true }{true AUC value}

}
\references{
  Rosner B, Tworoger S, Qiu W (2015)
  Correcting AUC for Measurement Error.
  J Biom Biostat 6:270. doi:10.4172/2155-6180.1000270
}
\author{
Bernard Rosner <stbar@channing.harvard.edu>,
Shelley Tworoger <nhsst@channing.harvard.edu>,
Weiliang Qiu <stwxq@channing.harvard.edu>
}
%\note{
%}
%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
  set.seed(1234567)
  tt=genSimDataReiser(
       nX = 100, 
       nY = 100, 
       sigma.X2 = 1, 
       mu.X = 0.25, 
       sigma.Y2 = 1, 
       mu.Y = 0, 
       sigma.epsilon2 = 0.5, 
       sigma.eta2 = 0.5) 

  print(dim(tt$datFrame))
  print(tt$datFrame[1:2,1:3])
  print(tt$theta2)
  print(tt$mu.true)
  print(tt$AUC.true)

}
\keyword{ method }
