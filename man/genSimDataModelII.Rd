\name{genSimDataModelII}
\alias{genSimDataModelII}
\title{
Generate one simulated data set based on Model II in Rosner et al's (2015) manuscript
}
\description{
Generate one simulated data set based on Model II in Rosner et al's (2015) manuscript.
}
\usage{
genSimDataModelII(
  nX, 
  nY, 
  mu,
  lambda,
  sigma.X2, 
  sigma.Y2, 
  sigma.e.X, 
  sigma.e.Y) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nX}{
integer. number of cases.
}
  \item{nY}{
integer. number of controls.
}
  \item{mu}{
difference of means between the case distribution and control distribution.
}
  \item{lambda}{
mean for controls.
}
  \item{sigma.X2}{
variance of the true value for cases.
}
  \item{sigma.Y2}{
variance of the true value for controls.
}
  \item{sigma.e.X}{
variance of the random error term for cases.
}
  \item{sigma.e.Y}{
variance of the random error term for controls.
}
}
\details{
The Model II in Rosner et al.'s (2005) manuscript:
\deqn{
X_{ik, obs}=X_{i,true}+\epsilon_{ik},\\
   \log\left(X_{i, true}\right) \sim N\left(\lambda+\mu, \sigma_X^2\right),\\
   \epsilon_{ik} \sim N\left(0, \sigma_{\epsilon}^2\right),\\
   i=1,\ldots, n_X, k=1, 2
}{%
 X_{ik, obs}=X_{i,true}+\epsilon_{ik},
   log(X_{i, true}) ~ N(lambda+mu, sigma.X2),
   \epsilon_{ik} ~ N(0, sigma.epsilon2),
   i=1,\ldots, nX, k=1, 2
}
\deqn{
 Y_{jl, obs}=Y_{j,true}+\xi_{jl},\\
   \log\left(Y_{j, true}\right) \sim N\left(\lambda, \sigma_Y^2\right),\\
   \xi_{jl} \sim N(0, \sigma_{\eta}^2),\\
   j=1,\ldots, n_Y, l=1, 2
}{%
 Y_{jl, obs}=Y_{j,true}+\xi{jl},
   log(Y_{j, true}) ~ N(lambda, sigma.Y2),
   \xi_{jl} ~ N(0, sigma.eta2),
   j=1,\ldots, nY, l=1, 2
}
}
\value{
A list of 2 elements: 
\item{datFrame}{A data frame with 4 elements: 
    \code{y}: observations;
    \code{subjID}: subject ID;
    \code{grp}: group indicator;
    \code{myrep}: replication indicator.
}
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
  tt=genSimDataModelII(
    nX = 100, 
    nY = 100, 
    mu = 0.25,
    lambda = 0,
    sigma.X2 = 1, 
    sigma.Y2 = 1, 
    sigma.e.X = 1, 
    sigma.e.Y = 1) 

  print(dim(tt$datFrame))
  print(tt$datFrame[1:2,1:3])
  print(tt$AUC.true)

}
\keyword{ method }
