\name{AUCest.Rosner}
\alias{AUCest.Rosner}
\title{
Calculate AUC.c for measurement error based on probit-shift model
}
\description{
Calculate AUC.c for measurement error based on probit-shift model.
}
\usage{
AUCest.Rosner(
  datFrame, 
  sidVar = "subjID", 
  obsVar = "y", 
  grpVar = "grp", 
  repVar = "myrep", 
  alpha = 0.05)
}
\arguments{
  \item{datFrame}{
a data frame with at least the following columns:
    \code{y}: numerical vector of observations;
    \code{subjID}: vector of subject ids;
    \code{grp}: group indicator: 1 means case; and 0 means conrol;
    \code{myrep}: integer vector indicating replication. should
be consecutive positive integer starting from 1.
}
  \item{sidVar}{character. variable name for subject id in the data frame \code{datFrame}.}
  \item{obsVar}{character. variable name for observations in the data frame \code{datFrame}.}
  \item{grpVar}{character. variable name for group indictor in the data frame \code{datFrame}.}
  \item{repVar}{character. variable name for replication indictor in the data frame \code{datFrame}.}
  \item{alpha}{
confidence interval level \eqn{100(1-\alpha}\%}
}
%\details{
%
%}
\value{
A list of  9 elements:
\item{AUC.obs}{AUC estimated based on the Mann-Whitney statistic.}
\item{AUC.c }{AUC corrected for measurement error based on
the probit-shift model.}
\item{ICC.x}{intra-class correlation for cases.}
\item{ICC.y}{intra-class correlation for controls}
\item{mu.mle}{maximum likelihood estimate of \eqn{\mu} (i.e., the 
shift between the case distribution and the control distribution 
after probit transformation)}
\item{AUC.obs.low}{lower bound of the \code{AUC.obs}.}
\item{AUC.obs.upp}{upper bound of the \code{AUC.obs}.}
\item{AUC.c.low}{lower bound of the \code{AUC.c}.}
\item{AUC.c.upp}{upper bound of the \code{AUC.c}.}
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
%%  ~~further notes~~
%}

%\seealso{
%}
\examples{
  set.seed(1234567)
  tt=genSimDataModelIII(
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

  res = AUCest.Rosner(
    datFrame = tt$datFrame, 
    sidVar = "subjID",
    obsVar = "y",
    grpVar = "grp",
    repVar = "myrep",
    alpha = 0.05)
  print(res)

}
\keyword{ method }
