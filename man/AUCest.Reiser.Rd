\name{AUCest.Reiser}
\alias{AUCest.Reiser}
\title{
Calculate AUC corrected for measurement error based on Reiser's (2000) method
}
\description{
Calculate AUC corrected for measurement error based on Reiser's (2000) method.
}
\usage{
AUCest.Reiser(
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
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
A list of 4 elements
\item{AUC.c }{AUC corrected for measurement error based on
Reiser's (2000) method.}
\item{sd.AUC.c }{standard error of the estimated 
AUC corrected for measurement error based on
Reiser's (2000) method.}
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

  res = AUCest.Reiser(
    datFrame = tt$datFrame, 
    sidVar = "subjID",
    obsVar = "y",
    grpVar = "grp",
    repVar = "myrep",
    alpha = 0.05)
  print(res)

}
\keyword{ method }
