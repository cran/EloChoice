\name{randompairs}
\alias{randompairs}
\title{
generate random data of pairwise preference ratings
}
\description{
generate random data of pairwise preference ratings
}
\usage{
randompairs(nstim=10, nint=100, reverse=0.1, skew=FALSE)
}
\arguments{
  \item{nstim}{number of stimuli, must be less than 2,602}
  \item{nint}{number of paired ratings to be created}
  \item{reverse}{proportion of ratings that go against the default preference, see below for details}
  \item{skew}{logical, by default \code{FALSE}, see below for details}
}
\details{
The default preference for a given pair is given by their alphanumerical order. E.g. \emph{A} is preferred over \emph{M}, and \emph{kf} over \emph{kz}. The \code{reverse=} argument specifies the proportion of ratings that go against this default order.

The number of appearances of a given stimulus in the data set is by default determined by uniform sampling of individual stimuli, i.e. all stimuli will roughly appear equally often in a data set. If a somewhat more realistic (i.e. unbalanced) distribution is desired, the argument \code{skew=TRUE} will achieve sampling based on a negative binomial distribution.
}
\value{
data.frame with winner and loser column. An additional column (\code{index}) serves as an index for the sequence in which the trials occurred.
}
\author{
Christof Neumann
}
\examples{
# a relatively balanced data set
xdata <- randompairs(20, 500, skew=FALSE)
table(c(as.character(xdata$winner), as.character(xdata$loser)))
range(table(c(as.character(xdata$winner), as.character(xdata$loser))))
# and a less balanced data set
xdata <- randompairs(20, 500, skew=TRUE)
table(c(as.character(xdata$winner), as.character(xdata$loser)))
range(table(c(as.character(xdata$winner), as.character(xdata$loser))))
}
