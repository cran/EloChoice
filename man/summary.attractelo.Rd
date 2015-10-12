\name{summary.elochoice}
\alias{summary.elochoice}
\title{
summarize elochoice object
}
\description{
overview of elochoice object
}
\usage{
\method{summary}{elochoice}(object, \dots)
}
\arguments{
  \item{object}{
    an object of class \code{"elochoice"}, usually the result of a call to \code{\link{elochoice}}
}  \item{\dots}{
    further arguments passed to or from other methods. Nothing relevant in this case.
}
}
\author{
Christof Neumann
}
\examples{
xdata <- randompairs(nstim=10, nint=500)
x <- elochoice(xdata$winner, xdata$loser, runs=5)
summary(x)
}
