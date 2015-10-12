\name{reliability}
\alias{reliability}
\title{
calculate reliability-index
}
\description{
calculate reliability-index of Elo-ratings
}
\usage{
reliability(x)
}
\arguments{
  \item{x}{elochoice-object, the result of \code{\link{elochoice}}}
}

\value{
a data.frame with as many rows as randomizations were run in the original call to \code{elochoice()}. The first column represents the unweighted and the second the weighted reliability index (\emph{R} and \emph{R'}), which is followed by the total number of trials that contributed to the calculation of the index. Note that this number cannot reach the total number of trials in the data set because at least for the very first trial we did not have an expectation for the outcome of that trial (and such trials do not contribute to the calculation of the reliability index).
}
\references{
Clark, A. P., Neumann, C., Howard, K. L., Woods, A. T. & Penton-Voak, I. S. (in preparation). Why rate when you could compare? Using the ``EloChoice'' package to assess pairwise comparisons of perceived strength.
}
\author{Christof Neumann}
\examples{
# create data set and calculate ratings (with five randomizations)
xdata <- randompairs(12, 500)
x <- elochoice(xdata$winner, xdata$loser, runs=5)
# extract the reliability values
(u <- reliability(x))
# calculate average reliability index
mean(u$upset)
# and in its weighted form
mean(u$upset.wgt)
}
