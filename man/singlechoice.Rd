\name{singlechoice}
\alias{singlechoice}
\title{update stimulus ratings after one rating event}
\description{update stimulus ratings after one rating event}
\usage{
singlechoice(val1, val2, k)
}
\arguments{
  \item{val1}{rating of the preferred stimulus \emph{before} the rating event
}
  \item{val2}{rating of the unpreferred stimulus \emph{before} the rating event
}
  \item{k}{value of \emph{k}-constant, which determines the maximum change of ratings after a single rating event
}
}
\value{
vector with two values: updated ratings \emph{after} the rating event for preferred and unpreferred stimulus
}
\references{
Elo, A. E. 1978. The Rating of Chess Players, Past and Present. New York: Arco.
}
\author{Christof Neumann}
\seealso{
\code{\link[EloRating:e.single]{EloRating}}
}
\examples{
# little change because rating difference is large (positive), i.e. expectation is clear
singlechoice(1200, 500, 100)
# no change because rating difference is very large (positive), i.e. expectation is clear
singlechoice(1500, 500, 100)
# large change because rating difference is small (negative), i.e. expectation is clearly violated
singlechoice(500, 1500, 100)
}

