\name{elochoice}
\alias{elochoice}
\alias{eloint}
\alias{elointnorm}
\title{
Elo-ratings for pairwise comparisons of visual stimuli
}
\description{
Elo-ratings for pairwise comparisons of visual stimuli
}
\usage{
elochoice(winner, loser, kval=100, startvalue=0, runs=1, normprob=FALSE)
eloint(winner, loser, allids, kval, startvalues, runs)
elointnorm(winner, loser, allids, kval, startvalues, runs)
}
\arguments{
  \item{winner, loser}{character, vectors with the IDs of the winning (preferred) and losing (not preferred) stimuli}
  \item{kval}{numeric, k-value, which determines the maximum number of points a stimulus' rating can change after a single rating event, by default 100}
  \item{startvalue, startvalues}{numeric, start value around which ratings are centered, by default 0}
  \item{runs}{numeric, number of randomizations}
  \item{normprob}{logical, by default \code{FALSE}, which indicates a logistic approach is taken for calculating winning probabilities (see Elo 1978). Alternatively (\code{TRUE}), winning probabilities are calculated from a normal distribution}
  \item{allids}{internal, character of all stimulus IDs in the data set}
}
\details{
\code{elochoice()} is the workhorse function of the package, which wraps up all the calculations for obtaining Elo-ratings and the information for the reliability index

\code{eloint()} and \code{elointnorm()} are internal functions (which \code{elochoice()} makes use of) that do most of the calculations, but are usually not directly addressed by the user.
}
\value{
a list with the following items:
  \item{ratmat}{numeric matrix with final ratings for each stimulus, one row per randomization}
  \item{decmat}{logical matrix showing for each randomization (row) and each single rating event (column) whether or not there was an expectation for that trial, i.e. whether the two stimuli's ratings differed before the rating}
  \item{upsmat}{logical matrix showing for each randomization (row) and each single rating event (column) whether or not the outcome of a trial was in the direction of the expectation, i.e. whether or not the higher rated stimulus won}
  \item{wgtmat}{numeric matrix showing for each randomization (row) and each single rating event (column) the absolute difference in ratings before the rating event}
  \item{misc}{various information}
  \item{ov}{data set overview, i.e. in how many trials was a stimulus involved and how many trials did each stimulus win and lose}
  \item{ias}{character matrix, with the original sequence of rating events}
}
\references{
Elo, A. E. 1978. The Rating of Chess Players, Past and Present. New York: Arco.

Clark, A. P., Neumann, C., Howard, K. L., Woods, A. T. & Penton-Voak, I. S. (in preparation). Why rate when you could compare? Using the ``EloChoice'' package to assess pairwise comparisons of perceived strength.
}
\author{
Christof Neumann
}
\examples{
\dontrun{
ws <- c("b","b","b","b"); ls <- c("a", "c", "D", "e")
elochoice(ws, ls, 100, 1000)
elochoice(ws, ls, 200, 1000)}
}
