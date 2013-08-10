\name{plot.halflife}
\alias{plot.halflife}
\title{Plot regression lines used for half-life estimation}
\description{This method plots objects of S3 class \code{"halflife"} (\link{biexp} and \link{lee}). }

\usage{\method{plot}{halflife}(x, xlab='Time', ylab='Concentration', 
        main='Half-life Estimation', xlim=NULL, ylim=NULL, add=FALSE, ...)}

\arguments{
  \item{x}{ An object of S3 class \code{"halflife"} (\link{biexp} and \link{lee}).}
  \item{xlab}{ A label for the x axis.}
  \item{ylab}{ A label for the y axis.}
  \item{main}{ A main title for the plot.}
  \item{xlim}{ The x limits (min, max) of the plot. }
  \item{ylim}{ The y limits (min, max) of the plot. }
  \item{add}{ A logical value indicating whether to add plot to current plot (default=\code{FALSE}). }
  \item{\dots}{ Other parameters to be passed through to plotting functions.}
}

\author{Martin J. Wolfsegger}

\value{none}
\keyword{hplot}

