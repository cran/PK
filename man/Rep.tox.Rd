\name{Rep.tox}
\alias{Rep.tox}
\title{Plasma levels in rats following daily intravenous administration in a repeated dose toxicity study}
\description{
  The \code{Rep.tox} data frame has 28 rows and 4 columns and gives the plasma levels of 6 rats in a repeated dose toxicity study. 
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{id}{
      a numeric vector identifying the animal whose plasma level is measured.
    }
    \item{conc}{
      a numeric vector giving the measured plasma level (IU/mL).
    }
    \item{time}{
      a numeric vector giving the time since administration (hours).
    }
    \item{day}{
      a numeric vector indicating the day the data were collected. 
    }
  }
}
\details{
  In this study the compound is administered daily to 6 Rats over 14 days. Plasma levels are measured at 7 time points over 12 hours after the first administration and after the 14th administration.\cr
}
\source{
Wolfsegger M. J. and Jaki T. (2009) Assessing systemic drug exposure in repeated dose toxicity studies in the case of complete and incomplete sampling. \emph{Biometrical Journal}, 51(6):1017:1029.
}
%\examples{}
\keyword{datasets}
