\name{biexp}
\alias{biexp}
\title{Two-Phase Half-Life Estimation by Biexponential Model}
\description{Estimation of inital and terminal half-life by fitting a biexponential model.}

\usage{biexp(conc, time, prev=0, tol=1E-9, maxit=500)}

\arguments{
  \item{time}{ time points of concentration assessments. }
  \item{conc}{ levels of concentrations. }
  \item{prev}{ pre-dosing value. }
  \item{tol}{ relative error tolerance. } 
  \item{maxit}{ maximum number of iterations.}
}

\details{
Estimation of inital and terminal half-life using the biexponential \code{y=a1*exp(-b1*x)+a2*exp(-b2*x)} model with a parametrization to ensure b1 > b2 > 0 fitted by the least squares criteria with function \code{optim} of package \code{base} with \code{method} "Nelder-Mead". Curve peeling (Foss, 1969) is used get start values for nonlinear model fitting. When no adequate starting values are determined by curve peeling, a single exponential model is fitted with starting values obtained from an OLS regression on log transformed values with a parametrization to ensure a slope > 0. \cr \cr

If the pre-dosing value indicating the intrinsic level is greater than 0, the pre-dosing value is subtracted from all concentration levels before calculation of initial and terminal half-life.
}
 
\value{
A list of S3 class \code{"halflife"} containing the following components: \cr
  \item{parms}{ half-life and model estimates.}
  \item{time}{ time points of concentration assessments. }
  \item{conc}{ levels of concentrations. } 
  \item{method}{ "biexp". }
}

\note{Records including missing values and values below or equal to zero are omitted. }

\references{
Foss S. D. (1969). A Method for Obtaining Initial Estimates of the Parameters in Exponential Curve Fitting. Biometrics. 25:580-584 \cr \cr
Pinheiro J. C. and Bates D. M. (200). Mixed-Effects Models in S and S-PLUS. Springer, New York. \cr \cr
}

\author{Martin J. Wolfsegger and Thomas Jaki}

\examples{
## examples from Pinheiro J.C. and Bates D.M. (2000, page 279) 
time <- c(0.25, 0.5, 0.75, 1, 1.25, 2, 3, 4, 5, 6, 8, 0.25, 0.5, 0.75, 1, 1.25, 
2, 3, 4, 5, 6, 8, 0.25, 0.5, 0.75, 1, 1.25, 2, 3, 4, 5, 6, 8, 0.25, 0.5, 0.75, 1, 
1.25, 2, 3, 4, 5, 6, 8, 0.25, 0.5, 0.75, 1, 1.25, 2, 3, 4, 5, 6, 8, 0.25, 0.5, 
0.75, 1, 1.25, 2, 3, 4, 5, 6, 8)

conc <- c(1.5, 0.94, 0.78, 0.48, 0.37, 0.19, 0.12, 0.11, 0.08, 0.07, 0.05, 2.03, 
1.63, 0.71, 0.7, 0.64, 0.36, 0.32, 0.2, 0.25, 0.12, 0.08, 2.72, 1.49, 1.16, 0.8, 
0.8, 0.39, 0.22, 0.12, 0.11, 0.08, 0.08, 1.85, 1.39, 1.02, 0.89, 0.59, 0.4, 0.16, 
0.11, 0.1, 0.07, 0.07, 2.05, 1.04, 0.81, 0.39, 0.3, 0.23, 0.13, 0.11, 0.08, 0.1, 
0.06, 2.31, 1.44, 1.03, 0.84, 0.64, 0.42, 0.24, 0.17, 0.13, 0.1, 0.09)

result <- biexp(conc=conc, time=time)
print(result)
plot(result)

}
\keyword{misc}