\name{pcse-package}
\alias{pcse-package}
\alias{pcse}
\docType{package}
\title{pcse: A Package for Panel-Corrected Standard Error Estimation in R}
\description{
This package contains a function to estimate panel-corrected standard errors. Data may contain balanced or unbalanced panels. 
}
\details{
\tabular{ll}{
Package: \tab pcse\cr
Type: \tab Package\cr
Version: \tab 1.3\cr
Date: \tab 2007-07-27\cr
License: \tab \cr
}
}

\usage{pcse(object, groupN, groupT, pairwise=FALSE)}
\arguments{
\item{object}{A lm object containing the initial run of OLS.}
\item{groupN}{A vector containing the cross-sectional group identifier for each observation.}
\item{groupT}{A vector containing the time identifier for each observation.}
\item{pairwise}{An optional logical flag indicating whether the X's used to estimate the "middle" matrix should be chosen in a pairwise fashion or casewise fashion. If pairwise, the correlation between observations $i$ and $j$ is based on the time periods common to $i$ and $j$. If casewise, the correlation between observations i and j is based on the largest rectangular subset of the data, i.e., $T_i$ = $T_j$ = $T^*$ for all $i$ and $j$ if casewise is selected.}
}

\examples{
## see demo file.
}
\references{
Alvarez, R. Michael, Geoffrey Garrett and Peter Lange. (1991) Government Partisanship, Labor Organization, and Macroeconomic Performance. \emph{American Political Science Review} 85:539-56.

Beck, Nathaniel, Jonathan N. Katz, R. Michael Alvarez, Geoffrey Garrett and Peter Lange. (1993) Government Partisanship, Labor Organization, and Macroeconomic Performance: A Corrigendum. \emph{American Political Science Review} 87:945-948.

Beck, Nathaniel and Jonathan N. Katz. (1995) What to do (and not to do) with Time-Series Cross-Section Data. \emph{American Political Science Review} 89:634-647.

}

\author{
Delia Bailey <dbailey@wustl.edu> and Jonathan N. Katz <jkatz@caltech.edu>
Maintainer: Delia Bailey <dbailey@wustl.edu>
}

\keyword{ package }

