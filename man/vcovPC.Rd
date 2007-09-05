\name{vcovPC}
\alias{vcovPC}
\alias{vcovPC.default}
\title{Extract Panel-Corrected Variance Covariance Matrix}
\description{
The package pcse contains a function to estimate panel-corrected standard errors. Data may contain balanced or unbalanced panels. This function extracts the resulting variance covariance matrix.
}

\usage{vcovPC(x, \dots)

\method{vcovPC}{default}(x, groupN, groupT, pairwise=FALSE, \dots)
}
\arguments{
\item{x}{A lm object containing the initial run of OLS.}
\item{groupN}{A vector containing the cross-sectional group identifier for each observation.}
\item{groupT}{A vector containing the time identifier for each observation.}
\item{pairwise}{An optional logical flag indicating whether the X's used to estimate the "middle" matrix should be chosen in a pairwise fashion or casewise fashion. If pairwise, the correlation between observations $i$ and $j$ is based on the time periods common to $i$ and $j$. If casewise, the correlation between observations i and j is based on the largest rectangular subset of the data, i.e., $T_i$ = $T_j$ = $T^*$ for all $i$ and $j$ if casewise is selected.}
\item{\dots}{Further arguments passed to methods.}
}



\examples{
## see demo file.
}

\author{
Delia Bailey <delia.bailey@gmail.com> and Jonathan N. Katz <jkatz@caltech.edu>
Maintainer: Delia Bailey <delia.bailey@gmail.com>
}

\keyword{file}
