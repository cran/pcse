\name{agl}
\alias{agl}
\docType{data}
\title{agl}
\description{
Cross-National panel data on the OECD countries containing information about economic performance, government partisanship and labor organization.
}
\usage{data(agl)}
\format{
  A data frame with 240 observations on the following 10 variables.
  \describe{
    \item{\code{year}}{a numeric vector}
    \item{\code{country}}{a character vector}
    \item{\code{growth}}{the OECD growth rate}
    \item{\code{lagg1}}{An instrument for lagged growth rates constructed with an auxilary regression.}
    \item{\code{opengdp}}{weighted OECD demand}
    \item{\code{openex}}{weighted OECD export}
    \item{\code{openimp}}{weighted OECD import}
    \item{\code{leftc}}{"Left" cabinet composition}
    \item{\code{central}}{labor organization index}
    \item{\code{inter}}{interaction between leftc and central}
  }
}
\details{
}
\source{
Alvarez, R. Michael, Geoffrey Garrett and Peter Lange. (1991) Government Partisanship, Labor Organization, and Macroeconomic Performance. \emph{American Political Science Review} 85:539-56.}
\references{
Beck, Nathaniel, Jonathan N. Katz, R. Michael Alvarez, Geoffrey Garrett and Peter Lange. (1993) Government Partisanship, Labor Organization, and Macroeconomic Performance: A Corrigendum. \emph{American Political Science Review} 87:945-948.
}
\examples{
data(agl)
summary(agl)
}
\keyword{datasets}
