\name{aglUn}
\alias{aglUn}
\docType{data}
\title{aglUn}
\description{
Cross-National panel data on the OECD countries containing information about economic performance, government partisanship and labor organization.
}
\usage{data(aglUn)}
\format{
  A data frame with 230 observations on the following 10 variables.
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
This data frame differs from 'agl' only by the random omission of 10 rows of data. This is to created an unbalanced data version.
}
\source{
Alvarez, R. Michael, Geoffrey Garrett and Peter Lange. (1991) Government Partisanship, Labor Organization, and Macroeconomic Performance. \emph{American Political Science Review} 85:539-56.}
\references{
Beck, Nathaniel, Jonathan N. Katz, R. Michael Alvarez, Geoffrey Garrett and Peter Lange. (1993) Government Partisanship, Labor Organization, and Macroeconomic Performance: A Corrigendum. \emph{American Political Science Review} 87:945-948.
}
\examples{
data(aglUn)
summary(aglUn)
}
\keyword{datasets}
