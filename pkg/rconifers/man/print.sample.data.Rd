
\name{print.sample.data}
\alias{print.sample.data}
\title{Print a CONIFERS sample.data object}

\description{
   These functions print a few summary values from a sample.data object.
}

\usage{
   \method{print}{sample.data}( x,
      digits = max( 3, getOption("digits") - 1 ), \dots )
}

\arguments{
   \item{x}{an object of class \code{sample.data}.}
   \item{digits}{number of digits to print.}
   \item{\dots}{other arguments.}
}



\references{

Ritchie, M.W. 2008. User's Guide and Help System for CONIFERS: A Simulator for Young Conifer Plantations Version
4.10. See \url{http://www.fs.fed.us/psw/programs/ecology_of_western_forests/projects/conifers/}

}

\author{Jeff D. Hamann \email{jeff.hamann@forestinformatics.com},\cr
	     Martin W. Ritchie \email{mritchie@fs.fed.us} }


\seealso{
  \code{\link{calc.max.sdi}},
  \code{\link{impute}},
  % \code{\link{plants.cips}},
  \code{\link{plants.smc}},
  \code{\link{plants.swo}},
  \code{\link{plants.swohybrid}},
  \code{\link{plot.sample.data}},
  % \code{\link{plots.cips}},
  \code{\link{plots.smc}},
  \code{\link{plots.swo}},
  \code{\link{plots.swohybrid}},
  \code{\link{print.sample.data}},
  \code{\link{project}},
  \code{\link{rand.seed}},
  \code{\link{rconifers}},
  \code{\link{sample.data}},
  \code{\link{set.species.map}},
  \code{\link{set.variant}},
  % \code{\link{species.cips}},
  \code{\link{species.smc}},
  \code{\link{species.swo}},
  \code{\link{species.swohybrid}},
  \code{\link{sp.sums}},
  \code{\link{summary.sample.data}},
  \code{\link{thin}},
  \code{\link{variants}}
}


\examples{

## set the number of digits for the output to five (5)
options(digits=5)

## load the CONIFERS library
library( rconifers )

## set the variant to the SWO variant.
set.variant( 0 )

## load the species coefficients into R as a data.frame object and
## set the species map
data( species.swo )
print( species.swo )
set.species.map( species.swo )

## grow the smc example data
## load and display CONIFERS example swo plots
data( plots.swo )
print( plots.swo )

## load and display CONIFERS example swo plants
data( plants.swo )
print( plants.swo )

## create the sample.data list object
sample.swo.3 <- list( plots=plots.swo, plants=plants.swo, age=3, x0=0.0,
  n.years.projected=0 )
class(sample.swo.3)  <- "sample.data"
print( sample.swo.3)



}

\keyword{models}
