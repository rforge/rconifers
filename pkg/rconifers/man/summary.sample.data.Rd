\name{summary.sample.data}
\alias{summary.sample.data}
\title{Summary of a CONIFERS sample.data object}

\description{
   These functions create and print summary results for a CONIFERS
   sample.data object.
}

\usage{
   \method{summary}{sample.data}( object, ... )
}

\arguments{
   \item{object}{an object of class \code{\link{sample.data}}.}
   \item{...}{not used by user.}
}

\value{Applying \code{summary} on an object of class
   \code{\link{sample.data}} prints out a simple summary of the
   \code{\link{sample.data}}. A \code{\link{sample.data}} object is
   the basic data type used for most CONIFERS functions.
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
options(digits=10)

## load the CONIFERS library
library( rconifers )

## set the variant to the SWO variant.
## the set.variant command also assigns the species map to the default species map data.frame
set.variant( 0 )

## load the SWO species coefficients into R as a data.frame object
data( species.swo )
set.species.map( species.swo )

## grow the smc example data
## load and display CONIFERS example swo plots
data( plots.swo )
print( plots.swo )

## load and display CONIFERS example swo plants
data( plants.swo )
print( plants.swo )

## create the sample.data list object
sample.swo <- list( plots=plots.swo, plants=plants.swo, age=3, x0=0.0, n.years.projected=0 )
class(sample.swo)  <- "sample.data"

## print and plot the sample
print( sample.swo )
plot( sample.swo )

## print a summary of the projected sample
summary( sample.swo )

}

\keyword{models}



