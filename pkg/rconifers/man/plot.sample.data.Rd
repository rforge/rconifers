\name{plot.sample.data}
\alias{plot.sample.data}
\title{Summary plot of a CONIFERS sample.data object}

\description{
   These functions create and plot summary results for a CONIFERS
   sample.data object.
}

\usage{
   \method{plot}{sample.data}(x, digits = max(3, getOption("digits") - 1), ...)
}

\arguments{
   \item{x}{an object of class \code{\link{sample.data}}.}
   \item{digits}{the number of digits to display.}
   \item{...}{not used by user.}
}

\value{Applying \code{plot} on an object of class
   \code{\link{sample.data}} generates a simple summary of the
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

## load the CONIFERS library
library( rconifers )

## set the variant to the SWO variant.
set.variant( 0 )

# load the SWO species specific coefficients into R as a data.frame object
data( species.swo )
set.species.map( species.swo )

## load and print the plots for the SWO example dataset.
data( plots.swo )
print( plots.swo )

## load and print the plants for the SWO example dataset.
data( plants.swo )
print( plants.swo )

# create the sample.data list object
sample.swo.3 <- list( plots=plots.swo, plants=plants.swo, age=3, x0=0.0,n.years.grown=0 )
class(sample.swo.3)  <- "sample.data"

## generate a generic summary plot 
plot(sample.swo.3)


}

\keyword{models}



