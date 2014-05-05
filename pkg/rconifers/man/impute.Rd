
%%	$Id $	

\name{impute}
\alias{impute}

\title{Imputes missing values using the CONIFERS imputation functions}

\description{
  Imputes missing values for a CONIFERS sample.data object.
}

\usage{
impute( x, control=list(fpr=11.78, min.dbh=5.6, baf=40.0) )
}

\arguments{
  \item{x}{a CONIFERS \code{\link{sample.data}} object}
  \item{control}{a list object that controls how missing expansion
  factors are calculated (See details)}
}

\details{The impute function returns a list object of class
  \code{\link{sample.data}} that contains all the required data
  elements needed to call most CONIFERS functions. The inbound data
  frame objects for the plants and plots associated with the
  \code{\link{sample.data}} object must contain all the fields
  described in the \code{\link{plants.swo}} and \code{\link{plots.swo}}. The control argument
  needs to contain at least the following elements:

\describe{
  \item{fpr}{is the radius, in feet, for a fixed area plot which is used to
  sample all stems with a breast height diameter (DBH) that is less than min
  \code{min.dbh} argument.}
  \item{min.dbh}{is the minimum value, in inches, for stems to be sampled on a
  fixed area plot.}
  \item{baf}{is the constant, in $ft^2/acre$ that defines the maximum distance for
  which stems with a diameter at breast height over the \code{min.dbh}}
}
}

\value{	
	\code{impute} returns a \code{\link{sample.data}} object with
	no missing values that can then be used in any \code{rconifers} function.
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

## load the rconifers library
library( rconifers )

## set the variant to the SWO variant and set species map
set.variant(0)

## load the species coefficients into R as a data.frame object and
## set the species map
data( species.swo )
print( species.swo )
set.species.map( species.swo )

## grow the data that was originally swo in the smc variant
## load and display CONIFERS example plots
data( plots.swo )
print( plots.swo )

## load and display CONIFERS example plants
data( plants.swo )
print( plants.swo )

## view the current crown widths for first 7 obs
plants.swo[1:7,]

## remove crown.width observations from plants 3 to 6 to test
## the impute function
blanks <- c(3:6)
plants.swo[blanks,]$crown.width <- NA
plants.swo[1:7,]

## create the sample.data list object
sample.swo.3 <- list( plots=plots.swo, plants=plants.swo, age=3, x0=0.0,
   n.years.projected=0 )
class(sample.swo.3)  <- "sample.data"

## fill in missing values
sample.swo.3.imp <- impute( sample.swo.3 )

## view imputed values
sample.swo.3.imp$plants[1:7,]

}

\keyword{models}


