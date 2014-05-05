
%%	$Id $	

\name{calc.max.sdi}
\alias{calc.max.sdi}

\title{Calculate the maximum stand density index using the CONIFERS forest growth model}

\description{
Returns the maximum stand density index, given the current species
mapping and variant.  
}

\usage{
calc.max.sdi( x )
}
		   
\arguments{
  \item{x}{a \code{\link{sample.data}} object.}
}

\details{

  The \code{calc.max.sdi} function returns a numeric value that
  represents the maximum stand density index limit for the input
  plants/plots data. The value is used to compute the $x_0$ element in
  the \code{\link{sample.data}} object, which is then used to control
  competition induced mortality.

}


\references{

Hann, D.W. and C.H. Wang. 1990. Mortality equations for individual trees
      in southwest Oregon. Oregon State University, Forest Research
      Laboratory, Corvallis, Oregon. Research Bulletin 67. 17p.
  
Ritchie, M. and J. Hamann. 2006. Modeling dynamics of competing
vegetation in young conifer plantations of northern California and
southern Oregon, USA. Canadian Journal of Forest Research 36(10):
2523-2532.

Ritchie, M. and J. Hamann. 2008. Individual-tree height-, diameter- and
crown-width increment equations for young Douglas-fir plantations. New
Forests 35(2):173-186.

Ritchie, M.W. 2008. User's Guide and Help System for CONIFERS: A
Simulator for Young Conifer Plantations Version 4.10. See
\url{http://www.fs.fed.us/psw/programs/ecology_of_western_forests/projects/conifers/}

Vaughn, Nicholas. 2007. An individual-tree model to predict the annual
growth of young stands of Douglas-fir (Pseudotsuga menziesii (Mirbel)
Franco) in the Pacific northwest. M.S. Thesis, University of
Washington. 91p.

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

## set the variant to the SWO variant
set.variant(0)

## load the species coefficients into R as a data.frame object and
## set the species map
data( species.swo )
print( species.swo )
set.species.map( species.swo )

## load and display CONIFERS swo example plots
data( plots.swo )
print( plots.swo )

## load and display CONIFERS example plants
data( plants.swo )
print( plants.swo )

## create the sample.data list object
sample.swo <- list( plots=plots.swo, plants=plants.swo, age=3, x0=0.0, n.years.projected=0 )
class(sample.swo)  <- "sample.data"

## print the maximum stand density index for the current settings
print( calc.max.sdi( sample.swo ) )


}

\keyword{models}


