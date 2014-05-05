\name{rconifers}
\alias{rconifers}

%- Also NEED an `\alias' for EACH other topic documented here.

\title{Using the CONIFERS growth model from R}

\description{
  General information about the \code{rconifers} package.
}

%- maybe also `usage' for other objects documented here.
\details{

  \code{rconifers} provides an alternative interface to the standard CONIFERS
  graphical user interface (GUI). This package enables a user to
  define arbitrary silvicultural prescriptions, using standard R
  scripts, to predict near term (<20 years) future forest conditions
  under a variety of model configurations, management decisions, and
  policies. To learn more about the CONIFERS, follow this link
  \url{http://www.fs.fed.us/psw/programs/ecology_of_western_forests/projects/conifers/}.

  Type \code{library(help = rconifers)} to see a list of (public)
  \code{rconifers} functions.

}


\references{

Hann, D.W. and C.H. Wang. 1990. Mortality equations for individual
      trees in southwest Oregon. Oregon State University, Forest
      Research Laboratory, Corvallis, Oregon. Research Bulletin
      67. 17p.
  
Ritchie, M. and J. Hamann. 2006. Modeling dynamics of competing
vegetation in young conifer plantations of northern California and
southern Oregon, USA. Canadian Journal of Forest Research 36(10):
2523-2532.

Ritchie, M. and J. Hamann. 2008. Individual-tree height-, diameter-
and crown-width increment equations for young Douglas-fir
plantations. New Forests 35(2):173-186.

Ritchie, M.W. 2008. User's Guide and Help System for CONIFERS: A Simulator for Young Conifer Plantations Version
4.10. See \url{http://www.fs.fed.us/psw/programs/ecology_of_western_forests/projects/conifers/}

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

library( rconifers )

## set the variant to the SWO variant and set species map
##set.species.map( set.variant(0) )
set.variant(0)


## set the species map
data( species.swo )
set.species.map( species.swo )

## grow the data that was originally swo in the smc variant
# load and display CONIFERS example plots
data( plots.swo )
print( plots.swo )

# load and display CONIFERS example plants
data( plants.swo )
print( plants.swo )

# create the sample.data list object
sample.swo.3 <- list( plots=plots.swo, plants=plants.swo, age=3, x0=0.0,
 n.years.projected=0)
class(sample.swo.3)  <- "sample.data"

# print the maximum stand density index for the current settings
print( calc.max.sdi( sample.swo.3 ) )

# print a summary of the sample
print( sample.swo.3 )

# now, project the sample forward for 20 years
# with all of the options turned off
sample.swo.23 <- project( sample.swo.3,
	  20, 
	  control=list(rand.err=0,rand.seed=0,endemic.mort=0,sdi.mort=0))

## print the projected summaries
print( sample.swo.23 )

## thin the stand to capture mortality, redistribute growth,
## and possibly generate revenue
## Proportional thin for selected tree species, does not remove shrubs
sample.swo.23.t1 <- thin( sample.swo.23,
	     control=list(type=1, target=50.0, target.sp="DF" ) )
print( sample.swo.23.t1 )

## Proportional thin across all tree species
sample.swo.23.t2 <- thin( sample.swo.23,
	     control=list(type=2, target=50.0 ) )
print( sample.swo.23.t2 )

## Thin from below, by dbh, all species
sample.swo.23.t3 <- thin( sample.swo.23,
	     control=list(type=3, target=50.0 ) )
print( sample.swo.23.t3 )


## Thin from below, by dbh for species "PM"
sample.swo.23.t4 <- thin( sample.swo.23,
	     control=list(type=4, target=50.0, target.sp="PM" ) )
print( sample.swo.23.t4 )

}

\keyword{models}


