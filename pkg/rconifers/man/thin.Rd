
\name{thin}
\alias{thin}

\title{Thin a CONIFERS sample.data object with user-defined settings.}

\description{
  Reduces the expansion factors for tree records in \code{\link{sample.data}} list object, to the specification of the user.
}

\usage{
thin( x,control=list(type=2,target=50.0,target.sp=NULL ) )
}

\arguments{
  \item{x}{an of type \code{\link{sample.data}}.}
  \item{control}{a list that controls how the sample is thinned. See details.}
}

\details{

\describe{
\item{type}{integer indicator of thinning type.}
\tabular{rl}{
  0 \tab Thin trees to a target sdi, does not remove shrubs (not used directly)\cr
  1 \tab Proportional thin for selected tree species, does not remove shrubs\cr
  2 \tab Proportional thin across all tree species\cr
  3 \tab Thin from below, by dbh, all species\cr
  4 \tab Thin from below, by dbh, species identified with \code{sp.code}\cr
}

% \item{target}{target thinning level. The value is the targeted value for
%   the thinning operation. For an SDI thin, 0, the target is SDI in
%   English units. For remaining thinning operations, target is expressed
%   in stems per acre. (See \code{\link{species.swo}},
%   \code{\link{species.smc}}, \code{\link{species.cips}}, or
%   \code{\link{species.swohybrid}}).  For proportional thinning
%   (i.e. type=2 or 3), the species must be set to \code{NULL}
%   (i.e. sp=NULL)}.

% }


\item{target}{target thinning level. The value is the targeted value for
  the thinning operation. For an SDI thin, 0, the target is SDI in
  English units. For remaining thinning operations, target is expressed
  in stems per acre. (See \code{\link{species.swo}},
  \code{\link{species.smc}}, or \code{\link{species.swohybrid}}).
  For proportional thinning
  (i.e. type=2 or 3), the species must be set to \code{NULL}
  (i.e. sp=NULL)}.

}

The default thinning regime is to remove 50 percent of all trees
across all diameter classes. 

For thinning from below, the trees are ranked in increasing diameter,
and stems are removed (i.e. expf=0.0) for all stems until the target
percent removal has been achieved or all stems have been removed. 

The function prints out a diagnostic statement reporting the number of
stems removed and the total basal area removed during the operation. 

To reduce the shrubs, use a type 4 thinning, which thins from below,
by dbh, which here can be zero, for a single species (i.e. sp="COCO")

The sp element of the control list is examined only is the thinning
type is one (1) (\code{DO_EXPF_SP_THIN}) or four (4)
(\code{DO\_EXPF_SP_THIN_FROM_BELOW}), otherwise the argument is ignored and
assumed to be zero, which includes all species.

}


\value{
  The function returns a list object of type \code{\link{sample.data}}.
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

## Example 1. Basic functions and features of the package.
## load the CONIFERS library
library( rconifers )

## set the variant to the SWO variant.
## the set.variant command also assigns the species map to the default species map data.frame
set.variant( 0 )

# load the SWO species coefficients into R as a data.frame object
data( species.swo )
set.species.map( species.swo )

## grow the data that was originally in the swo variant
# load and display CONIFERS example plots
data( plots.swo )
print( plots.swo )

# load and display CONIFERS example plants
data( plants.swo )
print( plants.swo )

# create the sample.data list object
sample.swo.3 <- list( plots=plots.swo, plants=plants.swo, age=3, x0=0.0,
  n.years.projected=0 )
class(sample.swo.3)  <- "sample.data"

# now, project the sample forward for 20 years
# with all of the options turned off
sample.swo.23 <- project( sample.swo.3,
	  20, 
	  control=list(rand.err=0,rand.seed=0,endemic.mort=0,sdi.mort=0))

## print the projected summaries
print( sample.swo.23 )

## plot the diagnostics from the fit a linear dbh-tht model 
## before thinning
opar <- par( mfcol=c(2,2 ) )
plot( lm( sample.swo.23$plants$tht ~ sample.swo.23$plants$dbh ) )
par( opar )

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


