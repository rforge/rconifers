
\name{sp.sums}
\alias{sp.sums}
\title{A summary of a CONIFERS sample.data object by species}

\description{ This function returns a \code{\link{data.frame}} object
   that contain common summaries by species for
   \code{\link{sample.data}} object.  
}

\usage{
   sp.sums( x )
}

\arguments{
   \item{x}{an object of class \code{\link{sample.data}}.}
}

\format{

The \code{\link{data.frame}} object returned from sp.sums contains the
following columns for each species:

\describe{
\item{bhexpf}{The number of plants per acre that are above 4.5 feet (breast height).}
\item{texpf}{The total plants per acre including those below breast height (4.5 feet).}
\item{tht}{The mean total height.}
\item{ba}{The total basal area for all plants above breast height.}
\item{qmd}{The quadratic mean diameter for all trees above breast height.}
}
}

\details{ This function returns a \code{\link{data.frame}} object that
	contains species level summary information. The rownames of the 
	data.frame object contain the species codes defined in the 
	species.map. It is intended demonstration only and users are 
  	encouraged to examine the source code and modify.
}

\references{

Ritchie, M.W. 2008. User's Guide and Help System for CONIFERS: A Simulator for Young Conifer 
Plantations Version 4.10. 
See \url{http://www.fs.fed.us/psw/programs/ecology_of_western_forests/projects/conifers/}

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

\dontrun{
################################################################################
## Example 1 Using the CIPS variant

## set the variant to the CIPS variant.
## the set.variant command also assigns the species map to the default species map data.frame
set.variant( 3 )

# load the CIPS species coefficients into R as a data.frame object
data( species.cips )
set.species.map( species.cips )

## load the plots for the CIPS variant example
data( plots.cips )
print( plots.cips )

# load and display CONIFERS example cips plants
data( plants.cips )
print( plants.cips )

# create the sample.data list object
sample.cips <- list( plots=plots.cips, plants=plants.cips, age=3, x0=0.0, n.years.projected=0 )
class(sample.cips)  <- "sample.data"

## print and plot the sample using the default functions
print( sample.cips )
plot( sample.cips )

## print the species summaries, and extract the quadratic mean diameter
## for Douglas-fir (DF) only.
print( sp.sums( sample.cips ) )
smry<-sp.sums( sample.cips )
dfqmd<-smry["DF","qmd"]
}


################################################################################
## Example 2 Using SWO
## set the variant to the SWO variant
set.variant( 0 )

## load the SWO species coefficients into R as a data.frame object
data( species.swo )
set.species.map( species.swo )

## grow the swo data
## load and display CONIFERS example plots
data( plots.swo )
print( plots.swo )

## load and display CONIFERS example plants
data( plants.swo )
print( plants.swo )

## create the sample.data list object
sample.swo.3 <- list( plots=plots.swo, plants=plants.swo, age=3, x0=0.0,
n.years.projected=0 )
class(sample.swo.3)  <- "sample.data"

## print the maximum stand density index for the current settings
print( calc.max.sdi( sample.swo.3 ) )

## print a summary of the sample
print( sample.swo.3 )

## now, project the sample forward for 20 years
sample.swo.23 <- project( sample.swo.3,
                          20,
                          control=list(rand.err=0,
                                       rand.seed=0,
                                       endemic.mort=0,
                                       sdi.mort=1,
                                       genetic.gains=0 ) )


## print the projected summaries
print( sample.swo.23 )

## thin the stand
## Proportional thin for selected tree species (DF),
## it does not remove shrubs
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

## print the differences, by species
print( sp.sums( sample.swo.23.t4 ) - sp.sums( sample.swo.23 ) )

## generate a plot of the species summaries
opar <- par(mfrow=c(2,2))

##hist( sample.swo.23$plants$dbh, xlim=c(0,50), ylim=c(0,50), main="Diameter Distribution" )
##hist( sample.swo.23$plants$tht, xlim=c(0,100), ylim=c(0,110), main="Height Distribution" )

hist( sample.swo.23$plants$dbh, main="Diameter Distribution" )
hist( sample.swo.23$plants$tht, main="Height Distribution" )

ba <- sp.sums(sample.swo.23)$ba
names(ba) <- rownames( sp.sums(sample.swo.23) )
pie( ba[ba>0], main="Basal Area" )
plot( sample.swo.23$plants$tht ~ sample.swo.23$plants$dbh, ylim=c(0,220), 
                                                           xlim=c(0,40), 
                                                           main="Height vs. Basal Diameter" )

par(opar)


}

\keyword{models}



