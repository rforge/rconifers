
%%	$Id $	

\name{project}
\alias{project}

\title{Projects a sample.data object using the CONIFERS forest growth model}

\description{
  Projects a CONIFERS sample.data object.
}

\usage{
project( x,years=1,control=list(rand.err=0,
                        rand.seed=0,
                        endemic.mort=0,
                        sdi.mort=0,
			genetic.gains=0) )
}
		   
\arguments{
  \item{x}{a sample.data object.}
  \item{years}{number of years to project the sample.data}
  \item{control}{A list of control parameters. See *Details*}
}

\details{
  The project function returns a sample.data object that contains the
  projected plant list. The plot data are simply copied from the input
  sample.data object, and the x0 member is updated if the sample has
  been projected past the point where stand growth becomes limited by
  competition. 

  By default this function projects the sample.data one year, with the
  default control parameters. The defaults project the sample with a
  random error component added to the height increment, with endemic
  mortality turned on and the stand density index induced mortality
  turned on as well. 

  The 'control' argument is a list that can supply any of the following
  components:

  \describe{
  \item{rand.err}{Non-negative integer. If 0, no random error will be
    applied in the height growth function. If 1, then a random error will
    be added to each plant, in each year to the height growth. This is
    where you include some explanation regarding the random error.}

  \item{rand.seed}{Non-negative integer. If 0, the random number
    generator will be seeded with the time. If >0, the random number
    generator will be seeded with this value. To use random error for each
    projection, you only need to initialize the random number generator
    once. If you wish to repeat the randomness, you need to reseed with
    the same seed, in the same sequence.}
  
  \item{endemic.mort}{Non-negative integer. If 0, no endemic mortality
    will be applied to the sample.data. If 1, the default mortality in the
    mortality settings table will be applied to the sample.data for each
    growth cycle.} 
  
  \item{sdi.mort}{Non-negative integer. If 0, no stand density index mortality
    will be applied to the sample.data. If 1, stand density index
    mortality will be applied to the sample.data for each growth
    cycle. This operation is applied independently of the endemic
    mortality.}    

  \item{genetic.gains}{Non-negative integer. If 0, genetic gains will not be applied to 
    sample.data. If 1, genetic.gains will be applied to the sample.data for each growth
    cycle. This operation is applied independently of the other control variables. 
    This option only applies to Douglas-fir species in the SMC variant. }    

  \item{cips.first.season}{Non-negative integer. If 0, the CIPS model will not assume the plantation has been grown past the first growing season. 
    If 1, the model will assume that the plantation has been grown past the first growing seasons. 
    This is important for users wanting to modify the plant records (e.g. vegetation control or thinning) and must call project multiple times during the simulation.
    This option only applies to the CIPS variant. }    



 }

}

\value{returns a projected sample.data object.}


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

## load the CONIFERS library
library( rconifers )

## set the variant to the SWO variant.
set.variant( 0 )

# load the SWO species coefficients into R as a data.frame object
data( species.swo )
print( species.swo )
set.species.map( species.swo )

## grow the smc example data
# load and display CONIFERS example swo plots
data( plots.swo )
print( plots.swo )

# load and display CONIFERS example swo plants
data( plants.swo )
print( plants.swo )

# create the sample.data list object
sample.swo.3 <- list( plots=plots.swo, plants=plants.swo, age=3, x0=0.0,
  n.years.projected=0 )
class(sample.swo.3)  <- "sample.data"
print( sample.swo.3)

## set the seed for the random number generator to provide a repeatable
## result for testing
rand.seed( control=list(use.random.error=1,random.seed=615) )

# now, project the sample forward for 20 years
# with all of the options turned off
sample.swo.23 <- project( sample.swo.3, 20,
	  control=list(rand.err=1,rand.seed=615,endemic.mort=0,sdi.mort=0))

print( sample.swo.23 )

## print the projected summaries
print( sample.swo.23 )

\dontrun{
# now, project the sample forward for 20 years
# with random error on and duplicate run 100 times, finish by calculating
#summaries of DF basal area
reps<-100
dfba<-vector(mode="double", length=reps)
for ( i in 1:reps){
sample.swo.23wr <- project( sample.swo.3, 20,
	  control=list(rand.err=1,rand.seed=615,endemic.mort=1,sdi.mort=1))
	  
#do a sp.sums of each run
summi<-sp.sums(sample.swo.23wr)

# select the Douglas-fir basal area
dfba[i]<-summi["DF","ba"]
}
summary(dfba)
} ## end of the dontrun block


################################################################################
## now, reset the variant to SMC and grow with genetic gains option

## load the CONIFERS library
library( rconifers )

## set the variant to the SWC variant.
set.variant( 1 )

# load the SWC species coefficients into R as a data.frame object
data( species.smc )
print( species.smc )
set.species.map( species.smc )

## grow the smc example data
# load and display CONIFERS example swo plots
data( plots.smc )
print( plots.smc )

# load and display CONIFERS example swo plants
data( plants.smc )
print( plants.smc )

# create the sample.data list object
sample.smc.3 <- list( plots=plots.smc, plants=plants.smc, age=3, x0=0.0,
  n.years.projected=0 )
class(sample.smc.3)  <- "sample.data"
print( sample.smc.3)

# now, project the sample forward for 20 years
# with all of the options turned off
sample.smc.23.0 <- project( sample.smc.3, 20,
	  control=list(rand.err=0,rand.seed=0,endemic.mort=0,sdi.mort=0,genetic.gains=0))
print( sample.smc.23.0 )


## project with genetic gains turned on.
data( species.smc )
species.smc[species.smc$code == 'DF',]$genetic.worth.d <- 20
species.smc[species.smc$code == 'DF',]$genetic.worth.h <- 20
set.species.map( species.smc )


## now, project the sample forward for 20 years
## with all of the options turned off
## except for the genetic gains
sample.smc.23.1 <- project( sample.smc.3, 20, 
	  control=list(rand.err=0,rand.seed=0,endemic.mort=0,sdi.mort=0,genetic.gains=1))
print( sample.smc.23.1 )


}

\keyword{models}


