
%%	$Id $	

\name{rand.seed}
\alias{rand.seed}

\title{Initialize or reset the random number generator for the CONIFERS forest growth model}

\description{
Initialize or reset the random number generator for the CONIFERS forest growth model.  
}

\usage{
rand.seed( control )
}

\arguments{
  \item{control}{a list of control options. See 'Details'.}
}


\details{
  The \code{rand.seed} function is used to either initialize or reset
  the random behavior of the CONIFERS growth model. 

  The 'control' argument is a list that can supply any of the following
  components:

\describe{
  \item{use.random.error}{Non-negative integer. If 0, no random error will be
    applied in the height growth function. If 1, then a random error will
    be added to each plant, in each year to the height growth.}

  \item{random.seed}{Non-negative integer. If 0, the random number
    generator will be seeded with the time. If >0, the random number
    generator will be seeded with this value. To use random error for each
    projection, you only need to initialize the random number generator
    once. If you wish to repeat the randomness, you need to reseed with
    the same seed, in the same sequence.}
}
}



\references{

Ritchie, M.W. 2008 . User's Guide and Help System for CONIFERS: A Simulator for Young Conifer Plantations Version
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

library( rconifers )

# turn off the random error component
rand.seed( control=list(use.random.error=0 ) )

# reset the random error component and set the seed to 4326
rand.seed( control=list(use.random.error=1,random.seed=4326 ) )

# reset the random error component and set the seed to a random number between 0 and 10000
germ<-sample(1:9999,1)
rand.seed( control=list(use.random.error=1, random.seed=germ) )

}

\keyword{models}


