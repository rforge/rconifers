

\keyword{misc}
\name{variants}
\alias{variants}

\title{Print the CONIFERS growth model variants}

\description{Print the variant indexes and labels available in the CONIFERS forest growth model.}

\usage{
variants()
}

\arguments{ None }

\details{

  The \code{variants} function simply prints out the variant indexes and code labels for the variants available within the CONIFERS growth model.

}

\value{
  none.
}


\references{
  
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

Maguire, Doug. 2010. Center for Intensive Planted-Forest Silviculture
Annual Report. Oregon State University. 60 p.

}

\author{Jeff D. Hamann \email{jeff.hamann@forestinformatics.com},\cr
	     Martin W. Ritchie \email{mritchie@fs.fed.us},\cr 
	     Doug Maguire \email{doug.maguire@oregonstate.edu} 
        }


\seealso{
  \code{\link{calc.max.sdi}},
  \code{\link{impute}},
  \code{\link{plants.cips}},
  \code{\link{plants.smc}},
  \code{\link{plants.swo}},
  \code{\link{plants.swohybrid}},
  \code{\link{plot.sample.data}},
  \code{\link{plots.cips}},
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
  \code{\link{species.cips}},
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

## print the possible variants to use for your simulations.
variants()

}

\keyword{models}


