

\keyword{misc}
\name{set.variant}
\alias{set.variant}

\title{Set the CONIFERS growth model variant}

\description{Assigns the variant for the CONIFERS forest growth model.}

\usage{
set.variant( var=0 )
}

\arguments{

  \item{var}{is the CONIFERS variant. The only possible values for the
    variant at the current time are:

\tabular{lll}{
\emph{id} \tab \emph{code} \tab \emph{citation?}\cr
0 \tab CONIFERS_SWO \tab The Southwest Oregon/Northern California variant.\cr
1 \tab CONIFERS_SMC \tab Stand Management Cooperative.\cr
2 \tab CONIFERS_SWOHYBRID \tab Ritchie, ???, and Hamann. 2011?\cr
}

}
}


%% 3 \tab CONIFERS_CIPS \tab Mawaring, Magure. 2012?\cr

\details{

  % The \code{set.variant} function only accepts an argument of zero (0), one (1), two (2), and (3). The default variant is the (\code{SWO\_CONIFERS})
  % variant. This function returns a value of 0, 1, 2, 3 or an error
  % message. 

    The \code{set.variant} function only accepts an argument of zero
  (0), one (1), and two (2). The default variant is the (\code{SWO\_CONIFERS})
  variant. This function returns a value of 0, 1, 2 or an error
  message. 

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

% Maguire, Doug. 2010. Center for Intensive Planted-Forest Silviculture
% Annual Report. Oregon State University. 60 p.

}

\author{Jeff D. Hamann \email{jeff.hamann@forestinformatics.com},\cr
  Martin W. Ritchie \email{mritchie@fs.fed.us}}

	     % Doug Maguire <doug.maguire@oregonstate.edu>}


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

## load the RCONIFERS library
library( rconifers )

# Example 1 (SWO)
## set the variant to the SWO variant and set species map
set.variant(0)
data( species.swo )
set.species.map( species.swo )

# Example 2 (SMC)
## set the variant to the SMC variant
set.variant(1)
data( species.smc )
set.species.map( species.smc )

## and manually set the genetic worth values in the functional species coefficients
data( species.smc )
species.smc[species.smc$code == 'DF',]$genetic.worth.d <- 20
species.smc[species.smc$code == 'DF',]$genetic.worth.h <- 20
set.species.map( species.smc, verbose=TRUE )

# Example 3 (SWOHYBRID)
## set the variant to the SWOHYBRID variant
set.variant(2)
data( species.swohybrid )
set.species.map( species.swohybrid )

\dontrun{
# Example 4 (CIPS)
## set the variant to the CIPS variant
set.variant(3)
data( species.cips )
set.species.map( species.cips )
}

}

\keyword{models}


