
\keyword{methods}
\name{set.species.map}
\alias{set.species.map}

\title{Set and update the species mapping table in the CONIFERS simulator}

\description{ This function assigns the mapping of user defined
  species codes to the functional species within the CONIFERS growth
  model.}

\usage{set.species.map( x, verbose=FALSE )}

		   
\arguments{

  % \item{x}{ is a data.frame object with the same columns as
  % \code{\link{species.swo}}, \code{\link{species.smc}},
  % \code{\link{species.cips}}, or \code{\link{species.swohybrid}}}

  
  \item{x}{ is a data.frame object with the same columns as
    \code{\link{species.swo}}, \code{\link{species.smc}}, or
    \code{\link{species.swohybrid}}}

  \item{verbose}{ is switch for controlling the verbosity of the operation.}
}


%at least the zero-offset species
%    index in the species data frame, the functional species code (See
%    conifers manual, a user defined species code as a vector of
%    characters, a vector that contains the endemic mortality
%    coefficients, a vector that contains the maximum stand density
%    index value, a vector that contains the browse damage
%    coefficients, a vector that contains the mechanical damage
%    coefficients, a vector that contains genetic worth coefficients
%    for height growth, and a vector that contains genetic worth
%    coefficients for diameter growth.}

% \details{The set.species.map is how the user controls the behavior of
% the functional species that is used to project individual plants
% forward. The idx *must* be a zero-offset index vector which is assigned
% the element \code{SPECIES_RECORD[idx]} in the model (See
% \code{\link{species.swo}}, \code{\link{species.smc}},
% \code{\link{species.cips}}, or \code{\link{species.swohybrid}}). This
% function is used to control the mapping between the user species code
% and the functional species codes in the growth model.
% }


\details{The set.species.map is how the user controls the behavior of
the functional species that is used to project individual plants
forward. The idx *must* be a zero-offset index vector which is assigned
the element \code{SPECIES_RECORD[idx]} in the model (See
\code{\link{species.swo}}, \code{\link{species.smc}}, or
\code{\link{species.swohybrid}}). This
function is used to control the mapping between the user species code
and the functional species codes in the growth model.
}

\value{
	no value is returned.
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


