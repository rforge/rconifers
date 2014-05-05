
%%	$Id $	

\keyword{datasets}
\name{plants.smc}
\alias{plants.smc}
\non_function{}
\title{Simple plant data set for the SMC variant of the CONIFERS growth model.}
\usage{data(plants.smc)}

\description{
This a simple plants \code{\link{data.frame}} that contains basic
plant attributes used in the SMC variant of the CONIFERS growth
model. It contains the same columns as the \code{\link{plants.swo}} data.frame.

The \code{plants} \code{\link{data.frame}} has 60 rows and 9 columns.
}

\format{
  This data.frame contains the following columns:
  
  \describe{
    \item{plot}{The plot identifier.}
    % \item{sp.code}{species code (See \code{\link{species.swo}},
    % \code{\link{species.smc}}, \code{\link{species.cips}}, or \code{\link{species.swohybrid}})}
    \item{sp.code}{species code (See \code{\link{species.swo}},
    \code{\link{species.smc}}, or \code{\link{species.swohybrid}})}
    \item{d6}{Diameter at 6 inches above soil interface.}
    \item{dbh}{Diameter at breast height, in inches.}
    \item{tht}{Total height of the stem, in feet.}
    \item{cr}{Crown ratio (0 to 1).}
    \item{n.stems}{Number of stems this record represents.}
    \item{expf}{Expansion factor (weight) the tree record represents, in
      stems per acre.}
    \item{crown.width}{Width of the crown at the widest point, in feet.}
  }
}

\details{The example data represents a field sample that can be used
to generate future forest conditions using the SMC variant of the
CONIFERS forest growth model. The data are in imperial units. The
required values are the \code{plot}, \code{sp.code}, and \code{d6}
measurements. Values that have not been observed, or are missing, are
denoted with an \code{NA} and can be imputed using the
\code{\link{impute}} function.}

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

## set the variant to the CIPS variant.
set.variant( 1 )

## load and set the species mapping
data( species.smc )
print( species.smc )
set.species.map( species.smc )

## load and display the plants
data( plants.smc )
print( plants.smc )

}

\keyword{datasets}

