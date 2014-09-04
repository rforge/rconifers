

\keyword{datasets}
\name{plots.cips}
\alias{plots.cips}
\non_function{}
\title{Example plots data set for Center for Intensive Planted-forest Silviculture (CIPS) Variant of the CONIFERS growth model}
\usage{data(plots.cips)}

\description{ This a simple plot file \code{\link{data.frame}} that
contains basic plot attributes used in the CONIFERS growth model,
corresponding to the plants.cips dataframe.

The \code{plots.cips} data frame has ? rows and 6 columns.
}

\format{
This data frame contains the following columns:
\describe{
  \item{plot}{A numeric plot identifier.}
  \item{lat}{Latitude, in decimal degrees, of the plot.}
  \item{lon}{Longitude, in decimal degrees, of the plot.}
  \item{elevation}{Plot elevation, in feet.}
  \item{slope}{Plot slope, in percent.}
  \item{aspect}{Plot aspect, in degrees from North.}
  \item{whc}{Water holding capacity of the plot, in inches.}
  \item{map}{Mean annual precipitation, in inches.}

  %% growing season precip.
  %% check the preciup
  \item{gsp}{Growing season precipitation, in inches, for the growing season which normally is between April, May, and June.}

  \item{si30}{Site Index 30, in meters (optional). It must be present
  to project with variant=1 (CONIFERS_SMC) and variant=3 (CONIFERS_CIPS.}

  \item{mt1}{Monthly temperature, in C, for January (optional). It must be present to project with variant=2 (CONIFERS_SWOHYBRID)}
  \item{mt2}{Monthly temperature, in C, for February  (optional). It must be present to project with variant=2 (CONIFERS_SWOHYBRID)}
  \item{mt3}{Monthly temperature, in C, for March (optional). It must be present to project with variant=2 (CONIFERS_SWOHYBRID)}
  \item{mt4}{Monthly temperature, in C, for April (optional). It must be present to project with variant=2 (CONIFERS_SWOHYBRID)}
  \item{mt5}{Monthly temperature, in C, for May (optional). It must be present to project with variant=2 (CONIFERS_SWOHYBRID)}
  \item{mt6}{Monthly temperature, in C, for June (optional). It must be present to project with variant=2 (CONIFERS_SWOHYBRID)}
  \item{mt7}{Monthly temperature, in C, for July (optional). It must be present to project with variant=2 (CONIFERS_SWOHYBRID)}
  \item{mt8}{Monthly temperature, in C, for August (optional). It must be present to project with variant=2 (CONIFERS_SWOHYBRID)}
  \item{mt9}{Monthly temperature, in C, for September (optional). It must be present to project with variant=2 (CONIFERS_SWOHYBRID)}
  \item{mt10}{Monthly temperature, in C, for October (optional). It must be present to project with variant=2 (CONIFERS_SWOHYBRID)}
  \item{mt11}{Monthly temperature, in C, for November (optional). It must be present to project with variant=2 (CONIFERS_SWOHYBRID)}
  \item{mt12}{Monthly temperature, in C, for December (optional). It must be present to project with variant=2 (CONIFERS_SWOHYBRID)}

  \item{srad1}{Solar Monthly temperature, in Megajoules/M^2, for January (optional). It must be present to project with variant=2 (CONIFERS_SWOHYBRID)}
  \item{srad2}{Solar Monthly temperature, in Megajoules/M^2, for February (optional). It must be present to project with variant=2 (CONIFERS_SWOHYBRID)}
  \item{srad3}{Solar Monthly temperature, in Megajoules/M^2, for March (optional). It must be present to project with variant=2 (CONIFERS_SWOHYBRID)}
  \item{srad4}{Solar Monthly temperature, in Megajoules/M^2, for April (optional). It must be present to project with variant=2 (CONIFERS_SWOHYBRID)}
  \item{srad5}{Solar Monthly temperature, in Megajoules/M^2, for May (optional). It must be present to project with variant=2 (CONIFERS_SWOHYBRID)}
  \item{srad6}{Solar Monthly temperature, in Megajoules/M^2, for June (optional). It must be present to project with variant=2 (CONIFERS_SWOHYBRID)}
  \item{srad7}{Solar Monthly temperature, in Megajoules/M^2, for July (optional). It must be present to project with variant=2 (CONIFERS_SWOHYBRID)}
  \item{srad8}{Solar Monthly temperature, in Megajoules/M^2, for August (optional). It must be present to project with variant=2 (CONIFERS_SWOHYBRID)}
  \item{srad9}{Solar Monthly temperature, in Megajoules/M^2, for September (optional). It must be present to project with variant=2 (CONIFERS_SWOHYBRID)}
  \item{srad10}{Solar Monthly temperature, in Megajoules/M^2, for October (optional). It must be present to project with variant=2 (CONIFERS_SWOHYBRID)}
  \item{srad11}{Solar Monthly temperature, in Megajoules/M^2, for November (optional). It must be present to project with variant=2 (CONIFERS_SWOHYBRID)}
  \item{srad12}{Solar Monthly temperature, in Megajoules/M^2, for December (optional). It must be present to project with variant=2 (CONIFERS_SWOHYBRID)}

}
}

\details{The example plot data represents plot level observations for
a field sample that can be used to generate future forest conditions
using the CONIFERS forest growth model. The data are in imperial
units. The required values are the \code{plot}, \code{elevation}, and
\code{slope}, and \code{aspect} measurements. Values that have not
been observed, or are missing (e.g. whc and map), are denoted with an
\code{NA} and can be imputed using the \code{\link{impute}} function,
which the exception of the \code{si30}. The current version of the
library does not include the functions to compute \code{si30} from
tree observations in the \code{\link{impute}} function.}



\references{

Ritchie, M.W. 2008. User's Guide and Help System for CONIFERS: A Simulator for Young Conifer Plantations Version
4.10. See \url{http://www.fs.fed.us/psw/programs/ecology_of_western_forests/projects/conifers/}

Flewelling, J., R. Collier, B. Gonyea. Marshall D, and Turnblom
E. 2001. Height-Age Curves for Planted Stands of Douglas-fir, with
Adjustments for Density. (Working Paper 1). See
\url{http://www.cfr.washington.edu/research.smc/working_papers/smc_working_paper_1.pdf}.


}

\author{Jeff D. Hamann \email{jeff.hamann@forestinformatics.com},\cr
	     Doug Maguire \email{doug.maguire@oregonstate.edu},\cr
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

## load the rconifers library
library( rconifers )

## set the variant to the CIPS variant.
set.variant( 3 )

## load and set the species mapping
data( species.cips )
print( species.cips )
set.species.map( species.cips )

## load and display the plant data
data( plants.cips )
print( plants.cips )

## load and display the plot data
data( plots.cips )
print( plots.cips )

## create the sample.data list object
sample.cips <- list( plots=plots.cips, plants=plants.cips, age=3, x0=0.0, n.years.projected=0 )
class(sample.cips)  <- "sample.data"

## print and plot the sample
print( sample.cips )
plot( sample.cips )

}
