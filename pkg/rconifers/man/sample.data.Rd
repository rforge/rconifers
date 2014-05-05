

\name{sample.data}
\Rdversion{1.1}
\title{
rconifers Objects
}
\alias{sample.data}

\description{
  In \code{rconifers}, forest sample data sets are stored in an S3 class
called \code{sample.data}, which contains a variety of information in
addition to the plant observations themselves.  \code{sample.data} objects are created
by \code{\link{class}} or similar functions (no others currently
exist).

}

\section{Structure}{
The structure of a \code{sample.data} object is a list of 5 elements and an attribute as follows:

\tabular{lll}{
\emph{element} \tab \emph{type} \tab \emph{description}\cr
$freq \tab num \tab A common frequency (or wavelength) axis for all the spectra.\cr
$data \tab num \tab The intensities for the spectra.  A matrix of dimension \cr
\tab \tab no. samples x no. frequency points.\cr
$names \tab chr \tab The sample names for the spectra; length must be no. samples.\cr
$groups \tab Factor \tab The group classification of the samples; length must be no. samples.\cr
$colors \tab chr \tab The colors for each sample; length must be no. samples. \cr
\tab \tab Groups and colors correspond.\cr
$sym \tab integer \tab As for \code{colors}, but symbols for plotting (if b/w is desired).\cr
$alt.sym \tab chr \tab Lower-case letters as alternate symbols for plotting.\cr
$unit \tab chr \tab Two entries, the first giving the x axis unit, the second the y axis unit.\cr
$desc \tab chr \tab A character string describing the data set.  This appears on plots and therefore \cr
\tab \tab should probably be kept to 40 characters or less.\cr
- attr \tab chr "Spectra" \tab The S3 class designation.\cr
}
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
options(digits=10)

library( rconifers )

## Example for SWO variant
## set the variant to the SWO variant and set species map
##set.species.map( set.variant(0) )
set.variant(0)

## grow the data that was originally swo in the smc variant
# load and display CONIFERS example plots
data( plots.swo )
print( plots.swo )

# load and display CONIFERS example plants
data( plants.swo )
print( plants.swo )

# create the sample.data list object
sample.swo.3 <- list( plots=plots.swo, plants=plants.swo, age=3, x0=0.0,n.years.projected=0)
class(sample.swo.3)  <- "sample.data"

## print a summary 
print( sample.swo.3 )

## generate a basic summary plot
plot( sample.swo.3 )

}


\keyword{ classes }
