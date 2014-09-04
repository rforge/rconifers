
%%	$Id: species.cips.Rd 861 2012-02-22 21:25:32Z hamannj $	

\keyword{datasets}
\name{species.cips}
\alias{species.cips}
\non_function{}
\title{Species table for the Center for Intensive Planted-forest Silviculture (CIPS) variant of the CONIFERS growth model}
\usage{data(species.cips)}

\description{
These are species codes used for the Center for Intensive Planted-forest Silviculture (CIPS) variant of the CONIFERS growth model.

The \code{species.smc} data frame has 3 rows and 13 columns.
}

\format{
This data frame contains the following columns:
\describe{

\item{idx}{The index of the species code that is used within the
  CONIFERS source code.}
\item{code}{The alpha-numeric code that is used during data input,
  normally in the field data recorder. These are used in the
  summaries and outputs rather than the species names.}
\item{fsp}{The functional species code (found in the coefficients)
  that is used to project the plant record. Every species code is
  mapped to a single functional species code. See \code{\link{set.species.map}}.}
\item{name}{The long name for the plant species.}
\item{organon}{The code that should be used to export the species to
  the ORGANON growth and yield model. See \url{http://www.cof.orst.edu/cof/fr/research/organon/index.htm} for more details}
\item{cactos}{The code that should be used to export the species to
  the CACTOS growth and yield model. See \url{http://www.cnr.berkeley.edu/~wensel/cactos/cactoss.htm} for
  more details}
\item{fvs}{The code that should be used to export the species to
  the FVS growth and yield model. See \url{http://www.fs.fed.us/fmsc/index.php} for more
  details}
\item{endemic.mort}{ is a percent for the endemic mortality, where
(0.002) means that two tenths of a percent of the trees of that
functional species type will be "killed" for each year of projection.}
\item{max.sdi}{ is the maximum stand density index the functional
  species can obtain before competition induced mortality is applied.}
\item{browse.damage}{ is a percent for the browse damage,
  where (0.002) means that two tenths of a percent of the trees of that
  species type will be "killed" as a result of browse damage.}
\item{mechanical.damage}{ is a percent for the mechanical damage,
  where (0.002) means that two tenths of a percent of the trees of that
  species type will be "killed" as a result of mechanical damage.}
\item{genetic.worth.h}{ is the average expected level of gain,
expressed as percentage, for diameter at a given age. Values are
bounded by 0 and 20.}

\item{genetic.worth.d}{ is the average expected level of gain,
expressed as percentage, for diameter at a given age. Values are
bounded by 0 and 20.}

\item{min.temp}{ is the minimum temperature, in C, during which photosynthasis will take place for a given species.}
\item{max.temp}{ is the maximum temperature, in C, during which photosynthasis will take place for a given species.}
\item{opt.temp}{ is the optimal temperature, in C, during which photosynthasis will take place for a given species.}

}
}

\details{
  The dataset is similar to the smc.txt file that is distributed
  with the GUI version of CONIFERS. A species lookup table is nothing
  more than a data.frame with specific column names. 
}

\references{

Ritchie, M.W. 2008 . User's Guide and Help System for CONIFERS: A
Simulator for Young Conifer Plantations Version 4.10. See
\url{http://www.fs.fed.us/psw/programs/ecology_of_western_forests/projects/conifers/}

Maguire, Doug. 2010. Center for Intensive Planted-Forest Silviculture
Annual Report. Oregon State University. 60 p.

}

\author{Jeff D. Hamann \email{jeff.hamann@forestinformatics.com},\cr
	     Martin W. Ritchie \email{mritchie@fs.fed.us},\cr and
	     Doug Maguire \email{doug.maguire@oregonstate.edu} }


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

## set the variant to the CIPS variant.
## the set.variant command also assigns the species map to the default species map data.frame
set.variant( 3 )

## load the species coefficients into R as a data.frame object and
## set the species map
data( species.cips )
print( species.cips )
set.species.map( species.cips )

## load and display CONIFERS example CIPS plots
data( plots.cips )
print( plots.cips )

## load and display CONIFERS example CIPS plants
data( plants.cips )
print( plants.cips )

## create the sample.data list object
sample.cips <- list( plots=plots.cips, plants=plants.cips,
  age=3, x0=0.0, n.years.projected=0 )
class(sample.cips)  <- "sample.data"

## print and plot the sample data
print( sample.cips )
plot( sample.cips )

## plot a sample projected 5 years into the future
plot( project( sample.cips, 5, control=list(rand.err=1,
                     rand.seed=615,
                     endemic.mort=1,
                     sdi.mort=1,
                     genetic.gains=0) ) )


## create a loop and store the projected values
s0 <- sample.cips
sums <- as.data.frame( cbind( age=s0$age, sp.sums( s0 )["DF",] ) )
for( i in 1:20 ) {
     s1 <- project( s0, 1, 
        control=list(rand.err=0,
                     rand.seed=0,
                     endemic.mort=1,
                     sdi.mort=1,
                     genetic.gains=0))

     sums <- rbind( sums, as.data.frame( cbind( age=s1$age, sp.sums( s1 )["DF",] ) ) )
     s0 <- s1
}

rownames( sums ) <- 1:nrow(sums)


## create a plot with four figures.
opar <- par( mfrow=c(2,2))

## now, plot a few stand level variables.
plot( sums$age, sums$qmd, type="b" )
grid()

plot( sums$age, sums$tht, type="b" )
grid()

plot( sums$age, sums$ba, type="b" )
grid()

plot( sums$age, sums$texpf, type="b" )
grid()

par( opar )

title( "CONIFERS_CIPS Simulations" )


}

\keyword{datasets}


