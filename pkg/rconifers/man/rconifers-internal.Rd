
% 	$Id: rconifers-internal.Rd 850 2012-01-11 21:30:38Z hamannj $	

\name{rconifers-internal}
\alias{build.sample.data}

\title{Internal rconifers functions}
\description{
  Internal rconifers functions
}
\usage{
build.sample.data( x )
}

\arguments{
  \item{x}{a CONIFERS sample.data object}
}

\details{
  These functions are not the called by the user.
  The build.sample.data is the function that converts the data, passed
  from the shared library, back to R. It is called in most of the user
  interface functions.} 


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

Ritchie, M.W. 2008 . User's Guide and Help System for CONIFERS: A Simulator for Young Conifer Plantations Version
4.10. See \url{http://www.fs.fed.us/psw/programs/ecology_of_western_forests/projects/conifers/}

}

\author{Jeff D. Hamann \email{jeff.hamann@forestinformatics.com} and Martin W. Ritchie \email{mritchie@fs.fed.us} }




\keyword{internal}
