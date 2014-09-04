###
###	$Id: rconifers.r  2014-09-04 mritchie $	
###
###            R interface package for conifers growth model
###
### Copyright 2003-2008 Jeff D. Hamann <jeff.hamann@forestinformatics.com>
###
### This file is part of the rconifers library for R and related languages.
### It is made available under the terms of the GNU General Public
### License, version 2, or at your option, any later version,
### incorporated herein by reference.
###
### This program is distributed in the hope that it will be
### useful, but WITHOUT ANY WARRANTY; without even the implied
### warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
### PURPOSE.  See the GNU General Public License for more
### details.
###
### You should have received a copy of the GNU General Public
### License along with this program; if not, write to the Free
### Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
### MA 02111-1307, USA


## To submit to CRAN:
## $ cd to conifers directory
## $ R CMD check rconifers
## if check doesn't pass with flying colors, fix the errors and warnings until it does
## on unix:
## $ R CMD build rconifers
## on windows:
## $ R CMD INSTALL --build rconifers 
## $ ftp cran.r-project.org
## Connected to cran.wu-wien.ac.at.
## 220 Welcome to the CRAN FTP service.
## Name (cran.r-project.org:hamannj): anonymous
## 331 Please specify the password.
## Password: [hit return - no password needed]
## 230-Welcome, CRAN useR!
## 230-
## 230-If you have any unusual problems,
## 230-please report them via e-mail to <cran-sysadmin@statmath.wu-wien.ac.at>.
## 230-
## 230 Login successful.
## Remote system type is UNIX.
## Using binary mode to transfer files.
## ftp> cd incoming
## 250 Directory successfully changed.
## ftp> put rconifers_0.0-9.tar.gz
## local: rconifers_0.0-9.tar.gz remote: rconifers_0.0-9.tar.gz
## 229 Entering Extended Passive Mode (|||58381|)
## 150 Ok to send data.
## 100% |**************************************************************************************************************| 99741     146.79 MB/s    00:00 ETA
## 226 File receive OK.
## 99741 bytes sent in 00:04 (24.20 KB/s)
## ftp> quit
## 221 Goodbye.
## $
## send this email
#CRAN Team, 
#
#We have updated the rconifers young forest simulator package (i.e. rconifers_version_number.tar.gz). 
#
#Thanks,
#Jeff.

## to cran@r-project.org.

## to do this all at once,
## $ ftp


##.onAttach(libname, pkgname)


.onAttach <-
  function(lib, pkg)
{
  ## this package doesn't need any additional package...
  ## if it did, you'd modify the code below...
  ## but for now, simply display this little message
  ## or here?

##  library.dynam("rconifers", pkg, lib)

  ## put package checks here...
###     if(!require(ts, quietly = TRUE))
###         stop("Package ts is needed.  Stopping")
    
  ##mylib <- dirname(.path.package("rconifers"))
  mylib <- dirname(path.package("rconifers"))
  ver <- packageDescription("rconifers", lib.loc = mylib)["Version"]

  vertxt <- paste("\n`rconifers' version:", ver, "\n")

  introtxt <-
    paste("\n`rconfiers' is a package that provides an R interface\n",
          "to the CONIFERS forest growth model software.\n\n",
          "For more information about the CONIFERS forest growth model,\n",
          "paste this url into your web browser:\n",
          "http://www.fs.fed.us/psw/programs/ecology_of_western_forests/projects/conifers/\n\n",
          "The user is able to define arbitrary silvicultural prescriptions,\n",
          "using standard R scripts, to predict near term (<20 years) future\n",
          "forest conditions under a variety of model configurations and\n",
          "management decisions.\n\n",
          "The output can be used in standard graphics, file input and output,\n",
          "and data analysis functions using R.\n\n",
          "Jeff D. Hamann, Forest Informatics, Inc.,\n",
          "Martin W. Ritchie, USFS PSW Research Station, Redding Silviculture Lab, and\n",
          "Doug Maguire, Oregon State University, Corvallis, Oregon\n",
          "Doug Mainwaring, Oregon State University, Corvallis, Oregon\n\n",
	    ##"http://www.fsl.orst.edu/cips/index.htm\n",
          ##"Doug Maguire, Oregon State University, Corvallis, Oregon\n\n",
          "See `library (help=rconifers)' for details.\n\n",
          sep = "")

  ## if(interactive() || getOption("verbose")) {
  ##   cat(paste(vertxt, introtxt))
  ## }

  packageStartupMessage( paste(vertxt, introtxt) )
  
  if( .Call( "r_set_variant", 0, PACKAGE="rconifers" ) ) {
    stop("Rconifers Error: Package could not be loaded. Could not load default coefficients. Stopping...")
  }

  ## "species.swo"
  
  ## set the species mappings
  ##data( species.swo )

  ## sp=NULL
  ##data( species.swo, package="rconifers", envir=environment() )

## * checking R code for possible problems ... NOTE
## Found the following calls to data() loading into the global environment:
## File ‘rconifers/R/rconifers.r’:
##   data(species.swo)
## See section ‘Good practice’ in ‘?data’.

 ## 'data(..., envir =
 ##     environment())'.  However, two alternatives are usually
 ##     preferable, both described in the 'Writing R Extensions' manual.
 
  data( "species.swo", package="rconifers", envir=environment() )
  
  ##.localstuff <- new.env()
  sp.map <- list(idx=species.swo$idx,
                 fsp=species.swo$fsp,
                 code=as.character(species.swo$code),
                 em=species.swo$endemic.mort,
                 msdi=species.swo$max.sdi,
                 b=species.swo$browse.damage,
                 m=species.swo$mechanical.damage,
                 gwh=species.swo$genetic.worth.h,
                 gwd=species.swo$genetic.worth.d,
                 ## all variants will need to include these variables now.
                 mint=species.swo$min.temp,
                 maxt=species.swo$max.temp,
                 optt=species.swo$opt.temp
                 )


  ## sp.map <- list(idx=sp$idx,
  ##                fsp=sp$fsp,
  ##                code=as.character(sp$code),
  ##                em=sp$endemic.mort,
  ##                msdi=sp$max.sdi,
  ##                b=sp$browse.damage,
  ##                m=sp$mechanical.damage,
  ##                gwh=sp$genetic.worth.h,
  ##                gwd=sp$genetic.worth.d,
  ##                ## all variants will need to include these variables now.
  ##                mint=sp$min.temp,
  ##                maxt=sp$max.temp,
  ##                optt=sp$opt.temp
  ##                )

  ##      .localstuff <- new.env() , verbose=FALSE,PACKAGE="rconifers" ) ) {
  ## val <- .Call( "r_set_species_map", sp.map, verbose=FALSE, PACKAGE="rconifers" )


  
##  if( !.Call( "r_set_species_map", .localstuff <- new.env() , verbose=FALSE,PACKAGE="rconifers" ) ) {
  if( !.Call( "r_set_species_map", sp.map,verbose=FALSE, PACKAGE="rconifers" ) ) {
  stop("Rconifers Error: Package could not be loaded. Could not load swo species map. Stopping...")
  } else {
    introtxt <- paste("\nInitialized species map using data(swo)\nType help.start() for documentation\n", sep = "")

    packageStartupMessage( introtxt )

    ## if(interactive() || getOption("verbose")) {
    ##   cat( introtxt)
    ## }

  }  
  
  ## do I register my routines here?

}


.Last.lib <- function( lib ) {
 .Call( "exit_conifers", PACKAGE="rconifers" )
## library.dynam.unload( "rconifers" , lib)  
}


rand.seed <- function( control ) {
  val <- .Call( "r_reseed", control, PACKAGE="rconifers" )
  val
}



## this function will set the species codes (internal)
## sp.map is a list
set.species.map <- function( x, verbose=FALSE ) {

  if( class( x ) != "data.frame" ) {
    stop( "Rconifers Error: x is not a data.frame object." )
    return
  }
 
  sp.map <- list(idx=x$idx,
                 fsp=x$fsp,
                 code=as.character(x$code),
                 em=x$endemic.mort,
                 msdi=x$max.sdi,
                 b=x$browse.damage,
                 m=x$mechanical.damage,
                 gwh=x$genetic.worth.h,
                 gwd=x$genetic.worth.d,

                 ## added for the swo-hybrid variant
                 mint=x$min.temp,
                 maxt=x$max.temp,
                 optt=x$opt.temp

                 )
  
  ## verify the lengths match
  if( length( x$idx ) != length( x$fsp ) ) {
    stop( "Rconifers Error: The lengths of the index and functional species vectors do not match." )
    return
  }
  
  val <- .Call( "r_set_species_map", sp.map, verbose=FALSE, PACKAGE="rconifers" )
  ## no need to return the number of species assigned
}


set.variant <- function( var=0 ).Call( "r_set_variant", var, PACKAGE="rconifers" )

## temp hack to simply list the variants numbers and labels.
variants <- function() {
    ## print( "0=CONIFERS_SWO" )
    ## print( "1=CONIFERS_SMC" )
    ## print( "2=CONIFERS_SWOHYBRID" )
    ## print( "3=CONIFERS_CIPS" )

print( "#define CONIFERS_SWO            0" )
print( "#define CONIFERS_SMC            1" )
print( "#define CONIFERS_SWOHYBRID      2" )
print( "#define CONIFERS_CIPS           3" )
    
}

## todo: this code needs to be cleaned up so that users only need
## supply the function with the string label (e.g. set.variant( "conifers_swo" ) )
## set.variant <- function( var=0 ).Call( "r_set_variant", var, PACKAGE="rconifers" )
## and maybe that'l have to come from a variants data.frame data set?


## this function, which is used internally, to convert data from the c code
build.sample.data <- function( x ) {

  ## build the two data.frame objects and add them to the list
  ret.val <- vector( "list" )
  ret.val$x0 <- x[[1]]
  ret.val$age <- x[[2]]

  ## construct a data.frame for the plots

  ## todo: update the variables that get passed.
  plots <-   as.data.frame( cbind( plot=x[[3]][[1]],

                                  ## agh. must include latitude and longitude.
                                  lat=x[[3]][[2]],
                                  lon=x[[3]][[3]],
                                  
                                  elevation=x[[3]][[4]],
                                  slope=x[[3]][[5]],
                                  aspect=x[[3]][[6]],
                                  whc=x[[3]][[7]],
                                  map=x[[3]][[8]],
                                  si30=x[[3]][[9]],
                                  
                                  ## put the other plot data in here.
                                  gsp=x[[3]][[10]],

                                  ## include the mean monthly temps
                                  mt1=x[[3]][[11]],
                                  mt2=x[[3]][[12]],
                                  mt3=x[[3]][[13]],
                                  mt4=x[[3]][[14]],
                                  mt5=x[[3]][[15]],
                                  mt6=x[[3]][[16]],
                                  mt7=x[[3]][[17]],
                                  mt8=x[[3]][[18]],
                                  mt9=x[[3]][[19]],
                                  mt10=x[[3]][[20]],
                                  mt11=x[[3]][[21]],
                                  mt12=x[[3]][[22]],

                                  ## include the mean monthly temps
                                  srad1=x[[3]][[23]],
                                  srad2=x[[3]][[24]],
                                  srad3=x[[3]][[25]],
                                  srad4=x[[3]][[26]],
                                  srad5=x[[3]][[27]],
                                  srad6=x[[3]][[28]],
                                  srad7=x[[3]][[29]],
                                  srad8=x[[3]][[30]],
                                  srad9=x[[3]][[31]],
                                  srad10=x[[3]][[32]],
                                  srad11=x[[3]][[33]],
                                  srad12=x[[3]][[34]]

                                  )
                           )

  ## construct a data.frame for the plants
  plants <-   as.data.frame( cbind( plot=x[[4]][[1]] ,
                                   sp.code=x[[4]][[2]],
                                   d6=x[[4]][[3]],
                                   dbh=x[[4]][[4]],
                                   tht=x[[4]][[5]],
                                   cr=x[[4]][[6]],
                                   n.stems=x[[4]][[7]],
                                   expf=x[[4]][[8]],
                                   crown.width=x[[4]][[9]],
                                   errors=x[[4]][[10]]
                                   )
                            )

  ret.val$n.years.projected <- x[[5]]

  ## inserted the following line on 12/11/2009 to try and fix plot assignment
  plants$plot<- as.numeric( levels( plants$plot) )[plants$plot]
  plants$d6 <- as.numeric( levels( plants$d6 ) )[plants$d6]
  plants$dbh <- as.numeric( levels( plants$dbh ) )[plants$dbh]
  plants$tht <- as.numeric( levels( plants$tht ) )[plants$tht]
  plants$cr <- as.numeric( levels( plants$cr ) )[plants$cr]
  plants$n.stems <- as.numeric( levels( plants$n.stems ) )[plants$n.stems]
  plants$expf <- as.numeric( levels( plants$expf ) )[plants$expf]
  plants$crown.width <- as.numeric( levels( plants$crown.width))[plants$crown.width]
  plants$errors <- as.numeric( levels( plants$errors))[plants$errors]

  ret.val$plots <- plots
  ret.val$plants <- plants

  class(ret.val)  <- "sample.data"
  ret.val

}

################################################################################
## file i/o interfaces for C code
################################################################################
project <- function( x,
                    years=1,
                    control=list(rand.err=0,
                      rand.seed=0,
                      endemic.mort=0,
                      sdi.mort=0,
                      genetic.gains=0) )
{

  ## make sure the class of the object pass into the function is a "sample.data" object
  if( class( x ) != "sample.data" ) {
    stop( "Rconifers Error: x is not a sample.data object." )
    return
  }
  
  ## make sure the plant level variables are available
  ## todo: this might have to be modified for each variant
  if( sum( names(  x$plants ) %in% c("plot","sp.code","d6","dbh","tht","cr","n.stems","expf","crown.width" ) ) != 9 )
    {
    stop( "Rconifers Error: the plant list data.frame does not have all the required columns. See impute help (?impute)" )
    return
  }

  ## make sure the plot level variables are available
  ## todo: this might have to be modified for each variant
  if( sum( names(  x$plots ) %in% c("plot","elevation","slope","aspect","whc","map","si30" ) ) != 7 )
    {
      stop( "Rconifers Error: the plot data.frame does not have all the required columns. See impute help (?impute)" )
      return
    }

  x$plants$sp.code <- as.character( x$plants$sp.code )
  val <- build.sample.data( .Call( "r_project_sample",  x, years, control, PACKAGE="rconifers" ) )
  val
}


################################################################################
## file i/o interfaces for C code
################################################################################

thin <- function( x,
                 control=list(type=2,
                   target=50.0,
                   target.sp=NULL ) )
{

  if( class( x ) != "sample.data" ) {
    stop( "Rconifers Error: x is not a sample.data object." )
    return
  }

  if( control$type < 1 || control$type > 4 ) {
    stop( "Rconifers Error: control$type must be 1, 2, 3, or 4. See help." )
    return    
  }
  
  if( control$type %in% c(2,3) & !is.null( control$target.sp ) ) {
    stop( "Rconifers Error: If you want to thin all species, sp=NULL." )
    return    
  }

  if( control$type %in% c(1,4) & is.null( control$target.sp ) ) {
    stop( "Rconifers Error: If you want to thin a particular species, sp != NULL." )
    return    
  }

  ## make sure the sp.code values are not factors.
  x$plants$sp.code <- as.character( x$plants$sp.code )
  val <- build.sample.data( .Call( "r_thin_sample", x, control, PACKAGE="rconifers" ) )
  val
}

## fill in missing values wrapper function
impute <- function( x,
                   control=list(fpr=11.78,
                     min.dbh=5.6,
                     baf=40.0) )
{
  
  if( class( x ) != "sample.data" ) {
    stop( "Rconifers Error: x is not a sample.data object." )
    return
  }

  if( sum( names(  x$plants ) %in% c("plot","sp.code","d6","dbh","tht","cr","n.stems","expf","crown.width" ) ) != 9 )
    {
    stop( "Rconifers Error: the plant list data.frame does not have all the required columns. See impute help (?impute)" )
    return
  }

  if( sum( names(  x$plots ) %in% c("plot","elevation","slope","aspect","whc","map","si30" ) ) != 7 )
    {
      stop( "Rconifers Error: the plot data.frame does not have all the required columns. See impute help (?impute)" )
      return
    }
     
  x$plants$sp.code <- as.character( x$plants$sp.code )
  val <- build.sample.data( .Call( "r_impute_missing_values",
                                  x,
                                  control,
                                  PACKAGE="rconifers" ) )
  val
  
}


calc.max.sdi <- function( x )
{
  
  if( class( x ) != "sample.data" ) {
    stop( "Rconifers Error: x is not a sample.data object." )
    return
  }

    if( sum( names(  x$plants ) %in% c("plot","sp.code","d6","dbh","tht","cr","n.stems","expf","crown.width" ) ) != 9 )
    {
    stop( "Rconifers Error: the plant list data.frame does not have all the required columns. See impute help (?impute)" )
    return
  }

  if( sum( names(  x$plots ) %in% c("plot","elevation","slope","aspect","whc","map","si30" ) ) != 7 )
    {
      stop( "Rconifers Error: the plot data.frame does not have all the required columns. See impute help (?impute)" )
      return
    }

  x$plants$sp.code <- as.character( x$plants$sp.code )
  val <- .Call( "r_calc_max_sdi", x, PACKAGE="rconifers" )
  
  val
}


################################################################################
## print a few results of the whole system
print.sample.data <- function( x, digits = max( 3, getOption("digits") - 1 ),... ) {

    if( class( x ) != "sample.data" ) {
    stop( "Rconifers Error: x is not a sample.data object." )
    return
  }

      if( sum( names(  x$plants ) %in% c("plot","sp.code","d6","dbh","tht","cr","n.stems","expf","crown.width" ) ) != 9 )
    {
    stop( "Rconifers Error: the plant list data.frame does not have all the required columns. See impute help (?impute)" )
    return
  }

  if( sum( names(  x$plots ) %in% c("plot","elevation","slope","aspect","whc","map","si30" ) ) != 7 )
    {
      stop( "Rconifers Error: the plot data.frame does not have all the required columns. See impute help (?impute)" )
      return
    }

  save.digits <- unlist(options(digits=digits))
  on.exit(options(digits=save.digits))
  
  cat("\n")
  cat( "sample contains", nrow(x$plots), "plots records\n" )
  cat( "sample contains", nrow(x$plants), "plant records\n" )

  cat( "n.years.projected = ", x$n.years.projected, "\n")
  cat( "age = ", x$age, "\n")
  cat( "x0 = ", x$x0, "\n")
  cat( "max sdi = ", calc.max.sdi( x ), "\n")

  ## print out some species summaries...
  print( sp.sums( x ) )

  if( sum( x$plants$errors ) ) {

    ## iterate over the plant records and report how many of what kinds
    ## of errors there are in the plant list...
    ## the #defines in conifers.h contain the masks for the type
    ## of errors
    ## this also means that you have to have a dependent package (bitops)
    
    print( "plant records contain errors..." )
  }
  
  ##   print( stand.table( x ) )  
  invisible( x )
}


summary.sample.data <- function(object,...) {
  summary.sample.data <- object
  ##summary.sample.data$stand.table <- stand.table( object )
  summary.sample.data
}





## split and sum the plants by species
sp.sums <- function( x ) {

  if( class( x ) != "sample.data" ) {
    stop( "Rconifers Error: x is not a sample.data object." )
    return
  }

    if( sum( names(  x$plants ) %in% c("plot","sp.code","d6","dbh","tht","cr","n.stems","expf","crown.width" ) ) != 9 )
    {
    stop( "Rconifers Error: the plant list data.frame does not have all the required columns. See impute help (?impute)" )
    return
  }

  if( sum( names(  x$plots ) %in% c("plot","elevation","slope","aspect","whc","map","si30" ) ) != 7 )
    {
      stop( "Rconifers Error: the plot data.frame does not have all the required columns. See impute help (?impute)" )
      return
    }

  x.by.sp <- split( x$plants, f=list(x$plants$sp.code ) )

  ## where does npts come from???
  ## this might be a bug.
  ##npts <- x#npts
  ## didn't pass check. fixed.
  npts <- length( unique( x$plots$plot ) )
  
  ## this function computes the summaries each species.
  sp.sums.f <- function( x, npts ) {
    expf <- sum( x$expf * x$n.stems ) / npts
    tht <- sum( x$expf * x$n.stems * x$tht ) / sum( x$expf * x$n.stems )
    ba <- sum( x$expf * x$n.stems * x$dbh^2*0.0054541539 ) / npts
    expgtbh<-sum(x$n.stems[x$tht>4.5]*x$expf[x$tht>4.5])/npts
    qmd <- sqrt( ba / expgtbh / 0.0054541359 )
    sp.sums <- c(qmd,tht,ba,expgtbh,expf)
  }
  
  df <- as.data.frame( t(sapply( x.by.sp, sp.sums.f, nrow(x$plots) )) )
  names( df ) <- c("qmd","tht","ba","bhexpf","texpf")

  df
}


## todo: this needs a manual page
## this function generates a simple set of charts to visually represent the data
plot.sample.data <- function( x, digits = max( 3, getOption("digits") - 1 ),... ) {

    if( class( x ) != "sample.data" ) {
    stop( "Rconifers Error: x is not a sample.data object." )
    return
  }

    if( sum( names(  x$plants ) %in% c("plot","sp.code","d6","dbh","tht","cr","n.stems","expf","crown.width" ) ) != 9 )
      {
        stop( "Rconifers Error: the plant list data.frame does not have all the required columns. See impute help (?impute)" )
        return
      }
    
    if( sum( names(  x$plots ) %in% c("plot","elevation","slope","aspect","whc","map","si30" ) ) != 7 )
      {
        stop( "Rconifers Error: the plot data.frame does not have all the required columns. See impute help (?impute)" )
        return
      }
    
    save.digits <- unlist(options(digits=digits))
    on.exit(options(digits=save.digits))


    ## generate a plot of the species summaries
    opar <- par(mfrow=c(2,2))
    
    ## plot the diameter histogram
    hist( x$plants$dbh, main="Diameter Distribution", xlab="DBH, in Inches" )
    
    ## plot the height histogram
    hist( x$plants$tht, main="Height Distribution", xlab="Total Height, in Feet" )

    ## plot the ht vs. dbh plot
    plot( x$plants$tht ~ x$plants$dbh,
         ylab="Total Height, in Feet",
         xlab="Diameter at Breast Height, in Inches",     
         main="Height vs. Diameter" )
    grid()

    ## generate a pie chart of the basal area, if there is any
    ba <- sp.sums( x )$ba
    if( sum( ba ) > 0 ) {
      names(ba) <- rownames( sp.sums( x ) )
      pie( ba[ba>0], main="Basal Area" )
    }

    ## close out the plot.
    par(opar)
    
    invisible( x )
  }


