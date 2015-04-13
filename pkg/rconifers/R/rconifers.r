## R Conifers: R Interface to the CONIFERS Model

# Tasks to complete after loading the package (ex. library, require)
.onAttach <- function(lib, pkg){

	# Identify path to rconifers in the user library
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
				"Nate Osborne, Oregon State University, Corvallis, Oregon\n",
				"Doug Mainwaring, Oregon State University, Corvallis, Oregon\n\n",
				"See `library (help=rconifers)' for details.\n\n",
				sep = "")

	packageStartupMessage( paste(vertxt, introtxt) )
  
	# Check to make sure rconifers could be properly loaded: stop if else
	if( .Call( "r_set_variant", 0, PACKAGE="rconifers" ) ) {
		stop("Rconifers Error: Package could not be loaded. Could not load default coefficients. Stopping...")
	}
	
	# Load the species.swo on attach: this could be done using LazyLoad.
	data( "species.swo", package="rconifers", envir=environment() )
  
	# Create a species map object
	sp.map <- list(
			# Standard variables for all variants
			idx=species.swo$idx, 
			fsp=species.swo$fsp, 
			code=as.character(species.swo$code),
            em=species.swo$endemic.mort, 
			msdi=species.swo$max.sdi, 
			b=species.swo$browse.damage,
            m=species.swo$mechanical.damage, 
			gwh=species.swo$genetic.worth.h, 
			gwd=species.swo$genetic.worth.d,
            
			# New variables for the southwest Oregon variant
			mint=species.swo$min.temp, 
			maxt=species.swo$max.temp, 
			optt=species.swo$opt.temp
     )

	# Test if you can call to set the species map: if not stop otherwise provided introduction message
  	if( !.Call( "r_set_species_map", sp.map,verbose=FALSE, PACKAGE="rconifers" ) ) {
  		stop("Rconifers Error: Package could not be loaded. Could not load swo species map. Stopping...")
  	} else {
    	introtxt <- paste("\nInitialized species map using data(swo)\nType help.start() for documentation\n", sep = "")
    	packageStartupMessage( introtxt )
	}   
}

# Call to close down conifers at end of procedure
.Last.lib <- function( lib ) {
	.Call( "exit_conifers", PACKAGE="rconifers" )
}

# Make a random seed for the simulation
rand.seed <- function( control ) {
  val <- .Call( "r_reseed", control, PACKAGE="rconifers" )
  val
}

# This function essentially 'checks' the plot data object (internal): x is an object formatted like species.swo
# Ex. Need to set the species map before calling functions downsteam, set.species.map(species.swo)
set.species.map <- function( x, verbose=FALSE ) {

  if( class( x ) != "data.frame" ) {
    stop( "Rconifers Error: x is not a data.frame object." )
    return
  }
 
  sp.map <- list(
		  		# Standard variables for rconifers
		  		idx=x$idx,
          		fsp=x$fsp,
          		code=as.character(x$code),
				em=x$endemic.mort,
          		msdi=x$max.sdi,
          		b=x$browse.damage,
          		m=x$mechanical.damage,
          		gwh=x$genetic.worth.h,
          		gwd=x$genetic.worth.d,

          		# New variables for the southwest Oregon variant
         		mint=x$min.temp,
             	maxt=x$max.temp,
                optt=x$opt.temp,
		 
		 		# New variables for the CIPS variant
				yrst = rep(1,nrow(x)) # Set equal to one for now
	)
  
  	# Verify the lengths match
 	if( length( x$idx ) != length( x$fsp ) ) {
    	stop( "Rconifers Error: The lengths of the index and functional species vectors do not match." )
    	return
	}
  	
	val <- .Call( "r_set_species_map", sp.map, verbose=FALSE, PACKAGE="rconifers" )
}

# Set the variant of rconifers to use internal to c-code
set.variant <- function( var=0 ).Call( "r_set_variant", var, PACKAGE="rconifers" )

# Function will define the kinds of rconifers avaliable for use (better in the documentation called using: ?)
variants <- function() {
	print( "#define CONIFERS_SWO            0" ) # Southwest Oregon variant
	print( "#define CONIFERS_SMC            1" ) # Stand Management Cooperative variant
	print( "#define CONIFERS_SWOHYBRID      2" ) # Southwest Oregon variant
	print( "#define CONIFERS_CIPS           3" ) # CIPS variant
}

# This function just combines the plant and plot data into a special class of object: 
# For example: x = list( plots=plots.swo, plants=plants.swo, age=3, yrst=3, x0=0.0, n.years.projected=0 )
build.sample.data <- function( x ) {
  
  	# * MUST BE ORDERED EXPLICITY in LIST and OBJECTS WITHIN LIST * 
		
	# Make a plots data frame
	plots = as.data.frame(matrix(0,nrow=nrow(x$plots), ncol=34))
	names(plots) <- c("plot","lat","lon","elevation","slope","aspect","whc","map","si30","gsp",paste("mt",1:12,sep=""),paste("srad",1:12,sep=""))
	
	# Assign information as possible (input avaliable differs by variant!)
	sapply(1:ncol(plots),function(k){
				z = names(plots)[k]
				if(z %in% names(x$plots)){
					z = get(names(plots)[k],x$plots)
					plots[,k] <<- z
				}
			}
	)
	
	# Make a plants data frame
	plants = as.data.frame(matrix(0,nrow=nrow(x$plants), ncol=10))
	names(plants) <- c("plot","sp.code","d6","dbh","tht","cr","n.stems","expf","crown.width","errors")
					
	# Assign information as possible (input avaliable differs by variant!)
	sapply(1:ncol(plants),function(k){
				z = names(plants)[k]
				if(z %in% names(x$plants)){
					z = get(names(plants)[k],x$plants)
					plants[,k] <<- z
				}
			}
	)
	
	# Make a blank list and then populate it
	ret.val = list(x0=NA,age=NA,plots=NA,plants=NA,n.years.projected=NA,yrst=NA)
	ret.val$x0 = x$x0 
	ret.val$age = x$age
	ret.val$plots = plots
	ret.val$plants = plants 
	ret.val$n.years.projected = x$n.years.projected
  	ret.val$yrst <- x$yrst
  	class(ret.val)  <- "sample.data"
  	
	return(ret.val)
}

# Make the output information from called to rconifer.dll clean 
process.output.data <- function( x ) {
	
	## construct a data.frame for the plots
	
	## todo: update the variables that get passed.
	plots <-   as.data.frame( 
	
			cbind( plot=x[[3]][[1]],
					
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
	
	ret.val = list(x0=NA,age=NA,plots=NA,plants=NA,n.years.projected=NA,yrst=NA)
	ret.val$x0 <- x[[1]]
	ret.val$age <- x[[2]]
	ret.val$plots <- droplevels(plots)
	ret.val$plants <- droplevels(plants)
	ret.val$n.years.projected <- x[[5]]
	ret.val$yrst <- x[[6]]
	
	class(ret.val)  <- "sample.data"

	return(ret.val)
}

# Project the plant list into the future
# Default conditions: one year of growth, no random error, seed, endemic mort.. etc.
project <- function(x,years=1,control=list(rand.err=0,rand.seed=0,endemic.mort=0,sdi.mort=0,genetic.gains=0)){
	
	  	# Make sure the class of the object passed into the function is a "sample.data" object
	  	if( class( x ) != "sample.data" ) {
			stop( "rconifers Error: x is not a sample.data object." )
			return
		}
	  
	  # Make sure the plant level variables are available
	  # To Do!: this might have to be modified for each variant 
	  if( sum( names(x$plants) %in% c("plot","sp.code","d6","dbh","tht","cr","n.stems","expf","crown.width" ) ) != 9 ){
	    	stop( "Rconifers Error: the plant list data.frame does not have all the required columns. See impute help (?impute)" )
	    	return
	  }

	  # Make sure the plot level variables are available
	  ## To Do: this might have to be modified for each variant
	  if( sum( names(  x$plots ) %in% c("plot","elevation","slope","aspect","whc","map","si30" ) ) != 7 ){
			stop( "Rconifers Error: the plot data.frame does not have all the required columns. See impute help (?impute)" )
			return
	  }

	  x$plants$sp.code <- as.character( x$plants$sp.code )
	  val <- process.output.data( .Call( "r_project_sample",  x, years, control, PACKAGE="rconifers" ) )
	  val
}

# Thinning function for rconifers
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
  val <- process.output.data( .Call( "r_thin_sample", x, control, PACKAGE="rconifers" ) )
  val
}

# Fill in missing values wrapper function 
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
  val <- process.sample.data( .Call( "r_impute_missing_values",
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
   
    print( "plant records contain errors..." )
  }
  
  invisible( x )
}

summary.sample.data <- function(object,...) {
  summary.sample.data <- object
  summary.sample.data
}

# Split and sum the plants by species
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
  
  df$ba[is.na(df$ba)] <- 0
  df$qmd[is.na(df$qmd)] <- 0
    
  return(df)
}

## To Do!: This needs a manual page
# This function generates a simple set of charts to visually represent the data
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
    if ( sum( ba ) > 0 ) {
      names(ba) <- rownames( sp.sums( x ) )
      pie( ba[ba>0], main="Basal Area" )
    }
    else{
      plot.new()
      legend("center", c("BA Equal Zero"), cex=3, bty='n')
    }
    
    ## close out the plot.
    par(opar)
    
    invisible( x )
}

# Function to treat stands competing vegetation provided the sample data and residual vegetation cover target (percent)
vegman <- function(sample,target){
	
	veg = subset(sample$plants,sp.code=="CV") # Isolate just the vegetation entry
	
	# Numerically converge on expansion factor yielding desired canopy cover level
	for(int in c(1000,100,10,1,0.001)) {
		repeat{
			cc = with(veg,(3.1416*(crown.width/2)^2*expf)/43560*100) # Calculate current canopy cover
			if(cc<target) {veg$expf=veg$expf+int; break} # Determine if canopy cover is at the target level
			veg$expf=veg$expf-int # Subtract vegetation expansion factor by an increasing reduction intension
		}
	}
	
	sample$plants$expf[sample$plant$sp.code=="CV"] <- veg$expf # Input new expansion factor for vegetation density
	
	return(sample) # Return the treated sample information
}

