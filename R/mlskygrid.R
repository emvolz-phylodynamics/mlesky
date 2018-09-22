invisible('#
DONE
x cross validation
x maxheight option 
x print and plot methods 

TODO 
x ?bootstrap CIs 
x add covars for 
	- size
	- logsize 
	- growth 
')


# derive timeseries of coalescent and ltt along appropriate time axis 
.tre2df <- function( tre, res, maxHeight = Inf, minLTT = NULL ){
	n <- Ntip( tre )
	if (is.null( minLTT)) 
	  minLTT <- floor( n / 5 ) 
	
	D <- ape::node.depth.edgelength( tre )
	rh <- max( D[1:n] )
	sts <- D[1:n]
	maxHeight <- min( rh, maxHeight )
	
	ne_haxis <- seq( maxHeight/res ,maxHeight, le = res )
	shs <- rh - sts 
	inhs <- rh - D[ (n+1):(n + tre$Nnode) ]
	u_shs <- unique( shs ) 
	u_inhs <- unique( inhs )
	
	#< h , event, ltt(descending), intervallength, nco, likterm, ne_bin >
	tredat <- data.frame( h= c( u_shs, u_inhs, ne_haxis) 
	 , type = c( rep('sample', length( u_shs))
	           , rep('node', length(u_inhs))
	           , rep('neswitch', length(ne_haxis))
	          )
	)
	tredat <- tredat[ tredat$h <= maxHeight , ]
	
	tredat$ne_bin <- sapply( tredat$h, function(x) sum( ne_haxis  < x ) + 1)
	
	ltt.h <- function(h) sum( shs < h ) - sum( inhs < h )
	tredat$ltt <- sapply( tredat$h, ltt.h )
	
	tredat$nco <- 0
	tredat$nco[ tredat$type=='node' ] <- sapply( tredat$h[ tredat$type=='node'  ], function(x) sum( x == inhs ))
	
	tredat <- tredat[ order( tredat$h ), ]
	tredat$intervalLength <- c( 0, diff( tredat$h ))
	tredat$ltt_terms <- tredat$ltt * (tredat$ltt-1) / 2
	
	tredat$dh <- ne_haxis[2] -  ne_haxis[1] 
	
	# combine ne bins:
	done <- ( minLTT == 1 )
	while(!done){
		done <- TRUE
		bin2maxltt <- sapply( 1:res, function( bin ){
			i = which(  tredat$ne_bin == bin )
			if (length(i) > 0)
				return( max( tredat$ltt[ i ] ) )
			-Inf
		})
		for (i in 1:nrow(tredat)){
			if ( bin2maxltt[ tredat$ne_bin[i] ] < minLTT ){
				tredat$ne_bin[i] <- max( 1, tredat$ne_bin[i] - 1 )
				done <- FALSE
			}
		}
	}
	tredat	
}



#' Objective function for cross validation; computes out-sample-iog likelihood and takes mean of all crosses 
.mlskygrid_oos <- function( tau 
  , tredat
  , ne0 
  , res = 50 
  , maxHeight = Inf 
  , quiet = FALSE
  , control = NULL
  , ncross = 5
  , ncpu = 1
){
	if ( ncross < 2 ) stop('*ncross* must be > 1')
	
	ne <- rlnorm( res , log( ne0 ), .2 ) # add some jitter
	
	dh <- tredat$rh[1] / res 
	
	cvsets = lapply( 1:ncross, function(icross) seq( icross, nrow(tredat), by = ncross ) )
	
	roughness_penalty <- function(logne){
		sum( dnorm( diff(diff( logne)), 0, sd = sqrt(dh/tau), log = TRUE) )
	}
	
	lterms <- function(logne)
	{
		ne <- exp( logne )
		sterms <- with(tredat, {
			-intervalLength * ltt_terms / ne [ ne_bin ] 
		})
		coterms <- with(tredat, {
			nco * ( log( ltt_terms ) - logne[ ne_bin ]  )
		})
		coterms[ is.na(coterms)] <- 0
		coterms + sterms 
	}
	
	of.cv.oos <- function(logne, icross  ){
		sum( lterms( logne )[ cvsets[[icross]] ] ) + roughness_penalty( logne )
	}
	
	of.cv.ws <- function(logne, icross){
		i = setdiff(1:nrow(tredat) , cvsets[[icross]]  )
		sum( lterms( logne )[ i ] ) + roughness_penalty( logne )
	}
	
	fits <- parallel::mclapply( 1:ncross, function(icross){
		optim( par = log(ne), fn = of.cv.ws
		  , method = 'BFGS'
		  , control = list( trace = ifelse( quiet ,0, 1), fnscale  = -1 )
		  , hessian = FALSE 
		  , icross = icross 
		)
	}, mc.cores = ncpu )
	oos_perfs <- sapply( 1:ncross, function(i){
		logne = fits[[i]]$par 
		of.cv.oos( logne, i )
	})
	
	sum( oos_perfs )
}



#' Maximum likelihood non-parametric estimation of effective population size through time
#'
#' 
#'
#' @param tre A dated phylogeny in ape::phylo format (see documentation for ape)
#' @param res Length of time axis over which to estimate Ne(t) (integer)
#' @param tau Precision parameter. Larger values generate smoother trajectories of Ne(t). If NULL, will optimize using cross-validation.
#' @param tau_lower Lower bound for precision parameter if estimating
#' @param tau_upper Upper bound for precision parameter if estimating
#' @param tau_tol Optimization tolerance when optimizing tau by cross-validation
#' @param ncross Number of folds in cross-validation
#' @param ncpu If doing cross-validation, each fold will be handled in parallel if ncpu > 1 (see parallel package)
#' @param quiet Provide verbose output from optimizer? 
#' @param NeStartTimeBeforePresent If <Inf, will only estimate Ne(t) between the most recent sample and this time before the most recent sample
#' @return A fitted model including effective size through time
#' @export
# @examples
# library(mlskygrid)
# tree <- read.tree( system.file( package='mlskygrid', 'mrsa.nwk' , mustWork=TRUE) ) 
# print( (fit <- mlskygrid( tree, tau = 10, NeStartTimeBeforePresent = 15) ))
# plot( fit , logy = FALSE)
mlskygrid <- function(tre
  , res = 25 
  , tau = 1
  , tau_lower = NULL 
  , tau_upper = NULL
  , tau_tol = 1e-3 
  , ncross = 5
  , ncpu = 1
  , quiet = FALSE
  , NeStartTimeBeforePresent = Inf 
  , minLTT = 1 
  , ne0 = NULL
){
	tredat <- .tre2df( tre = tre, res = res , maxHeight= NeStartTimeBeforePresent, minLTT = minLTT )
	if ( is.null( tau  ) ) {
		if ( is.null(tau_lower) | is.null(tau_upper))
		 stop('If *tau* is not specified, boundaries *tau_lower* and *tau_upper* must be specified.')
	}
	
	if ( is.null( ne0 )){
		coint <- ape::coalescent.intervals( tre )
		with( coint , {
			abs(interval.length) * ( lineages * (lineages-1) / 2) 
		}) -> .ne
		.ne[ .ne == 0 ] <- NA
		ne0 <- median( .ne, na.rm=T)
		ne <- rlnorm( res , log( ne0 ), .2 ) # add some jitter
	} else{
		if ( length(ne0)==1){
			ne0 <- rep( ne0, res )
		}
		ne = ne0
	}
	
	dh<- tredat$dh[1]
	
	
	# estimate tau 
	tauof <- function(tau){
		.mlskygrid_oos( tau, tredat, ne0, res =res, maxHeight = NeStartTimeBeforePresent, ncross = ncross, ncpu = ncpu)
	}
	if (is.null(tau )){
		cat('Precision parameter *tau* not provided. Computing now....\n')
		taustar <- optimize( tauof, lower = tau_lower, upper = tau_upper, maximum = TRUE , tol = tau_tol)
		tau = taustar$maximum 
		cat( paste( 'Precision parameter tau = ', tau , '\n') )
	}
	#/estimate tau 
	
	
	roughness_penalty <- function(logne){
		sum( dnorm( diff(diff( logne)), 0, sd = sqrt(dh/tau), log = TRUE) )
	}
	
	lterms <- function(logne)
	{
		ne <- exp( logne )
		sterms <- with(tredat, {
			-intervalLength * ltt_terms / ne [ ne_bin ] 
		})
		coterms <- with(tredat, {
			nco * ( log( ltt_terms ) - logne[ ne_bin ]  )
		})
		coterms[ is.na(coterms)] <- 0
#~ browser()
		coterms + sterms 
	}
	
	of <- function(logne ){
		sum( lterms( logne )) + roughness_penalty( logne )
	}
	
	cat( ' Estimating Ne(t)...\n')
	optim( par = log(ne), fn = of 
	  , method = 'BFGS'
	  , control = list( trace = ifelse(quiet, 0, 1), fnscale  = -1 )
	  , hessian = TRUE 
	) -> fit
	fi <- tryCatch( solve( -fit$hessian), error = function(e) {
		warning('Hessian could not be computed. Will not compute CIs.')
		NA
	})
	fsigma <- if (!any(is.na(fi))) {
		sqrt( diag( fi ) )
	} else{
		NA
	}
	
	#output 
	# note reverse axis 
	fsigma <- rev( fsigma )
	logne <- rev( fit$par )
	ne <- rev( exp( fit$par ))
	nelb <- exp( logne - fsigma*1.96 )
	neub <- exp( logne + fsigma*1.96 )
	ne_ci <- cbind( nelb, ne, neub )
	
	growthrate <-  c( diff ( ne  ) / (ne[-res] ) / dh , NA)
	loglik = of ( fit$par ) - roughness_penalty ( fit$par )
	rv <- list( 
		ne =  ne
	  , ne_ci = ne_ci  
	  , growthrate =  growthrate
	  , tau = tau
	  , time = sort( -tredat$h[ tredat$type == 'neswitch' ] )
	  , tredat = tredat
	  , tre = tre	
	  , sigma = fsigma 
	  , optim = fit 
	  , loglik = loglik
	  , rp = roughness_penalty( fit$par )
	)
	
	class(rv) <- 'mlskygrid'
	rv
}




##############
.neplot <- function( fit, ggplot=TRUE, logy = TRUE , ... )
{
	stopifnot(inherits(fit, "mlskygrid"))
	ne <- fit$ne_ci
	if ( 'ggplot2' %in% installed.packages()  & ggplot)
	{
		pldf <- data.frame( t = fit$time, nelb = ne[,1], nemed = ne[,2], neub = ne[,3] )
		pl <- ggplot2::ggplot( pldf, ggplot2::aes( x = t, y = nemed), ... ) + ggplot2::geom_line() + ggplot2::geom_ribbon( ggplot2::aes( ymin = nelb, ymax = neub), fill = 'blue', alpha = .2) + ggplot2::ylab('Effective population size') + ggplot2::xlab('Time before most recent sample')
		if (logy) pl <- pl + ggplot2::scale_y_log10()
		return(pl)
	} else{
		if (logy)
			plot( fit$time, ne[,2], ylim=range(ne[,1:3],na.rm=T),lwd =2, col = 'black', type = 'l', log='y',xlab='Time', ylab='Effective population size', ...)
		else
			plot( fit$time, ne[,2], ylim=range(ne[,1:3],na.rm=T),lwd =2, col = 'black', type = 'l',xlab='Time', ylab='Effective population size', ...)
		lines( fit$time, ne[,1] , lty=3)
		lines( fit$time, ne[,3] , lty=3)
		invisible(fit)
	}
}

.growthplot  <- function( fit , ggplot=TRUE, logy=FALSE, ...)
{
	stopifnot(inherits(fit, "mlskygrid"))
	if ( 'ggplot2' %in% installed.packages()  & ggplot)
	{
		pldf <- data.frame( t = fit$time, gr = fit$growth)
		pl <- ggplot2::ggplot( pldf, ggplot2::aes( x = t, y = gr), ... ) + ggplot2::geom_line() + ggplot2::ylab('Growth rate') + ggplot2::xlab('Time before most recent sample')
		if (logy) pl <- pl + ggplot2::scale_y_log10() 
		return(pl)
	} else{
		if (logy)
			plot( fit$time, fit$growth, lwd =2, col = 'black', type = 'l', log='y', xlab='Time', ylab='Growth rate',...)
		else
			plot( fit$time, fit$growth, lwd =2, col = 'black', type = 'l', xlab='Time', ylab='Growth rate', ...)
		
		invisible(fit)
	}
}


#' Plot mleskygrid 
#'
#' @param fit A fitted object
#' @param growth If TRUE will plot estimated growth rate instead of Ne(t) 
#' @param logy  If TRUE, the plot is returned with logarithmic y-axis
#' @param ggplot  If TRUE, returns a ggplot2 figure
#' @param ... Additional parameters are passed to ggplot or the base plotting function
#' @return Plotted object 
#' @export
plot.mlskygrid <- function(fit, growth=FALSE, ggplot=FALSE, logy=TRUE, ... ){
	if (growth) {
	  return(.growthplot(fit, ggplot, logy, ... ))
	} else{
		return(.neplot(fit, ggplot, logy, ... ))
	}
}

#' Print fitted mleskygrid 
#'
#' @param fit Fitted mleskygrid object 
#' @export 
print.mlskygrid <- function( fit ){
	stopifnot(inherits(fit, "mlskygrid"))
	d <- as.data.frame( fit$ne_ci )
	d <- cbind( fit$time, d )
	colnames( d ) <- c( 'Time before most recent sample', '2.5%', 'MLE', '97.5%' )
	cat(paste( 'mlskygrid fit
	Smoothing parameter tau =', fit$tau, '\n\n'))
	
	cat( 'Estimated Ne(t): \n')
	print ( d )
	invisible( fit )
}
