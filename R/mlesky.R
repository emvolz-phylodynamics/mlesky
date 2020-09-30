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
.tre2df <- function( apephylo, tre, res, maxHeight = Inf, minLTT = 1, adapt_time_axis = TRUE, sampleTimes = NULL ){
	n <- ape::Ntip( apephylo )
	if (is.null( minLTT)) 
	  minLTT <- floor( n / 5 ) 
	
	
	if ( inherits( tre, c('multiPhylo','list') ) ){
		phys <- lapply( tre, function(tr) {class(tr) <- 'phylo'; tr } )
		stopifnot( !is.null( sampleTimes ))
		rts <- sapply( phys, function( phy ) {
			mdepth <- max( node.depth.edgelength( phy )  ) 
			mst <- max( sampleTimes [ phy$tip.label  ] )
			mst - mdepth  
		})
		mst = max( sampleTimes )
		sts <- sampleTimes - min( rts )
		rhs = mst - rts 
		rh = max( rhs ) 
		maxHeight <- min( rh, maxHeight )
		shs = rh - sts 
		inhs_list <- lapply( 1:length(phys), function(k){
			phy <- phys[[k]] 
			rh = rhs[k] 
			ndel <- ape::node.depth.edgelength( phy )
			sort( rh - ndel[ (Ntip(phy)+1):(Ntip(phy) + phy$Nnode) ] )
		})
		inhs <- sort( do.call( c, inhs_list ) )
	}else{
		stopifnot( inherits (tre, c('phylo','treedater') ))
		D <- ape::node.depth.edgelength( apephylo )
		rh <- max( D[1:n] )
		rhs = rh # for compatibility with multitree version 
		sts <- D[1:n]
		maxHeight <- min( rh, maxHeight )
		
		shs <- rh - sts 
		inhs <- sort( rh - D[ (n+1):(n + apephylo$Nnode) ] )
	}
	
	u_shs <- unique( shs ) 
	u_inhs <- unique( inhs )
	nnode <- sum( inhs  <= maxHeight)
	
	ne_haxis <- seq( maxHeight/res ,maxHeight*(1-1/res), le = res-1 )
	if (adapt_time_axis){
		ne_haxis <- approx(  seq(0,1,length.out=nnode), inhs[inhs <= maxHeight], xout = seq(1/res, 1-1/res, length.out = res-1 ) )$y
	}
	dh_ne <- diff( c(0, ne_haxis, maxHeight ) )
	
	#< h , event, ltt(descending), intervallength, nco, likterm, ne_bin >
	tredat <- data.frame( h= c( u_shs, u_inhs, ne_haxis) 
	 , type = c( rep('sample', length( u_shs))
	           , rep('node', length(u_inhs))
	           , rep('neswitch', length(ne_haxis))
	          )
	)
	tredat <- tredat[ tredat$h <= maxHeight , ]
	
	tredat$ne_bin <- sapply( tredat$h, function(x) sum( ne_haxis  < x ) + 1)
	
	ltt.h <- function(h) max(1, sum( shs < h ) - sum( inhs < h ) - sum( rhs < h ) )
	tredat$ltt <- sapply( tredat$h, ltt.h )
	
	tredat$nco <- 0
	tredat$nco[ tredat$type=='node' ] <- sapply( tredat$h[ tredat$type=='node'  ], function(x) sum( x == inhs ))
	
	tredat <- tredat[ order( tredat$h ), ]
	tredat$intervalLength <- c( 0, diff( tredat$h ))
	tredat$ltt_terms <- tredat$ltt * (tredat$ltt-1) / 2
	
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
	
	tredat$dh <- dh_ne[ tredat$ne_bin ] 
	tredat	
}


#' Optimize the skygrid time axis resolution using AIC criterion
#' 
#' @param tree A dated phylogeny in ape::phylo format
#' @param res A vector of time axis resolution parameters to test 
#' @param ncpu Integer number of cores to use with parallel processing 
#' @param ... Remaining parameters are passed to mlskygrid 
#' @export 
optim_res_aic <- function(tree, res = c(3, seq(10, 100, by = 10)),  ncpu = 1, ... )
{ 
	res2aic <- function(r){
		ll1 <- mlskygrid( tree, res = r, ncpu =ncpu,  ...)$loglik
		 2 * r - 2 * ll1
	}
	aics <- unlist( parallel::mclapply( res,  res2aic, mc.cores = ncpu ) )
	res[ which.min( aics )]
}

#' Optimize the skygrid time axis resolution using BIC criterion
#' 
#' @param tree A dated phylogeny in ape::phylo format
#' @param res A vector of time axis resolution parameters to test 
#' @param ncpu Integer number of cores to use with parallel processing 
#' @param ... Remaining parameters are passed to mlskygrid 
#' @export 
optim_res_bic <- function(tree, res = c(3, seq(10, 100, by = 10)),  ncpu = 1, ... )
{ 
  res2bic <- function(r){
    ll1 <- mlskygrid( tree, res = r, ncpu =ncpu,  ...)$loglik
    r * log(tree$Nnode) - 2 * ll1
  }
  bics <- unlist( parallel::mclapply( res,  res2bic, mc.cores = ncpu ) )
  res[ which.min( bics )]
}


.mlskygrid_oos <- function( tau, tredat, ne0, res = 50, maxHeight = Inf, quiet = TRUE, control = NULL, ncross = 5, ncpu = 1){
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

.bind_tres <- function(tres, sts){
	phys <- lapply( tres, function(tr) {class(tr) <- 'phylo'; tr } )
	rts <- sapply( phys, function( phy ) {
		mdepth <- max( node.depth.edgelength( phy )  ) 
		mst <- max( sts [ phy$tip.label  ] )
		mst - mdepth  
	})
	minrt <- min( rts ) 
	maxrt <- max( rts ) 
	rels <- rts - minrt + 1e-3 
	phys <- lapply( 1:length(phys), function(k) {
		phy <- phys[[k]]
		phy$root.edge <- rels[k]; phy 
	})
	
	.phy <- phys[[1]]
	for ( k in 2:length( phys )){
		#.phy <- bind.tree( .phy , phys[[k]] )
		.phy <- .phy + phys[[k]] 
	}
	multi2di( .phy )
}


#' Maximum likelihood non-parametric estimation of effective population size through time
#'
#' 
#'
#' @param tre A dated phylogeny in ape::phylo or treedater format (see documentation for ape). This can also be a multiPhylo or list of trees, in which case each is treated as a clade sampled from within the same population. In this case the sampleTimes vector should be supplied so that clades can be aligned in time. 
#' @param sampleTimes An optional named vector of sample times for each taxon. Names should correspond to tip labels in trees. This is required if providing a list of trees. 
#' @param res Length of time axis over which to estimate Ne(t) (integer). If NULL, will search for a good value (see *optim_res_aic*)
#' @param tau Precision parameter. Larger values generate smoother trajectories of Ne(t). If NULL, will optimize using cross-validation.
#' @param tau_lower Lower bound for precision parameter if estimating
#' @param tau_upper Upper bound for precision parameter if estimating
#' @param tau_tol Optimization tolerance when optimizing tau by cross-validation
#' @param ncross Number of folds in cross-validation
#' @param ncpu If doing cross-validation, each fold will be handled in parallel if ncpu > 1 (see parallel package)
#' @param quiet Provide verbose output from optimizer? 
#' @param NeStartTimeBeforePresent If <Inf, will only estimate Ne(t) between the most recent sample and this time before the most recent sample
#' @param adapt_time_axis If TRUE will choose Ne(t) change points in periods with high frequency of phylogenetic branching
#' @param ne0 Vector of length *res* giving starting conditions of Ne(t) for optimization
#' @return A fitted model including effective size through time
#' @export
# @examples
# library(mlesky)
# tree <- read.tree( system.file( package='mlskygrid', 'mrsa.nwk' , mustWork=TRUE) ) 
# print( (fit <- mlskygrid( tree, tau = 10, NeStartTimeBeforePresent = 15) ))
# plot( fit , logy = FALSE)
mlskygrid <- function(tre
  , sampleTimes = NULL
  , res = 25 
  , tau = 1
  , tau_lower = NULL 
  , tau_upper = NULL
  , tau_tol = 1e-3 
  , ncross = 5
  , ncpu = 1
  , quiet = TRUE
  , NeStartTimeBeforePresent = Inf 
  , ne0 = NULL
  , adapt_time_axis = TRUE 
){
	apephylo <- tre
	if ( inherits( tre, c('multiPhylo','list') ) ){
		apephylo <- .bind_tres( tre , sampleTimes )
		class(apephylo) <- 'phylo'
		stopifnot( !is.null( sampleTimes ))
	}else{
		stopifnot( inherits (tre, c('phylo','treedater') ))
	}
	class( apephylo ) <- 'phylo'
	
	if (!is.null( sampleTimes )){
		stopifnot( is.numeric( sampleTimes ))
		sampleTimes <- sampleTimes[apephylo$tip]
	}
	
	if ( is.null( res )){
		res = optim_res_aic( tre, tau = tau, tau_lower = tau_lower, tau_upper = tau_upper, tau_tol = tau_tol, ncross = ncross, ncpu = ncpu , quiet = quiet, NeStartTimeBeforePresent = NeStartTimeBeforePresent, ne0 = ne0, adapt_time_axis = adapt_time_axis, sampleTimes=sampleTimes )
	}
	if ( res < 3) 
	  stop('The minimum allowable *res* value is 3.')
	
	tredat <- .tre2df(apephylo = apephylo,  tre = tre, res = res , maxHeight= NeStartTimeBeforePresent, adapt_time_axis = adapt_time_axis , sampleTimes = sampleTimes )
	if ( is.null( tau  ) ) {
		if ( is.null(tau_lower) | is.null(tau_upper))
		 stop('If *tau* is not specified, boundaries *tau_lower* and *tau_upper* must be specified.')
	}
	
	if ( is.null( ne0 )){
		coint <- ape::coalescent.intervals( apephylo )
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
	
	
	dh <- sapply( 1:res, function(i) tredat$dh[ which(tredat$ne_bin==i)[1]] )
	dh2 <- dh[ -c(1, length(dh)) ]
	
	# estimate tau 
	tauof <- function(tau){
		.mlskygrid_oos( tau, tredat, ne0, res =res, maxHeight = NeStartTimeBeforePresent, ncross = ncross, ncpu = ncpu,quiet=quiet)
	}
	if (is.null(tau )){
		cat('Precision parameter *tau* not provided. Computing now....\n')
		taustar <- optimize( tauof, lower = tau_lower, upper = tau_upper, maximum = TRUE , tol = tau_tol)
		tau = taustar$maximum 
		cat( paste( 'Precision parameter tau = ', tau , '\n') )
	}
	#/estimate tau 
	
	rp_terms <- function(logne){
		dnorm( diff(diff( logne)), 0, sd = sqrt(dh2/tau), log = TRUE)
	}
	
	roughness_penalty <- function(logne){
		sum( na.omit(rp_terms(logne)) )
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
	
	of <- function(logne ){
		sum( lterms( logne )) + roughness_penalty( logne )
	}
	
	#cat( ' Estimating Ne(t)...\n')
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
	
	growthrate <-  c( diff ( ne  ) / (ne[-res] ) / dh[-res] , NA)
	loglik = of ( fit$par ) - roughness_penalty ( fit$par )
	
	mst = ifelse( is.null(sampleTimes), 0, max(sampleTimes) )
	
	h2 <- -c( sort( -tredat$h[ tredat$type == 'neswitch' ] ) , 0)
	time <- h2  - diff( c(max(tredat$h), h2) )/2
	time <- mst - time 
	rv <- list( 
		ne =  ne
	  , ne_ci = ne_ci  
	  , growthrate =  growthrate
	  , tau = tau
	  , time = time 
	  , tredat = tredat
	  , tre = tre	
	  , sigma = fsigma 
	  , optim = fit 
	  , loglik = loglik
	  , rp = roughness_penalty( fit$par )
	  , rpterms = rp_terms( fit$par ) 
	  , lterms = lterms( fit$par )
	  , sampleTimes = sampleTimes
	)
	class(rv) <- 'mlskygrid'
	rv
}




##############
.neplot <- function( fit, ggplot=TRUE, logy = TRUE , ... )
{
  nemed <- nelb <- neub <- NULL
	stopifnot(inherits(fit, "mlskygrid"))
  if (!is.null(fit$tre$root.time)) dateLastSample=fit$tre$root.time+max(dist.nodes(fit$tre)[Ntip(fit$tre)+1,]) else dateLastSample=0
	ne <- fit$ne_ci
	if ( 'ggplot2' %in% installed.packages()  & ggplot)
	{
		pldf <- data.frame( t = dateLastSample+fit$time, nelb = ne[,1], nemed = ne[,2], neub = ne[,3] )
		pl <- ggplot2::ggplot( pldf, ggplot2::aes( x = t, y = nemed), ... ) + ggplot2::geom_line() + ggplot2::geom_ribbon( ggplot2::aes( ymin = nelb, ymax = neub), fill = 'blue', alpha = .2) + ggplot2::ylab('Effective population size') + ggplot2::xlab('Time before most recent sample')
		if (logy) pl <- pl + ggplot2::scale_y_log10()
		return(pl)
	} else{
		if (logy)
			plot( dateLastSample+fit$time, ne[,2], ylim=range(ne[,1:3],na.rm=T),lwd =2, col = 'black', type = 'l', log='y',xlab='Time', ylab='Effective population size', ...)
		else
			plot( dateLastSample+fit$time, ne[,2], ylim=range(ne[,1:3],na.rm=T),lwd =2, col = 'black', type = 'l',xlab='Time', ylab='Effective population size', ...)
		lines( dateLastSample+fit$time, ne[,1] , lty=3)
		lines( dateLastSample+fit$time, ne[,3] , lty=3)
		invisible(fit)
	}
}

.growthplot  <- function( fit , ggplot=TRUE, logy=FALSE, ...)
{
  gr<-NULL
	stopifnot(inherits(fit, "mlskygrid"))
	if (!is.null(fit$tre$root.time)) dateLastSample=fit$tre$root.time+max(dist.nodes(fit$tre)[Ntip(fit$tre)+1,]) else dateLastSample=0
	if ( 'ggplot2' %in% installed.packages()  & ggplot)
	{
		pldf <- data.frame( t = dateLastSample+fit$time, gr = fit$growth)
		pl <- ggplot2::ggplot( pldf, ggplot2::aes( x = t, y = gr), ... ) + ggplot2::geom_line() + ggplot2::ylab('Growth rate') + ggplot2::xlab('Time before most recent sample')
		if (logy) pl <- pl + ggplot2::scale_y_log10() 
		return(pl)
	} else{
		if (logy)
			plot( dateLastSample+fit$time, fit$growth, lwd =2, col = 'black', type = 'l', log='y', xlab='Time', ylab='Growth rate',...)
		else
			plot( dateLastSample+fit$time, fit$growth, lwd =2, col = 'black', type = 'l', xlab='Time', ylab='Growth rate', ...)
		
		invisible(fit)
	}
}


#' Plot mlskygrid 
#'
#' @param x A fitted object
#' @param growth If TRUE will plot estimated growth rate instead of Ne(t) 
#' @param logy  If TRUE, the plot is returned with logarithmic y-axis
#' @param ggplot  If TRUE, returns a ggplot2 figure
#' @param ... Additional parameters are passed to ggplot or the base plotting function
#' @return Plotted object 
#' @export
plot.mlskygrid <- function(x, growth=FALSE, ggplot=FALSE, logy=TRUE, ... ){
	if (growth) {
	  return(.growthplot(x, ggplot, logy, ... ))
	} else{
		return(.neplot(x, ggplot, logy, ... ))
	}
}

#' Print fitted mlskygrid 
#'
#' @param x Fitted mlskygrid object 
#' @param ... Additional parameters are passed on
#' @export 
print.mlskygrid <- function( x,... ){
	stopifnot(inherits(x, "mlskygrid"))
	d <- as.data.frame( x$ne_ci )
	d <- cbind( x$time, d )
	if ( is.null ( x$sampleTimes )){
		colnames( d ) <- c( 'Time before most recent sample', '2.5%', 'MLE', '97.5%' )
	} else {
		colnames( d ) <- c( 'Time', '2.5%', 'MLE', '97.5%' )
	}
	cat(paste( 'mlskygrid fit
	Smoothing parameter tau =', x$tau, '\n\n'))
	
	cat( 'Estimated Ne(t): \n')
	print ( d,... )
	invisible( x )
}
