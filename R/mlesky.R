# derive timeseries of coalescent and ltt along appropriate time axis 
.tre2df <- function( apephylo, tre, res, maxHeight = Inf, minLTT = 1, adapt_time_axis = F, sampleTimes = NULL ){
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

#' Heuristic selection of grid resolution. Selects res based on  
#' reduction in MSE of coalescent times relative to nearest grid point. 
suggest_res <- function(tree, th = .001 ) 
{
	#tree = ape::read.tree( system.file('mrsa.nwk', package = 'mlesky') )
	#tree = ape::read.tree( system.file('sn02ag2.0.nwk', package = 'mlesky') )
	ct = tail( ape::node.depth.edgelength( tree ) , ape::Nnode( tree ))
	#ct <- rnorm( 399, sd = 100 )
	#~ ct =  rnorm( 1000, mean = 1:100, sd = .25 )
	#~ ct = rep( 0, 399 )
	rct = range(ct) 
	k <- 1 
	x <- seq( rct[1], rct[2], length = k + 2 )
	l0 <-  sum( apply( sapply( ct, function(y) (y - x)^2 ), MAR = 2, FUN = min ) )
	ll <- l0 
	repeat{ 
		k <- k + 1 
		x <- seq( rct[1], rct[2], length = k + 2 )
		l <- sum( apply( sapply( ct, function(y) (y - x)^2 ), MAR = 2, FUN = min ) )
		print( c( k, l, x ))
		if ( (l < (th*l0)) | (ll<=l) ) break
		#if ( ((ll - l ) / l) < th ) break
		ll <- l
	}
	k <- k - 1 
	k
}


#' Optimize the skygrid time axis resolution using AIC criterion
#' 
#' @param tree A dated phylogeny in ape::phylo format
#' @param res A vector of time axis resolution parameters to test 
#' @param ncpu Integer number of cores to use with parallel processing 
#' @param ... Remaining parameters are passed to mlskygrid 
#~ optim_res_aic <- function(tree, res = c(1:5, seq(10, 100, by = 10)),  ncpu = 1, ... )
#~ { 
#~ 	res2aic <- function(r){
#~ 		ll1 <- mlskygrid( tree, res = r, ncpu =ncpu,  ...)$loglik
#~ 		 2 * r - 2 * ll1
#~ 	}
#~ 	aics <- unlist( parallel::mclapply( res,  res2aic, mc.cores = ncpu ) )
#~ 	res[ which.min( aics )]
#~ }

#' Optimize the skygrid time axis resolution using BIC criterion
#' 
#' @param tree A dated phylogeny in ape::phylo format
#' @param res A vector of time axis resolution parameters to test 
#' @param ncpu Integer number of cores to use with parallel processing 
#' @param ... Remaining parameters are passed to mlskygrid 
#~ optim_res_bic <- function(tree, res = c(1:5, seq(10, 100, by = 10)),  ncpu = 1, ... )
#~ { 
#~   res2bic <- function(r){
#~     ll1 <- mlskygrid( tree, res = r, ncpu =ncpu,  ...)$loglik
#~     r * log(tree$Nnode) - 2 * ll1
#~   }
#~   bics <- unlist( parallel::mclapply( res,  res2bic, mc.cores = ncpu ) )
#~   res[ which.min( bics )]
#~ }

roughness_penalty <- function(logne,dh,tau,b=NULL,model=1){
  y=0
  if (!is.null(b)) y=b
  if (model==1) {#skysigma model
    dh2 <- dh[ -c(1, length(dh)) ]
    rp_terms=dnorm(diff(diff( logne)), y, sd = sqrt(dh2/tau), log = TRUE)
  }
  if (model==2) {#skygrid model
    dh2 <- dh[ -1 ]
    rp_terms=dnorm(diff(logne), y, sd = sqrt(dh2/tau), log = TRUE)
  }
  if (model==3) {#skygrowth model
    dh2 <- dh[ -c(1, length(dh)) ]
    rhos=diff(exp(logne))/exp(logne[-length(logne)])
    rp_terms=dnorm( rhos, y, sd = sqrt(dh2/tau), log = TRUE)
  }
  sum(rp_terms)
}

#This function calculates the cross-validation score for a given tau
.mlskygrid_oos <- function( tau, tredat, ne0, res = 50, maxHeight = Inf, quiet = TRUE, control = NULL, ncross = 5, ncpu = 1,model=1, cvtype = 'interweaved' ){
	if ( ncross < 2 ) stop('*ncross* must be at least two')
	
	ne=ne0 #ne <- rlnorm( res , log( ne0 ), .2 ) # add some jitter
	
	if ( cvtype == 'interweaved' ){
		cvsets = lapply( 1:ncross, function(icross) seq( icross, nrow(tredat), by = ncross ) )
	} else if ( cvtype == 'segmented' ) {
		cvbounds <- cbind( 
			seq( 1/nrow(tredat), 1-1/ncross , length = ncross )
			, seq( 1/ncross, 1 , length = ncross )
		) * nrow(tredat)
		cvbounds[,1] <- ceiling( cvbounds[,1] )
		cvbounds[,2] <- floor( cvbounds[,2] )
		cvsets <- lapply( 1:ncross, function(icross){
			cvbounds[icross,1]:cvbounds[icross,2]
		})
	} else{
		stop('Incorrect cvtype')
	}
	
	dh <- sapply( 1:res, function(i) tredat$dh[ which(tredat$ne_bin==i)[1]] )
	
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
		sum( lterms( logne )[ cvsets[[icross]] ] )
	}
	
	of.cv.ws <- function(logne, icross){
		i = setdiff(1:nrow(tredat) , cvsets[[icross]]  )
		sum( lterms( logne )[ i ] ) + roughness_penalty( logne,dh,tau,model )
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
#' @param ne0 Vector of length *res* giving starting conditions of Ne(t) for optimization, or a single value which will be used as rep(ne0,res)
#' @param adapt_time_axis If TRUE will choose Ne(t) change points in periods with high frequency of phylogenetic branching
#' @param model Model to use, can be 1 (=skysigma model), 2 (=skygrid model) or 3 (=skygrowth model)
#' @param formula Formula for use of covariates, should not have left hand side
#' @param formula_order For use of covariates, lhs is 0'th 1st or 2nd deriv of ne wrt time 
#' @param data For use of covariates, data.frame must include 'time' 
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
  , adapt_time_axis = FALSE 
  , model=1
  , formula = NULL 
  , formula_order = c( 0, 1, 2) 
  , data = NULL 
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
		res = optim_res_aic( tre, tau = tau, tau_lower = tau_lower, tau_upper = tau_upper, tau_tol = tau_tol, ncross = ncross, ncpu = ncpu , quiet = quiet, NeStartTimeBeforePresent = NeStartTimeBeforePresent, ne0 = ne0, adapt_time_axis = adapt_time_axis, sampleTimes=sampleTimes,model=model )
	}
	if ( res < 1) 
	  stop('The minimum allowable *res* value is 1.')
	
	tredat <- .tre2df(apephylo = apephylo,  tre = tre, res = res , maxHeight= NeStartTimeBeforePresent, adapt_time_axis = adapt_time_axis , sampleTimes = sampleTimes )
	if ( is.null( tau  ) ) {
		if ( is.null(tau_lower) | is.null(tau_upper))
		 stop('If *tau* is not specified, boundaries *tau_lower* and *tau_upper* must be specified.')
	}
	if (res<3) tau=NULL
	
	if ( is.null( ne0 )){
		coint <- ape::coalescent.intervals( apephylo )
		with( coint , {
			abs(interval.length) * ( lineages * (lineages-1) / 2) 
		}) -> .ne
		.ne[ .ne == 0 ] <- NA
		ne0 <- median( .ne, na.rm=T)
		ne0<-rep(ne0,res);ne<-ne0#ne <- rlnorm( res , log( ne0 ), .2 ) # add some jitter
	} else{
		if ( length(ne0)==1){
			ne0 <- rep( ne0, res )
		}
		ne = ne0
	}
	
	
	# if covariates 
	ncovar = 0 
	covar.df = NA 
	if (!is.null( formula )){
		formula_order = formula_order[1]
		stopifnot( !is.null(sampleTimes)  ) 
		stopifnot( !is.null(data)  ) 
		
		# do an initial fit without covars to serve as initial condition
		fit0 <- mlskygrid(tre
		  , sampleTimes = sampleTimes
		  , res = res 
		  , tau = tau
		  , tau_lower = tau_lower 
		  , tau_upper = tau_upper
		  , tau_tol = tau_tol 
		  , ncross = ncross
		  , ncpu = ncpu
		  , quiet = quiet
		  , NeStartTimeBeforePresent = NeStartTimeBeforePresent 
		  , ne0 = ne0
		  , adapt_time_axis = adapt_time_axis 
		  , formula = NULL 
		  , data = NULL 
		)
		
		X0 <- as.data.frame( model.matrix(  formula , data ) )
		betanames <- colnames( X0 )[-1]
		ncovar = length( betanames )
		X0 <- cbind( time = data$time 
		 , X0 )
		
		covar.df <- data.frame( time =fit0$time[-c(1,length(fit0$time))], y = diff(diff(log(fit0$ne)))
		, logne = log(fit0$ne)[-c(1,length(fit0$time))]  )
		for ( bn in betanames ){
			itime <- setdiff( order( X0$time ), which(is.na( X0[[bn]] )) )
			
			covar.df[[bn]] <- approx( X0$time[itime] , X0[[bn]][itime], xout = covar.df$time,rule = 2)$y
			if ( formula_order == 0 ){
				covar.df[[bn]] <- c( NA, NA,  diff( diff( covar.df[[bn]] )  ) )
				covar.df[[bn]][1:2] <- covar.df[[bn]][3]
			}else if ( formula_order ==1  ){
				covar.df[[bn]] <- c( NA,  diff( covar.df[[bn]] )  )
				covar.df[[bn]][1] <- covar.df[[bn]][2]
			}else if ( formula_order == 2){
				# do nothing 
			} else{
				stop('formula_order > 2 not supported')
			}
		}
		covar.df <- covar.df[ order( covar.df$time) , ] 
		lmfit <- lm ( formula(
			   paste0( paste( 'y', paste( betanames , collapse = '+'), sep='~' ) , '-1' )
			 )
			 , data = covar.df 
			)
		
		beta0 <- coef(lmfit)
		tau = fit0$tau 
		ne0 = fit0$ne 
		beta2zxb <- function( beta ){
			if ( length( betanames ) > 1 ){
				zedCrossBeta <- as.vector( as.matrix(covar.df[, betanames]) %*% beta )
			} else{
				zedCrossBeta <- covar.df[, betanames] * beta 
			}
			rev( zedCrossBeta ) #  b/c reverse time axis of Ne in optimizer
		}
	}
	
	dh <- sapply( 1:res, function(i) tredat$dh[ which(tredat$ne_bin==i)[1]] )

	# estimate tau 
	tauof <- function(tau){
		.mlskygrid_oos( tau, tredat, ne0, res =res, maxHeight = NeStartTimeBeforePresent, ncross = ncross, ncpu = ncpu,quiet=quiet,model=model)
	}
	if (is.null(tau) && res>=3){
		cat('Precision parameter *tau* not provided. Computing now....\n')
		taustar <- optimize( tauof, lower = tau_lower, upper = tau_upper, maximum = TRUE , tol = tau_tol)
		tau = taustar$maximum 
		cat( paste( 'Precision parameter tau = ', tau , '\n') )
	}
	#/estimate tau 
	
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
	
	of <- function( theta ){ 
		if (ncovar == 0) b = NULL else b = tail( theta , ncovar ) 
		logne = head( theta, length(theta) - ncovar )
		lt = sum(lterms( logne)) 
		.beta2zxb = 0
		if ( ncovar > 0 ){
			.beta2zxb = beta2zxb(b)
		}
		rp = roughness_penalty(logne,dh,tau,.beta2zxb,model=model)
		lt + rp 
	}

	theta0 <- log( ne ) 
	if (ncovar > 0 ){
		theta0 <- c( theta0, beta0 )
	}
	
	parscale = rep(1, length( theta0 ))
	if ( ncovar > 0 ){
		parscale = c( rep(1, length(ne)), abs(median(log(ne)) / beta0  ) )
	}
	optim( par = theta0, fn = of 
	  , method = 'BFGS'
	  , control = list( trace = ifelse(quiet, 0, 1)
		, fnscale  = -1 
		, parscale = parscale
		)
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
	theta <- fit$par 
	logne = theta 
	beta = NULL 
	if ( ncovar > 0 ){
		logne <- head( theta , length( theta ) - ncovar )
		beta <- tail( theta, ncovar )
	}
	# reverse 
	fsigma_ne <- rev( head( fsigma, length( fsigma) - ncovar ) )
	logne <- rev( logne )
	ne <- exp( logne )
	nelb <- exp( logne - fsigma_ne*1.96 )
	neub <- exp( logne + fsigma_ne*1.96 )
	ne_ci <- cbind( nelb, ne, neub )
	
	loglik = of ( fit$par ) - roughness_penalty ( fit$par,dh,tau,model=model )
	
	mst = ifelse( is.null(sampleTimes), 0, max(sampleTimes) )
	
	h2 <- -c( sort( -tredat$h[ tredat$type == 'neswitch' ] ) , 0)
	time <- h2  - diff( c(max(tredat$h), h2) )/2
	time <- mst - time 
	growthrate  = c( diff( ne ) / ( ne[-length(ne)]* diff(time)  ) , NA) 
	rv <- list( 
		ne =  ne
	  , ne_ci = ne_ci  
	  , growthrate =  growthrate
	  , tau = tau
	  , time = time 
	  , tredat = tredat
	  , tre = tre	
	  , sigma = fsigma_ne 
	  , fsigma = fsigma 
	  , optim = fit 
	  , loglik = loglik
	  , lterms = lterms( theta )
	  , sampleTimes = sampleTimes
	  , beta = beta 
	  , covar.df  = covar.df
	)
	class(rv) <- 'mlskygrid'
	rv
}




##############
.neplot <- function( fit, ggplot=TRUE, logy = TRUE , ylim,... )
{
  nemed <- nelb <- neub <- NULL
	stopifnot(inherits(fit, "mlskygrid"))

	if (length(fit$ne)==1) {#This is in case res=1
	  fit$ne=rep(fit$ne,2)
	  fit$ne_ci=rbind(fit$ne_ci,fit$ne_ci)
	  fit$time=rep(fit$time,2)
	}
	
	if (is.null(fit$sampleTimes)) {
	  fit$time[1]=-sum(fit$tredat$intervalLength)#Force leftmost point to the root date
	  fit$time[length(fit$time)]=0#Force rightmost point to the most recent leaf date
    if (!is.null(fit$tre$root.time)) dateLastSample=fit$tre$root.time+max(dist.nodes(fit$tre)[Ntip(fit$tre)+1,]) else dateLastSample=0
	} else {
	  fit$time[1]=max(fit$sampleTimes)-sum(fit$tredat$intervalLength)#Force leftmost point to the root date
	  fit$time[length(fit$time)]=max(fit$sampleTimes)#Force rightmost point to the most recent leaf date
	  dateLastSample=0
	}
	
	ne <- fit$ne_ci
	if ( 'ggplot2' %in% installed.packages()  & ggplot)
	{
		pldf <- data.frame( t = dateLastSample+fit$time, nelb = ne[,1], nemed = ne[,2], neub = ne[,3] )
		pl <- ggplot2::ggplot( pldf, ggplot2::aes( x = t, y = nemed), ... ) + ggplot2::geom_line() + ggplot2::geom_ribbon( ggplot2::aes( ymin = nelb, ymax = neub), fill = 'blue', alpha = .2) + ggplot2::ylab('Effective population size') + ggplot2::xlab('Time before most recent sample')
		if (logy) pl <- pl + ggplot2::scale_y_log10()
		return(pl)
	} else{
	  if (!hasArg('ylim')) ylim=range(ne[,1:3],na.rm=T)
		if (logy)
			plot( dateLastSample+fit$time, ne[,2], ylim=ylim,lwd =2, col = 'black', type = 'l', log='y',xlab='Time', ylab='Effective population size', ...)
		else
			plot( dateLastSample+fit$time, ne[,2], ylim=ylim,lwd =2, col = 'black', type = 'l',xlab='Time', ylab='Effective population size', ...)
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
#' @param logy  If TRUE, the plot is returned with logarithmic y-axis (default is TRUE for Ne plot and FALSE for growth plot)
#' @param ggplot  If TRUE, returns a ggplot2 figure
#' @param ... Additional parameters are passed to ggplot or the base plotting function
#' @return Plotted object 
#' @export
plot.mlskygrid <- function(x, growth=FALSE, ggplot=FALSE, logy, ... ){
	if (growth) {
	  if (missing(logy)) logy=F
	  return(.growthplot(x, ggplot, logy, ... ))
	} else{
	  if (missing(logy)) logy=T
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
