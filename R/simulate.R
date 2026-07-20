#' Simulation of coalescent dated phylogeny
#' @param dates Sampling dates
#' @param alphaFun Population size function Ne(t)
#' @param alphaMin Minimum value of alphaFun
#' @return A simulated dated phylogeny
#' @export
simCoal = function(dates=1990:2010,alphaFun=function(x){return(10)},alphaMin=NA) {
	if (is.na(alphaMin)) alphaMin=optimize(alphaFun,c(-1e5,max(dates)))$objective
	ind <- order( dates, decreasing=TRUE )
	tim <- dates [ ind  ] #s <- sort(dates,decreasing=TRUE,index.return = TRUE)
	n <- length(tim)
	nodes <- cbind(-Inf,ind[1],-Inf)#Start with one node at time -Inf and with the first isolate connected to it
	i <- 2
	while (i <= n) {#Graft branches one by one
		curt <- tim[i]#Current time:start with date of isolate and go back in time until coalescence happens
		accept=F
		while (accept==F) {
		  r <- -log(runif(1)) * alphaMin
		  fi <- which( nodes[ ,1] < curt )[1]
		  if (fi<=nrow(nodes)) for (j in (fi:nrow(nodes)))  {
			if (r > (curt-nodes[j,1]) * (i-j))  {
			  r <- r-(curt-nodes[j,1]) * (i-j)
			  curt <- nodes[j,1]
			} else {
			  curt <- curt-r/(i-j)#Proposed coalescent time
			  break
			}
		  }
		  #Filtering the Poisson process to obtain non-homogeneous Poisson
		  if (runif(1)<alphaMin/alphaFun(curt)) accept=T
		}
		#Create new node
		a <- nodes[ ,2:3]
		a[a >= j + n] <- a[a >= j + n] + 1
		nodes[ ,2:3] <- a#Renumbering according to table insertion in next line
		nodes2=c(curt,ind[i],0)
		if (1<=j-1) nodes2=rbind(nodes[1:(j-1), ],nodes2)
		if (j<=nrow(nodes)) nodes2=rbind(nodes2,nodes[j:nrow(nodes),])
		nodes=unname(nodes2)
		#Now choose on which branch to coalesce among the branches alive at date curt
		no <- j
		side <- 2
		w <- 1 + floor(runif(1) * (nrow(nodes)-j))
		while (w > 0)  {
		  no <- no + side-1
		  side <- 3-side
		  if (nodes[no,side + 1] <= n ||(nodes[no,side + 1] > n && nodes[nodes[no,side + 1]-n,1] > curt))  {
			w <- w-1
		  }
		}
		nodes[j,3] <- nodes[no,side + 1]
		nodes[no,side + 1] <- n + j
		i <- i + 1
	}
	v=nrow(nodes)-1
	if (nrow(nodes)>2) v=c(v,1:(nrow(nodes)-2))
	nodes=nodes[v,,drop=F]
	m=nodes[,2:3]
	m[which(m>n)]=m[which(m>n)]+1
	nodes[,2:3]=m
	nodes <- rbind(matrix(0, nrow = n, ncol = 3),nodes)
	nodes[1:n,1] <- dates

	#Convert into phylo object from package ape
	t=list()
	t$Nnode=n-1
	t$tip.label=as.character(1:n)
	if (!is.null( names(dates )))
		t$tip.label <- names(dates)
	t$edge=matrix(NA,2*n-2,2)
	t$edge.length=rep(NA,n*2-2)
	t$root.time=nodes[n+1,1]
	c=1
	for (i in (n+1):nrow(nodes)) for (j in 2:3) {
		t$edge[c,1]=i
		t$edge[c,2]=nodes[i,j]
		t$edge.length[c]=nodes[nodes[i,j],1]-nodes[i,1]
		c=c+1
	}
	class(t)='phylo'
	return(t)
}

.solve_A <- function(dates=1990:2010
 , alphaFun=function(x){return(10)}
 , guessRootTime = 1950
 , res = 1e3)
{
	mst <- max(dates)
	stopifnot( guessRootTime <= min(dates) )
	dco <- function( t, y, parms ){
		cumco <- y['cumco'] 
		A <- max( 1 , y['A']  )
		x = mst - t # t is on reverse time axis 
		ne <- alphaFun( x )
		stopifnot( ne > 0 )
		list( c(
			cumco = unname( A * (A-1) / ne / 2 )
			, A = unname(  -A * (A-1) / ne / 2 )
			# rescaled coalescent time: dtau = dt/(2 Ne) makes dA/dtau = -A(A-1),
			# so (A-1) decays like exp(-tau) whatever the shape of alphaFun
			, tau = unname( 1 / ( 2 * ne ) )
		) )
	}
	eventdf <- data.frame(
		var = 'A'
		, time = sort( mst - dates )
		, value = 1
		, method = 'add'
	)
	y0 <- c( cumco = 0 , A = 0 , tau = 0 )
	taxis <- seq( 0, mst - guessRootTime , length = res )
	o = suppressWarnings( deSolve::ode( dco, times = taxis, y = y0 , parms = NULL, method = 'lsoda', events = list( data = eventdf ) ))
	o
}


# Extend the time axis until A reaches 1 (within A_tol), working in rescaled coalescent
# time tau = int dt/(2 Ne). Because dA/dtau = -A(A-1), the extra tau needed to reach
# A_tol is known in closed form, which both sizes the step and identifies the case where
# alphaFun admits no MRCA: if tau converges to a finite limit as t -> -Inf then A stalls
# above 1 and the tree has no root. That happens whenever Ne grows fast enough into the
# past, which the previous heuristic could not detect and looped on forever.
.horizon_tau <- function(dates, alphaFun, guessRootTime, res, A_tol, maxit)
{
	mst <- max( dates )
	stalls <- 0
	prevTau <- NA
	prevNeed <- NA
	for ( it in 1:maxit ){
		o <- suppressWarnings( .solve_A( dates = dates, alphaFun = alphaFun, guessRootTime = guessRootTime, res = res ) )
		A1 <- tail( o[, 'A'], 1 )
		h1 <- tail( o[, 'time'], 1 )
		tau1 <- tail( o[, 'tau'], 1 )
		if ( (A1 - 1) <= A_tol )
			return( o )
		need <- log( (A1 - 1) / A_tol ) # additional rescaled time required
		if ( !is.na( prevTau ) ){
			# far less tau bought than was asked for => the integral is converging
			if ( (tau1 - prevTau) < 0.25 * prevNeed ) stalls <- stalls + 1 else stalls <- 0
			if ( stalls >= 2 )
				stop( sprintf( paste0( 'ddSimCoal: alphaFun implies Ne grows without bound into the past, so the total rescaled coalescent time int dt/(2 Ne) converges (tau -> %.3g) and the number of extant lineages stalls at %.3g > 1; this population size function has no most recent common ancestor. Provide an alphaFun that is bounded going back in time.' ), tau1, A1 ) )
		}
		prevTau <- tau1
		prevNeed <- need
		ne <- alphaFun( mst - h1 )
		step <- if ( is.finite(ne) && ne > 0 ) 2 * ne * need else Inf # dh = 2 Ne dtau
		step <- min( max( step, 0.1 * h1 ), 9 * h1 ) # bound growth to [1.1x, 10x] per iteration
		guessRootTime <- mst - ( h1 + step )
	}
	stop( 'ddSimCoal: root time not reached in ', maxit, ' iterations; try increasing maxit or A_tol' )
}


# Previous root time search, retained for comparison. Kept as-is apart from the maxit
# backstop, without which it does not terminate for some alphaFun (see .horizon_tau).
.horizon_legacy <- function(dates, alphaFun, guessRootTime, res, A_tol, maxit)
{
	mst <- max( dates )
	for ( it in 1:maxit ){
		o <- suppressWarnings( .solve_A( dates = dates, alphaFun = alphaFun, guessRootTime = guessRootTime, res = res ) )
		A1 <- tail( o[, 'A'], 1 )
		h1 <- tail( o[, 'time'], 1 )
		if ( (A1 - 1) <= A_tol )
			return( o )
		if ( A1 < 3 )
			h2 <- h1 + A1 * alphaFun( mst - h1 )
		else
			h2 <- h1 + 2 * diff( range(dates) )
		guessRootTime <- mst - h2
	}
	stop( 'ddSimCoal: root time not reached in ', maxit, ' iterations; try increasing maxit or A_tol' )
}


#' Simulation of coalescent dated phylogeny using deterministic/'sideways' distribution of coalescent times
#' @param dates Sampling dates
#' @param alphaFun Population size function Ne(t)
#' @param res Time axis resolution
#' @param guessRootTime Optional (but recommended) guess for when root may be reached. 
#' @param A_tol Tolerance for reaching time of root >0
#' @param horizon Strategy for extending the time axis until the root is reached. 'tau' (the default) works in rescaled coalescent time and reports an error for population size functions which admit no most recent common ancestor; 'legacy' is the previous heuristic.
#' @param maxit Maximum number of iterations of the root time search
#' @return A simulated dated phylogeny
#' @export
ddSimCoal <- function(dates=1990:2010,alphaFun=function(x){return(10)}, guessRootTime = NA, res = 1e3, A_tol = 1e-2, horizon = c('tau','legacy'), maxit = 100) {
	horizon <- match.arg( horizon )
	mst <- max(dates )
	if ( is.na( guessRootTime )){
		guessRootTime = min(dates) - 2 * diff(range(dates))#alphaFun( min( dates ))
	}
	o <- if ( horizon == 'tau' )
		.horizon_tau( dates, alphaFun, guessRootTime, res, A_tol, maxit )
	else
		.horizon_legacy( dates, alphaFun, guessRootTime, res, A_tol, maxit )
	n <- length(dates)
	treedf1 <- data.frame( 
			node = 1:n
			, type = 'sample'
			, height = mst - dates
			, dA = 1
			, stringsAsFactors = FALSE 
			
	)
	tcofun <- suppressWarnings( approxfun( 
	 x = o[,'cumco'] 
	 , y = o[, 'time'] ) )
	tco <- suppressWarnings( tcofun ( runif(n-1 , 0, tail(o[,'cumco'],1) ) )  )
	treedf2 <- data.frame(
		node = (n+1):(n+n-1)
		, type = 'coalescent'
		, height = sort( tco , decreasing=TRUE )
		, dA = -1 
		, stringsAsFactors = FALSE
		
	)
	treedf <- rbind( treedf1, treedf2 )
	treedf <- treedf[ order( treedf$height ) , ]
	treedf$A <- cumsum( treedf$dA )
	# resampling to fix any negative branch lengths
	# note A==1 at the root is valid, so only A<1 indicates a coalescent event with too few extant lineages
	wresample <- which( treedf$A < 1 & treedf$dA < 0  )
	while( length( wresample ) > 0 ){
		treedf$height[wresample] <- suppressWarnings( tcofun ( runif( length( wresample ) , 0, tail(o[,'cumco'],1) ) )  )
		treedf <- treedf[ order( treedf$height ) , ]
		treedf$A <- cumsum( treedf$dA )
		wresample <- which( treedf$A < 1 & treedf$dA < 0  )
	}
	# resampling reorders the coalescent events, so node numbers must be re-derived from the
	# final height order: ape::write.tree assumes the root is node n+1 and will otherwise
	# misalign the internal branch lengths
	wco <- which( treedf$type == 'coalescent' )
	treedf$node[ wco ] <- rev( (n+1):(2*n-1) )
	nonbl <- Inf

	edge <- matrix( NA, nrow = 2*n-2 , ncol = 2 )
	edge.length <- rep( NA, 2*n-2)
	tip.label = names(dates); if (is.null( names(dates))) tip.label <- as.character(1:n)
	Nnode = n -1 
	 
	k <- 1
	extant <- c() 
	for ( i in 1:nrow( treedf )){
		if ( treedf$type[i] == 'sample' ){
			extant <- c( extant , treedf$node[i] )
		} else if ( treedf$type[i]  == 'coalescent' ){
			uv = sample( extant, size = 2, replace=FALSE)
			u = uv[1]; v = uv[2] 
			
			edge[k,1] <- treedf$node[i] 
			edge[k,2] <- u
			edge.length[k] <- treedf$height[i] - treedf$height[treedf$node==u]
			k <- k + 1
			edge[k,1] <- treedf$node[i] 
			edge[k,2] <- v
			edge.length[k] <- treedf$height[i] - treedf$height[treedf$node==v]
			k <- k +1 
			
			extant <- c( extant, treedf$node[i] )
			extant <- setdiff( extant, uv ) 
		}
	}
	.trd = list( edge = edge, edge.length = edge.length, Nnode = Nnode
		, tip.label = tip.label
		, treedf = treedf)
	tr = structure( .trd
	  , class = 'phylo' )
	ape::read.tree( text = ape::write.tree( tr ) )
}


# Tabulate the rescaling tau(h) = int_0^h dt/(2 Ne) on a grid of heights h below mst.
# alphaFun is evaluated scalar-wise since it need not be vectorised.
# The grid is broken at hsamp (the oldest sample) and given full resolution on either
# side. Extending H then leaves the grid over the sampling period, and hence tau there,
# bit-identical, so sample positions computed before an extension stay consistent with
# coalescent times inverted after one. A single uniform grid does not have that property:
# heterochronous sampling with a deep tree loses resolution over the sampling period as H
# grows, and the two maps can disagree enough to order a coalescence below its own child.
.tabulate_tau <- function(mst, alphaFun, H, res, hsamp = 0)
{
	h <- if ( hsamp > 0 && H > hsamp )
		c( seq( 0, hsamp, length.out = res ), seq( hsamp, H, length.out = res )[-1] )
	else
		seq( 0, H, length.out = res )
	ne <- vapply( mst - h, function(x) as.numeric( alphaFun(x) ), numeric(1) )
	if ( any( !is.finite(ne) ) || any( ne <= 0 ) )
		return( NULL )
	r <- 1 / ( 2 * ne )
	list( h = h
		, tau = c( 0, cumsum( (r[-1] + r[-length(r)]) / 2 * diff(h) ) ) # cumulative trapezoid
	)
}


#' Simulation of coalescent dated phylogeny by time rescaling
#'
#' Samples exactly from the coalescent with time varying effective population size.
#' Under the rescaling tau = int dt/(2 Ne(t)) 
#' 
#' @param dates Sampling dates
#' @param alphaFun Population size function Ne(t)
#' @param res Resolution of the grid used to tabulate the rescaling
#' @param maxit Maximum number of iterations when extending the time axis to reach the root
#' @return A simulated dated phylogeny
#' @export
tauSimCoal <- function(dates=1990:2010, alphaFun=function(x){return(10)}, res = 1e3, maxit = 100)
{
	n <- length( dates )
	stopifnot( n >= 2 )
	mst <- max( dates )
	hs <- mst - dates                   # sample heights
	hsamp <- max( hs )
	H <- hsamp
	if ( H <= 0 ) H <- 2 * as.numeric( alphaFun( mst ) )  # isochronous: use the coalescent timescale
	g <- .tabulate_tau( mst, alphaFun, H, res, hsamp )
	if ( is.null(g) ) stop( 'tauSimCoal: alphaFun must be finite and positive over the sampling period' )

	# sample positions in rescaled time; these are fixed and do not depend on how far
	# back tau is subsequently tabulated, so the coalescent can be simulated up front
	tau_s <- stats::approx( g$h, g$tau, xout = hs, rule = 2 )$y
	ord <- order( tau_s )
	# record the event sequence in rescaled time: this ordering is exact, whereas the
	# heights it maps to are interpolated, so the tree is assembled in this order
	nev <- 2L*n - 1L
	ev_type <- integer( nev ); ev_tau <- numeric( nev ); ev_tip <- integer( nev )
	nxt <- 1L; A <- 0L; tcur <- 0; nco <- 0L; m <- 0L
	while ( nco < n - 1L ){
		ts <- if ( nxt <= n ) tau_s[ ord[nxt] ] else Inf
		if ( A >= 2L ){
			w <- stats::rexp( 1, A * (A - 1) )
			if ( tcur + w < ts ){
				tcur <- tcur + w
				m <- m + 1L; ev_type[m] <- 2L; ev_tau[m] <- tcur
				nco <- nco + 1L; A <- A - 1L
				next
			}
		}
		if ( !is.finite(ts) ) stop( 'tauSimCoal: exhausted samples with fewer than two extant lineages' )
		tcur <- ts
		m <- m + 1L; ev_type[m] <- 1L; ev_tau[m] <- ts; ev_tip[m] <- ord[nxt]
		A <- A + 1L; nxt <- nxt + 1L
	}
	tau_root <- tcur

	# extend the height axis until it spans tau_root so that tau can be inverted;
	# step sized by the local Ne, as dh = 2 Ne dtau
	stalls <- 0
	for ( it in 1:maxit ){
		if ( max( g$tau ) >= tau_root ) break
		deficit <- tau_root - max( g$tau )
		ne <- as.numeric( alphaFun( mst - H ) )
		step <- if ( is.finite(ne) && ne > 0 ) 2 * ne * deficit else Inf
		H <- H + min( max( step, 0.1 * H ), 9 * H )
		g2 <- .tabulate_tau( mst, alphaFun, H, res, hsamp )
		# Ne overflowing while extending, and tau gaining far less than asked for, are both
		# signatures of Ne diverging into the past: int dt/(2 Ne) then converges short of
		# the root and no MRCA exists
		if ( is.null(g2) ) stalls <- 2
		else {
			if ( max(g2$tau) - max(g$tau) < 0.01 * deficit ) stalls <- stalls + 1 else stalls <- 0
			g <- g2
		}
		if ( stalls >= 2 )
			stop( sprintf( 'tauSimCoal: alphaFun implies Ne grows without bound into the past (Ne is no longer finite and positive by height %.4g), so the total rescaled coalescent time int dt/(2 Ne) converges (tau -> %.3g) short of the root (tau %.3g); this population size function has no most recent common ancestor.', H, max(g$tau), tau_root ) )
		if ( it == maxit )
			stop( 'tauSimCoal: root not reached in ', maxit, ' iterations; try increasing maxit or res' )
	}
	# heights: exact for samples, inverted from rescaled time for coalescent events
	hev <- numeric( nev )
	isco <- ev_type == 2L
	hev[ !isco ] <- hs[ ev_tip[ !isco ] ]
	hev[ isco ] <- stats::approx( g$tau, g$h, xout = ev_tau[ isco ], rule = 2 )$y

	# tree assembly, in rescaled time order; heights are indexed by node id so no scanning
	nodeh <- numeric( 2*n - 1L )
	nodeh[ 1:n ] <- hs
	ci <- 2L*n - 1L                                    # coalescent ids assigned downwards
	k <- 1L; extant <- integer(0)
	edge <- matrix( NA_integer_, 2*n - 2L, 2L )
	edge.length <- numeric( 2*n - 2L )
	for ( i in 1:nev ){
		if ( ev_type[i] == 1L ){
			extant <- c( extant, ev_tip[i] )
		} else {
			u <- sample.int( length(extant), 2L )
			uv <- extant[u]
			# the rescaled time ordering is exact, the interpolated height is not, so
			# never let interpolation error push a parent below one of its children
			hi <- max( hev[i], nodeh[ uv[1] ], nodeh[ uv[2] ] )
			nodeh[ci] <- hi
			edge[k,] <- c( ci, uv[1] ); edge.length[k] <- hi - nodeh[ uv[1] ]; k <- k + 1L
			edge[k,] <- c( ci, uv[2] ); edge.length[k] <- hi - nodeh[ uv[2] ]; k <- k + 1L
			extant <- c( extant[-u], ci )
			ci <- ci - 1L                              # so the root ends up as node n+1
		}
	}
	tip.label = names(dates); if (is.null( names(dates))) tip.label <- as.character(1:n)
	tr <- structure( list( edge = edge, edge.length = edge.length, Nnode = n - 1L
		, tip.label = tip.label )
	  , class = 'phylo' )
	ape::read.tree( text = ape::write.tree( tr ) )
}



#' Plot a dated phylogeny and alpha function
#' @param tree Dated phylogeny
#' @param alphaFun Population size function Ne(t)
#' @return Figure
#' @export
plotBoth = function(tree,alphaFun) {
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  par(mfrow=c(2,1),mar=c(4,4,1,4))
  plot(tree,show.tip.label = F)
  axisPhylo(1,backward = F)
  if (!is.null(tree$root.time)) from=tree$root.time else from=-max(dist.nodes(tree)[Ntip(tree)+1,])
  to=from+max(dist.nodes(tree)[Ntip(tree)+1,])
  xs=seq(from=from,to=to,length.out=100)
  ys=xs
  for (i in 1:length(ys)) ys[i]=alphaFun(ys[i])
  plot(xs,ys,type='l',xlab='',ylab='Effective population size', bty='l',ylim=c(0,1.05*max(ys)))
}
