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
		) )
	}
	eventdf <- data.frame( 
		var = 'A'
		, time = sort( mst - dates )
		, value = 1
		, method = 'add' 
	)
	y0 <- c( cumco = 0 , A = 0 )
	taxis <- seq( 0, mst - guessRootTime , length = res )
	o = suppressWarnings( deSolve::ode( dco, times = taxis, y = y0 , parms = NULL, method = 'lsoda', events = list( data = eventdf ) )) 
	o
}


#' Simulation of coalescent dated phylogeny using deterministic/'sideways' distribution of coalescent times
#' @param dates Sampling dates
#' @param alphaFun Population size function Ne(t)
#' @param res Time axis resolution
#' @param guessRootTime Optional (but recommended) guess for when root may be reached. 
#' @param A_tol Tolerance for reaching time of root >0
#' @return A simulated dated phylogeny
#' @export
ddSimCoal <- function(dates=1990:2010,alphaFun=function(x){return(10)}, guessRootTime = NA, res = 1e3, A_tol = 1e-2) {
	mst <- max(dates )
	if ( is.na( guessRootTime )){
		guessRootTime = min(dates) - 2 * diff(range(dates))#alphaFun( min( dates ))
	}
	A1 <- Inf 
	while( (A1 - 1) > A_tol  ){
		o <- suppressWarnings( .solve_A(dates=dates, alphaFun = alphaFun, guessRootTime = guessRootTime, res = res)  )
		A1 <- tail( o[, 'A'], 1)
		h1 <- tail( o[, 'time'], 1) 
		if ( A1 < 3 )
			h2 <- h1 + A1 * alphaFun( mst - h1 )
		else
			h2 <- h1 + 2 * diff( range(dates))
		guessRootTime <- mst - h2 
	}
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
	wresample <- which( treedf$A <= 1 & treedf$dA < 0  )
	while( length( wresample ) > 0 ){
		treedf$height[wresample] <- suppressWarnings( tcofun ( runif( length( wresample ) , 0, tail(o[,'cumco'],1) ) )  )
		treedf <- treedf[ order( treedf$height ) , ]
		treedf$A <- cumsum( treedf$dA )
		wresample <- which( treedf$A <= 0 & treedf$dA < 0  )
	}
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
  plot(xs,ys,type='l',xlab='',ylab='Population size', bty='l',ylim=c(0,1.05*max(ys)))
}
