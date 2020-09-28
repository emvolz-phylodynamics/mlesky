#' Simulation of coalescent dated phylogeny
#' @param dates Sampling dates
#' @param alphaFun Population size function Ne(t)
#' @param alphaMin Minimum value of alphaFun
#' @return A simulated dated phylogeny
#' @export
simCoal = function(dates=1990:2010,alphaFun=function(x){return(10)},alphaMin=NA) {
  if (is.na(alphaMin)) alphaMin=optimize(alphaFun,c(-1e5,max(dates)))$objective
  s <- sort(dates,decreasing=TRUE,index.return = TRUE)
  tim <- s$x
  ind <- s$ix
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
  plot(xs,ys,type='l',xlab='',ylab='Population size', bty='l')
}