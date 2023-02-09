match<-function(z,dist,dat,p=rep(1,length(z)),exact=NULL,fine=rep(1,length(z)),ncontrol=1,penalty=ifelse(is.matrix(dist),round(max(dist)*1000),round(max(dist$d)*1000)),s.cost=100,subX=NULL){
  #Check input
  stopifnot(is.data.frame(dat))
  stopifnot(is.vector(z))
  if (!is.numeric(fine)){
    if (!is.factor(fine)) fine=as.factor(fine)
    levels(fine)<-1:nlevels(fine)
    fine<-as.numeric(fine)
  }

  stopifnot(all((z==1)|(z==0)))
  stopifnot((ncontrol==round(ncontrol))&(ncontrol>=1))
  n<-length(z)
  ntreat<-sum(z)
  ncontr<-sum(1-z)
  stopifnot(ncontr>=(ncontrol*ntreat))
  stopifnot(length(z)==length(fine))

  stopifnot(length(z)==(dim(dat)[1]))

  if (is.matrix(dist)){
    distance<-t(dist)
    dim(distance)<-c(1,ntreat*ncontr)
    distance<-as.vector(distance)
    start<-rep(1:ntreat,each=ncontr)
    end<-rep((ntreat+1):n,ntreat)
    d0<-distance
    distance<-distance[which(d0<Inf)]
    start<-start[which(d0<Inf)]
    end<-end[which(d0<Inf)]
    dist=list(d=distance,start=start,end=end)
  }

  if (!is.null(subX)){
    if (is.factor(subX)){
      levels(subX)<-1:nlevels(subX)
      subX<-as.integer(subX)
    }
    stopifnot(is.vector(subX))
  }

  #sort input
  if (is.null(exact)){
    o<-order(1-p)
  }else{
    o<-order(exact,1-p)
    exact<-exact[o]
  }

  z<-z[o]
  p<-p[o]
  fine<-fine[o]
  dat<-dat[o,]

  #Must have treated first
  if(!(min(z[1:(n-1)]-z[2:n])>=0)){
    o<-order(1-z)
    z<-z[o]
    dat<-dat[o,]
    fine<-fine[o]
    if (!is.null(subX)) subX<-subX[o]
  }

  #do match
  net<-net(z,dist,ncontrol,fine,penalty,s.cost,subX)
  if (any(net$cost==Inf)) net$cost[net$cost==Inf]<-2*max(net$cost[net$cost!=Inf])

  callrelax <- function (net, solver = 'rlemon'){
    startn <- net$startn
    endn <- net$endn
    ucap <- net$ucap
    b <- net$b
    cost <- net$cost
    stopifnot(length(startn) == length(endn))
    stopifnot(length(startn) == length(ucap))
    stopifnot(length(startn) == length(cost))
    stopifnot(min(c(startn, endn)) >= 1)
    stopifnot(max(c(startn, endn)) <= length(b))
    stopifnot(all(startn != endn))

    nnodes <- length(b)
    if(solver == 'rrelaxiv'){
      if(requireNamespace('rrelaxiv', quietly = TRUE)) {
        rout <- rrelaxiv::RELAX_IV(startnodes = as.integer(startn),
                                   endnodes = as.integer(endn),
                                   arccosts = as.integer(cost),
                                   arccapacity = as.integer(ucap),
                                   supply = as.integer(b))
        return.obj <- list(crash = 0,
                           feasible = !all(rout == 0),
                           x = rout)
        return(return.obj)
      } else {
        solver = 'rlemon'
        warning('Package rrelaxiv not available, using rlemon instead.')
      }
    }
    if(solver == 'rlemon'){
      lout <- rlemon::MinCostFlow(arcSources = as.integer(startn),
                                  arcTargets = as.integer(endn),
                                  arcCapacities = as.integer(ucap),
                                  arcCosts = as.integer(cost),
                                  nodeSupplies = as.integer(b),
                                  numNodes = max(c(startn, endn)),
                                  algorithm = 'CycleCancelling')
      return.obj <- list(crash = 0,
                         feasible = !all(lout[[1]] == 0),
                         x = lout[[1]])
      return(return.obj)
    }else{
      stop(
        'Argument to solver not recognized: please use one of rlemon and rrelaxiv')
    }
  }

  output<-callrelax(net)

  if (output$feasible!=1){
    warning("Match is infeasible.  Change dist or ncontrol to obtain a feasible match.")
    m<-list(feasible=output$feasible,d=NULL)
  }else{
    x<-output$x[1:net$tcarcs]
    treated<-net$startn[1:net$tcarcs]
    control<-net$endn[1:net$tcarcs]
    treated=treated[which(x==1)]
    control=control[which(x==1)]
    match.df=data.frame('treat'=treated,'control'=control)
    match.df$treat<-as.factor(as.character(match.df$treat))
    matches<-as.matrix(plyr::daply(match.df, plyr::.(match.df$treat),
                                   function(treat.edges) treat.edges$control,.drop_o=FALSE))
    id1<-(1:n)[z==1]
    id0<-(1:n)[z==0]
    matchid<-matrix(c(id1[as.numeric(row.names(matches))], id0[as.vector((matches-sum(z)))]),ncol=ncontrol+1)
    matchid<-as.vector(t(matchid))
    dat1<-dat[matchid,]
    zm<-z[matchid]
    mset<-rep(1:nrow(matches),each=ncontrol+1)
    dat1$mset<-mset
    #dat1<-cbind(mset,dat1)
    m<-list(feasible=output$feasible,data=dat1,x=x)
  }

  if(m[[1]]==0) {
    warning("The match you requested is infeasible, reconsider caliper or ncontrol or exact for distance.")
  }
  m
}
