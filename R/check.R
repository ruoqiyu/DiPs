check<-function(fdata,mdata,fz,mz){
  stopifnot(dim(fdata)[2]==dim(mdata)[2])
  stopifnot(colnames(fdata)==colnames(mdata))
  stopifnot(all(fz%in%c(0,1)))
  stopifnot(all(mz%in%c(0,1)))

  if (is.vector(fdata)) fdata<-matrix(fdata,ncol=1)
  if (is.vector(mdata)) mdata<-matrix(mdata,ncol=1)
  fd0<-fdata[which(fz==0),]
  fd1<-fdata[which(fz==1),]
  md0<-mdata[which(mz==0),]
  md1<-mdata[which(mz==1),]
  if (is.vector(fd0)) fd0<-matrix(fd0,length(fd0),1)
  if (is.vector(md0)) md0<-matrix(md0,length(md0),1)
  if (is.vector(fd1)) fd1<-matrix(fd1,length(fd1),1)
  if (is.vector(md1)) md1<-matrix(md1,length(md1),1)

  varf0<-apply(fd0,2,var,na.rm=TRUE)
  varf1<-apply(fd1,2,var,na.rm=TRUE)
  meanf0<-apply(fd0,2,mean,na.rm=TRUE)
  meanf1<-apply(fd1,2,mean,na.rm=TRUE)
  meanm0<-apply(md0,2,mean,na.rm=TRUE)
  meanm1<-apply(md1,2,mean,na.rm=TRUE)
  smdf<-(meanf1-meanf0)/sqrt((varf0+varf1)/2)
  smdm<-(meanm1-meanm0)/sqrt((varf0+varf1)/2)

  r<-cbind(meanf1,meanm0,meanf0,smdm,smdf)
  rownames(r)<-colnames(fdata)
  colnames(r)<-c('Treated Mean','Control Match Mean','Control All Mean',
                 'Control Match SMD','Control All SMD')
  r
}
