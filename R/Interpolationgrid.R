



Interpolationgrid <- function(Geodata,loca,data,size=21,intensity=0.85){

  loca <- as.matrix(loca)

  alphaM <- Geodata[[1]]
  betaM <- Geodata[[2]]
  Mb <-Geodata[[3]]
  Mv <- Geodata[[4]]
  MPsi <- Geodata[[6]]
  MW <- Geodata[[5]]


  rxinf=min(loca[,1])
  rxsup=max(loca[,1])

  ryinf=min(loca[,2])
  rysup=max(loca[,2])

  rrea<-size

  seqtemp1=seq(rxinf, rxsup, length.out = rrea)
  seq1=seqtemp1

  seqtemp2=seq(ryinf,rysup,length.out=rrea)
  seq2=seqtemp2


  ort=length(seq1)

  DNOxtemp=array(NA,dim=c(ort,2))
  DNOx=NULL


  for(u in 1:ort){
    for(hj in 1:ort){
      DNOxtemp[hj,1]=seq1[u]
      DNOxtemp[hj,2]=seq2[hj]
    }

    DNOx=rbind(DNOx,DNOxtemp)

  }

  res=rbind(DNOx,loca)




  tt=nrow(res)
  X1=as.matrix(rep(1,tt))
  X2=as.matrix(res[,1])
  X3=as.matrix(res[,2])
  Xr=cbind(X1,X2,X3)



  iter<-length(alphaM)
  rra=seq(1,iter,100)



  PsiM=as.matrix(MPsi[rra,])
  bM=as.matrix(Mb[rra])
  WM=as.matrix(MW[rra,])
  vM=as.matrix(Mv[rra])
  betaM=as.matrix(betaM[rra])
  alphaM=as.matrix(alphaM[rra])


  #Sz é loca

  rxinf=min(loca[,1])
  rxsup=max(loca[,1])

  ryinf=min(loca[,2])
  rysup=max(loca[,2])

  rrea<-size
  seqtemp1=seq(rxinf, rxsup, length.out = rrea)
  seq1=seqtemp1

  seqtemp2=seq(ryinf, rysup, length.out = rrea)
  seq2=seqtemp2

  ort=length(seq1)



  DNOtemp=array(NA,dim=c(ort,2))
  DNO=NULL
  for(u in 1:ort){
    for(hj in 1:ort){
      DNOtemp[hj,1]=seq1[u]
      DNOtemp[hj,2]=seq2[hj]
    }

    DNO=rbind(DNO,DNOtemp)

  }


  DO=loca
  jj=nrow(DNO)

  iter1=nrow(PsiM)

  MWNO=array(NA,dim=c(iter1,jj))

  for(h in 1:iter1){
    sigma=gSigma(bM[h,1],vM[h,1],res)
    Mpro=Xr%*%as.matrix(PsiM[h,])

    A1=as.matrix(Mpro[(jj+1):tt,])
    A2=as.matrix(Mpro[1:jj,])

    SSA1=sigma[(jj+1):tt,(jj+1):tt]
    SSA2=sigma[1:jj,1:jj]
    SSA12=sigma[1:jj,(jj+1):tt]

    A2est=A2+SSA12%*%solve(SSA1)%*%(as.matrix(WM[h,])-A1)
    SSA2est=SSA2-SSA12%*%solve(SSA1)%*%t(SSA12)
    WNO=as.matrix(mvrnorm(1,A2est,SSA2est))
    MWNO[h,]=t(WNO)
  }



  tinter=nrow(data)

  #deixar número de dias do dadinho
  #por default, o máximo dos dados transformados


  ITT=as.matrix(rep(1,ncol(MWNO)))


  tres1=exp(MWNO)

  tres2=-betaM%*%t(ITT)



  tres3=alphaM%*%t(ITT)


  tres4=tinter^tres3


  minter=tres1*(1-exp(tres2*tres4))

  minterest=apply(minter,2,mean,na.rm=T)

  Mminterest=matrix(minterest,byrow=T,ncol=ort)




  dimen<-max((rxsup-rxinf),(rysup-ryinf))

  tam<-dimen/(rrea-1)


  df <- expand.grid(x = seq(rxinf+tam,rxinf+rrea*tam,tam), y = seq(ryinf+tam,ryinf+rrea*tam,tam))

  df$z <- as.vector(Mminterest)



  interpolationmap<-c(left=rxinf,bottom=ryinf,right=(rxinf+tam*rrea),top=(ryinf+tam*rrea))





  zoomdim <- 2
  if(dimen < 50 ){
    zoomdim <- 3
  }
  if(dimen < 20 ){
    zoomdim <- 6
  }
  if(dimen < 5 ){
    zoomdim <- 8
  }
  if(dimen < 0.5){
    zoom <- 11
  }

  if(dimen < 0.02 ){
    zoomdim <- 14
  }




  p <-  get_stamenmap(interpolationmap, zoom = zoomdim, maptype = "toner-lite")


  station <- matrix(c(1:nrow(loca)),nrow=29)
  station<-as.character(station)

  locad<-data.frame(loca,station)

  local<- cbind(rep(seq1,size),rep(seq2,each=size),cbind(as.vector(Mminterest)))

  colnames(local) <- c("seq1","seq2","numero")
  local<- as.data.frame(local)

  ggmap(p)+stat_contour(data=local,aes(x=seq1,y=seq2,z=numero))+
    geom_tile(data = local, aes(x = seq1, y = seq2, z = numero, fill = numero), alpha = intensity)+
    scale_fill_continuous(name = "Número de eventos anômalos",
                          low = "brown", high = "blue")+
    geom_point(data=locad,aes(x=x,y=y))+geom_text(data = locad, aes(x = x, y = y, label = station),
                                                  size = 3, vjust = 0, hjust = -0.5)


}

