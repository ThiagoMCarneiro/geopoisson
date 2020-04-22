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







##################################################################333333############################################################333




GeoPoisson <- function(data,limuser,c1,d1,iter,burn_in,X,loca){
  
  #A matriz de W para as colunas
  data <- matrixtransf(data,limuser)
  
  # Chute inicial de W
  
  W <- matrix( c(rep(1,ncol(data))),ncol=1,nrow=ncol(data))
  
  loca <- data.matrix(loca,rownames.force = NA)
  
  X <- data.matrix(X,rownames.force = NA)
  
  
  #HiperParam
  
  
  
  ##?? 
  V <- diag(1000,ncol(X))
  
  M <- c(rep(0,ncol(X)))
  
  ##??
  
  
  ab <- 0.1
  bb <- 0.1
  
  #Chute inicial 
  
  alpha <-1
  beta <- 0.1
  b <- 5
  v <- 5
  Psi <- c(5,5,5)
  
  ###
  
  masoqeisso <- c()
  
  ###
  
  #Beta conhecido proposto (para teste apenas, no código final, receberá um chute inicial)
  
  
  n <- ncol(data)
  
  #Limpeza
  
  a <- c()
  d <- c()
  
  l <- c()
  m <- c()
  
  b_i <- c()
  b_ii <- c()
  
  v_i<-c()
  v_ii<-c()
  
  psy <- c()
  psycho <- Psi
  psycho1 <- matrix()
  bochecha<-c()
  
  dablio <- matrix()
  dab <- matrix()
  dabliz <- t(W)
  
  aceitaalpha <- NULL
  aceitabeta <- NULL
  aceitabe <- NULL
  aceitave <- NULL
  aceitaW <- NULL
  aceitaPsi <- NULL
  
  nj <- colSums(!is.na(data))
  
  #Arrays talvez sejam substitu?dos
  
  N=array(NA,dim=c(1,n))
  
  Tt=array(NA,dim=c(1,n))
  
  ff <- c(5,5,5,5,0.005,5)
  
  #Cpp ? amor
  
  for(y in 1:n){
    Tt[1,y]=data[nj[y],y]
  }
  
  Tt=t(Tt)
  
  
  for(j in 1:iter){
    
    theta=exp(W)
    
    
    #Fazer em Cpp. S?o 100.000 vezes fazendo de 1:22(ou n?mero de colunas);
    
    for(e in 1:n){
      #ultimo tempo de cada coluna
      
      t1<-data[nj[e],e]
      #até 22 - recebendo o vetor N
      
      N[1,e]<-rpois(1,theta[e,1]*(1-pgamma(beta*t1^alpha,1,1)))
    }
    
    g<-unlist(samplealpha(alpha,beta,N,Tt,nj,c1,d1,data,ff[1]),use.names=FALSE)
    
    #a[j] <- g[1]
    d[j]<- g[2]
    alpha <- g[1]
    
    g<-unlist(samplebeta(alpha,beta,N,Tt,nj,c1,d1,data,ff[2]),use.names=FALSE)
    
    #l[j] <- g[1]
    m[j]<- g[2]
    beta <- g[1]
    
    g <- unlist(sampleb(W,v,b,loca,ab,bb,X,Psi,ff[3]),use.names=FALSE)
    
    #b_i[j] <- g[1]
    b_ii[j] <-g[2]
    b <- g[1]
    
    g <- unlist(samplev(W,v,b,loca,ab,bb,X,Psi,ff[4]),use.names=FALSE)
    
    #v_i[j] <- g[1]
    v_ii[j] <- g[2]
    v <- g[1]
    
    g <- sampleW(W,loca,X,Psi,b,v,nj,N,ff[5])
    
    W <- as.matrix(g[[1]])
    #dabliz <- rbind(dabliz,t(W))
    dab[j] <- g[[2]]
    
    psy <- unlist(samplePsi(W,X,Psi,V,M,ff[6],loca,b,v))
    
    psycho1 <- psy[1:3]
    #psycho <- rbind(psycho,psycho1)
    bochecha[j] <- psy[4]
    Psi <- psy[1:3]
    
    #masoqeisso <- c(masoqeisso,ff[5])
    
    
    if((j%%50==0)&&(j<=burn_in)){
      
      aceitaalpha[j/50] <- (sum(d[(j-49):j])/50)
      aceitabeta[j/50] <-(sum(m[(j-49):j])/50)
      aceitabe[j/50] <-(sum(b_ii[(j-49):j])/50)
      aceitave[j/50] <-(sum(v_ii[(j-49):j])/50)
      aceitaW[j/50] <- (sum(dab[(j-49):j])/50)
      aceitaPsi[j/50] <-(sum(bochecha[(j-49):j])/50)
      
      
      if(sum(d[(j-49):j])/50 > 0.44){
        ff[1] <-  ff[1]/exp(min(0.01,((j/50)^-0.5)))
      } 
      
      if(sum(d[(j-49):j])/50 < 0.44){
        ff[1] <-  ff[1]*exp(min(0.01,((j/50)^-0.5)))
      }
      
      if(sum(m[(j-49):j])/50 > 0.44){
        ff[2] <-  ff[2]/exp(min(0.01,((j/50)^-0.5)))
      } 
      
      if(sum(m[(j-49):j])/50 < 0.44){
        ff[2] <-  ff[2]*exp(min(0.01,((j/50)^(-0.5))))
      }
      if(sum(b_ii[(j-49):j])/50 > 0.44){
        ff[3] <-  ff[3]/exp(min(0.01,((j/50)^-0.5)))
      } 
      
      if(sum(b_ii[(j-49):j])/50 < 0.44){
        ff[3] <-  ff[3]*exp(min(0.01,((j/50)^(-0.5))))
      }
      if(sum(v_ii[(j-49):j])/50 > 0.44){
        ff[4] <-  ff[4]/exp(min(0.01,((j/50)^-0.5)))
      } 
      
      if(sum(v_ii[(j-49):j])/50 < 0.44){
        ff[4] <-  ff[4]*exp(min(0.01,((j/50)^(-0.5))))
      }
      
      if(sum(dab[(j-49):j])/50 < 0.24){
        ff[5] <-  ff[5]/exp(min(0.01,((j/50)^-0.5)))
      } 
      
      if(sum(dab[(j-49):j])/50 > 0.24){
        ff[5] <-  ff[5]*exp(min(0.01,((j/50)^-0.5)))
      }
      
      if(sum(bochecha[(j-49):j])/50 < 0.24){
        ff[6] <-  ff[6]/exp(min(0.01,((j/50)^-0.5)))
      } 
      
      if(sum(bochecha[(j-49):j])/50 > 0.24){
        ff[6] <-  ff[6]*exp(min(0.01,((j/50)^(-0.5))))
      }
      
    }
    
    if(j>burn_in){
      a[j-burn_in] <- alpha
      l[j-burn_in] <- beta
      b_i[j-burn_in] <- b
      v_i[j-burn_in] <- v
      dabliz <- rbind(dabliz,t(W))
      psycho <- rbind(psycho,psycho1)
    }
    
    
  }
  #return(list(a,l,b_i,v_i,dabliz,psycho,aceitaalpha,aceitabeta,aceitabe,aceitave,aceitaW,aceitaPsi,masoqeisso))
  return(list(a,l,b_i,v_i,dabliz,psycho))
}

####################################################################################################


samplealpha <- function(alpha,beta,N,Tt,nj,c1,d1,x,ff){
  
  
  n=ncol(x)
  
  #FF é feito por nós; Devemos alcançar uma taxa de aceitação dee 40%;
  
  
  alphaprop<-rgamma(1,alpha*ff, rate=ff)
  
  
  ## Temp e temp1 parecem desnecessários
  
  #temp=array(NA,dim=c(n,1))
  #temp1=array(NA,dim=c(n,1))
  
  #Tt é auxiliar;
  #c1, d1 são hiperparâmetros;
  
  # dados <- na.ommit*(drop(data))
  #Transformar x em data
  
  palpha=(sum(nj)+c1-1)*log(alpha)+alpha*sum(log(x),na.rm = T)-beta*( sum(x^alpha,na.rm = T)+sum(t(N)*(Tt^alpha)) ) -d1*alpha
  palphaprop=(sum(nj)+c1-1)*log(alphaprop)+alphaprop*sum(log(x),na.rm = T)-beta*( sum(x^alphaprop,na.rm = T)+sum(t(N)*(Tt^alphaprop)) )-d1*alphaprop
  
  #Tirar o 0.0??01
  
  logprob<-palphaprop+log(dgamma(alpha,alphaprop*ff, rate=ff)+0.0000001)-(palpha + log(dgamma(alphaprop,alpha*ff, rate=ff)+0.0000001))
  
  probac<-min(c(1,exp(logprob)))
  
  u<-runif(1)
  
  if(u<probac){
    res<-alphaprop
    rejei=1
    
    
  }
  else{
    res<-alpha
    rejei=0
  }
  
  res=list(res,rejei)
  res
  
}


###################################################################################################



samplebeta <- function(alpha,beta,N,Tt,nj,c,d,x,ff)
{
  
  n=ncol(x)
  
  betaprop<-rgamma(1,shape=beta*ff, rate=ff)
  
  
  
  temp=array(NA,dim=c(n,1))
  temp1=array(NA,dim=c(n,1))
  
  
  pbeta<-(sum(nj)+c-1)*log(beta)-beta*( sum(x^alpha,na.rm = T)+d+sum(t(N)*(Tt^alpha)) )
  pbetaprop<-(sum(nj)+c-1)*log(betaprop)-betaprop*( sum(x^alpha,na.rm = T)+d+sum(t(N)*(Tt^alpha)) )
  
  
  logprob<-pbetaprop+log(dgamma(beta,shape=betaprop*ff, rate=ff)+0.0000001)-(pbeta + log(dgamma(betaprop,shape=beta*ff, rate=ff)+0.0000001))
  
  probac<-min(c(1,exp(logprob)))
  
  u<-runif(1)
  
  if(u<probac){
    res<-betaprop
    rejei=1
    
    
  }
  else{
    res<-beta
    rejei=0
  }
  
  res=list(res,rejei)
  
}

####################################################################################################


##################################################################

#def ? loca

###################################################################################################
#usu?rio colocaa ab e bb

####################################################################################################

###################################################################################################

####################################################################################################


matrixtransf <- function(matrixuser,limuser){
  
  matrixuser <- data.matrix(matrixuser,rownames.force = NA)
  
  maior <- apply(matrixuser,2,function(x) length(x[x>limuser]))
  maior <- max(maior)
  
  
  ajustado <- c()
  
  
  for(i in 1:ncol(matrixuser)){
    linhaprov <- c()
    
    linha <- matrixuser[,i]
    
    ## pensar em outro apply
    for(j in 1:nrow(matrixuser)){
      if(linha[j]>limuser){
        
        linhaprov <- c(linhaprov,j)
        
      }
      
    }
    linhaprov
    
    if(length(linhaprov)<maior){
      
      linhaprov <- c(linhaprov,rep(0,(maior-length(linhaprov))))
      
      ajustado<- rbind(ajustado,linhaprov)
      
    }else{
      
      ajustado <- rbind(ajustado,linhaprov)
    }
    
  }
  
  ajustado <- t(ajustado)
  ajustado[ajustado==0] <- NA
  
  return(ajustado)
}

