
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
  
  alpha <- c(quantile(a,0.025),quantile(a,0.5),quantile(a,0.975),mean(a))
  beta <- c(quantile(l,0.025),quantile(l,0.5),quantile(l,0.975),mean(l))
  b <-   c(quantile(b_i,0.025),quantile(b_i,0.5),quantile(b_i,0.975),mean(b_i))
  v <-  c(quantile(v_i,0.025),quantile(v_i,0.5),quantile(v_i,0.975),mean(v_i))
  dada <- rbind(alpha,beta,b,v)
  colnames(dada) <-  c("2.5%","50%","97.5%","Mean")
  
  print(dada)
  
  return(list(a,l,b_i,v_i,dabliz,psycho))
}
