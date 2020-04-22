
samplePsi<-function(W,X,Psi,V,M,u1,loca,b,v){
  
  n=nrow(Psi)
  Psiprop=mvrnorm(1,Psi,u1*solve(t(X)%*%X))
  SSig=gSigma(b,v,loca)
  
  postPsi=-0.5*t(Psi-M)%*%solve(V)%*%(Psi-M)-0.5*t(W-X%*%Psi)%*%solve(SSig)%*%(W-X%*%Psi)
  
  postPsiprop=-0.5*t(Psiprop-M)%*%solve(V)%*%(Psiprop-M)-0.5*t(W-X%*%Psiprop)%*%solve(SSig)%*%(W-X%*%Psiprop)
  
  
  prob=min(exp((postPsiprop)-(postPsi)),1)
  
  
  u=runif(1,0,1)
  
  if(u<prob){
    
    Psiprox=Psiprop
    
    rejei=1
    
  }else{
    
    Psiprox=Psi		
    rejei=0
  }
  
  
  
  
  
  res=as.matrix(Psiprox)
  res=list(res,rejei)
  res
  
}



