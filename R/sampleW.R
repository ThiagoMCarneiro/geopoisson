
sampleW <-function(W,loca,X,Psi,b,v,nj,N,u1){
  
  n=nrow(W)
  Wprop=mvrnorm(1,W,u1*diag(1,n))
  SSig=gSigma(b,v,loca)
  
  postW=sum(as.matrix(nj)*W)-sum(exp(W))+sum(t(N)*W)-0.5*t(W-X%*%Psi)%*%solve(SSig)%*%(W-X%*%Psi)
  postWprop=sum(as.matrix(nj)*Wprop)-sum(exp(Wprop))+sum(t(N)*Wprop)-0.5*t(Wprop-X%*%Psi)%*%solve(SSig)%*%(Wprop-X%*%Psi)
  
  prob=min(exp((postWprop)-(postW)),1)
  
  
  u=runif(1,0,1)
  
  if(u<prob){
    
    Wprox=Wprop
    
    rejei=1
    
    
  }else{
    
    Wprox=W		
    rejei=0
  }
  
  res=as.matrix(Wprox)
  res=list(Wprox,rejei)
  res
  
}
