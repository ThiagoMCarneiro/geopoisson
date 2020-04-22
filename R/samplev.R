
samplev <- function(W,v,b,loca,ab,bb,X,Psi,u1){
  
  vprop=rgamma(1,shape=v*u1, rate = u1)
  
  SSigprop=gSigma(b,vprop,loca)
  
  if((det(SSigprop)==0)|(vprop< 0.005)){
    return(list(v,0))
  }
  
  SSig=gSigma(b,v,loca)
  logp=-0.5*t(W-X%*%Psi)%*%solve(SSig)%*%(W-X%*%Psi)-0.5*log(det(SSig))+(ab-1)*log(v)-bb*v
  
  SSigprop=gSigma(b,vprop,loca)
  logpprop=-0.5*t(W-X%*%Psi)%*%solve(SSigprop)%*%(W-X%*%Psi)-0.5*log(det(SSigprop))+(ab-1)*log(vprop)-bb*vprop
  
  logprob=logpprop+log(dgamma(v,shape=vprop*u1,rate=u1)+1e-17)-(logp+log(dgamma(vprop,shape=v*u1,rate=u1)+1e-17))
  prob<-min(c(1,exp(logprob)))
  
  u=runif(1,0,1)
  
  if(u<prob){
    
    vprox=vprop
    
    rejei=1
    
    
  }else{
    
    vprox=v
    
    rejei=0;
    
  }
  
  
  
  res=list(vprox,rejei)
  res
}
