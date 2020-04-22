
sampleb <- function(W,v,b,loca,ab,bb,X,Psi,u1){
  
  bprop=rgamma(1,shape=b*u1, rate = u1)
  
  SSigprop=gSigma(bprop,v,loca)
  
  if((det(SSigprop)==0)|(bprop< 0.005)){
    return(list(b,0))
  }
  
  SSig=gSigma(b,v,loca)
  logp=-0.5*t(W-X%*%Psi)%*%solve(SSig)%*%(W-X%*%Psi)-0.5*log(det(SSig))+(ab-1)*log(b)-bb*b
  
  SSigprop=gSigma(bprop,v,loca)
  logpprop=-0.5*t(W-X%*%Psi)%*%solve(SSigprop)%*%(W-X%*%Psi)-0.5*log(det(SSigprop))+(ab-1)*log(bprop)-bb*bprop
  
  logprob=logpprop+log(dgamma(b,shape=bprop*u1,rate=u1)+1e-17)-(logp+log(dgamma(bprop,shape=b*u1,rate=u1)+1e-17))
  prob<-min(c(1,exp(logprob)))
  
  u=runif(1,0,1)
  
  if(u<prob){
    
    bprox=bprop
    
    rejei=1
    
    
  }else{
    
    bprox=b
    
    rejei=0;
    
  }
  
  
  
  res=list(bprox,rejei)
  res
}