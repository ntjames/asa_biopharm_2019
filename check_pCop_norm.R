# check pCop_norm

# Stan code
if (0){
  real pCop_norm(real u1, real u2, real theta){
    if (u1==0 || u2==0){ // grounded
      return 0;
    }
    else if (u1==1){ // uniform margin
      return u2;
    }
    else if (u2==1){ // uniform margin
      return u1;
    } else {
      return binormal_cdf(inv_Phi(u1), inv_Phi(u2), theta);
    }
  }
}

# R code
library(mvtnorm)

my_pCop_norm<-function(u1,u2,theta){
  if (u1==0|u2==0){ return(0) }
  else if (u1==1){ return(u2) }
  else if (u2==1){ return(u1) }
  else {  
    out<-pmvnorm(lower=c(-Inf,-Inf), upper=c(qnorm(u1),qnorm(u2)), corr=diag(2)+(1-diag(2))*theta)
    return(out[1])
    }
}


p1<-0.00001
p2<-0.7
thet<-0.4
y1<-1
y2<-1

# pseudo-obs from bernoulli dist
u1<-pbinom(y1,1,p1)
u1_prime <-pbinom(y1-1,1,p1)
u2<-pbinom(y2,1,p2)
u2_prime <-pbinom(y2-1,1,p2)

# likelihood for joint dist. using differences
# real cP0 = pCop_norm(U1, U2, theta) -
# pCop_norm(U1_prime, U2, theta) -
# pCop_norm(U1, U2_prime, theta) +
# pCop_norm(U1_prime, U2_prime, theta);

a<-my_pCop_norm(u1,u2,thet) 
b<-my_pCop_norm(u1_prime,u2,thet)
c<-my_pCop_norm(u1,u2_prime,thet)
d<-my_pCop_norm(u1_prime,u2_prime,thet)

a-b-c+d
