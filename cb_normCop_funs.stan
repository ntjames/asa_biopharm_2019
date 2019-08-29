// functions for bivariate Normal copula with normal and binary margins

// copula likelihood for normal and binary margins
real binorm_cop_lp(real y1, real mu, real sigma, real y2, real p, real theta) {
real targ = normal_lpdf(y1|mu,sigma);
if (y2==0) {
targ = targ + normal_lcdf((inv_Phi(1-p)-theta*inv_Phi(normal_cdf(y1,mu,sigma)))/sqrt(1-theta^2)|0,1);
} else {
targ = targ + normal_lccdf((inv_Phi(1-p)-theta*inv_Phi(normal_cdf(y1,mu,sigma)))/sqrt(1-theta^2)|0,1);
}
return targ;
}
