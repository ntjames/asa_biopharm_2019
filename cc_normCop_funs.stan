// functions for bivariate Normal copula with normal margins

// normal copula likelihood for 2 continuous normal margins
real bi_cop_lp(real y1, real mu1, real sigma1, real y2, real mu2, real sigma2, real rho) {
real Z1;
real Z2;
real lcP;

// since margins are normal can just make Z1, Z2 directly
// no need for pseudo-obs U1, U2 to plug into copula 
//real U1 = Phi((y1 - mu1)/sigma1); 
//real U2 = Phi((y2 - mu2)/sigma2); 
//Z = [inv_Phi(U1), inv_Phi(U2)]';

Z1 = (y1 - mu1)/sigma1;
Z2 = (y2 - mu2)/sigma2;

// log likelihood for joint dist. from Joe text p. 226
// use copula density from Meyer - bivariate normal copula paper
lcP = -0.5*log1m(rho^2) + ((2*rho*Z1*Z2 - rho^2*(Z1^2+Z2^2))/(2*(1-rho^2))) + normal_lpdf(y1|mu1, sigma1) + normal_lpdf(y2|mu2, sigma2);

// log-likelihood
return lcP;
}
