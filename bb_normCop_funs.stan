// functions for bivariate Normal copula with binary margins

// bivariate normal cdf (from Stan manual)
real binormal_cdf(real z1, real z2, real rho){
if (z1!=0 || z2 !=0){
real denom = fabs(rho) < 1.0 ? sqrt((1+rho)*(1+rho)): not_a_number();
real a1 = (z2/z1 - rho)/denom;
real a2 = (z1/z2 - rho)/denom;
real product = z1*z2;
real delt = product < 0 || (product==0 && (z1+z2)<0);
return 0.5*(Phi(z1)+Phi(z2)-delt)-owens_t(z1,a1)-owens_t(z2,a2);
}
return 0.25 + asin(rho)/(2*pi());
}

// bivariate normal copula distribution function 
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

// copula likelihood for 2 binary margins
real bi_cop_lp(int y1, real p1, int y2, real p2, real theta) {

// pseudo-obs from bernoulli dist
real U1 = bernoulli_cdf(y1, p1); 
real U1_prime = bernoulli_cdf(y1-1, p1);
real U2 = bernoulli_cdf(y2, p2); 
real U2_prime = bernoulli_cdf(y2-1, p2);

// likelihood for joint dist. using differences
// replace with sum() function?
real cP0 = pCop_norm(U1, U2, theta) -
pCop_norm(U1_prime, U2, theta) -
pCop_norm(U1, U2_prime, theta) +
pCop_norm(U1_prime, U2_prime, theta);

// put lower bound on cP to avoid log(0)
// is there a better way to do this??
real cP = cP0 < 1e-200 ? 1e-200: cP0;

// log-likelihood
return log(cP);
}