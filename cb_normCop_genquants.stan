// generated mean, risk difference, polyserial and pearson correlation quantities
vector[2] mu_e;
vector[2] p_s;
vector[2] theta_;
vector[2] rho;
real mu_e_diff;
real p_s_diff;

mu_e = beta_e;
p_s = Phi(beta_s);
theta_ = omega_;

rho[1] = theta_[1]*exp(normal_lpdf(inv_Phi(p_s[1])|0,1))/sqrt(p_s[1]*(1-p_s[1]));
rho[2] = theta_[2]*exp(normal_lpdf(inv_Phi(p_s[2])|0,1))/sqrt(p_s[2]*(1-p_s[2]));

mu_e_diff = mu_e[2]-mu_e[1]; 
p_s_diff = p_s[2]-p_s[1];
