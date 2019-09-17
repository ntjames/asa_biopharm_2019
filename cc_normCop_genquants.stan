// generated mean, mean difference, and pearson correlation quantities
vector[2] mu_e;
vector[2] mu_s;
vector[2] theta_;
//vector[2] rho;
real mu_e_diff;
real mu_s_diff;

mu_e = beta_e;
mu_s = beta_s;
theta_ = omega_;

//rho[1] = theta_[1]*exp(normal_lpdf(inv_Phi(p_s[1])|0,1))/sqrt(p_s[1]*(1-p_s[1]));
//rho[2] = theta_[2]*exp(normal_lpdf(inv_Phi(p_s[2])|0,1))/sqrt(p_s[2]*(1-p_s[2]));

mu_e_diff = mu_e[2]-mu_e[1]; 
mu_s_diff = mu_s[2]-mu_s[1]; 
