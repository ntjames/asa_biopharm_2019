// generated risk difference quantities
vector[2] p_e;
vector[2] p_s;
real p_e_diff;
real p_s_diff;

p_e = Phi(beta_e);
p_s = Phi(beta_s);
p_e_diff = p_e[2]-p_e[1];
p_s_diff = p_s[2]-p_s[1];