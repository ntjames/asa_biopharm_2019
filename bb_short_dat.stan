// data statement for 2 binary outcomes, covar x, 
// and freq in each of 2^3 = 8 outcome/covar/
int<lower=1> N;
matrix[N, 2] x;
int<lower=0, upper=1> y_e[N];
int<lower=0, upper=1> y_s[N];
int<lower=0> freq[N];