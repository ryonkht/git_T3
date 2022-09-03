data {
  int N;
  real A[N]; //
  real S[N];
  int N_new;
  real A_new[N_new];
}

parameters {
  real<lower=0> a;
  real<lower=0, upper=1> b;
  real<lower=0> c;
  // real<lower=0> sigma;
  real<lower=0> beta;
}

transformed parameters {
  real S_base[N];
  real alpha[N];
  for  (n in 1:N) {
    S_base[n] = a*(1-exp(-b*A[n]))^c;
    // S_base[n] = a*A[n] + b;
    alpha[n] = beta * S_base[n];
  }
}

model {
  for (n in 1:N) {
    S[n] ~ gamma(alpha[n], beta);
  }
}

generated quantities {
  real S_base_new[N_new];
  real S_new[N_new];
  real alpha_new[N_new];
  for (n in 1:N_new) {
    S_base_new[n] = a*(1-exp(-b*A_new[n]))^c;
    // S_base_new[n] = a*A_new[n] + b;
    alpha_new[n] = beta * S_base_new[n];
    S_new[n] = gamma_rng(alpha_new[n], beta);
  }
}

