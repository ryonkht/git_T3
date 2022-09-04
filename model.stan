data {
  int N;
  real A[N]; //
  real S[N];
  int P[N];
  int N_P;
  int N_new;
  real A_new[N_new];
}

parameters {
  real<lower=0> a;
  real<lower=0> a_p[N_P];
  real<lower=0> b_p[N_P];
  real<lower=0> c_p[N_P];
  real<lower=0, upper=1> b;
  real<lower=0> c;
  real<lower=0> s_a;
  real<lower=0> s_b;
  real<lower=0> s_c;
  real<lower=0> beta;
}

transformed parameters {
  real S_pref[N];
  real alpha[N];
  for  (n in 1:N) {
    S_pref[n] = a_p[P[n]]*(1-exp(-b_p[P[n]]*A[n]))^c_p[P[n]];
    alpha[n] = beta * S_pref[n];
  }
}

model {
  for (n in 1:N_P) {
    a_p[n] ~ normal(a, s_a);
    b_p[n] ~ normal(b, s_b);
    c_p[n] ~ normal(c, s_c);
  }
  for (n in 1:N) {
    S[n] ~ gamma(alpha[n], beta);
  }
}

generated quantities {
  real S_base_new[N_new];
  real S_pref_new[N_new,N_P];
  real S_new[N_new,N_P];
  real alpha_new[N_new,N_P];
  for (n in 1:N_new) {
    S_base_new[n] = a*(1-exp(-b*A_new[n]))^c;
    for (p in 1:N_P) {
      S_pref_new[n,p] = a_p[p]*(1-exp(-b_p[p]*A_new[n]))^c_p[p];
      alpha_new[n,p] = beta * S_pref_new[n,p];
      S_new[n,p] = gamma_rng(alpha_new[n,p], beta);
    }
  }
}

