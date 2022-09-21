// hierarchical model temp a[pref] b[pref] c[pref]
data {
  int N;
  real A[N]; //
  real S[N];
  real T[N];
  int P[N];
  int N_P;
  int N_Anew;
  real A_new[N_Anew];
  int N_Tnew;
  real T_new[N_Tnew];
}

parameters {
  real a;
  real a_p[N_P];
  real<lower=0> b_p[N_P];
  real<lower=0> c_p[N_P];
  real<lower=0, upper=1> b;
  real<lower=0> c;
  real<lower=0> d;
  real<lower=0> s_a;
  real<lower=0> s_b;
  real<lower=0> s_c;
  real<lower=0> beta;
}

transformed parameters {
  real S_pref[N];
  real alpha[N];
  for  (n in 1:N) {
    S_pref[n] = (a_p[P[n]]*T[n] + d) * (1-exp(-b_p[P[n]]*A[n]))^c_p[P[n]];
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
  real S_base_new[N_Anew];
  real S_baseL_new[N_Anew];
  real S_baseH_new[N_Anew];
  real S_pref_new[N_Anew,N_P];
  real S_new[N_Anew,N_P];
  real alpha_new[N_Anew,N_P];
  real A_base_new[N_Tnew];
  real A_area_new[N_Tnew,N_P];
  for (n in 1:N_Anew) {
    S_base_new[n] = (a*mean(T) + d) * (1-exp(-b*A_new[n]))^c;
    S_baseL_new[n] = (a*min(T) + d) * (1-exp(-b*A_new[n]))^c;
    S_baseH_new[n] = (a*max(T) + d) * (1-exp(-b*A_new[n]))^c;
    for (p in 1:N_P) {
      S_pref_new[n,p] = (a_p[p]*mean(T) + d) * (1-exp(-b_p[p]*A_new[n]))^c_p[p];
      alpha_new[n,p] = beta * S_pref_new[n,p];
      S_new[n,p] = gamma_rng(alpha_new[n,p], beta);
    }
  }
  for (n in 1:N_Tnew) {
    A_base_new[n] = a*T_new[n] + d;
    for (p in 1:N_P) {
      A_area_new[n,p] = a_p[p]*T_new[n] + d;
    }
  }
}

