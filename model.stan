// hierarchical model temp a[pref] b[pref] c[pref]
data {
  int N;
  real A[N]; //
  real S[N];
  real T[N];
  // vector[N] T;
  int R[N];
  int N_R;
  int N_Anew;
  real A_new[N_Anew];
  int N_Tnew;
  real T_new[N_Tnew];
}

parameters {
  real a;
  real a_r[N_R];
  // real<lower=-5, upper=5> a;
  // row_vector<lower=-2.5, upper=2.5>[N_R] a_r;
  real<lower=0> b_r[N_R];
  real<lower=0> c_r[N_R];
  real<lower=0, upper=1> b;
  real<lower=0> c;
  // real<lower=max(-(T*a_r))> d;
  real<lower=0> d;
  real<lower=0> s_a;
  real<lower=0> s_b;
  real<lower=0> s_c;
  real<lower=0> beta;
}

transformed parameters {
  real S_rand[N];
  real alpha[N];
  for  (n in 1:N) {
    S_rand[n] = (a_r[R[n]]*T[n] + d) * (1-exp(-b_r[R[n]]*A[n]))^c_r[R[n]];
    alpha[n] = beta * S_rand[n];
  }
}

model {
  for (n in 1:N_R) {
    a_r[n] ~ normal(a, s_a);
    b_r[n] ~ normal(b, s_b);
    c_r[n] ~ normal(c, s_c);
  }
  for (n in 1:N) {
    S[n] ~ gamma(alpha[n], beta);
  }
}

generated quantities {
  real S_base_new_aveT[N_Anew]; // basic S with A at mean T
  real S_rand_new_aveT[N_Anew,N_R]; // basic S with A by rand at mean T
  real S_base_new[N_Anew,N_Tnew];
  real S_rand_new[N_Anew,N_Tnew,N_R];
  real S_new[N_Anew,N_Tnew,N_R];
  real alpha_new[N_Anew,N_Tnew,N_R];
  real Smax_base_new[N_Tnew];
  real Smax_rand_new[N_Tnew,N_R];
  for (n in 1:N_Anew) {
    S_base_new_aveT[n] = (a*mean(T) + d) * (1-exp(-b*A_new[n]))^c;
    for (r in 1:N_R) {
      S_rand_new_aveT[n,r] = (a_r[r]*mean(T) + d) * (1-exp(-b_r[r]*A_new[n]))^c_r[r];
    }
    for (t in 1:N_Tnew) {
      S_base_new[n,t] = (a*T_new[t] + d) * (1-exp(-b*A_new[n]))^c;
      for (r in 1:N_R) {
        S_rand_new[n,t,r] = (a_r[r]*T_new[t] + d) * (1-exp(-b_r[r]*A_new[n]))^c_r[r];
        alpha_new[n,t,r] = beta * S_rand_new[n,t,r];
        S_new[n,t,r] = gamma_rng(alpha_new[n,t,r], beta);
      }
    }
  }
  for (n in 1:N_Tnew) {
    Smax_base_new[n] = a*T_new[n] + d;
    for (r in 1:N_R) {
      Smax_rand_new[n,r] = a_r[r]*T_new[n] + d;
    }
  }
}

