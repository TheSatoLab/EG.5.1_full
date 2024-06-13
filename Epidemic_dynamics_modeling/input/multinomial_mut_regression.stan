

data {
  int K;
  int N;
  int D;
  real S_w;
  matrix[N,2] TS;
  matrix[K,D] COEF;
  int Y[N,K];
  real generation_time;
  int bin_size;
  int Y_sum[N];
}

transformed data {
  real Max_date;
  matrix[N,2] TS_norm;
  Max_date = TS[N,2];
  for (n in 1:N) {
    TS_norm[n,1] = 1;
    TS_norm[n,2] = TS[n,2] / Max_date;
  }


  for (n in 1:N) {
    TS_norm[n] = TS[n] / Max_date;
  }
}

parameters {
  vector[K-1] b0_raw;
  vector[K] b1_raw;
  vector[D] w;
  real<lower=0> S_b1;
}

transformed parameters {
  vector[1] Zero;
  vector[K] b0;
  vector[K] b1;
  matrix[K,2] b;
  matrix[K,N] mu;
  Zero[1] = 0;
  b0 = append_row(Zero,b0_raw);
  b1 = b1_raw * S_b1;
  b = append_col(b0,b1);
  mu = b * TS_norm'; 
}

model {
  b1_raw ~ student_t(5, COEF * w,1);
  w ~ double_exponential(0, S_w);
  for (n in 1:N)
    Y[n,] ~ multinomial_logit(mu[,n]);
  S_b1 ~ student_t(5, 0, 10);

}



generated quantities {
  vector[K] growth_rate;
  vector[D] growth_gain;
  matrix[N,K] theta;
  int Y_predict[N,K];
  
  for(k in 1:K){
      growth_rate[k] = exp(((b1[k] / Max_date) / bin_size)  * generation_time);
  }
  for(d in 1:D){
      growth_gain[d] = exp(((w[d] / Max_date) / bin_size)  * generation_time);
  }

  for(n in 1:N){
    theta[n,] = softmax(mu[,n])';
    Y_predict[n,] = multinomial_rng(softmax(mu[,n]),Y_sum[n]);
  }

}







