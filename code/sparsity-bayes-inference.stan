data {
    int<lower=0> N; //number of individuals
    int<lower=0> K; //number of SNPs
    vector[K] z; //yes or no for effect size
    matrix[N, K] x; //scaled genotype
    vector[N] y; //expression
}

parameters {
    real alpha; //intercept
    vector[K] beta_raw; //beta before mixing
    real<lower=0> sigma; //error
}

transformed parameters {
    vector[K] beta;
    for (k in 1:K) {
        beta[k] = z[k] * beta_raw[k];  // Apply sparsity
    }
}

model {
    // generate:
    alpha ~ normal(0, 1);
    sigma ~ normal(0, 1);
    beta_raw ~ normal(0, 1);
    y ~ normal(x * beta + alpha, sigma); //likelihood
}
