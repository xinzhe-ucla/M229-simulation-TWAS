data {
    int<lower=0> N; //number of individuals
    int<lower=0> K; //number of SNPs
    matrix[N, K] x; //scaled genotype
    vector[N] y; //expression
}

parameters {
    real alpha; //intercept
    vector[K] beta; //slope
    real<lower=0> sigma; //error
}

model {
    //priors
    alpha ~ normal(0, 1);
    beta ~ normal(0, 1);
    sigma ~ normal(0, 1);

    y ~ normal(x * beta + alpha, sigma); //likelihood
}
