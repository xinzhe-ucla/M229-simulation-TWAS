data {
    int<lower=0> N; //number of individuals
    int<lower=0> K; //number of SNPs
    real<lower=0> b; //background variance
    real<lower=0> a; //additional variance
    real<lower=0, upper=1> p; //probability
    matrix[N, K] x; //scaled genotype
    vector[N] y; //expression
}

parameters {
    real alpha; //intercept
    vector[K] background;
    vector[K] additional;
    real<lower=0> sigma; //error
}

transformed parameters {
    vector[K] beta;
    beta = p * background + (1 - p) * additional;
}


model {
    //component
    background ~ normal(0, b);
    additional ~ normal(0, a + b);

    // generate:
    alpha ~ normal(0, 1);
    sigma ~ normal(0, 1);
    y ~ normal(x * beta + alpha, sigma); //likelihood
}
