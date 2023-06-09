### PREAMBLE ######################################################################################
library(rstan);

# load in the true causal effect:
simulated.dir <- '/Users/tardigrade/Documents/me/Grad_school/sriram-class/expression-simulation-1000G/';
stan.dir <- '/Users/tardigrade/Documents/me/Grad_school/sriram-class/simulation-M229/code/';
causal.effect <- read.table(
    file = paste0(simulated.dir, 'simulated-causal-effects.txt'),
    sep = '\t',
    stringsAsFactors = FALSE,
    header = TRUE
    );

# load in the expression:
simulated.expr <- read.table(
    file = paste0(simulated.dir, 'simulated-expression-unrelated.txt'),
    sep = '\t',
    stringsAsFactors = FALSE,
    header = TRUE, 
    row.names = 1
    );

# directory of simulated genotype:
genotype.dir <- '/Users/tardigrade/Documents/me/Grad_school/sriram-class/genotype-simulations-1000G/';
result.dir <- '/Users/tardigrade/Documents/me/Grad_school/sriram-class/simulation-M229/results/'
# find out the number of siulated genotypes:
genotype.simulations <- list.files(path = genotype.dir, full.names = FALSE);

# load in one instance of genotype simulation to get simulation parameters:
first.instance <- read.table(
    file = paste0(genotype.dir, genotype.simulations[1]),
    sep = '\t',
    stringsAsFactors = FALSE,
    row.names = 1,
    header = TRUE
    );

number.individuals <- nrow(first.instance);
number.variants <- ncol(first.instance);
number.simulation = length(genotype.simulations);

# define the function for progress bar:
progress_bar <- function(current, total) {
    width = options("width")$width;
    percent <- current / total;
    progress <- round(percent * width);
    left <- width - progress;

    cat(
        '\r', '[', paste(rep('*', progress), collapse = ''), '*',
        paste(rep(' ', left), collapse = ''), '] ',
        sprintf('%3.0f%%', percent * 100),
        sep = ''
        )
  
    if (current == total) {
        cat('\n Complete \n')
        }
    }

### MODEL #########################################################################################
# for each of the instance of simulation, we will call the inference:
fit.coef <- vector('list', length = number.simulation);
bernoulli.vector <- rbinom(n = number.variants, size = 1, prob = 0.1);
non.convergence.indication <- fit.coef;

for (instance in seq(1, number.simulation)) {
    # load in the genotype:
    genotype.instance <- read.table(
        file = paste0(genotype.dir, genotype.simulations[instance]),
        sep = '\t',
        stringsAsFactors = FALSE,
        row.names = 1,
        header = TRUE
        );

    # scale the genotype instance:
    genotype.z <- scale(genotype.instance);

    # call model:
    data <- list(
        N = number.individuals,
        K = number.variants,
        x = genotype.z,
        y = simulated.expr[, instance],
        z = bernoulli.vector
        );
    invisible({capture.output({
        start <- Sys.time();
        stan.fit <- stan(
            file = paste0(stan.dir, 'sparsity-bayes-inference.stan'),
            data = data,
            chains = 6,
            iter = 1000,
            cores = 6,
            verbose = FALSE,
            open_progress = FALSE,
            control = list(adapt_delta = 0.99, max_treedepth = 15)
            );
        end <- Sys.time()
        cat('time of code chunk: ', end - start, '\n')
        })})
    # see divergence issue with rhat > 1.1
    problematic.coef <- sum((summary(stan.fit))$summary[,'Rhat'] > 1.1);
    non.convergence.indication[[instance]] <- problematic.coef;

    # store results:
    fit.coef[[instance]] <- as.data.frame(stan.fit);

    # print progress:
    progress_bar(current = instance, total = number.simulation);
    }
saveRDS(fit.coef, paste0(result.dir, 'sparsity-assumption-coefs.rds'));
