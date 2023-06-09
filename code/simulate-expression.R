### simulate-expression.R #########################################################################
# purpose: simualte the expression from the simulated genotype data:

### PREAMBLE ######################################################################################
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

# set seed for reproducbility:
set.seed(103)

# define some parameters for the simulation:
causal.number <- 30;
causal.sd <- 1;
scaling.factor <- 1;

# directory of simulated genotype:
genotype.dir <- '/Users/tardigrade/Documents/me/Grad_school/sriram-class/genotype-simulations-1000G/';
output.dir <- '/Users/tardigrade/Documents/me/Grad_school/sriram-class/expression-simulation-1000G/';

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

### SIMULATE EXPRESSION ###########################################################################
# initiate some placeholder matrix:
expression.instances <- matrix(NA, nrow = number.individuals, ncol = number.simulation);

# randomly determine the effect sizes:
causal.variants <- sample(number.variants, size = causal.number, replace = FALSE);
causal.effect <- rnorm(n = causal.number, sd = causal.sd);
names(causal.effect) <- causal.variants;

cat('COMMENCING SIMULATION \n');
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

    # compute the genetic component toward expression:
    genetic.component <- genotype.z[, causal.variants] %*% causal.effect;

    # now compute the amount of variance from the genetic component:
    genetic.variance <- var(genetic.component);

    # next simulate the environmental component:
    environment.variance <- scaling.factor * genetic.variance;
    residual <- rnorm(n = nrow(genotype.z), sd = sqrt(environment.variance));

    # add the two component to get the simulated expression value:
    expression.data <- genetic.component + residual;

    # store the simulated expression:
    expression.instances[, instance] <- expression.data;

    # print progress:
    progress_bar(current = instance, total = number.simulation);
    }

### DATA SAVING ###################################################################################
# save the expression data:
write.table(
    expression.instances,
    file = paste0(output.dir, 'simulated-expression-unrelated.txt'),
    sep = '\t',
    quote = FALSE,
    row.names = TRUE,
    col.names = TRUE
    );

# save the causal snp and its effect size:
causality.matrix <- data.frame(
    snp = names(causal.effect),
    beta = causal.effect
    );

write.table(
    causality.matrix,
    file = paste0(output.dir, 'simulated-causal-effects.txt'),
    sep = '\t',
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
    );
