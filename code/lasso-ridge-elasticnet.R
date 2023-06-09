### lasso-ridge-elasticnet.R ######################################################################
# purpose: using lasso ridge and elastic net for inference:

### PREAMBLE ######################################################################################
# load in the glmnet
library(glmnet);

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
result.dir <- '/Users/tardigrade/Documents/me/Grad_school/sriram-class/results/'

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

### APPLY GLMNET FOR INFERENCE ####################################################################
# for each of the instance of simulation, we will call the inference:
glm.collection <- vector('list', length = number.simulation);
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

    # deploy glmnet with l1 penalty:
    lasso <- cv.glmnet(
        x = genotype.z,
        y = simulated.expr[, instance],
        family = 'gaussian',
        alpha = 1
        );

    # deploy glmnet with l2 penalty:
    ridge <- cv.glmnet(
        x = genotype.z,
        y = simulated.expr[, instance],
        family = 'gaussian',
        alpha = 0
        );

    # deploy elastic net:
    elastic.net <- cv.glmnet(
        x = genotype.z,
        y = simulated.expr[, instance],
        family = 'gaussian',
        alpha = 0.5
        );  

    # store the coefs:
    glm.models <- list(
        lasso,
        ridge,
        elastic.net
        )
    names(glm.models) <- c('lasso', 'ridge', 'elastic');
    glm.collection[[instance]] <- glm.models;

    # print progress:
    progress_bar(current = instance, total = number.simulation);
    }
saveRDS(glm.collection, paste0(result.dir, 'glm-coefs.rds'));
