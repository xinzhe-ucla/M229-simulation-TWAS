### PREAMBLE ######################################################################################
library(sim1000G)

# set seed to ensure reproducibility in simulation:
set.seed(103)

# load in the vcf from 1000 genome:
examples_dir <- system.file("examples", package = "sim1000G");
vcf_file <- file.path(examples_dir, "region.vcf.gz");
vcf <- readVCF(
    vcf_file,
    maxNumberOfVariants = 500,
    min_maf = 0.02,
    max_maf = NA,
    );

# downloadGeneticMap for chromosome 4:
genetic_map_of_region = system.file(
    "examples",
    "chr4-geneticmap.txt",
    package = "sim1000G"
    );
readGeneticMapFromFile(genetic_map_of_region)

# set up global parameters:
number.individuals <- 100;
simulation.runs <- 100;
output.dir <- '/Users/tardigrade/Documents/me/Grad_school/sriram-class/genotype-simulations-1000G/';

# write a function that prints out progress bar:
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
        cat('\n Prgress Complete! \n')
        }
    }

### SIMULATE GENOTYPE #############################################################################
for(simulation in seq(1, simulation.runs)) {
    # initiate simulation:
    invisible({capture.output({
        startSimulation(
            vcf = vcf,
            totalNumberOfIndividuals = number.individuals,
            );

        # simulate unrelated individuals:
        id = c()
        for(individual in 1:number.individuals) id[individual] <- SIM$addUnrelatedIndividual();

        # obtain unphased simulated genotype:
        genotypes <- SIM$gt1[1:number.individuals, ] + SIM$gt2[1:number.individuals, ];
        rownames(genotypes) <- paste0('individual', seq(1, number.individuals));
        colnames(genotypes) <- paste0('snp', seq(1, ncol(genotypes)));

        # save the genotypes that is imputed:
        write.table(
            genotypes,
            file = paste0(output.dir, 'seed', simulation, '-simulated-genotype-1000G-unrelated.txt'),
            sep = '\t',
            quote = FALSE,
            row.names = TRUE,
            col.names = TRUE
            );

        # reset the stage of simulation:
        SIM$reset();
        })})

    # print out progress bar during simulation:
    progress_bar(current = simulation, total = simulation.runs);
    }
