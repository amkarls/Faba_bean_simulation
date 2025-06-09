
# Load required library
library(AlphaSimR)

#---------------------------------------------
# STEP 1: Simulate founder genomes
#---------------------------------------------
founderGenomes <- runMacs2(
  nInd = 187,           # Number of individuals in the diversity panel
  nChr = 6,             # Number of chromosomes
  segSites = 1100,      # Approximate number of segregating sites per chromosome
  Ne = 1104,            # Effective population size (estimated using SNeP)
  bp = 2166666666,      # Total base pairs per genome
  genLen = 2.16,        # Genetic length in Morgans (1 cM per 1M bp)
  mutRate = 2.5e-08,    # Mutation rate
  histNe = c(2500, 5000, 10000),  # Historical effective population sizes
  histGen = c(100, 1000, 2000),   # Generations corresponding to histNe
  inbred = FALSE,
  ploidy = 2
)

#---------------------------------------------
# STEP 2: Set up simulation parameters and trait architecture
#---------------------------------------------
SP <- SimParam$new(founderGenomes)

# Add a quantitative trait with additive and dominance effects
# meanDD specifies the average dominance deviation:
#   0   = no dominance
#   0.5 = partial dominance
#   1   = complete dominance
SP$addTraitADG(
  nQtlPerChr = 500,   # Number of QTL per chromosome
  mean = 1,           # Trait mean
  var = 1,            # Trait variance
  meanDD = 1,         # Mean dominance deviation (0, 0.5, or 1)
  varDD = 0.2,        # Variance in dominance effects
  varGxE = 0.2        # Genotype × environment interaction variance
)

# Add a SNP chip with minor allele frequency filtering
SP$addSnpChip(nSnpPerChr = 500, minSnpFreq = 0.01)

# Define heritability
H2 <- 0.2

#---------------------------------------------
# STEP 3: Create and phenotype diversity panel
#---------------------------------------------
Divpanel_D <- newPop(founderGenomes)
Divpanel_D <- setPheno(Divpanel_D, H2 = H2)

#---------------------------------------------
# STEP 4: Save output
#---------------------------------------------
save(Divpanel_D, file = "Divpanel_ADG.Rdata")
save(SP, file = "SP_ADG.Rdata")