## =============================================================================
## ARGUMENTS / INTERACTIVE MODE
## If you're in R, define rep_id & workdir first; otherwise read commandArgs.
## =============================================================================

if (exists("rep_id", inherits = FALSE) && exists("workdir", inherits = FALSE)) {
  # Interactive: use the pre-defined variables
  rep_id  <- as.integer(rep_id)
  workdir <- as.character(workdir)
} else {
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) < 1) {
    stop("Usage: Linebreeding_phenotypic_DD1.R <rep_id> [workdir]")
  }
  
  rep_id  <- as.integer(args[1])
  workdir <- if (length(args) >= 2) args[2] else getwd()
}

if (!dir.exists(workdir)) {
  stop("Workdir does not exist: ", workdir)
}

setwd(workdir)

suppressPackageStartupMessages({
  library(AlphaSimR)
  library(dplyr)
  library(tibble)
})

## Always keep BLAS single-threaded inside an array task
Sys.setenv(
  OMP_NUM_THREADS       = "1",
  MKL_NUM_THREADS       = "1",
  OPENBLAS_NUM_THREADS  = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS   = "1"
)

## =============================================================================
## USE WORKDIR FOR ALL I/O
## =============================================================================

load("SP_ADG_DD1.RData")
load(sprintf("workspace_ADG_DD1_%d_2.RData", rep_id))

## =============================================================================
## SIMULATION FUNCTION
## =============================================================================

model.mse <- function(rep_id) {
  
  ##### function
  
  # ---------------------------------------------------------------------------
  # Simulation Parameters
  # ---------------------------------------------------------------------------
  
  # Replication at each stage
  varE  = 4   # Error variance for field trials
  repF2 = 1 / 50
  repF3 = 1 / 4
  repF4 = 1
  repF5 = 2
  repF6 = 6
  repF7 = 9
  
  # Number of crosses per year
  nCrosses = 50
  
  # Individuals produced per cross
  nind = 16
  
  # Degree of selfing when accession allowed to freely open pollinate
  # in isolation from other accessions
  selfing = 0.5
  
  # Number of individuals selected per selection stage
  nF1sel = 256  # Random plants selected for seed multiplication
  nF2sel = 20   # Individual plants selected per cross in greenhouse
  nF3fam = 300  # Number of families selected in F3
  nF4sel = 100  # Lines selected in F4
  nF5sel = 25   # Lines selected in F5
  nF6sel = 10   # Lines selected in F6
  nF7sel = 2    # Lines selected in F7
  
  # ---------------------------------------------------------------------------
  # Helper: pick top indices by mean phenotype
  # ---------------------------------------------------------------------------
  
  top_idx <- function(lst, k) {
    vals <- vapply(lst, function(p) mean(p@pheno), numeric(1))
    order(vals, decreasing = TRUE)[seq_len(min(k, length(vals)))]
  }
  
  # ---------------------------------------------------------------------------
  # Create output structure
  # ---------------------------------------------------------------------------
  
  nYears = 20
  output = data.frame(year = 1:nYears)
  
  for (year in 1:nYears) {
    
    # =========================================================================
    # Year 7: Create variety seed (F8)
    # =========================================================================
    
    variety = vector("list", length(select_F7))  # Keep lines separate
    
    for (i in 1:length(select_F7)) {
      variety_i = selectOP(
        select_F7[[i]],
        nInd = 200,
        nSeeds = 10,
        probSelf = selfing, ###
        pollenControl = FALSE,
        use = "rand"
      )
      
      variety[[i]] = selectInd(variety_i, nInd = 2000, use = "rand")
    }
    
    outputvariety = mergePops(variety)
    
    # =========================================================================
    # Year 6: Phenotype F7 in Svalöv/Sweden and select two lines that will be
    # developed into varieties. Create parents for crossing.
    # =========================================================================
    
    # Multiply seeds for F7 and phenotype in multilocation trials
    F7 = vector("list", length(select_F6))  # Keep lines separate
    
    for (i in 1:length(select_F6)) {
      F7_i = selectOP(
        select_F6[[i]],
        nInd = 200,
        nSeeds = 10,
        probSelf = selfing, ###
        pollenControl = FALSE,
        use = "rand"
      )
      
      F7_i = setPheno(
        F7_i,
        reps = repF7,
        varE = varE
      ) ############Refers to phenotyping in Svalöv/Sweden at F7 stage
      
      F7[[i]] = selectInd(F7_i, nInd = 2000, use = "rand")
    }
    
    output_F7 = mergePops(F7)
    
    # Select top 2 lines
    select_F7 <- F7[top_idx(F7, nF7sel)]
    
    # Create parents
    Parents = c(F6, F7)  ### F7 are seeds produced by F6 in off-season location
    
    # Select one parental individual from each F5 and F6 accession randomly
    Parentssel = vector("list", length(Parents))  # Keep crosses seperate
    
    for (i in 1:length(Parents)) {
      Parentssel[[i]] <- selectInd(Parents[[i]], nInd = 1, use = "rand")
    }
    
    Parentssel = mergePops(Parentssel)
    
    # =========================================================================
    # Year 5: Multiply seeds from lines selected in F5 and trial F6 in Svalöv
    # =========================================================================
    
    # Multiply F5 in off-season location to create F6
    F6 = vector("list", length(select_F5))  # Keep lines separate
    
    for (i in 1:length(select_F5)) {
      F6_i = selectOP(
        select_F5[[i]],
        nInd = 200,
        nSeeds = 10,
        probSelf = selfing, ###
        pollenControl = FALSE,
        use = "rand"
      )
      
      F6_i = setPheno(
        F6_i,
        reps = repF6,
        varE = varE
      ) ############Refers to phenotyping in Svalöv/Sweden at F6 stage
      
      F6[[i]] = selectInd(F6_i, nInd = 2000, use = "rand")
    }
    
    output_F6 = mergePops(F6)
    
    # Select 10 lines from F6 based on trial in Svalov
    select_F6 <- F6[top_idx(F6, nF6sel)]
    
    # =========================================================================
    # Year 4: Multiplication of seeds from F4 in Chile and F5 trial in
    # Svalöv/Sweden
    # =========================================================================
    
    # Multiply F4 seeds in off-season location to create F5
    F5 = vector("list", length(select_F4))  # Keep lines separate
    
    for (i in 1:length(select_F4)) {
      F5_i = selectOP(
        select_F4[[i]],
        nInd = 200,
        nSeeds = 10,
        probSelf = selfing, ###
        pollenControl = FALSE,
        use = "rand"
      )
      
      F5_i = setPheno(
        F5_i,
        reps = repF5,
        varE = varE
      ) ##########Refers to phenotyping in Svalöv/Sweden of F5s
      
      F5[[i]] = selectInd(F5_i, nInd = 2000, use = "rand")
    }
    
    output_F5 = mergePops(F5)
    
    # Select top 25 lines
    select_F5 <- F5[top_idx(F5, nF5sel)]
    
    # =========================================================================
    # Year 3: F3 selection in Chile and trial of F4 in Svalöv
    # =========================================================================
    
    # Select among F3 to create the F4 to be trialled in Svalov.
    F3_merged = mergePops(F3)
    
    F3_select = selectFam(
      F3_merged,
      nFam = 300,
      trait = 1,
      use = "pheno",
      famType = "F"
    )
    
    # Split population by half-sib and collect within-family open-pollinated
    # seeds for each "line" in Chile
    motherSval_F4 = unique(F3_select@mother)
    F4_split = vector("list", length(motherSval_F4))
    
    for (i in 1:length(motherSval_F4)) {  # Loop over maternal parent
      motherpop = subset(F3_select, F3_select@mother == paste(motherSval_F4[i]))
      
      F4_split_i = selectOP(
        pop = motherpop,
        nInd = 40,
        nSeeds = 40,
        probSelf = 0.95, ###
        pollenControl = FALSE,
        use = "rand"
      )
      
      F4_split_i = setPheno(
        F4_split_i,
        reps = repF4,
        varE = varE
      ) ###Corresponds to the phenotype of the F4 in Svalov
      
      F4_split[[i]] = selectInd(F4_split_i, nInd = 1600)
    }
    
    output_F4 = mergePops(F4_split)
    
    # Select in F4 for seed multiplication in off-season location
    select_F4 <- F4_split[top_idx(F4_split, nF4sel)]
    
    # =========================================================================
    # Year 2: F2 planted to create seed for F3
    # =========================================================================
    
    # Create F3
    F3 = vector("list", nCrosses)  # Keep crosses separate
    
    for (i in 1:nCrosses) {  # Loop over crosses
      F3_i = selectOP(
        F2[[i]],
        nInd = 20,
        nSeeds = 40,
        probSelf = 0.95, ###
        pollenControl = FALSE,
        use = "pheno"
      )
      
      F3_i = setPheno(F3_i, reps = repF3, varE = varE)
      F3[[i]] = selectInd(F3_i, nInd = 800, use = "rand")
    }
    
    # =========================================================================
    # Year 1: Crosses, F1 and create seed for F2
    # =========================================================================
    
    F1 = vector("list", nCrosses)  # Keep crosses seperate
    
    for (i in 1:nCrosses) {  # Loop over crosses
      F1[[i]] = randCross(Parentssel, nCrosses = 1, nProgeny = nind)
    }
    
    F2 = vector("list", nCrosses)  # Keep crosses separate
    
    for (i in 1:nCrosses) {  # Loop over crosses
      F2_i = selectOP(
        F1[[i]],
        nInd = nind,
        nSeeds = 16,
        probSelf = 0.95, ###
        pollenControl = FALSE,
        use = "rand"
      )
      
      F2_i = setPheno(F2_i, reps = repF2, varE = varE)
      F2[[i]] = selectInd(F2_i, nInd = 200, use = "rand")
    }
    
    outputF2 = mergePops(F2)
    
    # -------------------------------------------------------------------------
    # Report results
    # -------------------------------------------------------------------------
    
    output$meanGF3[year]      = meanG(F3_merged)
    output$varGF3[year]       = varG(F3_merged)
    
    output$meanGvariety[year] = meanG(outputvariety)
    output$varGvariety[year]  = varG(outputvariety)
    
    output$meanGParents[year] = meanG(Parentssel)
    output$varGParents[year]  = varG(Parentssel)
    output$varAParents[year]  = varA(Parentssel, simParam = SP)
    output$varDParents[year]  = varD(Parentssel, simParam = SP)
    
    output$accF2[year] = cor(outputF2@gv, outputF2@pheno)
    output$accF3[year] = cor(F3_merged@gv, F3_merged@pheno)
    output$accF4[year] = cor(output_F4@gv, output_F4@pheno)
    output$accF5[year] = cor(output_F5@gv, output_F5@pheno)
    output$accF6[year] = cor(output_F6@gv, output_F6@pheno)
    output$accF7[year] = cor(output_F7@gv, output_F7@pheno)
  }
  
  output
}

## =============================================================================
## RUN MODEL
## =============================================================================

set.seed(1000 + rep_id)  # reproducible per replicate
out <- model.mse(rep_id)

saveRDS(
  out,
  file = file.path(
    workdir,
    sprintf("Linebreeding_phenotypic_DD1_rep%03d.rds", rep_id)
  )
)