## =============================================================================
## ARGUMENTS / INTERACTIVE MODE
## =============================================================================

## If running interactively in R, define rep_id and workdir first.
## Otherwise, read them from commandArgs().

if (exists("rep_id", inherits = FALSE) && exists("workdir", inherits = FALSE)) {
  # Interactive mode: use pre-defined variables
  rep_id  <- as.integer(rep_id)
  workdir <- as.character(workdir)
} else {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 1) {
    stop("Usage: Linebreeding_GS_1stage_DD1.R <rep_id> [workdir]")
  }
  rep_id  <- as.integer(args[1])
  workdir <- if (length(args) >= 2) args[2] else getwd()
}

if (!dir.exists(workdir)) {
  stop("Workdir does not exist: ", workdir)
}

setwd(workdir)


## =============================================================================
## LIBRARIES
## =============================================================================

suppressPackageStartupMessages({
  library(AlphaSimR)
  library(dplyr)
  library(tibble)
  library(BGLR)
})


## =============================================================================
## THREADING
## =============================================================================

## Always keep BLAS single-threaded inside an array task
Sys.setenv(
  OMP_NUM_THREADS        = "1",
  MKL_NUM_THREADS        = "1",
  OPENBLAS_NUM_THREADS   = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS    = "1"
)


## =============================================================================
## INPUT FILES
## =============================================================================

load("SP_ADG_DD1.RData")
load("Divpanel_ADG_DD1.RData")
load(sprintf("workspace_ADG_DD1_%d_2.RData", rep_id))


## =============================================================================
## MAIN SIMULATION FUNCTION
## =============================================================================

model.mse <- function(rep_id) {
  
  ## ---------------------------------------------------------------------------
  ## Initial parent population
  ## ---------------------------------------------------------------------------
  
  Parents   = c(Parents_17, Parents_18, Parents_19)
  newtrain2 = Parentssel
  
  
  ## ---------------------------------------------------------------------------
  ## Simulation parameters
  ## ---------------------------------------------------------------------------
  
  # Error variance for field trials
  varE  = 4
  repF6 = 9
  repF7 = 9
  
  # Number of individuals selected per stage
  nF2sel = 1000
  nF3fam = 100
  nF6sel = 10
  nF7sel = 2
  
  # Number of crosses per year
  nCrosses = 50
  
  # Number of individuals produced per cross
  nind = 16
  
  # Degree of selfing
  selfing = 0.5
  
  
  ## ---------------------------------------------------------------------------
  ## Initial training population
  ## ---------------------------------------------------------------------------
  
  # Create one DH line from each accession of the parent set and diversity panel
  training = makeDH(mergePops(list(Parents, Divpanel)), nDH = 1)
  training = setPheno(training, varE = varE, reps = repF7)
  
  # Add 100 random accessions from F4 (not phenotyped until year 1)
  add_training = vector("list", length(F4_split))
  for (i in 1:length(F4_split)) {
    add_training[[i]] <- selectInd(F4_split[[i]], nInd = 1, use = "rand")
  }
  add_training = mergePops(add_training)
  add_training = selectInd(add_training, nInd = 100, use = "rand")
  
  
  ## ---------------------------------------------------------------------------
  ## Training genotypes and phenotypes
  ## ---------------------------------------------------------------------------
  
  M_tr <- pullSnpGeno(training)   # rows = individuals, cols = SNPs, coded 0/1/2
  y    <- as.numeric(training@pheno)
  
  
  ## ---------------------------------------------------------------------------
  ## Additive and dominance design matrices
  ## ---------------------------------------------------------------------------
  
  # Center using training allele frequencies
  p_tr  <- colMeans(M_tr) / 2
  Z_A   <- sweep(M_tr, 2, 2 * p_tr, "-")
  
  D_raw <- 1 * (M_tr == 1)
  Z_D   <- sweep(D_raw, 2, 2 * p_tr * (1 - p_tr), "-")
  
  
  ## ---------------------------------------------------------------------------
  ## Fit BGLR model (additive only at initialization)
  ## ---------------------------------------------------------------------------
  
  ETA <- list(
    list(X = Z_A, model = "BRR")
    # list(X = Z_D, model = "BRR")  # dominance not used; all training are DH
  )
  
  fit <- BGLR(
    y       = y,
    ETA     = ETA,
    nIter   = 10000,
    burnIn  = 3000,
    verbose = FALSE
  )
  
  # Marker effects
  bA <- fit$ETA[[1]]$b
  # bD <- fit$ETA[[2]]$b
  mu <- as.numeric(fit$mu)
  
  
  ## ---------------------------------------------------------------------------
  ## NA-safe cleanup
  ## ---------------------------------------------------------------------------
  
  nNA_A <- sum(is.na(bA))
  # nNA_D <- sum(is.na(bD))
  
  if (nNA_A > 0) {
    message("Warning: ", nNA_A, " NA in bA and setting them to 0")
    bA[is.na(bA)] <- 0
    # bD[is.na(bD)] <- 0
  }
  
  if (is.na(mu)) {
    message("Warning: mu is NA; setting mu = 0")
    mu <- 0
  }
  
  
  ## ---------------------------------------------------------------------------
  ## Prepare breeding cycle
  ## ---------------------------------------------------------------------------
  
  # Modify input files so breeding cycle can be run as 20 consecutive sequences
  F1 = mergePops(F1)
  
  nYears = 20
  output = data.frame(year = 1:nYears)
  
  
  ## ===========================================================================
  ## Main breeding loop
  ## ===========================================================================
  
  for (year in 1:nYears) {
    
    print("Breeding cycle year:", year)
    
    if (year == 1) {
      
      ## -----------------------------------------------------------------------
      ## Year 6: Create variety seed from selected F7
      ## -----------------------------------------------------------------------
      
      variety = vector("list", length(select_F7))
      for (i in 1:length(select_F7)) {
        variety_i = selectOP(
          select_F7[[i]],
          nInd          = 200,
          nSeeds        = 10,
          probSelf      = selfing,
          pollenControl = FALSE,
          use           = "rand"
        )
        variety[[i]] = selectInd(variety_i, nInd = 2000, use = "rand")
      }
      
      outputvariety = mergePops(variety)
      
      
      ## -----------------------------------------------------------------------
      ## Year 5: Select F7 families
      ## -----------------------------------------------------------------------
      
      F7_fampheno = data.frame(1:length(F7), NA)
      for (i in 1:length(F7)) {
        F7_fampheno_i   = mean(F7[[i]]@pheno)
        F7_fampheno[i,2] = F7_fampheno_i
      }
      
      colnames(F7_fampheno) <- c("id", "mean_pheno")
      topx_F7   <- dplyr::top_n(F7_fampheno, nF7sel, wt = mean_pheno)
      select_F7 <- F7[topx_F7$id]
      
      
      ## -----------------------------------------------------------------------
      ## Year 4: Select F6 families and advance to F7
      ## -----------------------------------------------------------------------
      
      F6_fampheno = data.frame(1:length(F6), NA)
      for (i in 1:length(F6)) {
        F6_fampheno_i   = mean(F6[[i]]@pheno)
        F6_fampheno[i,2] = F6_fampheno_i
      }
      
      colnames(F6_fampheno) <- c("id", "mean_pheno")
      topx_F6   <- dplyr::top_n(F6_fampheno, nF6sel, wt = mean_pheno)
      select_F6 <- F6[topx_F6$id]
      
      F7 = vector("list", length(select_F6))
      for (i in 1:length(select_F6)) {
        F7_i = selectOP(
          select_F6[[i]],
          nInd          = 200,
          nSeeds        = 10,
          probSelf      = selfing,
          pollenControl = FALSE,
          use           = "rand"
        )
        F7_i    = setPheno(F7_i, varE = varE, reps = repF7)
        F7[[i]] = selectInd(F7_i, nInd = 2000, use = "rand")
      }
      
      output_F7 = mergePops(F7)
      
      
      ## -----------------------------------------------------------------------
      ## Year 3: Advance F4 to F5 and F6
      ## -----------------------------------------------------------------------
      
      F5 = vector("list", length(F4_split))
      for (i in 1:length(F4_split)) {
        F5_i = selectOP(
          F4_split[[i]],
          nInd          = 200,
          nSeeds        = 10,
          probSelf      = selfing,
          pollenControl = FALSE,
          use           = "rand"
        )
        F5[[i]] = selectInd(F5_i, nInd = 2000, use = "rand")
      }
      
      F6 = vector("list", length(F5))
      for (i in 1:length(F5)) {
        F6_i = selectOP(
          F5[[i]],
          nInd          = 200,
          nSeeds        = 10,
          probSelf      = selfing,
          pollenControl = FALSE,
          use           = "rand"
        )
        F6_i    = setPheno(F6_i, varE = varE, reps = repF6)
        F6[[i]] = selectInd(F6_i, nInd = 2000, use = "rand")
      }
      
      output_F6 = mergePops(F6)
      
      
      ## -----------------------------------------------------------------------
      ## Year 2: Create F2, predict EBVs, select, and create F3/F4
      ## -----------------------------------------------------------------------
      
      # Create seed for F2
      F2 = selectOP(
        F1,
        nInd          = 800,
        nSeeds        = 16,
        probSelf      = 0.99,
        pollenControl = TRUE,
        use           = "rand"
      )
      
      # Predict EBV for F2
      M_F2 <- pullSnpGeno(F2)
      
      if (anyNA(M_F2)) {
        for (j in seq_len(ncol(M_F2))) {
          nas <- is.na(M_F2[, j])
          if (any(nas)) M_F2[nas, j] <- 2 * p_tr[j]
        }
      }
      
      Z_A_F2 <- sweep(M_F2, 2, 2 * p_tr, "-")
      D_F2   <- 1 * (M_F2 == 1)
      Z_D_F2 <- sweep(D_F2, 2, 2 * p_tr * (1 - p_tr), "-")
      
      gF2 <- as.numeric(Z_A_F2 %*% bA) + mu
      
      F2@ebv <- as.matrix(gF2)
      
      # Select and create F3
      F3 = selectOP(
        F2,
        nInd          = nF2sel,
        nSeeds        = 10,
        probSelf      = 0.99,
        pollenControl = TRUE,
        use           = "ebv"
      )
      
      # Predict EBV for F3
      M_F3 <- pullSnpGeno(F3)
      
      if (anyNA(M_F3)) {
        for (j in seq_len(ncol(M_F3))) {
          nas <- is.na(M_F3[, j])
          if (any(nas)) M_F3[nas, j] <- 2 * p_tr[j]
        }
      }
      
      Z_A_F3 <- sweep(M_F3, 2, 2 * p_tr, "-")
      D_F3   <- 1 * (M_F3 == 1)
      Z_D_F3 <- sweep(D_F3, 2, 2 * p_tr * (1 - p_tr), "-")
      
      gF3 <- as.numeric(Z_A_F3 %*% bA) + mu
      
      F3@ebv <- as.matrix(gF3)
      
      # Select families from F3
      F3_select = selectFam(F3, nF3fam, use = "ebv", famType = "F")
      
      # Split by mother and generate F4 seed
      motherSval_F4 = unique(F3_select@mother)
      F4_split      = vector("list", length(motherSval_F4))
      
      for (i in 1:length(motherSval_F4)) {
        motherpop = subset(F3_select, F3_select@mother == paste(motherSval_F4[i]))
        F4_split_i = selectOP(
          pop            = motherpop,
          nInd           = 10,
          nSeeds         = 40,
          probSelf       = 0.95,
          pollenControl  = FALSE,
          use            = "rand"
        )
        F4_split[[i]] = selectInd(F4_split_i, nInd = 400, use = "rand")
      }
      
      
      ## -----------------------------------------------------------------------
      ## Year 1: Create new F1 and update training population
      ## -----------------------------------------------------------------------
      
      # Parents come from last year's F6 and F7 populations
      F1 = randCross(newtrain2, nCrosses = nCrosses, nProgeny = nind)
      
      # Phenotype F4 additions from burn-in year 20
      add_training = setPheno(add_training, reps = repF7, varE = varE)
      
      # Add new accessions from F6 and F7 trials
      newtrain = c(F6, F7)
      
      newtrain2 = vector("list", length(newtrain))
      for (i in 1:length(newtrain)) {
        newtrain2[[i]] <- selectInd(newtrain[[i]], nInd = 1, use = "rand")
      }
      newtrain2 = mergePops(newtrain2)
      
      # Update training population
      training <- mergePops(list(training, newtrain2, add_training))
      
      M_tr <- pullSnpGeno(training)
      y    <- as.numeric(training@pheno)
      
      if (anyNA(M_tr)) {
        cm <- colMeans(M_tr, na.rm = TRUE)
        for (j in seq_len(ncol(M_tr))) {
          nas <- is.na(M_tr[, j])
          if (any(nas)) M_tr[nas, j] <- cm[j]
        }
      }
      
      p_tr  <- colMeans(M_tr) / 2
      Z_A   <- sweep(M_tr, 2, 2 * p_tr, "-")
      
      D_raw <- 1 * (M_tr == 1)
      Z_D   <- sweep(D_raw, 2, 2 * p_tr * (1 - p_tr), "-")
      
      ETA <- list(
        list(X = Z_A, model = "BRR"),
        list(X = Z_D, model = "BRR")
      )
      
      set.seed(1)
      fit <- BGLR(
        y       = y,
        ETA     = ETA,
        nIter   = 10000,
        burnIn  = 3000,
        verbose = FALSE
      )
      
      bA <- fit$ETA[[1]]$b
      bD <- fit$ETA[[2]]$b
      mu <- as.numeric(fit$mu)
      
      nNA_A <- sum(is.na(bA))
      nNA_D <- sum(is.na(bD))
      if (nNA_A > 0 || nNA_D > 0) {
        message("Warning: ", nNA_A, " NA in bA and ", nNA_D, " NA in bD; setting them to 0")
        bA[is.na(bA)] <- 0
        bD[is.na(bD)] <- 0
      }
      
      if (is.na(mu)) {
        message("Warning: mu is NA; setting mu = 0")
        mu <- 0
      }
      
      
      ## -----------------------------------------------------------------------
      ## Report results
      ## -----------------------------------------------------------------------
      
      output$meanGF3[year]      = meanG(F3)
      output$varGF3[year]       = varG(F3)
      output$meanGvariety[year] = meanG(outputvariety)
      output$varGvariety[year]  = varG(outputvariety)
      output$meanGParents[year] = meanG(newtrain2)
      output$varGParents[year]  = varG(newtrain2)
      output$varAParents[year]  = varA(newtrain2)
      output$varDParents[year]  = varD(newtrain2)
      output$accF2[year]        = cor(F2@gv, F2@ebv)
      output$accF3[year]        = cor(F3@gv, F3@ebv)
      output$accF6[year]        = cor(output_F6@gv, output_F6@pheno)
      output$accF7[year]        = cor(output_F7@gv, output_F7@pheno)
    }
    
    if (year > 1) {
      
      ## -----------------------------------------------------------------------
      ## Year 6: Create variety seed from selected F7
      ## -----------------------------------------------------------------------
      
      variety = vector("list", length(select_F7))
      for (i in 1:length(select_F7)) {
        variety_i = selectOP(
          select_F7[[i]],
          nInd          = 200,
          nSeeds        = 10,
          probSelf      = selfing,
          pollenControl = FALSE,
          use           = "rand"
        )
        variety[[i]] = selectInd(variety_i, nInd = 2000, use = "rand")
      }
      
      outputvariety = mergePops(variety)
      
      
      ## -----------------------------------------------------------------------
      ## Year 5: Select F7 families
      ## -----------------------------------------------------------------------
      
      F7_fampheno = data.frame(1:length(F7), NA)
      for (i in 1:length(F7)) {
        F7_fampheno_i   = mean(F7[[i]]@pheno)
        F7_fampheno[i,2] = F7_fampheno_i
      }
      
      colnames(F7_fampheno) <- c("id", "mean_pheno")
      topx_F7   <- dplyr::top_n(F7_fampheno, nF7sel, wt = mean_pheno)
      select_F7 <- F7[topx_F7$id]
      
      
      ## -----------------------------------------------------------------------
      ## Year 4: Select F6 families and advance to F7
      ## -----------------------------------------------------------------------
      
      F6_fampheno = data.frame(1:length(F6), NA)
      for (i in 1:length(F6)) {
        F6_fampheno_i   = mean(F6[[i]]@pheno)
        F6_fampheno[i,2] = F6_fampheno_i
      }
      
      colnames(F6_fampheno) <- c("id", "mean_pheno")
      topx_F6   <- dplyr::top_n(F6_fampheno, nF6sel, wt = mean_pheno)
      select_F6 <- F6[topx_F6$id]
      
      F7 = vector("list", length(select_F6))
      for (i in 1:length(select_F6)) {
        F7_i = selectOP(
          select_F6[[i]],
          nInd          = 200,
          nSeeds        = 10,
          probSelf      = selfing,
          pollenControl = FALSE,
          use           = "rand"
        )
        F7_i    = setPheno(F7_i, varE = varE, reps = repF7)
        F7[[i]] = selectInd(F7_i, nInd = 2000, use = "rand")
      }
      
      output_F7 = mergePops(F7)
      
      
      ## -----------------------------------------------------------------------
      ## Year 3: Advance F4 to F5 and F6
      ## -----------------------------------------------------------------------
      
      F5 = vector("list", length(F4_split))
      for (i in 1:length(F4_split)) {
        F5_i = selectOP(
          F4_split[[i]],
          nInd          = 200,
          nSeeds        = 10,
          probSelf      = selfing,
          pollenControl = FALSE,
          use           = "rand"
        )
        F5[[i]] = selectInd(F5_i, nInd = 2000, use = "rand")
      }
      
      F6 = vector("list", length(F5))
      for (i in 1:length(F5)) {
        F6_i = selectOP(
          F5[[i]],
          nInd          = 200,
          nSeeds        = 10,
          probSelf      = selfing,
          pollenControl = FALSE,
          use           = "rand"
        )
        F6_i    = setPheno(F6_i, varE = varE, reps = repF6)
        F6[[i]] = selectInd(F6_i, nInd = 2000, use = "rand")
      }
      
      output_F6 = mergePops(F6)
      
      
      ## -----------------------------------------------------------------------
      ## Year 2: Create F2, predict EBVs, select, and create F3/F4
      ## -----------------------------------------------------------------------
      
      F2 = selectOP(
        F1,
        nInd          = 800,
        nSeeds        = 16,
        probSelf      = 0.99,
        pollenControl = TRUE,
        use           = "rand"
      )
      
      M_F2 <- pullSnpGeno(F2)
      
      if (anyNA(M_F2)) {
        for (j in seq_len(ncol(M_F2))) {
          nas <- is.na(M_F2[, j])
          if (any(nas)) M_F2[nas, j] <- 2 * p_tr[j]
        }
      }
      
      Z_A_F2 <- sweep(M_F2, 2, 2 * p_tr, "-")
      D_F2   <- 1 * (M_F2 == 1)
      Z_D_F2 <- sweep(D_F2, 2, 2 * p_tr * (1 - p_tr), "-")
      
      gF2 <- as.numeric(Z_A_F2 %*% bA + Z_D_F2 %*% bD) + mu
      
      F2@ebv <- as.matrix(gF2)
      
      F3 = selectOP(
        F2,
        nInd          = nF2sel,
        nSeeds        = 10,
        probSelf      = 0.99,
        pollenControl = TRUE,
        use           = "ebv"
      )
      
      M_F3 <- pullSnpGeno(F3)
      
      if (anyNA(M_F3)) {
        for (j in seq_len(ncol(M_F3))) {
          nas <- is.na(M_F3[, j])
          if (any(nas)) M_F3[nas, j] <- 2 * p_tr[j]
        }
      }
      
      Z_A_F3 <- sweep(M_F3, 2, 2 * p_tr, "-")
      D_F3   <- 1 * (M_F3 == 1)
      Z_D_F3 <- sweep(D_F3, 2, 2 * p_tr * (1 - p_tr), "-")
      
      gF3 <- as.numeric(Z_A_F3 %*% bA + Z_D_F3 %*% bD) + mu
      
      F3@ebv <- as.matrix(gF3)
      
      F3_select = selectFam(F3, nF3fam, use = "ebv", famType = "F")
      
      motherSval_F4 = unique(F3_select@mother)
      F4_split      = vector("list", length(motherSval_F4))
      
      for (i in 1:length(motherSval_F4)) {
        motherpop = subset(F3_select, F3_select@mother == paste(motherSval_F4[i]))
        F4_split_i = selectOP(
          pop            = motherpop,
          nInd           = 10,
          nSeeds         = 40,
          probSelf       = 0.95,
          pollenControl  = FALSE,
          use            = "rand"
        )
        F4_split[[i]] = selectInd(F4_split_i, nInd = 400, use = "rand")
      }
      
      
      ## -----------------------------------------------------------------------
      ## Year 1: Create new F1 and update training population
      ## -----------------------------------------------------------------------
      
      F1 = randCross(newtrain2, nCrosses = nCrosses, nProgeny = nind)
      
      newtrain = c(F6, F7)
      
      newtrain2 = vector("list", length(newtrain))
      for (i in 1:length(newtrain)) {
        newtrain2[[i]] <- selectInd(newtrain[[i]], nInd = 1, use = "rand")
      }
      newtrain2 = mergePops(newtrain2)
      
      training <- mergePops(list(training, newtrain2))
      
      M_tr <- pullSnpGeno(training)
      y    <- as.numeric(training@pheno)
      
      if (anyNA(M_tr)) {
        cm <- colMeans(M_tr, na.rm = TRUE)
        for (j in seq_len(ncol(M_tr))) {
          nas <- is.na(M_tr[, j])
          if (any(nas)) M_tr[nas, j] <- cm[j]
        }
      }
      
      p_tr  <- colMeans(M_tr) / 2
      Z_A   <- sweep(M_tr, 2, 2 * p_tr, "-")
      
      D_raw <- 1 * (M_tr == 1)
      Z_D   <- sweep(D_raw, 2, 2 * p_tr * (1 - p_tr), "-")
      
      ETA <- list(
        list(X = Z_A, model = "BRR"),
        list(X = Z_D, model = "BRR")
      )
      
      set.seed(1)
      fit <- BGLR(
        y       = y,
        ETA     = ETA,
        nIter   = 10000,
        burnIn  = 3000,
        verbose = FALSE
      )
      
      bA <- fit$ETA[[1]]$b
      bD <- fit$ETA[[2]]$b
      mu <- as.numeric(fit$mu)
      
      nNA_A <- sum(is.na(bA))
      nNA_D <- sum(is.na(bD))
      if (nNA_A > 0 || nNA_D > 0) {
        message("Warning: ", nNA_A, " NA in bA and ", nNA_D, " NA in bD; setting them to 0")
        bA[is.na(bA)] <- 0
        bD[is.na(bD)] <- 0
      }
      
      if (is.na(mu)) {
        message("Warning: mu is NA; setting mu = 0")
        mu <- 0
      }
      
      
      ## -----------------------------------------------------------------------
      ## Report results
      ## -----------------------------------------------------------------------
      
      output$meanGF3[year]      = meanG(F3)
      output$varGF3[year]       = varG(F3)
      output$meanGvariety[year] = meanG(outputvariety)
      output$varGvariety[year]  = varG(outputvariety)
      output$meanGParents[year] = meanG(newtrain2)
      output$varGParents[year]  = varG(newtrain2)
      output$varAParents[year]  = varA(newtrain2)
      output$varDParents[year]  = varD(newtrain2)
      output$accF2[year]        = cor(F2@gv, F2@ebv)
      output$accF3[year]        = cor(F3@gv, F3@ebv)
      output$accF6[year]        = cor(output_F6@gv, output_F6@pheno)
      output$accF7[year]        = cor(output_F7@gv, output_F7@pheno)
    }
  }
  
  output
}


## =============================================================================
## RUN MODEL AND SAVE OUTPUT
## =============================================================================

set.seed(1000 + rep_id)  # reproducible per replicate
out <- model.mse(rep_id)

saveRDS(
  out,
  file = file.path(
    workdir,
    sprintf("Linebreeding_GS_1stage_BGLR_DD1_rep%03d.rds", rep_id)
  )
)