############ ------------------------------------------------------------
############ Training population data based on progeny test
############ followed by synthetic performance
############ ------------------------------------------------------------

## Workflow
## 1. Split heterotic groups by markers
## 2. Train prediction models on performance of synthetics
##    derived from each accession in the training population
## 3. Separate the training population by heterotic group


## ======================================================================
## ARGUMENTS / INTERACTIVE MODE
## ======================================================================
## If running interactively in R, define rep_id and workdir first.
## Otherwise, read from commandArgs().

if (exists("rep_id", inherits = FALSE) && exists("workdir", inherits = FALSE)) {
  rep_id  <- as.integer(rep_id)
  workdir <- as.character(workdir)
} else {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 1) {
    stop("Usage: GS_synthetic_Twopart_DD1.R <rep_id> [workdir]")
  }
  rep_id  <- as.integer(args[1])
  workdir <- if (length(args) >= 2) args[2] else getwd()
}

if (!dir.exists(workdir)) {
  stop("Workdir does not exist: ", workdir)
}
setwd(workdir)


## ======================================================================
## PACKAGES
## ======================================================================
suppressPackageStartupMessages({
  library(AlphaSimR)
  library(dplyr)
  library(rrBLUP)
  library(tibble)
  library(stats)
  library(cluster)
  library(BGLR)
})

## Keep BLAS single-threaded inside an array task
Sys.setenv(
  OMP_NUM_THREADS       = "1",
  MKL_NUM_THREADS       = "1",
  OPENBLAS_NUM_THREADS  = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS   = "1"
)


## ======================================================================
## INPUT FILES
## ======================================================================
load("SP_ADG_DD1.RData")
load("Divpanel_ADG_DD1.RData")
load(sprintf("workspace_ADG_DD1_%d_2.RData", rep_id))


## ======================================================================
## SIMULATION FUNCTION
## ======================================================================
model.mse <- function(rep_id) {
  
  ## --------------------------------------------------------------------
  ## PARAMETERS
  ## --------------------------------------------------------------------
  
  # Parents
  nCyclesPI <- 2
  nParents  <- 20
  nOP       <- 8
  nElite    <- 1
  
  nSelf     <- 16   # number of seeds per selfed F1
  nself1    <- 1
  nself2    <- 10
  
  s2sel     <- 200  # number of individuals selected from S2 per heterotic group
  s3sel     <- 20   # number of lines selected from each heterotic group
  
  nsel_syn1 <- 4    # number of lines per heterotic group
  nsel_syn2 <- 2    # number of varieties selected
  
  # Degree of selfing
  selfing   <- 0.5
  
  # Trial accuracy
  varE                <- 4   # error variance for field trials
  reptraining         <- 3
  repsyn1             <- 3
  repsyn2             <- 9
  reporiginaltraining <- 6
  
  
  ## --------------------------------------------------------------------
  ## HELPER FUNCTIONS
  ## --------------------------------------------------------------------
  
  # Predict EBV using additive and dominance effects estimated by BGLR
  predictEBV_BGLR <- function(pop, tester_pop, bA, bD) {
    
    # Genotypes: parents and testers
    geno_par  <- pullSnpGeno(pop)        # nInd x m, 0/1/2
    geno_test <- pullSnpGeno(tester_pop) # nTester x m, 0/1/2
    
    # Safety checks
    m_par <- ncol(geno_par)
    m_bA  <- length(bA)
    m_bD  <- length(bD)
    
    if (!identical(m_par, m_bA) || !identical(m_par, m_bD)) {
      stop("Length of bA/bD does not match number of markers in pop.")
    }
    
    # 1. Allele frequencies
    p_parent  <- geno_par / 2
    p_testers <- colMeans(geno_test / 2)
    
    # 2. Expected synthetic allele frequency
    p_syn <- (sweep(p_parent, 2, p_testers, "+")) / 2
    
    # 3. Additive and dominance genotype codes
    genoA <- 2 * (p_syn - 0.5)
    genoD <- 2 * p_syn * (1 - p_syn)
    
    # 4. EBV = A*bA + D*bD
    ebv <- as.vector(genoA %*% bA + genoD %*% bD)
    
    # 5. Store and return
    pop@ebv <- matrix(ebv)
    pop
  }
  
  
  # Compute expected additive and dominance genotype codes per synthetic
  build_design_for_year <- function(tr_pop, tst_pop) {
    
    # Genotypes: training parents and testers
    geno_par  <- pullSnpGeno(tr_pop)   # nTrain x m, 0/1/2
    geno_test <- pullSnpGeno(tst_pop)  # nTest x m, 0/1/2
    
    nTrain <- nrow(geno_par)
    m      <- ncol(geno_par)
    
    # 1. Allele frequencies
    p_parent  <- geno_par / 2
    p_testers <- colMeans(geno_test / 2)
    
    # 2. Expected synthetic allele frequencies
    p_syn <- (sweep(p_parent, 2, p_testers, "+")) / 2
    
    # 3. Additive and dominance codes
    A_blk <- 2 * (p_syn - 0.5)
    D_blk <- 2 * p_syn * (1 - p_syn)
    
    # Keep marker names
    colnames(A_blk) <- colnames(geno_par)
    colnames(D_blk) <- colnames(geno_par)
    
    list(A = A_blk, D = D_blk)
  }
  
  
  ## --------------------------------------------------------------------
  ## TRAINING POPULATION
  ## --------------------------------------------------------------------
  
  # Add 100 random accessions from F4_split to boost the training population
  add_training <- vector("list", length(F4_split))
  
  for (i in seq_along(F4_split)) {
    add_training[[i]] <- selectInd(F4_split[[i]], nInd = 1, use = "rand")
  }
  
  add_training <- mergePops(add_training)
  add_training <- selectInd(add_training, nInd = 100, use = "rand")
  # These will be phenotyped in year 1
  
  
  ## --------------------------------------------------------------------
  ## CREATE HETEROTIC GROUPS
  ## --------------------------------------------------------------------
  
  nGroups <- 5
  
  # Starting population
  start_pop <- c(Parents_17, Parents_18, Parents_19)
  
  # 1. Extract marker matrix
  geno <- pullSnpGeno(start_pop)
  
  # 2. Genetic distance matrix
  G <- A.mat(geno, min.MAF = 0.05, impute.method = "mean")
  D2 <- outer(diag(G), diag(G), "+") - 2 * G
  D2[D2 < 0] <- 0
  D <- as.dist(sqrt(D2))
  
  # 3. Hierarchical clustering
  hc <- hclust(D, method = "ward.D2")
  
  # 4. Cut dendrogram into clusters
  clusters <- cutree(hc, k = nGroups)
  
  # 5. Assign cluster labels
  start_pop@misc$cluster <- clusters
  
  rm(geno, G, D, hc, clusters)
  
  
  ## --------------------------------------------------------------------
  ## CREATE TESTERS FOR TRAINING POPULATION GCA EVALUATION
  ## --------------------------------------------------------------------
  
  tester_list <- vector("list", nGroups)
  
  for (k in 1:nGroups) {
    cluster_inds <- which(start_pop@misc$cluster == k)
    selected_ind <- sample(cluster_inds, 1)
    tester_list[[k]] <- start_pop[selected_ind]
  }
  
  # One inbred line developed per tester
  tester_pop <- makeDH(mergePops(tester_list))
  
  
  ## --------------------------------------------------------------------
  ## CREATE TRAINING POPULATION
  ## --------------------------------------------------------------------
  
  # Create one inbred line from each accession of diversity panel
  # (start_pop should already be largely homozygous)
  training <- makeDH(mergePops(list(start_pop, Divpanel)), nDH = 1)
  
  
  ## --------------------------------------------------------------------
  ## SYNTHETIC PERFORMANCE AS TRAINING DATA
  ## --------------------------------------------------------------------
  
  # Simulate open pollination to produce one synthetic line per training individual
  for (individual in 1:length(training)) {
    training_ind <- training[individual]
    train_syn    <- mergePops(list(training_ind, tester_pop))
    
    train_syn <- selectOP(
      pop           = train_syn,
      nInd          = 6,
      nSeeds        = 200,
      probSelf      = selfing,
      pollenControl = FALSE,
      use           = "rand"
    )
    
    train_syn <- setPheno(train_syn, varE = varE, reps = reporiginaltraining)
    training@pheno[individual] <- meanP(train_syn)
  }
  
  
  ## --------------------------------------------------------------------
  ## STORE TRAINING AND TESTER POPULATIONS BY YEAR
  ## --------------------------------------------------------------------
  
  training_by_year <- list()
  tester_by_year   <- list()
  
  training_by_year[["Y0"]] <- training
  tester_by_year[["Y0"]]   <- tester_pop
  
  year_key <- paste0("Y", 0)
  training_by_year[[year_key]] <- training
  tester_by_year[[year_key]]   <- tester_pop
  
  
  ## --------------------------------------------------------------------
  ## BUILD DESIGN MATRICES BY YEAR
  ## --------------------------------------------------------------------
  
  genoA_list <- vector("list", length(training_by_year))
  genoD_list <- vector("list", length(training_by_year))
  y_list     <- vector("list", length(training_by_year))
  
  for (k in seq_along(training_by_year)) {
    tr_pop  <- training_by_year[[k]]
    tst_pop <- tester_by_year[[k]]
    
    y_list[[k]] <- as.numeric(tr_pop@pheno)
    
    dmats <- build_design_for_year(tr_pop, tst_pop)
    genoA_list[[k]] <- dmats$A
    genoD_list[[k]] <- dmats$D
  }
  
  genoA <- do.call(rbind, genoA_list)
  genoD <- do.call(rbind, genoD_list)
  y     <- unlist(y_list)
  
  
  ## --------------------------------------------------------------------
  ## FIT INITIAL GS MODEL
  ## --------------------------------------------------------------------
  
  ETA <- list(
    list(X = genoA, model = "BRR"),  # additive
    list(X = genoD, model = "BRR")   # dominance
  )
  
  fit <- BGLR(
    y       = y,
    ETA     = ETA,
    nIter   = 10000,
    burnIn  = 3000,
    verbose = FALSE
  )
  
  bA <- fit$ETA[[1]]$b
  bD <- fit$ETA[[2]]$b
  
  # NA-safe cleanup
  nNA_A <- sum(is.na(bA))
  nNA_D <- sum(is.na(bD))
  
  if (nNA_A > 0 || nNA_D > 0) {
    message(
      "Warning: ", nNA_A, " NA in bA and ", nNA_D,
      " NA in bD; setting them to 0"
    )
    bA[is.na(bA)] <- 0
    bD[is.na(bD)] <- 0
  }
  
  
  ## --------------------------------------------------------------------
  ## POPULATIONS FOR FILLING THE PIPELINE
  ## --------------------------------------------------------------------
  
  parent_groups <- lapply(1:nGroups, function(k) {
    subset(start_pop, start_pop@misc$cluster == k)
  })
  
  benchmarking <- lapply(
    parent_groups,
    function(pop) selectInd(pop, nInd = 3, use = "rand")
  )
  
  elite_groups <- lapply(parent_groups, function(group) {
    selectInd(group, nElite, use = "rand")
  })
  
  
  ## --------------------------------------------------------------------
  ## PRODUCE ELITE LINE GROUPS
  ## --------------------------------------------------------------------
  
  G1 <- c(elite_groups[[2]], elite_groups[[3]], elite_groups[[4]], elite_groups[[5]])
  G2 <- c(elite_groups[[1]], elite_groups[[3]], elite_groups[[4]], elite_groups[[5]])
  G3 <- c(elite_groups[[4]], elite_groups[[1]], elite_groups[[2]], elite_groups[[5]])
  G4 <- c(elite_groups[[1]], elite_groups[[2]], elite_groups[[3]], elite_groups[[5]])
  G5 <- c(elite_groups[[1]], elite_groups[[2]], elite_groups[[3]], elite_groups[[4]])
  
  synthetics <- list(G1, G2, G3, G4, G5)
  
  
  ## --------------------------------------------------------------------
  ## START BREEDING CYCLE
  ## --------------------------------------------------------------------
  
  nYears <- 20
  output <- data.frame(year = 1:nYears)
  
  for (year in 1:nYears) {
    
    cat("Cycle year:", year, "of 20\n")
    
    if (year == 1) {
      
      ## ================================================================
      ## YEAR 1
      ## ================================================================
      
      ## ---------------------------------------------------------------
      ## Create first synthetic varieties
      ## ---------------------------------------------------------------
      
      # SYN0
      syn0 <- replicate(n = 5, expr = { list() }, simplify = FALSE)
      
      for (group in seq_along(parent_groups)) {
        for (i in 1:length(parent_groups[[group]])) {
          open1 <- mergePops(list(parent_groups[[group]][i], synthetics[[group]]))
          syn0[[group]][i] <- selectOP(
            open1,
            nInd          = 5,
            nSeeds        = 500,
            probSelf      = selfing,
            pollenControl = FALSE,
            use           = "rand"
          )
        }
      }
      
      # SYN1
      syn1 <- replicate(n = 5, expr = { list() }, simplify = FALSE)
      
      for (group in seq_along(syn1)) {
        for (i in 1:length(syn0[[group]])) {
          syn1[[group]][[i]] <- selectOP(
            syn0[[group]][[i]],
            nInd     = 2000,
            nSeeds   = 1,
            use      = "rand",
            probSelf = selfing
          )
          syn1[[group]][[i]] <- setPheno(
            syn1[[group]][[i]],
            varE = varE,
            reps = repsyn1
          )
        }
      }
      
      
      ## ---------------------------------------------------------------
      ## Use synthetic performance as training data
      ## ---------------------------------------------------------------
      
      for (group in 1:length(parent_groups)) {
        for (individual in 1:length(parent_groups[[group]])) {
          training_ind <- parent_groups[[group]][individual]
          train_syn    <- mergePops(list(training_ind, tester_pop))
          
          train_syn <- selectOP(
            pop           = train_syn,
            nInd          = 6,
            nSeeds        = 200,
            probSelf      = selfing,
            pollenControl = FALSE,
            use           = "rand"
          )
          
          train_syn <- setPheno(train_syn, varE = varE, reps = reporiginaltraining)
          syn1[[group]][[individual]]@pheno[, 1] <- meanP(train_syn)
        }
      }
      
      
      ## ---------------------------------------------------------------
      ## Select best SYN1 subpopulations
      ## ---------------------------------------------------------------
      
      Syn1_mean <- vector("list", length(syn1))
      
      for (group in seq_along(syn1)) {
        group_means <- sapply(
          syn1[[group]],
          function(pop) mean(pop@pheno, na.rm = TRUE)
        )
        Syn1_mean[[group]] <- data.frame(
          index      = seq_along(group_means),
          mean_pheno = group_means
        )
      }
      
      Syn1_sel <- vector("list", length(syn1))
      
      for (group in seq_along(syn1)) {
        mean_df     <- Syn1_mean[[group]]
        top_indices <- order(mean_df$mean_pheno, decreasing = TRUE)[1:nsel_syn1]
        Syn1_sel[[group]] <- syn1[[group]][top_indices]
      }
      
      Syn1_sel <- unlist(Syn1_sel, recursive = FALSE)
      
      
      ## ---------------------------------------------------------------
      ## Produce SYN2 and select SYN2
      ## ---------------------------------------------------------------
      
      syn2 <- lapply(Syn1_sel, function(line) {
        pop <- selectOP(
          line,
          nInd          = 500,
          nSeeds        = 4,
          probSelf      = selfing,
          pollenControl = FALSE,
          use           = "rand"
        )
        setPheno(pop, varE = varE, reps = repsyn2)
      })
      
      Syn2_sel <- sample(syn2, size = nsel_syn2, replace = FALSE)
      
      
      ## ---------------------------------------------------------------
      ## Create SYN3
      ## ---------------------------------------------------------------
      
      syn3 <- vector("list", length(Syn2_sel))
      
      for (group in seq_along(Syn2_sel)) {
        syn3[[group]] <- selectOP(
          pop           = Syn2_sel[[group]],
          nInd          = 200,
          nSeeds        = 10,
          probSelf      = selfing,
          pollenControl = FALSE,
          use           = "rand"
        )
      }
      
      synthetic_var <- mergePops(syn3)
      
      
      ## ---------------------------------------------------------------
      ## YEAR 3: Split F3 into heterotic groups and select
      ## ---------------------------------------------------------------
      
      F3      <- mergePops(F3)
      maf_min <- 0.05
      
      G_ind   <- pullSnpGeno(F3)
      AF_grpL <- lapply(parent_groups, \(p) colMeans(pullSnpGeno(p), na.rm = TRUE) / 2)
      AF_grp  <- do.call(rbind, AF_grpL)
      
      p_all <- colMeans(rbind(G_ind / 2, AF_grp), na.rm = TRUE)
      
      keep  <- which(p_all > maf_min & p_all < 1 - maf_min)
      G_ind <- G_ind[, keep, drop = FALSE]
      AF_grp <- AF_grp[, keep, drop = FALSE]
      p_all <- p_all[keep]
      m     <- length(keep)
      
      std <- sqrt(2 * p_all * (1 - p_all))
      std[std == 0] <- 1
      
      Z_ind <- sweep(G_ind, 2, 2 * p_all, "-")
      Z_ind <- sweep(Z_ind, 2, std, "/")
      Z_ind[!is.finite(Z_ind)] <- 0
      
      Z_grp <- sweep(2 * AF_grp, 2, 2 * p_all, "-")
      Z_grp <- sweep(Z_grp, 2, std, "/")
      Z_grp[!is.finite(Z_grp)] <- 0
      
      S <- Z_ind %*% t(Z_grp) / m
      group_assign <- max.col(S, ties.method = "random")
      
      F3@misc$cluster <- group_assign
      S3 <- lapply(
        split(seq_len(nInd(F3)), group_assign),
        function(idx) F3[idx]
      )
      
      S3 <- lapply(S3, function(F3) {
        predictEBV_BGLR(F3, tester_pop, bA, bD)
      })
      
      S3s_sel <- lapply(S3, function(F3) {
        selectFam(F3, nFam = s3sel, use = "ebv", famType = "F")
      })
      
      F3_split_list <- vector("list", length(S3s_sel))
      
      for (group in seq_along(S3s_sel)) {
        group_pop  <- S3s_sel[[group]]
        mother_ids <- unique(group_pop@mother)
        F3_split   <- vector("list", length(mother_ids))
        
        for (i in seq_along(mother_ids)) {
          F3_split[[i]] <- subset(group_pop, group_pop@mother == mother_ids[i])
        }
        
        names(F3_split) <- as.character(mother_ids)
        F3_split_list[[group]] <- F3_split
      }
      
      S3s_sel <- F3_split_list
      
      
      ## ---------------------------------------------------------------
      ## Produce S4 and S5
      ## ---------------------------------------------------------------
      
      S4 <- vector("list", length(S3s_sel))
      for (group in seq_along(S3s_sel)) {
        for (line in 1:length(S3s_sel[[group]])) {
          S4[[group]][[line]] <- self(S3s_sel[[group]][[line]], nProgeny = 200)
        }
      }
      
      S5 <- vector("list", length(S4))
      for (group in seq_along(S4)) {
        for (line in 1:length(S4[[group]])) {
          S5[[group]][[line]] <- self(S4[[group]][[line]], nProgeny = 1)
        }
      }
      
      
      ## ---------------------------------------------------------------
      ## YEAR 2: Split F2 into heterotic groups and select
      ## ---------------------------------------------------------------
      
      F2_i <- mergePops(F2)
      
      G_ind   <- pullSnpGeno(F2_i)
      AF_grpL <- lapply(parent_groups, \(p) colMeans(pullSnpGeno(p), na.rm = TRUE) / 2)
      AF_grp  <- do.call(rbind, AF_grpL)
      
      p_all <- colMeans(rbind(G_ind / 2, AF_grp), na.rm = TRUE)
      
      keep  <- which(p_all > maf_min & p_all < 1 - maf_min)
      G_ind <- G_ind[, keep, drop = FALSE]
      AF_grp <- AF_grp[, keep, drop = FALSE]
      p_all <- p_all[keep]
      m     <- length(keep)
      
      std <- sqrt(2 * p_all * (1 - p_all))
      std[std == 0] <- 1
      
      Z_ind <- sweep(G_ind, 2, 2 * p_all, "-")
      Z_ind <- sweep(Z_ind, 2, std, "/")
      Z_ind[!is.finite(Z_ind)] <- 0
      
      Z_grp <- sweep(2 * AF_grp, 2, 2 * p_all, "-")
      Z_grp <- sweep(Z_grp, 2, std, "/")
      Z_grp[!is.finite(Z_grp)] <- 0
      
      S <- Z_ind %*% t(Z_grp) / m
      group_assign <- max.col(S, ties.method = "random")
      
      F2_i@misc$cluster <- group_assign
      assigned_groups <- lapply(
        split(seq_len(nInd(F2_i)), group_assign),
        function(idx) F2_i[idx]
      )
      
      S2 <- lapply(assigned_groups, function(group) {
        predictEBV_BGLR(group, tester_pop, bA, bD)
      })
      
      S2s_sel <- lapply(S2, function(group) {
        selectInd(group, nInd = s2sel, use = "ebv")
      })
      
      
      ## ---------------------------------------------------------------
      ## YEAR 1: Two cycles of population improvement
      ## ---------------------------------------------------------------
      
      for (cycle in 1:nCyclesPI) {
        
        if (cycle == 1) {
          
          parent_groups <- lapply(parent_groups, function(group) {
            predictEBV_BGLR(group, tester_pop, bA, bD)
          })
          
          Parents_output <- mergePops(parent_groups)
          
          output$meanGParents[year] <- meanG(Parents_output)
          output$varGParents[year]  <- varG(Parents_output)
          output$varAParents[year]  <- varA(Parents_output)
          output$varDParents[year]  <- varD(Parents_output)
          
          parent_groups <- lapply(parent_groups, function(group) {
            selectOP(
              group,
              nInd          = length(group),
              nSeeds        = nOP,
              probSelf      = 0.2,
              pollenControl = FALSE,
              use           = "ebv"
            )
          })
          
          parent_groups <- lapply(parent_groups, function(group) {
            predictEBV_BGLR(group, tester_pop, bA, bD)
          })
          
        } else {
          
          parent_groups <- lapply(parent_groups, function(group) {
            selectOP(
              group,
              nInd          = nParents,
              nSeeds        = nOP,
              probSelf      = 0.2,
              pollenControl = FALSE,
              use           = "ebv"
            )
          })
        }
      }
      
      S1 <- lapply(parent_groups, function(group) {
        self(group, nProgeny = nSelf, keepParents = FALSE)
      })
      
      S2 <- lapply(S1, function(group) {
        self(group, nProgeny = nself1, keepParents = FALSE)
      })
      
      
      ## ---------------------------------------------------------------
      ## Prepare new training data
      ## ---------------------------------------------------------------
      
      S5_training <- replicate(n = 5, expr = { list() }, simplify = FALSE)
      
      for (group in seq_along(S4)) {
        for (line in 1:length(S4[[group]])) {
          S5_training[[group]][[line]] <- selectInd(
            S4[[group]][[line]],
            nInd = 1,
            use  = "rand"
          )
        }
      }
      
      for (group in seq_along(S5_training)) {
        S5_training[[group]] <- mergePops(unlist(S5_training[[group]]))
      }
      
      for (group in 1:length(S5_training)) {
        for (individual in 1:length(S5_training[[group]])) {
          training_ind <- S5_training[[group]][individual]
          train_syn    <- mergePops(list(training_ind, tester_pop))
          
          train_syn <- selectOP(
            pop           = train_syn,
            nInd          = 6,
            nSeeds        = 200,
            probSelf      = selfing,
            pollenControl = FALSE,
            use           = "rand"
          )
          
          train_syn <- setPheno(train_syn, varE = varE, reps = reptraining)
          S5_training[[group]]@pheno[individual] <- meanP(train_syn)
        }
      }
      
      for (group in 1:length(S5_training)) {
        for (individual in 1:length(S5_training[[group]])) {
          S4[[group]][[individual]]@pheno[, 1] <- S5_training[[group]]@pheno[individual]
        }
      }
      
      S5_training <- mergePops(S5_training)
      
      # Additional training population
      for (individual in 1:length(add_training)) {
        training_ind <- add_training[individual]
        train_syn    <- mergePops(list(training_ind, tester_pop))
        
        train_syn <- selectOP(
          pop           = train_syn,
          nInd          = 6,
          nSeeds        = 200,
          probSelf      = selfing,
          pollenControl = FALSE,
          use           = "rand"
        )
        
        train_syn <- setPheno(train_syn, varE = varE, reps = reptraining)
        add_training@pheno[individual] <- meanP(train_syn)
      }
      
      year_key <- paste0("Y", year)
      training_by_year[[year_key]] <- mergePops(list(S5_training, add_training))
      tester_by_year[[year_key]]   <- tester_pop
      
      
      ## ---------------------------------------------------------------
      ## Rebuild design matrices and refit model
      ## ---------------------------------------------------------------
      
      genoA_list <- vector("list", length(training_by_year))
      genoD_list <- vector("list", length(training_by_year))
      y_list     <- vector("list", length(training_by_year))
      
      for (k in seq_along(training_by_year)) {
        tr_pop  <- training_by_year[[k]]
        tst_pop <- tester_by_year[[k]]
        
        y_list[[k]] <- as.numeric(tr_pop@pheno)
        
        dmats <- build_design_for_year(tr_pop, tst_pop)
        genoA_list[[k]] <- dmats$A
        genoD_list[[k]] <- dmats$D
      }
      
      genoA <- do.call(rbind, genoA_list)
      genoD <- do.call(rbind, genoD_list)
      y     <- unlist(y_list)
      
      ETA <- list(
        list(X = genoA, model = "BRR"),
        list(X = genoD, model = "BRR")
      )
      
      fit <- BGLR(
        y       = y,
        ETA     = ETA,
        nIter   = 10000,
        burnIn  = 3000,
        verbose = FALSE
      )
      
      bA <- fit$ETA[[1]]$b
      bD <- fit$ETA[[2]]$b
      
      nNA_A <- sum(is.na(bA))
      nNA_D <- sum(is.na(bD))
      
      if (nNA_A > 0 || nNA_D > 0) {
        message(
          "Warning: ", nNA_A, " NA in bA and ", nNA_D,
          " NA in bD; setting them to 0"
        )
        bA[is.na(bA)] <- 0
        bD[is.na(bD)] <- 0
      }
      
      
      ## ---------------------------------------------------------------
      ## Benchmark accuracy
      ## ---------------------------------------------------------------
      
      S2 <- lapply(S2, function(group) {
        predictEBV_BGLR(group, tester_pop, bA, bD)
      })
      
      S2_dat <- mergePops(S2)
      S3_dat <- mergePops(S3)
      
      S2_dat <- setPhenoProgTest(S2_dat, tester_pop, nMatePerInd = 100, use = "gv", H2 = 1)
      S3_dat <- setPhenoProgTest(S3_dat, tester_pop, nMatePerInd = 100, use = "gv", H2 = 1)
      
      output$acc_progenytest_S2[year] <- cor(S2_dat@ebv, S2_dat@pheno, use = "complete.obs")
      output$acc_progenytest_S3[year] <- cor(S3_dat@ebv, S3_dat@pheno, use = "complete.obs")
      
      S2_dat <- setPhenoGCA(S2_dat, testers = tester_pop, H2 = 1.0)
      output$acc_GCA_S2[year] <- cor(S2_dat@ebv, S2_dat@pheno, use = "complete.obs")
      
      S3_dat <- setPhenoGCA(S3_dat, testers = tester_pop, H2 = 1.0)
      output$acc_GCA_S3[year] <- cor(S3_dat@ebv, S3_dat@pheno, use = "complete.obs")
      
      
      ## ---------------------------------------------------------------
      ## Create new elites / tester lines
      ## ---------------------------------------------------------------
      
      elite <- vector("list", length(S4))
      
      elite <- lapply(S4, function(sublist) {
        means    <- sapply(sublist, function(pop) mean(pop@pheno))
        best_idx <- which.max(means)
        sublist[[best_idx]]
      })
      
      for (group in seq_along(elite)) {
        elite[[group]] <- selectInd(elite[[group]], nInd = 200, use = "rand")
      }
      
      G1 <- c(elite[[2]], elite[[3]], elite[[4]], elite[[5]])
      G2 <- c(elite[[1]], elite[[3]], elite[[4]], elite[[5]])
      G3 <- c(elite[[4]], elite[[1]], elite[[2]], elite[[5]])
      G4 <- c(elite[[1]], elite[[2]], elite[[3]], elite[[5]])
      G5 <- c(elite[[1]], elite[[2]], elite[[3]], elite[[4]])
      
      synthetics <- list(G1, G2, G3, G4, G5)
      
      for (group in seq_along(elite)) {
        elite[[group]] <- selectInd(elite[[group]], nInd = 1, use = "rand")
      }
      
      S5twoyearsago <- S5
      S5lastyear    <- S5
      
      
      ## ---------------------------------------------------------------
      ## Output for year 1
      ## ---------------------------------------------------------------
      
      syn0_dat <- mergePops(unlist(syn0))
      syn1_dat <- mergePops(unlist(syn1))
      syn2_dat <- mergePops(unlist(syn2))
      
      output$meanSyn0[year]      <- meanG(syn0_dat)
      output$varGSyn0[year]      <- varG(syn0_dat)
      output$varDSyn0[year]      <- varD(syn0_dat)
      output$varASyn0[year]      <- varA(syn0_dat)
      output$meanSyn2[year]      <- meanG(syn2_dat)
      output$varGSyn2[year]      <- varG(syn2_dat)
      output$varASyn2[year]      <- varA(syn2_dat)
      output$varDSyn2[year]      <- varD(syn2_dat)
      output$varGS2[year]        <- varG(S2_dat)
      output$varAS2[year]        <- varA(S2_dat)
      output$varDS2[year]        <- varD(S2_dat)
      output$meanGvariety[year]  <- meanG(synthetic_var)
      output$varGvariety[year]   <- varG(synthetic_var)
      output$accSyn1[year]       <- cor(syn1_dat@gv, syn1_dat@pheno)
      output$accSyn2[year]       <- cor(syn2_dat@gv, syn2_dat@pheno)
      
      output$varGS2_1[year] <- varG(S2[[1]])
      output$varGS2_2[year] <- varG(S2[[2]])
      output$varGS2_3[year] <- varG(S2[[3]])
      output$varGS2_4[year] <- varG(S2[[4]])
      output$varGS2_5[year] <- varG(S2[[5]])
      
      output$varDS2_1[year] <- varD(S2[[1]])
      output$varDS2_2[year] <- varD(S2[[2]])
      output$varDS2_3[year] <- varD(S2[[3]])
      output$varDS2_4[year] <- varD(S2[[4]])
      output$varDS2_5[year] <- varD(S2[[5]])
      
    } else {
      
      ## ================================================================
      ## YEARS 2-20
      ## ================================================================
      
      ## ---------------------------------------------------------------
      ## Year 7: create SYN3 from selected SYN2 lines
      ## ---------------------------------------------------------------
      
      syn3 <- vector("list", length(Syn2_sel))
      
      for (group in seq_along(Syn2_sel)) {
        syn3[[group]] <- selectOP(
          pop           = Syn2_sel[[group]],
          nInd          = 500,
          nSeeds        = 4,
          probSelf      = selfing,
          pollenControl = FALSE,
          use           = "rand"
        )
      }
      
      synthetic_var <- mergePops(syn3)
      
      
      ## ---------------------------------------------------------------
      ## Year 6: selection from SYN1
      ## ---------------------------------------------------------------
      
      Syn1_mean <- vector("list", length(syn1))
      
      for (group in seq_along(syn1)) {
        group_means <- sapply(
          syn1[[group]],
          function(pop) mean(pop@pheno, na.rm = TRUE)
        )
        Syn1_mean[[group]] <- data.frame(
          index      = seq_along(group_means),
          mean_pheno = group_means
        )
      }
      
      Syn1_sel <- vector("list", length(syn1))
      
      for (group in seq_along(syn1)) {
        mean_df     <- Syn1_mean[[group]]
        top_indices <- order(mean_df$mean_pheno, decreasing = TRUE)[1:nsel_syn1]
        Syn1_sel[[group]] <- syn1[[group]][top_indices]
      }
      
      Syn1_sel <- unlist(Syn1_sel, recursive = FALSE)
      
      syn2 <- vector("list", length(Syn1_sel))
      
      for (group in seq_along(Syn1_sel)) {
        syn2[[group]] <- selectOP(
          pop           = Syn1_sel[[group]],
          nInd          = 500,
          nSeeds        = 4,
          probSelf      = selfing,
          pollenControl = FALSE,
          use           = "rand"
        )
        
        syn2[[group]] <- setPheno(syn2[[group]], varE = varE, reps = repsyn2)
      }
      
      Syn2_mean <- data.frame(
        index      = seq_along(syn2),
        mean_pheno = sapply(syn2, function(pop) mean(pop@pheno, na.rm = TRUE))
      )
      
      top_indices <- order(Syn2_mean$mean_pheno, decreasing = TRUE)[1:nsel_syn2]
      Syn2_sel    <- syn2[top_indices]
      
      
      ## ---------------------------------------------------------------
      ## Year 5: create SYN1 from SYN0
      ## ---------------------------------------------------------------
      
      syn1 <- vector("list", length(syn0))
      
      for (group in seq_along(syn0)) {
        syn0_group <- syn0[[group]]
        syn1_group <- vector("list", length(syn0_group))
        
        for (i in seq_along(syn0_group)) {
          pop_syn1 <- selectOP(
            pop      = syn0_group[[i]],
            nInd     = 2000,
            nSeeds   = 1,
            use      = "rand",
            probSelf = selfing
          )
          
          pop_syn1 <- setPheno(pop_syn1, varE = varE, reps = repsyn1)
          syn1_group[[i]] <- pop_syn1
        }
        
        syn1[[group]] <- syn1_group
      }
      
      
      ## ---------------------------------------------------------------
      ## Year 4: create SYN0
      ## ---------------------------------------------------------------
      
      S5_syn0 <- vector("list", length(S5))
      
      for (group in seq_along(S5)) {
        for (individual in 1:length(S5[[group]])) {
          S5_syn0[[group]][[individual]] <- selectInd(
            S5[[group]][[individual]],
            nInd = 200,
            use  = "rand"
          )
        }
      }
      
      syn0 <- replicate(n = 5, expr = { list() }, simplify = FALSE)
      
      for (group in seq_along(S5_syn0)) {
        for (i in 1:length(S5_syn0[[group]])) {
          open1 <- mergePops(list(S5_syn0[[group]][[i]], synthetics[[group]]))
          syn0[[group]][[i]] <- selectOP(
            open1,
            nInd          = 1000,
            nSeeds        = 2,
            probSelf      = selfing,
            pollenControl = FALSE,
            use           = "rand"
          )
        }
      }
      
      
      ## ---------------------------------------------------------------
      ## Year 3: create S3, select families, produce S4 and S5
      ## ---------------------------------------------------------------
      
      S3 <- lapply(S2s_sel, function(s2) {
        self(s2, nProgeny = nself2)
      })
      
      S3 <- lapply(S3, function(group) {
        predictEBV_BGLR(group, tester_pop, bA, bD)
      })
      
      S3s_sel <- lapply(S3, function(group) {
        selectFam(group, nFam = s3sel, use = "ebv", famType = "F")
      })
      
      S3_split_list <- vector("list", length(S3s_sel))
      
      for (group in seq_along(S3s_sel)) {
        group_pop  <- S3s_sel[[group]]
        mother_ids <- unique(group_pop@mother)
        S3_split   <- vector("list", length(mother_ids))
        
        for (i in seq_along(mother_ids)) {
          S3_split[[i]] <- subset(group_pop, group_pop@mother == mother_ids[i])
        }
        
        names(S3_split) <- as.character(mother_ids)
        S3_split_list[[group]] <- S3_split
      }
      
      S3s_sel <- S3_split_list
      
      S4 <- vector("list", length(S3s_sel))
      for (group in seq_along(S3s_sel)) {
        for (line in 1:length(S3s_sel[[group]])) {
          S4[[group]][[line]] <- self(S3s_sel[[group]][[line]], nProgeny = 200)
        }
      }
      
      S5 <- vector("list", length(S4))
      for (group in seq_along(S4)) {
        for (line in 1:length(S4[[group]])) {
          S5[[group]][[line]] <- self(S4[[group]][[line]], nProgeny = 1)
        }
      }
      
      
      ## ---------------------------------------------------------------
      ## Year 2: select S2
      ## ---------------------------------------------------------------
      
      S2s_sel <- lapply(S2, function(s2) {
        selectInd(s2, nInd = s2sel, use = "ebv")
      })
      
      
      ## ---------------------------------------------------------------
      ## Year 1: two cycles of population improvement
      ## ---------------------------------------------------------------
      
      for (cycle in 1:nCyclesPI) {
        
        if (cycle == 1) {
          
          parent_groups <- lapply(parent_groups, function(group) {
            predictEBV_BGLR(group, tester_pop, bA, bD)
          })
          
          Parents_output <- mergePops(parent_groups)
          
          output$meanGParents[year] <- meanG(Parents_output)
          output$varGParents[year]  <- varG(Parents_output)
          output$varAParents[year]  <- varA(Parents_output)
          output$varDParents[year]  <- varD(Parents_output)
          
          parent_groups <- lapply(parent_groups, function(group) {
            selectOP(
              group,
              nInd          = nParents,
              nSeeds        = nOP,
              probSelf      = 0.2,
              pollenControl = FALSE,
              use           = "ebv"
            )
          })
          
          parent_groups <- lapply(parent_groups, function(group) {
            predictEBV_BGLR(group, tester_pop, bA, bD)
          })
          
        } else {
          
          parent_groups <- lapply(parent_groups, function(group) {
            selectOP(
              group,
              nInd          = nParents,
              nSeeds        = nOP,
              probSelf      = 0.2,
              pollenControl = FALSE,
              use           = "ebv"
            )
          })
        }
      }
      
      S1 <- lapply(parent_groups, function(group) {
        self(group, nProgeny = nSelf, keepParents = FALSE)
      })
      
      S2 <- lapply(S1, function(group) {
        self(group, nProgeny = nself1, keepParents = FALSE)
      })
      
      S2 <- lapply(S2, function(group) {
        predictEBV_BGLR(group, tester_pop, bA, bD)
      })
      
      
      ## ---------------------------------------------------------------
      ## Prepare new training data
      ## ---------------------------------------------------------------
      
      S5_training <- replicate(n = 5, expr = { list() }, simplify = FALSE)
      
      for (group in seq_along(S5twoyearsago)) {
        for (line in 1:length(S5twoyearsago[[group]])) {
          S5_training[[group]][[line]] <- selectInd(
            S5twoyearsago[[group]][[line]],
            nInd = 1,
            use  = "rand"
          )
        }
      }
      
      for (group in seq_along(S5_training)) {
        S5_training[[group]] <- mergePops(unlist(S5_training[[group]]))
      }
      
      for (group in 1:length(S5_training)) {
        for (individual in 1:length(S5_training[[group]])) {
          training_ind <- S5_training[[group]][individual]
          train_syn    <- mergePops(list(training_ind, tester_pop))
          
          train_syn <- selectOP(
            pop           = train_syn,
            nInd          = 6,
            nSeeds        = 200,
            probSelf      = selfing,
            pollenControl = FALSE,
            use           = "rand"
          )
          
          train_syn <- setPheno(train_syn, varE = varE, reps = reptraining)
          S5_training[[group]]@pheno[individual] <- meanP(train_syn)
        }
      }
      
      for (group in 1:length(S5_training)) {
        for (individual in 1:length(S5_training[[group]])) {
          S5twoyearsago[[group]][[individual]]@pheno[, 1] <- S5_training[[group]]@pheno[individual]
        }
      }
      
      S5_training <- mergePops(S5_training)
      
      year_key <- paste0("Y", year)
      training_by_year[[year_key]] <- S5_training
      tester_by_year[[year_key]]   <- tester_pop
      
      
      ## ---------------------------------------------------------------
      ## Refit GS model
      ## ---------------------------------------------------------------
      
      genoA_list <- vector("list", length(training_by_year))
      genoD_list <- vector("list", length(training_by_year))
      y_list     <- vector("list", length(training_by_year))
      
      for (k in seq_along(training_by_year)) {
        tr_pop  <- training_by_year[[k]]
        tst_pop <- tester_by_year[[k]]
        
        y_list[[k]] <- as.numeric(tr_pop@pheno)
        
        dmats <- build_design_for_year(tr_pop, tst_pop)
        genoA_list[[k]] <- dmats$A
        genoD_list[[k]] <- dmats$D
      }
      
      genoA <- do.call(rbind, genoA_list)
      genoD <- do.call(rbind, genoD_list)
      y     <- unlist(y_list)
      
      ETA <- list(
        list(X = genoA, model = "BRR"),
        list(X = genoD, model = "BRR")
      )
      
      fit <- BGLR(
        y       = y,
        ETA     = ETA,
        nIter   = 10000,
        burnIn  = 3000,
        verbose = FALSE
      )
      
      bA <- fit$ETA[[1]]$b
      bD <- fit$ETA[[2]]$b
      
      nNA_A <- sum(is.na(bA))
      nNA_D <- sum(is.na(bD))
      
      if (nNA_A > 0 || nNA_D > 0) {
        message(
          "Warning: ", nNA_A, " NA in bA and ", nNA_D,
          " NA in bD; setting them to 0"
        )
        bA[is.na(bA)] <- 0
        bD[is.na(bD)] <- 0
      }
      
      
      ## ---------------------------------------------------------------
      ## Benchmark model accuracy
      ## ---------------------------------------------------------------
      
      S2_dat <- mergePops(S2)
      S3_dat <- mergePops(S3)
      
      S2_dat <- setPhenoProgTest(S2_dat, tester_pop, nMatePerInd = 100, use = "gv", H2 = 1)
      S3_dat <- setPhenoProgTest(S3_dat, tester_pop, nMatePerInd = 100, use = "gv", H2 = 1)
      
      output$acc_progenytest_S2[year] <- cor(S2_dat@ebv, S2_dat@pheno, use = "complete.obs")
      output$acc_progenytest_S3[year] <- cor(S3_dat@ebv, S3_dat@pheno, use = "complete.obs")
      
      S2_dat <- setPhenoGCA(S2_dat, testers = tester_pop, H2 = 1.0)
      output$acc_GCA_S2[year] <- cor(S2_dat@ebv, S2_dat@pheno, use = "complete.obs")
      
      S3_dat <- setPhenoGCA(S3_dat, testers = tester_pop, H2 = 1.0)
      output$acc_GCA_S3[year] <- cor(S3_dat@ebv, S3_dat@pheno, use = "complete.obs")
      
      
      ## ---------------------------------------------------------------
      ## Create new elite / tester lines
      ## ---------------------------------------------------------------
      
      elite <- vector("list", length(S5twoyearsago))
      
      elite <- lapply(S5twoyearsago, function(sublist) {
        means    <- sapply(sublist, function(pop) mean(pop@pheno))
        best_idx <- which.max(means)
        sublist[[best_idx]]
      })
      
      for (group in seq_along(elite)) {
        elite[[group]] <- selectInd(elite[[group]], nInd = 200, use = "rand")
      }
      
      G1 <- c(elite[[2]], elite[[3]], elite[[4]], elite[[5]])
      G2 <- c(elite[[1]], elite[[3]], elite[[4]], elite[[5]])
      G3 <- c(elite[[4]], elite[[1]], elite[[2]], elite[[5]])
      G4 <- c(elite[[1]], elite[[2]], elite[[3]], elite[[5]])
      G5 <- c(elite[[1]], elite[[2]], elite[[3]], elite[[4]])
      
      synthetics <- list(G1, G2, G3, G4, G5)
      
      for (group in seq_along(elite)) {
        elite[[group]] <- selectInd(elite[[group]], nInd = 1, use = "rand")
      }
      
      S5twoyearsago <- S5lastyear
      S5lastyear    <- S5
      
      
      ## ---------------------------------------------------------------
      ## Report results
      ## ---------------------------------------------------------------
      
      syn0_dat <- mergePops(unlist(syn0))
      syn1_dat <- mergePops(unlist(syn1))
      syn2_dat <- mergePops(unlist(syn2))
      
      output$meanSyn0[year]      <- meanG(syn0_dat)
      output$varGSyn0[year]      <- varG(syn0_dat)
      output$varDSyn0[year]      <- varD(syn0_dat)
      output$varASyn0[year]      <- varA(syn0_dat)
      output$meanSyn2[year]      <- meanG(syn2_dat)
      output$varGSyn2[year]      <- varG(syn2_dat)
      output$varASyn2[year]      <- varA(syn2_dat)
      output$varDSyn2[year]      <- varD(syn2_dat)
      output$varGS2[year]        <- varG(S2_dat)
      output$varAS2[year]        <- varA(S2_dat)
      output$varDS2[year]        <- varD(S2_dat)
      output$meanGvariety[year]  <- meanG(synthetic_var)
      output$varGvariety[year]   <- varG(synthetic_var)
      output$accSyn1[year]       <- cor(syn1_dat@gv, syn1_dat@pheno)
      output$accSyn2[year]       <- cor(syn2_dat@gv, syn2_dat@pheno)
      
      output$varGS2_1[year] <- varG(S2[[1]])
      output$varGS2_2[year] <- varG(S2[[2]])
      output$varGS2_3[year] <- varG(S2[[3]])
      output$varGS2_4[year] <- varG(S2[[4]])
      output$varGS2_5[year] <- varG(S2[[5]])
      
      output$varDS2_1[year] <- varD(S2[[1]])
      output$varDS2_2[year] <- varD(S2[[2]])
      output$varDS2_3[year] <- varD(S2[[3]])
      output$varDS2_4[year] <- varD(S2[[4]])
      output$varDS2_5[year] <- varD(S2[[5]])
    }
  }
  
  output
}


## ======================================================================
## RUN MODEL AND SAVE OUTPUT
## ======================================================================
out <- model.mse(rep_id)

saveRDS(
  out,
  file = file.path(
    workdir,
    sprintf("GS_synthetic_2part_BGLR_DD1_V2_rep%03d.rds", rep_id)
  )
)