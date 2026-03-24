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
    stop("Usage: Synthetic_breeding_phenotypic_DD1.R <rep_id> [workdir]")
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
  OMP_NUM_THREADS        = "1",
  MKL_NUM_THREADS        = "1",
  OPENBLAS_NUM_THREADS   = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS    = "1"
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
  
  # ---------------------------------------------------------------------------
  # Parameters
  # ---------------------------------------------------------------------------
  
  nParents = 21
  nCrosses = 160
  nElite   = 1
  
  nSelf       = 16   ## number of seeds per selfed F1
  nself2      = 1    ## number of seeds per selfed S1
  s1sel       = 160  ## number selected from S1
  s2sel       = 16   # number of lines selected from S2
  nsel_test1  = 10   ## Number of lines selected from testcross 1 per heterotic group
  nsel_test2  = 4    ## Number of lines selected from testcross 2 per heterotic group
  nsel_syn1   = 2    ## number of lines per heterotic group
  nsel_syn2   = 2    ## Number of varieties selected
  
  # Degree of selfing
  selfing = 0.5
  
  # Accuracy of trials
  # Replication at each stage
  varE     = 4     # error variance for field trials
  repS2    = 1 / 50
  reptest1 = 1
  reptest2 = 2
  repsyn1  = 6
  repsyn2  = 9
  
  # ---------------------------------------------------------------------------
  # Helper: pick top indices by mean phenotype
  # ---------------------------------------------------------------------------
  
  top_idx <- function(lst, k) {
    vals <- vapply(lst, function(p) mean(p@pheno), numeric(1))
    order(vals, decreasing = TRUE)[seq_len(min(k, length(vals)))]
  }
  
  # ---------------------------------------------------------------------------
  # Populations for filling pipeline
  # ---------------------------------------------------------------------------
  
  # Split heterotic pools to form initial parents
  nGroups <- 5
  selected_ids <- character()
  parent_groups <- list()
  
  # Starting population
  start_pop = c(Parents_17, Parents_18, Parents_19)
  
  # Iteratively select parents from the remaining population
  for (i in 1:nGroups) {
    remaining_ids <- setdiff(start_pop@id, selected_ids)
    selected <- selectInd(
      start_pop,
      nInd = nParents,
      use = "rand",
      candidates = remaining_ids
    )
    parent_groups[[i]] <- selected
    selected_ids <- c(selected_ids, selected@id)
  }
  
  # Random selection of individuals that will be used to evaluate selection accuracy
  benchmarking <- lapply(
    parent_groups,
    function(pop) {
      selectInd(pop, nInd = 3, use = "rand")
    }
  )
  
  # ---------------------------------------------------------------------------
  # Create testers
  # ---------------------------------------------------------------------------
  
  # Select elites from each group based on genetic value
  elite_groups <- lapply(parent_groups, function(group) {
    selectInd(group, nElite, use = "gv")
  })
  
  # Self elite 4 generations
  for (i in 1:4) {
    elite_groups <- lapply(elite_groups, function(group) {
      self(group, nProgeny = 5)
    })
  }
  
  G  <- length(elite_groups)
  k1 <- 2  # number of other groups to use for Tester1 (random)
  
  # 1) Pick representative/random individual(s) from each group
  elite_groups <- lapply(elite_groups, function(g) {
    selectInd(g, nInd = 1, use = "rand")
  })
  
  # 2) Build testers per focal group
  testers1 <- vector("list", G)  # random subset of others
  testers2 <- vector("list", G)  # all other groups
  
  for (i in seq_len(G)) {
    others_idx <- setdiff(seq_len(G), i)
    
    # Testcross 1: random k1 groups from the others
    pick <- sample(others_idx, size = k1, replace = FALSE)
    testers1[[i]] <- do.call(c, elite_groups[pick])
    
    # Testcross 2: all the other groups
    testers2[[i]] <- do.call(c, elite_groups[others_idx])
  }
  
  # ---------------------------------------------------------------------------
  # Start running breeding cycle
  # ---------------------------------------------------------------------------
  
  nYears = 20
  output = data.frame(year = 1:nYears)
  
  for (year in 1:nYears) {
    
    cat("Cycle year:", year, "of 20\n")
    
    if (year == 1) {
      
      # =======================================================================
      # Creation of synthetic varieties
      # Year 7-10 run within one year (last few years of burnin) to enable
      # creation of synthetics
      # =======================================================================
      
      # First synthetic varieties created with parents selected in year 17, 18
      # and 19 of burnin. Filling last steps of pipeline with synthetics derived
      # from these lines
      
      # Starting population
      firstsyn = c(Parents_17, Parents_18, Parents_19)
      
      firstsyn2 <- list()
      selected_ids <- character()
      
      # Iteratively select parents from the remaining population
      for (i in 1:nGroups) {
        remaining_ids <- setdiff(firstsyn@id, selected_ids)
        selected <- selectInd(
          firstsyn,
          nInd = nsel_test2,
          use = "rand",
          candidates = remaining_ids
        )
        firstsyn2[[i]] <- selected
        selected_ids <- c(selected_ids, selected@id)
      }
      
      # -----------------------------------------------------------------------
      # Create SYN0 by crossing all the synthetic components
      # -----------------------------------------------------------------------
      
      syn0 = replicate(n = 5, expr = { list() }, simplify = FALSE)
      
      for (group in seq_along(firstsyn2)) {
        for (i in 1:length(firstsyn2[[group]])) {
          open1 = mergePops(list(firstsyn2[[group]][i], testers2[[group]]))
          syn0[[group]][i] = selectOP(
            open1,
            nInd = 4,
            nSeeds = 500,  # unrealistic but I think best way to simulate this open pollination
            probSelf = selfing, ###
            pollenControl = FALSE,
            use = "rand"
          )
        }
      }
      
      # -----------------------------------------------------------------------
      # Random open pollination of Syn0 to produce Syn1
      # -----------------------------------------------------------------------
      
      syn1 = replicate(n = 5, expr = { list() }, simplify = FALSE)  # Keep lines separate
      
      for (group in seq_along(syn1)) {
        for (i in 1:length(syn0[[group]])) {
          syn1[[group]][[i]] = selectOP(
            syn0[[group]][[i]],
            nInd = 2000,
            nSeeds = 1,
            use = "rand",
            probSelf = selfing
          )
          syn1[[group]][[i]] = setPheno(
            syn1[[group]][[i]],
            varE = varE,
            reps = repsyn1
          )
        }
      }
      
      # Random selection from syn1 as no phenotyping prior to breeding cycles
      Syn1_sel <- lapply(syn1, function(group) {
        n <- length(group)
        keep <- sample(
          seq_len(n),
          size = min(nsel_syn1, n),
          replace = FALSE
        )
        group[keep]
      })
      
      Syn1_sel = unlist(Syn1_sel)
      
      # -----------------------------------------------------------------------
      # Open pollination of Syn-2 to produce Syn-3 and phenotype
      # -----------------------------------------------------------------------
      
      syn2 <- lapply(Syn1_sel, function(line) {
        pop <- selectOP(
          line,
          nInd = 500,
          nSeeds = 4,
          probSelf = selfing,
          pollenControl = FALSE,
          use = "rand"
        )
        setPheno(pop, varE = varE, reps = repsyn2)
      })
      
      # Selection from Syn2: random selection as no phenotyping prior to breeding cycle
      Syn2_sel <- sample(syn2, size = nsel_syn2, replace = FALSE)
      
      # Creation of Syn3
      syn3 <- vector("list", length(Syn2_sel))  # Initialize list for SYN2
      
      for (group in seq_along(Syn2_sel)) {
        syn3[[group]] <- selectOP(
          pop = Syn2_sel[[group]],
          nInd = 200,         # Number of seeds per line for multiplication
          nSeeds = 10,
          probSelf = selfing, # Mostly selfed, some open pollination
          pollenControl = FALSE,
          use = "rand"
        )
      }
      
      synthetic_var = mergePops(syn3)
      
      # -----------------------------------------------------------------------
      # Year 6
      # -----------------------------------------------------------------------
      
      # Random selection of lines from F5
      nGroups <- 5
      
      # Randomly sample nsel_test2 * nGroups indices (without replacement)
      idx <- sample(length(F5), size = nsel_test2 * nGroups, replace = FALSE)
      
      # Split those sampled indices into 5 groups of nsel_test2 each
      split_idx <- split(idx, rep(1:nGroups, each = nsel_test2))
      
      # Extract the objects into 5 separate lists
      Testcross2_sel <- lapply(split_idx, function(i) F5[i])
      
      # Open-pollination of selected lines from F5 to multiply
      S5 = replicate(n = 5, expr = { list() }, simplify = FALSE)  # Keep lines separate
      
      for (group in seq_along(Testcross2_sel)) {
        for (i in 1:length(Testcross2_sel[[group]])) {
          S5[[group]][[i]] = selectOP(
            Testcross2_sel[[group]][[i]],
            nInd = 1,
            nSeeds = 2000,
            probSelf = selfing, ###
            pollenControl = FALSE,
            use = "rand"
          )
        }
      }
      
      # -----------------------------------------------------------------------
      # Year 5
      # Crossing for testcross 2 performed
      # -----------------------------------------------------------------------
      
      # Testcross-2: composite plot from 5 parents per F4 line
      nMixPerLine <- 5  # number of F4 individuals contributing to the composite plot
      
      # Split F4_split into heterotic groups
      split_4 <- split(F4_split, rep(1:nGroups, each = (length(F4_split) / nGroups)))
      
      testcross2 <- lapply(seq_along(split_4), function(g) {
        lapply(seq_along(split_4[[g]]), function(i) {
          line_pop <- split_4[[g]][[i]]
          
          # sample up to nMixPerLine individuals from the S4 family
          k <- min(nMixPerLine, nInd(line_pop))
          mix_parents <- selectInd(line_pop, nInd = k, use = "rand")
          
          # 1) Genetic expectation only (no residual): GCA vs testers for the k parents
          gca_true <- setPhenoGCA(
            mix_parents,
            testers = testers2[[g]],
            H2 = 1.0  # pure genetic expectation; no plot error
          )
          
          # equal-contribution composite (use weighted mean if needed)
          mu_g <- mean(gca_true@pheno, na.rm = TRUE)
          
          # 2) Single-plot residual for the composite plot
          eps <- rnorm(1, mean = 0, sd = sqrt(varE / reptest2))
          
          fam_pheno <- mu_g + eps
          
          # assign the composite-plot phenotype to every individual in the family
          line_pop@pheno[] <- fam_pheno
          
          line_pop
        })
      })
      
      # -----------------------------------------------------------------------
      # Year 4
      # -----------------------------------------------------------------------
      
      # select top based on per se performance of lines in burnin
      Testcross1_sel <- selectFam(
        mergePops(F3),
        nFam = nsel_test1 * nGroups,
        trait = 1,
        use = "pheno",
        famType = "F"
      )
      
      motherSval_F3 = unique(Testcross1_sel@mother)
      Testcross1_sel_split = vector("list", length(motherSval_F3))
      
      for (i in 1:length(motherSval_F3)) {
        motherpop = subset(
          Testcross1_sel,
          Testcross1_sel@mother == paste(motherSval_F3[i])
        )
        Testcross1_sel_split[[i]] = selectInd(motherpop, nInd = 40)
      }
      
      # Split Testcross1_sel_split into heterotic groups
      Testcross1_sel <- split(
        Testcross1_sel_split,
        rep(1:nGroups, each = (length(Testcross1_sel_split) / nGroups))
      )
      
      S4 <- lapply(seq_along(Testcross1_sel), function(group) {
        lapply(seq_along(Testcross1_sel[[group]]), function(i) {
          
          # Step 1: Generate full OP progeny
          S4_full <- selectOP(
            Testcross1_sel[[group]][[i]],
            nInd = 40,
            nSeeds = 50,
            probSelf = 0.95,
            pollenControl = FALSE,
            use = "rand"
          )
          
          return(S4_full)
        })
      })
      
      # -----------------------------------------------------------------------
      # Year 3
      # -----------------------------------------------------------------------
      
      # Split F2 into heterotic groups
      F2_split <- split(F2, rep(1:nGroups, each = (length(F2) / nGroups)))
      
      # Step 5: Prepare empty lists for S3 populations (one per mother)
      S3s <- lapply(seq_along(F2_split), function(group) {
        lapply(seq_along(F2_split[[group]]), function(i) {
          
          # Step 1: Generate full OP progeny (40 seeds)
          S3_full <- selectOP(
            pop = F2_split[[group]][[i]],
            nInd = 40,
            nSeeds = 10,
            probSelf = 0.95,
            pollenControl = FALSE,
            use = "rand"
          )
          
          return(S3_full)
        })
      })
      
      # Testcross-1: composite plot from 5 parents per S4 line
      nMixPerLine <- 5  # number of S4 individuals contributing to the composite plot
      
      testcross1 <- lapply(seq_along(S3s), function(g) {
        lapply(seq_along(S3s[[g]]), function(i) {
          line_pop <- S3s[[g]][[i]]
          
          # sample up to nMixPerLine individuals from the S4 family
          k <- min(nMixPerLine, nInd(line_pop))
          mix_parents <- selectInd(line_pop, nInd = k, use = "rand")
          
          # 1) Genetic expectation only (no residual): GCA vs testers for the k parents
          gca_true <- setPhenoGCA(
            mix_parents,
            testers = testers1[[g]],
            H2 = 1.0  # pure genetic expectation; no plot error
          )
          
          # equal-contribution composite (use weighted mean if needed)
          mu_g <- mean(gca_true@pheno, na.rm = TRUE)
          
          # 2) Single-plot residual for the composite plot
          eps <- rnorm(1, mean = 0, sd = sqrt(varE / reptest1))
          
          fam_pheno <- mu_g + eps
          
          # assign the composite-plot phenotype to every individual in the family
          line_pop@pheno[] <- fam_pheno
          
          line_pop
        })
      })
      
      # -----------------------------------------------------------------------
      # Year 2
      # -----------------------------------------------------------------------
      
      # Selfing of S1 to produce S2
      
      # Split F1 into heterotic groups
      F1_split <- split(F1, rep(1:nGroups, each = (length(F1) / nGroups)))
      
      for (i in seq_along(F1_split)) {
        F1_split[[i]] = mergePops(F1_split[[i]])
      }
      
      # Step 1: Self S1-selected lines to produce S2
      S2s <- lapply(F1_split, function(sel) {
        self(sel, nProgeny = nself2 * 16, keepParents = FALSE)
      })
      
      # Phenotyping of S2 lines in Chile
      S2s <- lapply(S2s, function(sel) {
        setPheno(sel, varE = varE, reps = repS2)
      })
      
      S2s_sel <- lapply(S2s, function(sel) {
        selectFam(sel, nFam = s2sel, use = "pheno")
      })
      
      # -----------------------------------------------------------------------
      # Year 1
      # -----------------------------------------------------------------------
      
      # Step 1: Cross each group
      F1s <- lapply(parent_groups, function(parents) {
        randCross(parents, nCrosses = nCrosses, nProgeny = 1)
      })
      
      # Step 2: Self the F1s
      S1s <- lapply(F1s, function(f1) {
        self(f1, nProgeny = nSelf, keepParents = FALSE)
      })
      
      # -----------------------------------------------------------------------
      # Accuracy testcrosses
      # Do accuracy before testers are updated
      # -----------------------------------------------------------------------
      
      # 1. GCA Benchmark
      s2_dat = S2s
      
      for (group in seq_along(S2s)) {
        s2_dat[[group]] <- setPhenoGCA(
          s2_dat[[group]],
          testers = mergePops(elite_groups),
          H2 = 1.0
        )
      }
      
      s3_dat = testcross1
      
      for (group in seq_along(testcross1)) {
        for (line in 1:length(testcross1[[group]])) {
          s3_dat[[group]][[line]] <- setPhenoGCA(
            s3_dat[[group]][[line]],
            testers = mergePops(elite_groups),
            H2 = 1.0
          )
        }
      }
      
      s4_dat = testcross2
      
      for (group in seq_along(testcross2)) {
        for (line in 1:length(testcross2[[group]])) {
          s4_dat[[group]][[line]] <- setPhenoGCA(
            s4_dat[[group]][[line]],
            testers = mergePops(elite_groups),
            H2 = 1.0
          )
        }
      }
      
      # Calculate mean pheno GCA for each object in the sublist
      mean_S3_pheno <- lapply(testcross1, function(sublist) {
        sapply(sublist, function(pop) mean(pop@pheno))
      })
      
      mean_S4_pheno <- lapply(testcross2, function(sublist) {
        sapply(sublist, function(pop) mean(pop@pheno))
      })
      
      # Calculate mean "True" GCA for each object in the sublist
      mean_GCA_S3 <- lapply(s3_dat, function(sublist) {
        sapply(sublist, function(pop) mean(pop@pheno))
      })
      
      mean_GCA_S4 <- lapply(s4_dat, function(sublist) {
        sapply(sublist, function(pop) mean(pop@pheno))
      })
      
      # Unlist
      s2_pheno    <- mergePops(S2s)
      s3_pheno    <- unlist(mean_S3_pheno, recursive = FALSE)
      s4_pheno    <- unlist(mean_S4_pheno, recursive = FALSE)
      mean_GCA_S2 <- mergePops(s2_dat)
      mean_GCA_S3 <- unlist(mean_GCA_S3, recursive = FALSE)
      mean_GCA_S4 <- unlist(mean_GCA_S4, recursive = FALSE)
      
      output$accGCA_S2[year]         <- cor(s2_pheno@pheno, mean_GCA_S2@pheno)
      output$accGCA_Testcross1[year] <- cor(s3_pheno, mean_GCA_S3)
      output$accGCA_Testcross2[year] <- cor(s4_pheno, mean_GCA_S4)
      
      # 2. Synthetic performance benchmark – progeny test
      for (group in seq_along(S2s)) {
        s2_dat[[group]] <- setPhenoProgTest(
          s2_dat[[group]],
          mergePops(elite_groups),
          nMatePerInd = 100,
          H2 = 1
        )
      }
      
      for (group in seq_along(S3s)) {
        for (line in 1:length(S3s[[1]])) {
          s3_dat[[group]][[line]] <- setPhenoProgTest(
            s3_dat[[group]][[line]],
            mergePops(elite_groups),
            nMatePerInd = 100,
            H2 = 1
          )
        }
      }
      
      for (group in seq_along(S4)) {
        for (line in 1:length(S4[[group]])) {
          s4_dat[[group]][[line]] <- setPhenoProgTest(
            s4_dat[[group]][[line]],
            mergePops(elite_groups),
            nMatePerInd = 100,
            H2 = 1
          )
        }
      }
      
      mean_prog_S3 <- lapply(s3_dat, function(sublist) {
        sapply(sublist, function(pop) mean(pop@pheno))
      })
      
      mean_prog_S4 <- lapply(s4_dat, function(sublist) {
        sapply(sublist, function(pop) mean(pop@pheno))
      })
      
      mean_prog_S2 <- mergePops(s2_dat)
      mean_prog_S3 <- unlist(mean_prog_S3, recursive = FALSE)
      mean_prog_S4 <- unlist(mean_prog_S4, recursive = FALSE)
      
      output$accprogeny_S2[year]         <- cor(s2_pheno@pheno, mean_prog_S2@pheno)
      output$accprogeny_Testcross1[year] <- cor(s3_pheno, mean_prog_S3)
      output$accprogeny_Testcross2[year] <- cor(s4_pheno, mean_prog_S4)
      
      # -----------------------------------------------------------------------
      # Prepare testers
      # -----------------------------------------------------------------------
      
      # Create list of testers/lines for synthetic variety from random selection
      # from each heterotic group as no previous phenotype data
      sel_elite <- lapply(S5, function(group) {
        if (length(group) == 0) return(group)
        keep <- sample.int(length(group), size = min(nElite, length(group)))
        group[keep]
      })
      
      elite = list()
      sel_elite = unlist(sel_elite)
      
      for (group in seq_along(sel_elite)) {
        elite[[group]] = selectOP(
          sel_elite[[group]],
          nInd = 200,
          nSeeds = 10,
          probSelf = selfing, ###
          pollenControl = FALSE,
          use = "rand"
        )
      }
      
      # Separate heterotic groups
      synthetics <- list(
        list(elite[[2]], elite[[3]], elite[[4]], elite[[5]]),
        list(elite[[1]], elite[[3]], elite[[4]], elite[[5]]),
        list(elite[[1]], elite[[2]], elite[[4]], elite[[5]]),
        list(elite[[1]], elite[[2]], elite[[3]], elite[[5]]),
        list(elite[[1]], elite[[2]], elite[[3]], elite[[4]])
      )
      
      testers2 <- lapply(seq_along(synthetics), function(group) {
        selected <- lapply(synthetics[[group]], function(pop) {
          selectInd(pop, nInd = 1, use = "rand")
        })
        mergePops(selected)
      })
      
      # Testers for testcross 1
      testers1 <- lapply(testers2, selectInd, nInd = 2, use = "rand")
      
      # Tester for benchmarking accuracy of selection
      tester_pop <- list(elite[[1]], elite[[2]], elite[[3]], elite[[4]], elite[[5]])
      
      for (group in seq_along(tester_pop)) {
        tester_pop[[group]] = selectInd(tester_pop[[group]], nInd = 1, use = "rand")
      }
      
      tester_pop = mergePops(tester_pop)
      
      # Create parents from previous years S4
      Parents_breed <- lapply(seq_along(S4), function(group) {
        lapply(S4[[group]], selectInd, nInd = 1, use = "rand")
      })
      Parents_breed <- lapply(Parents_breed, mergePops)
      
      S5twoyearsago = S5
      S5lastyear    = S5
      
      # -----------------------------------------------------------------------
      # Report results
      # -----------------------------------------------------------------------
      
      Parents_output = mergePops(parent_groups)
      
      output$meanSyn0[year]     = meanG(mergePops(unlist(syn0)))
      output$varSyn0[year]      = varG(mergePops(unlist(syn0)))
      output$meanGvariety[year] = meanG(synthetic_var)
      output$varGvariety[year]  = varG(synthetic_var)
      output$meanGParents[year] = meanG(Parents_output)
      output$meanEBVParents[year] = meanEBV(Parents_output)
      output$varGParents[year]  = varG(Parents_output)
      output$varAParents[year]  = varA(Parents_output)
      output$varDParents[year]  = varD(Parents_output)
      
      syn1_dat = mergePops(unlist(syn1))
      syn2_dat = mergePops(unlist(syn2))
      
      output$accSyn1[year] = cor(syn1_dat@gv, syn1_dat@pheno)
      output$accSyn2[year] = cor(syn2_dat@gv, syn2_dat@pheno)
      
      output$varGS2[year] = varG(s2_pheno)
      output$varDS2[year] = varD(s2_pheno)
      
      output$varGS2_1[year] = varG(S2s[[1]])
      output$varGS2_2[year] = varG(S2s[[2]])
      output$varGS2_3[year] = varG(S2s[[3]])
      output$varGS2_4[year] = varG(S2s[[4]])
      output$varGS2_5[year] = varG(S2s[[5]])
      
      output$varDS2_1[year] = varD(S2s[[1]])
      output$varDS2_2[year] = varD(S2s[[2]])
      output$varDS2_3[year] = varD(S2s[[3]])
      output$varDS2_4[year] = varD(S2s[[4]])
      output$varDS2_5[year] = varD(S2s[[5]])
      
    } else {
      
      # =======================================================================
      # Subsequent breeding cycles
      # =======================================================================
      
      # -----------------------------------------------------------------------
      # Year 10: Variety release
      # Open pollination of selected SYN2 lines to create SYN3
      # -----------------------------------------------------------------------
      
      syn3 <- vector("list", length(Syn2_sel))  # Initialize list for SYN3
      
      for (group in seq_along(Syn2_sel)) {
        syn3[[group]] <- selectOP(
          pop = Syn2_sel[[group]],
          nInd = 500,         # Number of seeds per line for multiplication
          nSeeds = 4,
          probSelf = selfing, # Mostly selfed, some open pollination
          pollenControl = FALSE,
          use = "rand"
        )
      }
      
      synthetic_var = mergePops(syn3)
      
      # -----------------------------------------------------------------------
      # Year 9
      # Open pollination of Syn-1 to produce Syn-2
      # -----------------------------------------------------------------------
      
      Syn1_sel = unlist(Syn1_sel)
      syn2 = list()
      
      for (group in seq_along(Syn1_sel)) {
        syn2[[group]] = selectOP(
          Syn1_sel[[group]],
          nInd = 500,
          nSeeds = 4,
          probSelf = selfing, ###
          pollenControl = FALSE,
          use = "rand"
        )
        syn2[[group]] = setPheno(syn2[[group]], varE = varE, reps = repsyn2)
      }
      
      # Selection from Syn2
      syn2_means <- vapply(syn2, function(pop) mean(pop@pheno, na.rm = TRUE), numeric(1))
      keep_idx <- head(
        order(syn2_means, decreasing = TRUE, na.last = NA),
        n = min(nsel_syn2, length(syn2_means))
      )
      Syn2_sel <- syn2[keep_idx]
      
      # -----------------------------------------------------------------------
      # Year 8
      # Random open pollination of Syn0 to produce Syn1
      # -----------------------------------------------------------------------
      
      syn1 = replicate(n = 5, expr = { list() }, simplify = FALSE)  # Keep lines separate
      
      for (group in seq_along(syn0)) {
        for (i in 1:length(syn0[[group]])) {
          syn1[[group]][[i]] = selectOP(
            syn0[[group]][[i]],
            nInd = 2000,
            nSeeds = 1,
            use = "rand",
            probSelf = selfing
          )
          syn1[[group]][[i]] = setPheno(
            syn1[[group]][[i]],
            varE = varE,
            reps = repsyn1
          )
        }
      }
      
      # Selection from Syn1
      Syn1_sel <- lapply(syn1, function(group) {
        means <- vapply(group, function(pop) mean(pop@pheno), numeric(1))
        keep  <- order(means, decreasing = TRUE)[seq_len(min(nsel_syn1, length(means)))]
        group[keep]
      })
      
      # -----------------------------------------------------------------------
      # Year 7
      # Create SYN0 by crossing all the synthetic components
      # -----------------------------------------------------------------------
      
      syn0 = replicate(n = 5, expr = { list() }, simplify = FALSE)
      
      for (group in seq_along(S5)) {
        for (i in 1:length(S5[[group]])) {
          open1 = mergePops(list(S5[[group]][[i]], mergePops(synthetics[[group]])))
          syn0[[group]][[i]] = selectOP(
            open1,
            nInd = 2000,
            nSeeds = 1,
            probSelf = selfing, ###
            pollenControl = FALSE,
            use = "rand"
          )
        }
      }
      
      open2 <- lapply(seq_along(S5), function(group) {
        pops <- lapply(seq_along(S5[[group]]), function(i) {
          mergePops(list(S5[[group]][[i]], mergePops(synthetics[[group]])))
        })
        mergePops(pops)
      })
      
      open2 <- mergePops(open2)
      output$meanGS5[year] <- meanG(open2)
      
      # -----------------------------------------------------------------------
      # Year 6
      # Select top nsel_test2 S4 objects (by mean @pheno) within each group of S4
      # -----------------------------------------------------------------------
      
      Testcross2_sel <- lapply(testcross2, function(group) {
        means <- vapply(group, function(pop) mean(pop@pheno), numeric(1))
        keep  <- order(means, decreasing = TRUE)[seq_len(min(nsel_test2, length(means)))]
        group[keep]
      })
      
      # Open-pollination of selected lines from S4
      S5 = replicate(n = 5, expr = { list() }, simplify = FALSE)  # Keep lines separate
      
      for (group in seq_along(Testcross2_sel)) {
        for (i in 1:length(Testcross2_sel[[group]])) {
          S5[[group]][[i]] = selectOP(
            Testcross2_sel[[group]][[i]],
            nInd = 200,
            nSeeds = 10,
            probSelf = 0.95, ###
            pollenControl = FALSE,
            use = "rand"
          )
        }
      }
      
      # -----------------------------------------------------------------------
      # Year 5
      # Crossing for testcross 2 performed
      # -----------------------------------------------------------------------
      
      nMixPerLine <- 5  # number of S4 individuals contributing to the composite plot
      
      testcross2 <- lapply(seq_along(S4), function(g) {
        lapply(seq_along(S4[[g]]), function(i) {
          line_pop <- S4[[g]][[i]]
          
          # sample up to nMixPerLine individuals from the S4 family
          k <- min(nMixPerLine, nInd(line_pop))
          mix_parents <- selectInd(line_pop, nInd = k, use = "rand")
          
          # 1) Genetic expectation only (no residual): GCA vs testers for the k parents
          gca_true <- setPhenoGCA(
            mix_parents,
            testers = testers2[[g]],
            H2 = 1.0  # pure genetic expectation; no plot error
          )
          
          # equal-contribution composite (use weighted mean if needed)
          mu_g <- mean(gca_true@pheno, na.rm = TRUE)
          
          # 2) Single-plot residual for the composite plot
          eps <- rnorm(1, mean = 0, sd = sqrt(varE / reptest2))
          
          fam_pheno <- mu_g + eps
          
          # assign the composite-plot phenotype to every individual in the family
          line_pop@pheno[] <- fam_pheno
          
          line_pop
        })
      })
      
      # -----------------------------------------------------------------------
      # Year 4
      # Selection from testcross 1
      # -----------------------------------------------------------------------
      
      Testcross1_sel <- lapply(testcross1, function(group) {
        means <- vapply(group, function(pop) mean(pop@pheno), numeric(1))
        keep  <- order(means, decreasing = TRUE)[seq_len(min(nsel_test1, length(means)))]
        group[keep]
      })
      
      S4 <- lapply(seq_along(Testcross1_sel), function(group) {
        lapply(seq_along(Testcross1_sel[[group]]), function(i) {
          
          # Step 1: Generate full OP progeny
          S4_full <- selectOP(
            Testcross1_sel[[group]][[i]],
            nInd = 160,
            nSeeds = 10,
            probSelf = 0.95,
            pollenControl = FALSE,
            use = "rand"
          )
          
          return(S4_full)
        })
      })
      
      # -----------------------------------------------------------------------
      # Year 3
      # Selection from testcross 1
      # Production of open pollinated seeds from field in Chile
      # -----------------------------------------------------------------------
      
      # Step 4: Identify unique mothers for open pollination
      mothers <- lapply(S2s_sel, function(sel) {
        unique(sel@mother)
      })
      
      # Step 5: Prepare empty lists for S3 populations (one per mother)
      S3s <- lapply(seq_along(mothers), function(group) {
        lapply(seq_along(mothers[[group]]), function(i) {
          
          mama   <- mothers[[group]]
          S2     <- S2s_sel[[group]]
          tester <- testers1[[group]]
          
          motherpop <- subset(S2, S2@mother == mama[i])
          
          # Step 1: Generate full OP progeny (40 seeds)
          S3_full <- selectOP(
            pop = motherpop,
            nInd = 16,
            nSeeds = 10,
            probSelf = 0.95,
            pollenControl = FALSE,
            use = "rand"
          )
          
          return(S3_full)
        })
      })
      
      # Testcross-1: composite plot from 5 parents per S4 line
      nMixPerLine <- 5  # number of S4 individuals contributing to the composite plot
      
      testcross1 <- lapply(seq_along(S3s), function(g) {
        lapply(seq_along(S3s[[g]]), function(i) {
          line_pop <- S3s[[g]][[i]]
          
          # sample up to nMixPerLine individuals from the S4 family
          k <- min(nMixPerLine, nInd(line_pop))
          mix_parents <- selectInd(line_pop, nInd = k, use = "rand")
          
          # 1) Genetic expectation only (no residual): GCA vs testers for the k parents
          gca_true <- setPhenoGCA(
            mix_parents,
            testers = testers1[[g]],
            H2 = 1.0  # pure genetic expectation; no plot error
          )
          
          # equal-contribution composite (use weighted mean if needed)
          mu_g <- mean(gca_true@pheno, na.rm = TRUE)
          
          # 2) Single-plot residual for the composite plot
          eps <- rnorm(1, mean = 0, sd = sqrt(varE / reptest1))
          
          fam_pheno <- mu_g + eps
          
          # assign the composite-plot phenotype to every individual in the family
          line_pop@pheno[] <- fam_pheno
          
          line_pop
        })
      })
      
      # -----------------------------------------------------------------------
      # Year 2
      # Selfing of S1 to produce S2
      # -----------------------------------------------------------------------
      
      S2s <- lapply(S1s, function(sel) {
        self(sel, nProgeny = nself2, keepParents = TRUE)
      })
      
      # Phenotyping of S2 lines in Chile
      S2s <- lapply(S2s, function(sel) {
        setPheno(sel, varE = varE, reps = repS2)
      })
      
      S2s_sel <- lapply(S2s, function(sel) {
        selectFam(sel, nFam = s2sel, use = "pheno")
      })
      
      # -----------------------------------------------------------------------
      # Year 1
      # -----------------------------------------------------------------------
      
      F1s <- lapply(Parents_breed, function(parents) {
        randCross(parents, nCrosses = nCrosses, nProgeny = 1)
      })
      
      S1s <- lapply(F1s, function(f1) {
        self(f1, nProgeny = nSelf, keepParents = FALSE)
      })
      
      # -----------------------------------------------------------------------
      # Create parents for next cycle
      # -----------------------------------------------------------------------
      
      Parents_breed = replicate(n = 5, expr = { list() }, simplify = FALSE)
      
      for (group in seq_along(S4)) {
        for (line in 1:length(S4[[group]])) {
          Parents_breed[[group]][[line]] = selectInd(
            S4[[group]][[line]],
            nInd = 1,
            use = "rand"
          )
        }
      }
      
      # Merge Parent populations
      for (group in seq_along(Parents_breed)) {
        Parents_breed[[group]] = mergePops(Parents_breed[[group]])
      }
      
      Parents_output = mergePops(Parents_breed)
      
      # -----------------------------------------------------------------------
      # Accuracy testcrosses
      # Do accuracy tests before tester pop is updated
      # -----------------------------------------------------------------------
      
      # 1. GCA Benchmark
      s2_dat = S2s
      
      for (group in seq_along(s2_dat)) {
        s2_dat[[group]] <- setPhenoGCA(
          s2_dat[[group]],
          testers = tester_pop,
          H2 = 1.0
        )
      }
      
      s3_dat = testcross1
      
      for (group in seq_along(testcross1)) {
        for (line in 1:length(testcross1[[group]])) {
          s3_dat[[group]][[line]] <- setPhenoGCA(
            s3_dat[[group]][[line]],
            testers = tester_pop,
            H2 = 1.0
          )
        }
      }
      
      s4_dat = testcross2
      
      for (group in seq_along(testcross2)) {
        for (line in 1:length(testcross2[[group]])) {
          s4_dat[[group]][[line]] <- setPhenoGCA(
            s4_dat[[group]][[line]],
            testers = tester_pop,
            H2 = 1.0
          )
        }
      }
      
      # Calculate mean pheno GCA for each object in the sublist
      mean_S3_pheno <- lapply(testcross1, function(sublist) {
        sapply(sublist, function(pop) mean(pop@pheno))
      })
      
      mean_S4_pheno <- lapply(testcross2, function(sublist) {
        sapply(sublist, function(pop) mean(pop@pheno))
      })
      
      # Calculate mean "True" GCA for each object in the sublist
      mean_GCA_S3 <- lapply(s3_dat, function(sublist) {
        sapply(sublist, function(pop) mean(pop@pheno))
      })
      
      mean_GCA_S4 <- lapply(s4_dat, function(sublist) {
        sapply(sublist, function(pop) mean(pop@pheno))
      })
      
      # Unlist
      s2_pheno    <- mergePops(S2s)
      s3_pheno    <- unlist(mean_S3_pheno, recursive = FALSE)
      s4_pheno    <- unlist(mean_S4_pheno, recursive = FALSE)
      mean_GCA_S2 <- mergePops(s2_dat)
      mean_GCA_S3 <- unlist(mean_GCA_S3, recursive = FALSE)
      mean_GCA_S4 <- unlist(mean_GCA_S4, recursive = FALSE)
      
      output$accGCA_S2[year]         <- cor(s2_pheno@pheno, mean_GCA_S2@pheno)
      output$accGCA_Testcross1[year] <- cor(s3_pheno, mean_GCA_S3)
      output$accGCA_Testcross2[year] <- cor(s4_pheno, mean_GCA_S4)
      
      # 2. Synthetic performance benchmark – progeny test
      for (group in seq_along(S2s)) {
        s2_dat[[group]] <- setPhenoProgTest(
          s2_dat[[group]],
          tester_pop,
          nMatePerInd = 100,
          H2 = 1
        )
      }
      
      for (group in seq_along(S3s)) {
        for (line in 1:length(S3s[[group]])) {
          s3_dat[[group]][[line]] <- setPhenoProgTest(
            s3_dat[[group]][[line]],
            tester_pop,
            nMatePerInd = 100,
            H2 = 1
          )
        }
      }
      
      for (group in seq_along(S4)) {
        for (line in 1:length(S4[[group]])) {
          s4_dat[[group]][[line]] <- setPhenoProgTest(
            s4_dat[[group]][[line]],
            tester_pop,
            nMatePerInd = 100,
            H2 = 1
          )
        }
      }
      
      mean_prog_S3 <- lapply(s3_dat, function(sublist) {
        sapply(sublist, function(pop) mean(pop@pheno))
      })
      
      mean_prog_S4 <- lapply(s4_dat, function(sublist) {
        sapply(sublist, function(pop) mean(pop@pheno))
      })
      
      mean_prog_S2 <- mergePops(s2_dat)
      mean_prog_S3 <- unlist(mean_prog_S3, recursive = FALSE)
      mean_prog_S4 <- unlist(mean_prog_S4, recursive = FALSE)
      
      output$accprogeny_S2[year]         <- cor(s2_pheno@pheno, mean_prog_S2@pheno)
      output$accprogeny_Testcross1[year] <- cor(s3_pheno, mean_prog_S3)
      output$accprogeny_Testcross2[year] <- cor(s4_pheno, mean_prog_S4)
      
      # -----------------------------------------------------------------------
      # Update synthetics and testers
      # -----------------------------------------------------------------------
      
      # Assign phenotype from syn1 to S5 object that was created two cycles ago
      # (S5twoyearsago)
      for (group in seq_along(S5twoyearsago)) {
        for (line in 1:length(S5twoyearsago[[group]])) {
          S5twoyearsago[[group]][[line]]@pheno[, 1] = meanP(syn1[[group]][[line]])
        }
      }
      
      # Create list of testers/lines for synthetic variety from best performers
      # from each heterotic group in last years S5. This is placed here as
      # object S5 updated in Year 6
      sel_elite <- lapply(S5twoyearsago, function(group) {
        means <- vapply(group, function(pop) mean(pop@pheno), numeric(1))
        keep  <- order(means, decreasing = TRUE)[seq_len(min(nElite, length(means)))]
        group[keep]
      })
      
      elite = list()
      sel_elite = unlist(sel_elite)
      
      for (group in seq_along(sel_elite)) {
        elite[[group]] = selectOP(
          sel_elite[[group]],
          nInd = 200,
          nSeeds = 10,
          probSelf = selfing, ###
          pollenControl = FALSE,
          use = "rand"
        )
      }
      
      # Separate heterotic groups
      synthetics <- list(
        list(elite[[2]], elite[[3]], elite[[4]], elite[[5]]),
        list(elite[[1]], elite[[3]], elite[[4]], elite[[5]]),
        list(elite[[1]], elite[[2]], elite[[4]], elite[[5]]),
        list(elite[[1]], elite[[2]], elite[[3]], elite[[5]]),
        list(elite[[1]], elite[[2]], elite[[3]], elite[[4]])
      )
      
      testers2 <- lapply(synthetics, function(group_list) {
        selected <- lapply(group_list, function(pop) {
          selectInd(pop, nInd = 1, use = "rand")
        })
        mergePops(selected)
      })
      
      # Testers for testcross 1
      testers1 <- lapply(testers2, selectInd, nInd = 2, use = "rand")
      
      # Tester for benchmarking accuracy of selection
      tester_pop <- list(elite[[1]], elite[[2]], elite[[3]], elite[[4]], elite[[5]])
      
      for (group in seq_along(tester_pop)) {
        tester_pop[[group]] = selectInd(tester_pop[[group]], nInd = 1, use = "rand")
      }
      
      tester_pop = mergePops(tester_pop)
      
      # Update S5lastyear for next cycle
      S5twoyearsago = S5lastyear
      S5lastyear    = S5
      
      # -----------------------------------------------------------------------
      # Report results
      # -----------------------------------------------------------------------
      
      output$meanSyn0[year]     = meanG(mergePops(unlist(syn0)))
      output$varSyn0[year]      = varG(mergePops(unlist(syn0)))
      output$meanGvariety[year] = meanG(synthetic_var)
      output$varGvariety[year]  = varG(synthetic_var)
      output$meanGParents[year] = meanG(Parents_output)
      output$meanEBVParents[year] = meanEBV(Parents_output)
      output$varGParents[year]  = varG(Parents_output)
      output$varAParents[year]  = varA(Parents_output)
      output$varDParents[year]  = varD(Parents_output)
      
      syn1_dat = mergePops(unlist(syn1))
      syn2_dat = mergePops(unlist(syn2))
      
      output$accSyn1[year] = cor(syn1_dat@gv, syn1_dat@pheno)
      output$accSyn2[year] = cor(syn2_dat@gv, syn2_dat@pheno)
      
      output$varGS2[year] = varG(s2_pheno)
      output$varDS2[year] = varD(s2_pheno)
      
      output$varGS2_1[year] = varG(S2s[[1]])
      output$varGS2_2[year] = varG(S2s[[2]])
      output$varGS2_3[year] = varG(S2s[[3]])
      output$varGS2_4[year] = varG(S2s[[4]])
      output$varGS2_5[year] = varG(S2s[[5]])
      
      output$varDS2_1[year] = varD(S2s[[1]])
      output$varDS2_2[year] = varD(S2s[[2]])
      output$varDS2_3[year] = varD(S2s[[3]])
      output$varDS2_4[year] = varD(S2s[[4]])
      output$varDS2_5[year] = varD(S2s[[5]])
    }
  }
  
  output
}

out <- model.mse(rep_id)

saveRDS(
  out,
  file = file.path(workdir, sprintf("SyntheticPheno_DD1_rep%03d.rds", rep_id))
)