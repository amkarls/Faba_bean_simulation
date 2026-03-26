
############------Training population data made based on per se performance of training population

####Below script uses the following steps for implementing GS
## Splits heterotic groups by markers
## Trains prediction models on per se performance of lines


## --- ARGUMENTS / INTERACTIVE MODE ---
## If you're in R, define rep_id & workdir first; otherwise read commandArgs.
if (exists("rep_id", inherits = FALSE) && exists("workdir", inherits = FALSE)) {
  # interactive: use the pre-defined variables
  rep_id  <- as.integer(rep_id)
  workdir <- as.character(workdir)
} else {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 1) stop("Usage: GS_synthetic_BGLR_perse_DD1.R <rep_id> [workdir]")
  rep_id  <- as.integer(args[1])
  workdir <- if (length(args) >= 2) args[2] else getwd()
}

if (!dir.exists(workdir)) stop("Workdir does not exist: ", workdir)
setwd(workdir)



suppressPackageStartupMessages({
  library(AlphaSimR)
  library(dplyr)
  library(rrBLUP)
  library(tibble)
  library(stats)
  library(cluster)
  library(BGLR)
})

## Always keep BLAS single-threaded inside an array task
Sys.setenv(OMP_NUM_THREADS="1", MKL_NUM_THREADS="1",
           OPENBLAS_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1",
           NUMEXPR_NUM_THREADS="1")

## --- Use workdir for ALL I/O ---
load("SP_ADG_DD1.RData")
load("Divpanel_ADG_DD1.RData")
load(sprintf("workspace_ADG_DD1_%d_2.RData", rep_id))



## ---- your simulation function (trimmed for brevity) ----
model.mse <- function(rep_id) {
  
  
  
  ##############Parameters
  
  
  
  #Parents
  nParents=21
  nCrosses=160
  nElite=1
  
  nSelf =16 ## number of seeds per selfed F1
  nself1=1
  nself2=10
  s2sel =200 #number of individuals selected from S2 per heterotic group
  s3sel=20 #number of lines selected from each heterotic group
  #parentsel=10 #number of lines selected from S2 to be selfed and used for crossing
  nsel_syn1=4 ##number of lines per heterotic group
  nsel_syn2=2##Number of varieties selected
  
  ###Degree of selfing
  selfing=0.5
  
  ####Accuracy of trials
  #Repliation at each stage
  varE=4 #error variance for field trials
  reptraining=3
  repsyn1=3
  repsyn2=9
  reporiginaltraining=6
  
  
  #--------------------------------------------Function for training------------
  predictEBV_BGLR <- function(pop, p_tr, bA, bD = NULL, mu) {
    # pop  : AlphaSimR Pop to predict
    # p_tr : allele freqs from training population (length m)
    # bA   : additive marker effects (length m)
    # bD   : dominance marker effects (length m) or NULL
    # mu   : intercept
    
    G <- pullSnpGeno(pop)          # nInd x m, 0/1/2
    G <- as.matrix(G)
    storage.mode(G) <- "numeric"
    
    stopifnot(
      ncol(G) == length(p_tr),
      length(bA) == length(p_tr)
    )
    
    q_tr <- 1 - p_tr
    
    # Additive code: centered on 2p
    A <- sweep(G, 2, 2 * p_tr, "-")
    gv <- as.numeric(A %*% bA)
    
    # Dominance contribution (only if bD exists)
    if (!is.null(bD) && length(bD) > 0) {
      stopifnot(length(bD) == length(p_tr))
      
      het <- (G == 1) * 1
      D   <- sweep(het, 2, 2 * p_tr * q_tr, "-")
      gv  <- gv + as.numeric(D %*% bD)
    }
    
    gv <- gv + mu
    pop@ebv <- matrix(gv)
    pop
  }
  
  
  #-------------------------------------------Training population---------------
  
  
  #Add 100 random accessions from F4 (F4_split) to boost training population
  add_training = vector("list",length(F4_split)) # Keep lines separate
  for(i in 1:length(F4_split)){ 
    add_training[[i]]<-selectInd(F4_split[[i]], nInd = 1, use = "rand")
  }
  add_training=mergePops(add_training)
  add_training=selectInd(add_training, nInd=100, use="rand")
  #These will be phenotyped in year 1
  
  
  #-------------------------Create heterotic groups-----------------------
  # Split heterotic pools to form initial parents
  # Initialize variables
  nGroups <- 5
  
  #Starting population
  start_pop = c(Parents_17, Parents_18, Parents_19)
  
  
  # 1. Extract genotype marker matrix (0/1/2 coding)
  geno <- pullSnpGeno(start_pop)  # rows = individuals, cols = markers
  
  
  # 2. Calculate genetic distance matrix
  # You can use Euclidean distance or a custom metric
  G=A.mat(geno, min.MAF = 0.05,impute.method = "mean")  
  D2 <- outer(diag(G), diag(G), "+") - 2*G
  D2[D2 < 0] <- 0                 # numerical safeguard
  D  <- as.dist(sqrt(D2))  
  
  # 3. Perform hierarchical clustering
  hc <- hclust(D, method = "ward.D2")
  
  # 4. Cut dendrogram into 5 clusters
  clusters <- cutree(hc, k = nGroups)
  
  # 5. Assign cluster labels to individuals
  start_pop@misc$cluster <- clusters
  
  rm(geno,G, D, hc, clusters)
  
  
  #Create testers for training population GCA evaluation
  # Initialize list to store tester individuals
  tester_list <- vector("list", nGroups)
  
  # Select one random individual from each of the clusters
  for (k in 1:nGroups) {
    cluster_inds <- which(start_pop@misc$cluster == k)
    selected_ind <- sample(cluster_inds, 1)
    tester_list[[k]] <- start_pop[selected_ind]
  }
  
  # Combine into one tester population from which one inbred line is developed per tester
  tester_pop <- makeDH(mergePops(tester_list))
  
  
  ##Create training population by creating one inbred line from each accession of diversity panel (start_pop should already be largely homozygous))
  training =makeDH(mergePops(list(start_pop, Divpanel)), nDH=1)
  
  
  #Per se performance as training data
  
  training <- setPheno(training, varE = varE, reps = reporiginaltraining)
  
  # === Step 2: initialize the containers ===
  training_by_year <- list()
  
  
  # === Step 3: assign the first-year training and tester ===
  training_by_year[["Y0"]] <- training
  
  
  # Append this year 
  year_key <- paste0("Y", 0)  
  training_by_year[[year_key]] <- training
  
  # Build combined training to get global allele frequencies
  training_all <- mergePops(training_by_year)
  G_tr <- pullSnpGeno(training_all)
  G_tr <- as.matrix(G_tr)
  storage.mode(G_tr) <- "numeric"
  p_tr <- colMeans(G_tr) / 2
  q_tr <- 1 - p_tr
  
  genoA_list <- list()
  genoD_list <- list()
  y_list     <- list()
  
  for (k in seq_along(training_by_year)) {
    tr_pop <- training_by_year[[k]]
    y_list[[k]] <- as.numeric(tr_pop@pheno)
    
    G <- pullSnpGeno(tr_pop)           # nTrain x m (0/1/2)
    G <- as.matrix(G)
    storage.mode(G) <- "numeric"
    
    stopifnot(ncol(G) == length(p_tr))
    
    # Use global p_tr (or you can use year-specific if you prefer)
    A_blk <- sweep(G, 2, 2 * p_tr, "-")
    het   <- (G == 1) * 1
    D_blk <- sweep(het, 2, 2 * p_tr * q_tr, "-")
    
    colnames(A_blk) <- colnames(G)
    colnames(D_blk) <- colnames(G)
    
    genoA_list[[k]] <- A_blk
    genoD_list[[k]] <- D_blk
  }
  
  genoA <- do.call(rbind, genoA_list)
  genoD <- do.call(rbind, genoD_list)
  y     <- unlist(y_list)
  
  ETA <- list(
    list(X = genoA, model = "BRR")  # additive
    # list(X = genoD, model = "BRR")   # no dominance in first training population
  )
  
  fit <- BGLR(y = y, ETA = ETA, nIter = 10000, burnIn = 3000, verbose = FALSE)
  bA <- fit$ETA[[1]]$b
  # bD <- fit$ETA[[2]]$b
  mu <- as.numeric(fit$mu)
  bD <- NULL
  
  
  
  ############Populations for filling pipeline
  
  # Split heterotic pools to form initial parents
  #  Split population into groups based on clusters
  parent_groups <- lapply(1:nGroups, function(k) {
    subset(start_pop, start_pop@misc$cluster == k)
  })
  
  ###Random selection of individuals that will be used to evaluate selection accuracy
  benchmarking<- lapply(
    parent_groups,
    function(pop) selectInd(pop, nInd = 3, use = "rand")  # or use = "random"
  )
  
  
  # Select random elites from each group based on genetic value
  elite_groups <- lapply(parent_groups, function(group) {
    selectInd(group, nElite, use = "rand")
  })
  
  
  
  ###Produce elite line groups
  G1 = c(elite_groups[[2]], elite_groups[[3]], elite_groups[[4]], elite_groups[[5]])
  G2 = c(elite_groups[[1]], elite_groups[[3]], elite_groups[[4]], elite_groups[[5]])
  G3= c(elite_groups[[4]], elite_groups[[1]], elite_groups[[2]], elite_groups[[5]])
  G4=c(elite_groups[[1]], elite_groups[[2]], elite_groups[[3]], elite_groups[[5]])
  G5=c(elite_groups[[1]], elite_groups[[2]], elite_groups[[3]], elite_groups[[4]])
  
  
  #############Elite lines to use in first synthetics
  synthetics=list(G1,G2,G3,G4,G5)
  
  tester_pop=c(elite_groups[[1]], elite_groups[[2]], elite_groups[[3]], elite_groups[[4]],elite_groups[[5]])
  
  
  #------------------------------------Start running breeding cycle-------------
  
  # Creating empty vectors to store genetic values
  
  nYears = 20
  
  
  
  output = data.frame(year = 1:nYears)
  
  
  for (year in 1:nYears) {
    
    cat("Cycle year:", year, "of 20\n")
    
    if(year == 1){
      
      #------------------------------------Creation of synthetic varieties------
      
      #------------------ These activities are run as if they have been performed during the last few years of burnin
      
      #First synthetic varieties created with selected parents from year 17,18,19 of burnin. 
      #Filling last steps of pipeline with synthetics derived from these lines
      
      #This activity is done during burnin
      ###Create SYN0 by crossing all the synthetic components
      syn0=replicate(n=5, expr={list()}, simplify=F)
      for(group in seq_along(parent_groups)){
        for (i in 1:length(parent_groups[[group]])){
          open1 =mergePops(list(parent_groups[[group]][i], synthetics[[group]]))
          syn0[[group]][i]= selectOP(
            open1,
            nInd=5, 
            nSeeds=500, #unrealistic but I think best way to simulate this open pollination
            probSelf = selfing,###
            pollenControl = FALSE,
            use = "rand")
        }}
      
      #This activity is done during burnin
      ####Random open pollination of Syn0 to produce Syn1
      syn1 = replicate(n=5, expr={list()}, simplify=F) # Keep lines separate
      for(group in seq_along(syn1)){
        for (i in 1:length(syn0[[group]])){
          syn1[[group]][[i]] =selectOP(syn0[[group]][[i]],nInd=2000,nSeeds=1, use="rand", probSelf=selfing)
          syn1[[group]][[i]]=setPheno(syn1[[group]][[i]], varE=varE, reps=repsyn1)
        }}
      
      
      
      # Step 1: Calculate mean phenotype per SYN1 subpopulation
      Syn1_mean <- vector("list", length(syn1))
      
      for (group in seq_along(syn1)) {
        group_means <- sapply(syn1[[group]], function(pop) mean(pop@pheno, na.rm = TRUE))
        Syn1_mean[[group]] <- data.frame(index = seq_along(group_means), mean_pheno = group_means)
      }
      
      # Step 2: Select top-performing subpopulations in each group
      Syn1_sel <- vector("list", length(syn1))
      
      for (group in seq_along(syn1)) {
        mean_df <- Syn1_mean[[group]]
        top_indices <- order(mean_df$mean_pheno, decreasing = TRUE)[1:nsel_syn1]
        Syn1_sel[[group]] <- syn1[[group]][top_indices]
      }
      
      # Step 3: Combine selected subpopulations from all groups
      Syn1_sel <- unlist(Syn1_sel, recursive = FALSE)
      
      
      #This step is performed during first year of breeding
      ### Open pollination of Syn-2 to produce Syn-3 and phenotype
      syn2 <- lapply(Syn1_sel, function(line) {
        pop <- selectOP(line,
                        nInd=500, 
                        nSeeds=4,
                        probSelf = selfing,
                        pollenControl = FALSE,
                        use = "rand"
        )
        setPheno(pop, varE=varE, reps=repsyn2)
      })
      
      ### Selection from Syn2: random selection as phenotype not collected until first year of breeding
      
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
        )}
      
      synthetic_var=mergePops(syn3)
      
      
      
      
      #---------------------------Year 3--------------------------------------------
      F3=mergePops(F3)
      
      maf_min=0.05
      
      #Split F3 population into heterotic groups aligned with parental groups
      # Genotypes
      G_ind   <- pullSnpGeno(F3)               # n × m (0/1/2), 
      AF_grpL <- lapply(parent_groups, \(p) colMeans(pullSnpGeno(p), na.rm=TRUE)/2)
      AF_grp  <- do.call(rbind, AF_grpL)            # 5 × m
      
      # Global allele frequency (across inds + group refs) for standardization
      p_all <- colMeans(rbind(G_ind/2, AF_grp), na.rm=TRUE)
      
      # MAF filter
      keep  <- which(p_all > maf_min & p_all < 1 - maf_min)
      G_ind <- G_ind[, keep, drop=FALSE]
      AF_grp<- AF_grp[, keep, drop=FALSE]
      p_all <- p_all[keep]
      m     <- length(keep)
      
      # Standardization denominator
      std <- sqrt(2 * p_all * (1 - p_all))
      std[std == 0] <- 1  # safeguard
      
      # Z for individuals: (X - 2p)/sqrt(2p(1-p))
      Z_ind <- sweep(G_ind, 2, 2*p_all, "-")
      Z_ind <- sweep(Z_ind, 2, std, "/")
      Z_ind[!is.finite(Z_ind)] <- 0
      
      # Z for group centers: (2*AF - 2p)/sqrt(2p(1-p))
      Z_grp <- sweep(2*AF_grp, 2, 2*p_all, "-")
      Z_grp <- sweep(Z_grp, 2, std, "/")
      Z_grp[!is.finite(Z_grp)] <- 0
      
      # Similarity of each individual to each group center: mean product over markers
      # 12800×m  %*%  m×5  ->  12800×5
      S <- Z_ind %*% t(Z_grp) / m
      
      # Assign group = argmax similarity
      group_assign <- max.col(S, ties.method = "random")  # values 1..5
      
      # Store on Pop and (optionally) split
      F3@misc$cluster <- group_assign
      S3 <- lapply(split(seq_len(nInd(F3)), group_assign),
                   function(idx) F3[idx])
      
      
      # Step 2: Apply genomic selection on F3, selecting best lines
      S3<- lapply(S3, function(F3) {
        predictEBV_BGLR(F3, p_tr, bA, bD, mu)
      })
      
      S3s_sel <- lapply(S3, function(F3) {
        selectFam(F3, nFam= s3sel, use = "ebv", famType="F")
      })
      
      # Split lines by mother ID
      F3_split_list <- vector("list", length(S3s_sel))  # will store the split-by-mother pops per group
      
      for (group in seq_along(S3s_sel)) {
        
        # Current group population
        group_pop <- S3s_sel[[group]]
        
        # Split this group into sub-pops by maternal ID
        mother_ids <- unique(group_pop@mother)
        F3_split <- vector("list", length(mother_ids))
        
        for (i in seq_along(mother_ids)) {
          F3_split[[i]] <- subset(group_pop, group_pop@mother == mother_ids[i])
        }
        
        # Optionally name each family by its maternal ID
        names(F3_split) <- as.character(mother_ids)
        
        # Store this list of maternal families in the corresponding heterotic group
        F3_split_list[[group]] <- F3_split
      }
      
      
      S3s_sel =F3_split_list
      
      #Collection of seeds from selected lines
      S4= vector("list", length(S3s_sel))
      for (group in seq_along(S3s_sel)){
        for (line in 1:length(S3s_sel[[group]])) {
          S4[[group]][[line]]=self(S3s_sel[[group]][[line]], nProgeny=200)
        }
      }
      
      #Selfing in counter-season location to produce S5
      S5= vector("list", length(S4))
      for (group in seq_along(S4)){
        for (line in 1:length(S4[[group]])) {
          S5[[group]][[line]]=self(S4[[group]][[line]], nProgeny=1)
        }
      }
      
      #---------------------------Year 2----------------------------------------
      
      #Selection of S2 (in this case F2) and creation of S3 
      
      #Split F2 into heterotic groups, disregarding line
      
      F2_i =mergePops(F2)
      
      
      # Genotypes
      G_ind   <- pullSnpGeno(F2_i)               # n × m (0/1/2), 
      AF_grpL <- lapply(parent_groups, \(p) colMeans(pullSnpGeno(p), na.rm=TRUE)/2)
      AF_grp  <- do.call(rbind, AF_grpL)            # 5 × m
      
      # Global allele frequency (across inds + group refs) for standardization
      p_all <- colMeans(rbind(G_ind/2, AF_grp), na.rm=TRUE)
      
      # MAF filter
      keep  <- which(p_all > maf_min & p_all < 1 - maf_min)
      G_ind <- G_ind[, keep, drop=FALSE]
      AF_grp<- AF_grp[, keep, drop=FALSE]
      p_all <- p_all[keep]
      m     <- length(keep)
      
      # Standardization denominator
      std <- sqrt(2 * p_all * (1 - p_all))
      std[std == 0] <- 1  # safeguard
      
      # Z for individuals: (X - 2p)/sqrt(2p(1-p))
      Z_ind <- sweep(G_ind, 2, 2*p_all, "-")
      Z_ind <- sweep(Z_ind, 2, std, "/")
      Z_ind[!is.finite(Z_ind)] <- 0
      
      # Z for group centers: (2*AF - 2p)/sqrt(2p(1-p))
      Z_grp <- sweep(2*AF_grp, 2, 2*p_all, "-")
      Z_grp <- sweep(Z_grp, 2, std, "/")
      Z_grp[!is.finite(Z_grp)] <- 0
      
      # Similarity of each individual to each group center: mean product over markers
      # 12800×m  %*%  m×5  ->  12800×5
      S <- Z_ind %*% t(Z_grp) / m
      
      # Assign group = argmax similarity
      group_assign <- max.col(S, ties.method = "random")  # values 1..5
      
      # Store on Pop and (optionally) split
      F2_i@misc$cluster <- group_assign
      assigned_groups <- lapply(split(seq_len(nInd(F2_i)), group_assign),
                                function(idx) F2_i[idx])
      
      
      
      # Step 2: Apply genomic selection on S2
      S2<- lapply(assigned_groups, function(group) {
        predictEBV_BGLR(group, p_tr, bA, bD, mu)
      })
      
      S2s_sel <- lapply(S2, function(group) {
        selectInd(group, nInd= s2sel, use = "ebv")
      })
      
      
      
      #---------------------------Year 1----------------------------------------
      
      
      # Step 1: Cross each group
      F1 <- lapply(parent_groups, function(parents) {
        randCross(parents, nCrosses = nCrosses, nProgeny = 1)
      })
      
      
      #-------------------------------------Prepare training data and parents-------
      
      
      #----------------------------New breeding parents-----------------------------
      
      ##Create parents and lines for synthetics from previous years S3, 
      #use tester crosses with synthetics for parental selection and training
      
      # parent_groups: list of populations, one per heterotic group
      # parentsel: number of divergent parents to select per group
      
      # For each heterotic group i:
      parent_sel= vector("list", length(S3s_sel))
      for (group in seq_along(S3s_sel)){
        for (line in 1:length(S3s_sel[[group]])) {
          parent_sel[[group]][[line]]=selectInd(S3s_sel[[group]][[line]], nInd = 1, use = "ebv")
        }
      }
      
      #merge individuals in each heterotic group
      parent_groups= vector("list", length(parent_sel))
      for (group in seq_along(parent_sel)){
        parent_groups[[group]]=mergePops(parent_sel[[group]])
      }
      
      
      
      
      #----------------------------New training data--------------------------------
      
      ################Use per se performance data from S4 (only in first year's cycle) for training
      
      S5_training=replicate(n=5, expr={list()}, simplify=F) ##Named S5_training to conform with later years training data
      for (group in seq_along(S4)){
        for (line in 1:length(S4[[group]])){
          S5_training[[group]][[line]]=selectInd(S4[[group]][[line]], nInd=1, use="rand")
        }}
      
      for (group in seq_along(S5_training)){
        S5_training[[group]]=mergePops(unlist(S5_training[[group]])) }
      
      # S5_training already holds 1 individual per S4 family per group
      for (group in seq_along(S5_training)) {
        S5_training[[group]] <- setPheno(S5_training[[group]], varE = varE, reps = reptraining)
      }
      
      
      ##Assign phenotypes to motherlines before merging S5 training
      for (group in 1:length(S5_training)){
        for (individual in 1:length(S5_training[[group]])){
          S4[[group]][[individual]]@pheno[,1]=S5_training[[group]]@pheno[individual]}}
      
      #merge S5 training
      S5_training=mergePops(S5_training)
      
      # Same idea for added training individuals from F4
      add_training <- setPheno(add_training, varE = varE, reps = reptraining)
      
      # Append this year 
      year_key <- paste0("Y", year)  
      training_by_year[[year_key]] <- mergePops(list(S5_training,add_training))
      
      # Build combined training to get global allele frequencies
      training_all <- mergePops(training_by_year)
      G_tr <- pullSnpGeno(training_all)
      G_tr <- as.matrix(G_tr)
      storage.mode(G_tr) <- "numeric"
      p_tr <- colMeans(G_tr) / 2
      q_tr <- 1 - p_tr
      
      genoA_list <- list()
      genoD_list <- list()
      y_list     <- list()
      
      for (k in seq_along(training_by_year)) {
        tr_pop <- training_by_year[[k]]
        y_list[[k]] <- as.numeric(tr_pop@pheno)
        
        G <- pullSnpGeno(tr_pop)           # nTrain x m (0/1/2)
        G <- as.matrix(G)
        storage.mode(G) <- "numeric"
        
        stopifnot(ncol(G) == length(p_tr))
        # Use global p_tr (or you can use year-specific if you prefer)
        A_blk <- sweep(G, 2, 2 * p_tr, "-")
        het   <- (G == 1) * 1
        D_blk <- sweep(het, 2, 2 * p_tr * q_tr, "-")
        
        colnames(A_blk) <- colnames(G)
        colnames(D_blk) <- colnames(G)
        
        genoA_list[[k]] <- A_blk
        genoD_list[[k]] <- D_blk
      }
      
      genoA <- do.call(rbind, genoA_list)
      genoD <- do.call(rbind, genoD_list)
      y     <- unlist(y_list)
      
      ETA <- list(
        list(X = genoA, model = "BRR"),  # additive
        list(X = genoD, model = "BRR")   # dominance
      )
      
      fit <- BGLR(y = y, ETA = ETA, nIter = 10000, burnIn = 3000, verbose = FALSE)
      bA <- fit$ETA[[1]]$b
      bD <- fit$ETA[[2]]$b
      mu <- as.numeric(fit$mu)
      
      
      
      
      #--------------------------Creation of new elites/tester lines---------------
      #--------Accuracy of GS models based on benchmarking progeny test------------- 
      #EValuate accuracy before updating testers
      #2. Progeny test benchmark
      
      S2_dat=mergePops(S2)
      S3_dat=mergePops(S3)
      
      
      S2_dat=setPhenoProgTest(S2_dat,tester_pop, nMatePerInd = 100, use="gv", H2=1)
      S3_dat=setPhenoProgTest(S3_dat,tester_pop, nMatePerInd = 100, use="gv",H2=1)
      
      output$acc_progenytest_S2[year] <- cor(S2_dat@ebv, S2_dat@pheno, use = "complete.obs")
      output$acc_progenytest_S3[year] <- cor(S3_dat@ebv, S3_dat@pheno, use = "complete.obs")
      
      ## 2. GCA Benchmark
      
      S2_dat <- setPhenoGCA(S2_dat, testers = tester_pop, H2 = 1.0)
      output$acc_GCA_S2[year] <- cor(S2_dat@ebv, S2_dat@pheno, use = "complete.obs")
      
      S3_dat <- setPhenoGCA(S3_dat, testers = tester_pop, H2 = 1.0)
      output$acc_GCA_S3[year] <- cor(S3_dat@ebv, S3_dat@pheno, use = "complete.obs")
      
      
      ##Create list of testers/lines for synthetic variety from best performers from each heterotic group.
      #For first cycle S4 will be used but in subsequent cycles S5lastyear
      
      
      # Step 1: Select one elite individual per group from S4 based on performance of synthetics
      
      
      elite <- vector("list", length(S4))
      
      
      elite <- lapply(S4, function(sublist) {
        # Calculate means for each object in the sublist
        means <- sapply(sublist, function(pop) mean(pop@pheno)) 
        
        # Get index of highest mean
        best_idx <- which.max(means)
        
        # Return the object with highest mean
        sublist[[best_idx]]
      })
      
      for (group in seq_along(elite)){
        elite[[group]]=selectInd(elite[[group]], nInd=200, use="rand")
      }
      
      ###Produce elite line groups
      G1 = c(elite[[2]], elite[[3]], elite[[4]], elite[[5]])
      G2 = c(elite[[1]], elite[[3]], elite[[4]], elite[[5]])
      G3= c(elite[[4]], elite[[1]], elite[[2]], elite[[5]])
      G4=c(elite[[1]], elite[[2]], elite[[3]], elite[[5]])
      G5=c(elite[[1]], elite[[2]], elite[[3]], elite[[4]])
      
      
      #############Elite lines to use in first synthetics
      synthetics=list(G1,G2,G3,G4,G5)
      
      for (group in seq_along(elite)){
        elite[[group]]=selectInd(elite[[group]], nInd=1, use="rand")
      }
      
      #Create new tester pop for training
      tester_pop=c(elite[[1]], elite[[2]], elite[[3]], elite[[4]], elite[[5]])
      
      S5lastyear=S5 
      
      
      
      
      #---------------------------------------Output file------------------------------
      syn0_dat=mergePops(unlist(syn0))
      syn1_dat=mergePops(unlist(syn1))
      syn2_dat=mergePops(unlist(syn2))
      Parents_output=mergePops(parent_groups)
      
      output$meanSyn0[year] = meanG(syn0_dat)
      output$varGSyn0[year] = varG(syn0_dat)
      output$varDSyn0[year] = varD(syn0_dat)
      output$varASyn0[year] = varA(syn0_dat)
      output$meanSyn2[year] = meanG(syn2_dat)
      output$varGSyn2[year] = varG(syn2_dat)
      output$varASyn2[year] = varA(syn2_dat)
      output$varDSyn2[year] = varD(syn2_dat)
      output$meanGvariety[year] = meanG(synthetic_var)
      output$varGvariety[year] = varG(synthetic_var)
      output$meanGParents[year]= meanG(Parents_output)
      output$varGParents[year] = varG(Parents_output)
      output$varAParents[year] = varA(Parents_output)
      output$varDParents[year] = varD(Parents_output)
      output$accSyn1[year] = cor(syn1_dat@gv, syn1_dat@pheno)
      output$accSyn2[year] = cor(syn2_dat@gv, syn2_dat@pheno)
      
      
    } else {
      
      #-------------------------------Breeding year 2-20----------------------------
      
      
      
      ###------------------------------Year 7---------------------------------------
      
      # Year 7
      # Open pollination of selected SYN2 lines to create SYN3
      
      syn3 <- vector("list", length(Syn2_sel))  # Initialize list for SYN2
      
      for (group in seq_along(Syn2_sel)) {
        syn3[[group]] <- selectOP(
          pop = Syn2_sel[[group]],
          nInd = 500,         # Number of seeds per line for multiplication
          nSeeds = 4,
          probSelf = selfing, # Mostly selfed, some open pollination
          pollenControl = FALSE,
          use = "rand"
        )}
      
      
      ##Variety release
      synthetic_var=mergePops(syn3)
      
      ###------------------------------Year 6---------------------------------------
      # Selection from SYN1 
      
      # Step 1: Calculate mean phenotype per SYN1 subpopulation
      Syn1_mean <- vector("list", length(syn1))
      
      for (group in seq_along(syn1)) {
        group_means <- sapply(syn1[[group]], function(pop) mean(pop@pheno, na.rm = TRUE))
        Syn1_mean[[group]] <- data.frame(index = seq_along(group_means), mean_pheno = group_means)
      }
      
      # Step 2: Select top-performing subpopulations in each group
      Syn1_sel <- vector("list", length(syn1))
      
      for (group in seq_along(syn1)) {
        mean_df <- Syn1_mean[[group]]
        top_indices <- order(mean_df$mean_pheno, decreasing = TRUE)[1:nsel_syn1]
        Syn1_sel[[group]] <- syn1[[group]][top_indices]
      }
      
      # Step 3: Combine selected subpopulations from all groups
      Syn1_sel <- unlist(Syn1_sel, recursive = FALSE)
      
      # Open pollination of selected SYN1 lines (Syn1_sel) to produce SYN2 (syn2)
      syn2 <- vector("list", length(Syn1_sel))  # Initialize list for SYN2
      
      for (group in seq_along(Syn1_sel)) {
        syn2[[group]] <- selectOP(
          pop = Syn1_sel[[group]],
          nInd = 500,         # Number of seeds per line for multiplication
          nSeeds = 4,
          probSelf = selfing, # Mostly selfed, some open pollination
          pollenControl = FALSE,
          use = "rand"
        )
        
        # Assign phenotypes to SYN2 population
        syn2[[group]] <- setPheno(syn2[[group]], varE=varE, reps=repsyn2)
      }
      
      
      ###Selection from Syn2  
      # Step 1: Compute mean phenotype for each SYN2 population
      Syn2_mean <- data.frame(
        index = seq_along(syn2),
        mean_pheno = sapply(syn2, function(pop) mean(pop@pheno, na.rm = TRUE))
      )
      
      # Step 2: Identify top populations
      top_indices <- order(Syn2_mean$mean_pheno, decreasing = TRUE)[1:nsel_syn2]
      
      # Step 3: Select top SYN2 populations
      Syn2_sel <- syn2[top_indices]
      
      ###------------------------------Year 5---------------------------------------
      
      # Random open pollination of SYN0 to produce SYN1
      syn1 <- vector("list", length(syn0))  # Keep lines separate by group
      
      for (group in seq_along(syn0)) {
        syn0_group <- syn0[[group]]
        syn1_group <- vector("list", length(syn0_group))
        
        for (i in seq_along(syn0_group)) {
          # Open-pollinate SYN0 to produce SYN1
          pop_syn1 <- selectOP(
            pop = syn0_group[[i]],
            nInd = 2000,
            nSeeds = 1,
            use = "rand",
            probSelf = selfing
          )
          
          # Assign phenotypes
          pop_syn1 <- setPheno(pop_syn1, varE=varE, reps=repsyn1)
          
          # Save
          syn1_group[[i]] <- pop_syn1
        }
        
        syn1[[group]] <- syn1_group
      }
      
      
      ###------------------------------Year 4---------------------------------------
      
      
      ###Create SYN0 by mixing synthetic components with testers and open-pollinating
      
      ###Create SYN0 by crossing all the synthetic components
      S5_syn0 <- vector("list", length(S5))  
      for(group in seq_along(S5)){
        for (individual in 1:length(S5[[group]])){
          S5_syn0[[group]][[individual]]=selectInd(S5[[group]][[individual]], nInd=200, use="rand")
        }}
      
      syn0=replicate(n=5, expr={list()}, simplify=F)
      for(group in seq_along(S5_syn0)){
        for (i in 1:length(S5_syn0[[group]])){
          open1 =mergePops(list(S5_syn0[[group]][[i]], synthetics[[group]]))
          syn0[[group]][[i]]= selectOP(
            open1,
            nInd=1000, 
            nSeeds=2, 
            probSelf = selfing,###
            pollenControl = FALSE,
            use = "rand")
        }
      }
      
      
      
      ###------------------------------Year 3---------------------------------------
      
      #  Self each line in isolation to create S3
      S3 <- vector("list", length(nGroups)) 
      S3<- lapply(S2s_sel, function(s2) {
        self(s2, nProgeny=nself2)
      })
      
      #  Apply genomic selection on S3, selecting best lines
      S3<- lapply(S3, function(group) {
        predictEBV_BGLR(group, p_tr, bA, bD, mu)
      })
      
      S3s_sel <- lapply(S3, function(group) {
        selectFam(group, nFam= s3sel, use = "ebv", famType="F")
      })
      
      # Split lines by mother ID
      S3_split_list <- vector("list", length(S3s_sel))  # will store the split-by-mother pops per group
      
      for (group in seq_along(S3s_sel)) {
        
        # Current group population
        group_pop <- S3s_sel[[group]]
        
        # Split this group into sub-pops by maternal ID
        mother_ids <- unique(group_pop@mother)
        S3_split <- vector("list", length(mother_ids))
        
        for (i in seq_along(mother_ids)) {
          S3_split[[i]] <- subset(group_pop, group_pop@mother == mother_ids[i])
        }
        
        # Optionally name each family by its maternal ID
        names(S3_split) <- as.character(mother_ids)
        
        # Store this list of maternal families in the corresponding heterotic group
        S3_split_list[[group]] <- S3_split
      }
      
      
      S3s_sel =S3_split_list
      
      #Collection of seeds from selected lines
      S4= vector("list", length(S3s_sel))
      for (group in seq_along(S3s_sel)){
        for (line in 1:length(S3s_sel[[group]])) {
          S4[[group]][[line]]=self(S3s_sel[[group]][[line]], nProgeny=200)
        }
      }
      
      #Selfing in counter-season location to produce S5
      S5= vector("list", length(S4))
      for (group in seq_along(S4)){
        for (line in 1:length(S4[[group]])) {
          S5[[group]][[line]]=self(S4[[group]][[line]], nProgeny=1)
        }
      }
      
      
      ##------------------------------Year 2---------------------------------------
      
      # Step 2: Self the F1s
      S1 <- lapply(F1, function(group) {
        self(group, nProgeny = nSelf, keepParents = FALSE)
      })
      
      # Step 2: Self the S1s
      S2 <- lapply(S1, function(group) {
        self(group, nProgeny = nself1, keepParents = FALSE)
      })
      
      
      # Step 2: Apply genomic selection on S2
      S2<- lapply(S2, function(group) {
        predictEBV_BGLR(group, p_tr, bA, bD, mu)
      })
      
      S2s_sel <- lapply( S2, function(s2) {
        selectInd(s2, nInd= s2sel, use = "ebv")
      })
      
      
      
      
      #---------------------------Year 1----------------------------------------
      
      
      # Step 1: Cross each group
      F1 <- lapply(parent_groups, function(parents) {
        randCross(parents, nCrosses = nCrosses, nProgeny = 1)
      })
      
      #-------------------------------------Prepare training data and parents-------
      
      
      #----------------------------New breeding parents-----------------------------
      
      ##Create parents and lines for synthetics from previous years S3, 
      #use tester crosses with synthetics for parental selection and training
      
      # parent_groups: list of populations, one per heterotic group
      # parentsel: number of divergent parents to select per group
      
      # For each heterotic group i:
      parent_sel= vector("list", length(S3s_sel))
      for (group in seq_along(S3s_sel)){
        for (line in 1:length(S3s_sel[[group]])) {
          parent_sel[[group]][[line]]=selectInd(S3s_sel[[group]][[line]], nInd = 1, use = "ebv")
        }
      }
      
      #merge individuals in each heterotic group
      parent_groups= vector("list", length(parent_sel))
      for (group in seq_along(parent_sel)){
        parent_groups[[group]]=mergePops(parent_sel[[group]])
      }
      
      
      
      #----------------------------New training data--------------------------------
      
      ################Use per se performance data from S5 last year  for training
      
      S5_training=replicate(n=5, expr={list()}, simplify=F) 
      for (group in seq_along(S5lastyear)){
        for (line in 1:length(S5lastyear[[group]])){
          S5_training[[group]][[line]]=selectInd(S5lastyear[[group]][[line]], nInd=1, use="rand")
        }}
      
      for (group in seq_along(S5_training)){
        S5_training[[group]]=mergePops(unlist(S5_training[[group]])) }
      
      # S5_training is 1 individual per S5lastyear family per group
      for (group in seq_along(S5_training)) {
        S5_training[[group]] <- setPheno(S5_training[[group]], varE = varE, reps = reptraining)
      }
      
      ##Assign phenotypes to motherlines before merging S5 training
      for (group in 1:length(S5_training)){
        for (individual in 1:length(S5_training[[group]])){
          S5lastyear[[group]][[individual]]@pheno[,1]=S5_training[[group]]@pheno[individual]}}
      
      #merge S5 training
      S5_training=mergePops(S5_training)
      
      # Append this year 
      year_key <- paste0("Y", year)  
      training_by_year[[year_key]] <- S5_training
      
      # Build combined training to get global allele frequencies
      training_all <- mergePops(training_by_year)
      G_tr <- pullSnpGeno(training_all)
      G_tr <- as.matrix(G_tr)
      storage.mode(G_tr) <- "numeric"
      p_tr <- colMeans(G_tr) / 2
      q_tr <- 1 - p_tr
      
      genoA_list <- list()
      genoD_list <- list()
      y_list     <- list()
      
      for (k in seq_along(training_by_year)) {
        tr_pop <- training_by_year[[k]]
        y_list[[k]] <- as.numeric(tr_pop@pheno)
        
        G <- pullSnpGeno(tr_pop)           # nTrain x m (0/1/2)
        G <- as.matrix(G)
        storage.mode(G) <- "numeric"
        
        stopifnot(ncol(G) == length(p_tr))
        # Use global p_tr (or you can use year-specific if you prefer)
        A_blk <- sweep(G, 2, 2 * p_tr, "-")
        het   <- (G == 1) * 1
        D_blk <- sweep(het, 2, 2 * p_tr * q_tr, "-")
        
        colnames(A_blk) <- colnames(G)
        colnames(D_blk) <- colnames(G)
        
        genoA_list[[k]] <- A_blk
        genoD_list[[k]] <- D_blk
      }
      
      genoA <- do.call(rbind, genoA_list)
      genoD <- do.call(rbind, genoD_list)
      y     <- unlist(y_list)
      
      ETA <- list(
        list(X = genoA, model = "BRR"),  # additive
        list(X = genoD, model = "BRR")   # dominance
      )
      
      fit <- BGLR(y = y, ETA = ETA, nIter = 10000, burnIn = 3000, verbose = FALSE)
      bA <- fit$ETA[[1]]$b
      bD <- fit$ETA[[2]]$b
      mu <- as.numeric(fit$mu)
      
      
      
      
      
      #--------------------------Creation of new elites/tester lines---------------
      
      #-------------------------------Accuracy of GS models based on benchmarking progeny test-------------------------------------------------- 
      ###Evaluate accuracy before testers are updated
      
      #2. Progeny test benchmark
      
      S2_dat=mergePops(S2)
      S3_dat=mergePops(S3)
      
      
      
      S2_dat=setPhenoProgTest(S2_dat,tester_pop, nMatePerInd = 100, use="gv", H2=1)
      S3_dat=setPhenoProgTest(S3_dat,tester_pop, nMatePerInd = 100, use="gv", H2=1)
      
      output$acc_progenytest_S2[year] <- cor(S2_dat@ebv, S2_dat@pheno, use = "complete.obs")
      output$acc_progenytest_S3[year] <- cor(S3_dat@ebv, S3_dat@pheno, use = "complete.obs")
      
      ## 2. GCA Benchmark
      
      S2_dat <- setPhenoGCA(S2_dat, testers = tester_pop, H2 = 1.0)
      output$acc_GCA_S2[year] <- cor(S2_dat@ebv, S2_dat@pheno, use = "complete.obs")
      
      S3_dat <- setPhenoGCA(S3_dat, testers = tester_pop, H2 = 1.0)
      output$acc_GCA_S3[year] <- cor(S3_dat@ebv, S3_dat@pheno, use = "complete.obs")
      
      ##Create list of testers/lines for synthetic variety from best performers from each heterotic group.
      #For first cycle S4 will be used but in subsequent cycles S5lastyear
      
      
      # Step 1: Select one elite individual per group from S4 based on performance of synthetics
      
      
      elite <- vector("list", length(S5lastyear))
      
      
      elite <- lapply(S5lastyear, function(sublist) {
        # Calculate means for each object in the sublist
        means <- sapply(sublist, function(pop) mean(pop@pheno)) 
        
        # Get index of highest mean
        best_idx <- which.max(means)
        
        # Return the object with highest mean
        sublist[[best_idx]]
      })
      
      for (group in seq_along(elite)){
        elite[[group]]=selectInd(elite[[group]], nInd=200, use="rand")
      }
      
      ###Produce elite line groups
      G1 = c(elite[[2]], elite[[3]], elite[[4]], elite[[5]])
      G2 = c(elite[[1]], elite[[3]], elite[[4]], elite[[5]])
      G3= c(elite[[4]], elite[[1]], elite[[2]], elite[[5]])
      G4=c(elite[[1]], elite[[2]], elite[[3]], elite[[5]])
      G5=c(elite[[1]], elite[[2]], elite[[3]], elite[[4]])
      
      
      #############Elite lines to use in first synthetics
      synthetics=list(G1,G2,G3,G4,G5)
      
      for (group in seq_along(elite)){
        elite[[group]]=selectInd(elite[[group]], nInd=1, use="rand")
      }
      
      #Create new tester pop for training
      tester_pop=c(elite[[1]], elite[[2]], elite[[3]], elite[[4]], elite[[5]])   
      
      
      S5lastyear=S5 
      
      
      
      #-----------------------------------End breeding cycle---------------------------------------
      
      
      #-------------------------------Report results--------------------------------------------------
      
      
      
      #---------------------------------------Output file------------------------------
      syn0_dat=mergePops(unlist(syn0))
      syn1_dat=mergePops(unlist(syn1))
      syn2_dat=mergePops(unlist(syn2))
      Parents_output=mergePops(parent_groups)
      
      output$meanSyn0[year] = meanG(syn0_dat)
      output$varGSyn0[year] = varG(syn0_dat)
      output$varDSyn0[year] = varD(syn0_dat)
      output$varASyn0[year] = varA(syn0_dat)
      output$meanSyn2[year] = meanG(syn2_dat)
      output$varGSyn2[year] = varG(syn2_dat)
      output$varASyn2[year] = varA(syn2_dat)
      output$varDSyn2[year] = varD(syn2_dat)
      output$meanGvariety[year] = meanG(synthetic_var)
      output$varGvariety[year] = varG(synthetic_var)
      output$meanGParents[year]= meanG(Parents_output)
      output$varGParents[year] = varG(Parents_output)
      output$varAParents[year] = varA(Parents_output)
      output$varDParents[year] = varD(Parents_output)
      output$accSyn1[year] = cor(syn1_dat@gv, syn1_dat@pheno)
      output$accSyn2[year] = cor(syn2_dat@gv, syn2_dat@pheno)
      
      
    }
  } 
  output
  
}

out <- model.mse(rep_id)

saveRDS(out, file = file.path(workdir, sprintf("GS_synthetic_BGLR_perse_DD1_rep%03d.rds", rep_id)))
