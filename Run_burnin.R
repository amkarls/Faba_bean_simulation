# ------------------------------------------------------------------------------
# Burn-in Simulation Script for Additive + Dominance Trait
# ------------------------------------------------------------------------------
# This script runs a 7-year burn-in phase to build generations (F1 → F7)
# followed by a 20-year simulation loop using phenotypic selection.
#
# INPUT FILES:
#   - Divpanel_ADG.RData : Founder population with additive + dominance trait
#   - SP_ADG.RData       : Simulation parameters (with meanDD for trait)
#
# OUTPUT:
#   - start_pop_ADG.RData     : Final population at end of burn-in
#   - Burninperiod_output_DD1.csv : Metrics recorded over 20 years
#
# NOTE:
# Run this script separately for each trait architecture (different meanDD values)
# ------------------------------------------------------------------------------


library(AlphaSimR)
library(dplyr)

# --- Load input data ----------------------------------------------------------
load("Divpanel_ADG.RData")   # Founder population
load("SP_ADG.RData")         # Trait architecture (with dominance model)

# --- Simulation Parameters -----------------------------------------------------

#Heritability for different stages of selection
H2F3=0.1 #Heritability in F3 single plot
H2F4=0.2 #Heritability in F4 field trials
H2F5=0.3 #Heritability in F5 field trials
H2F6=0.5 #Heritability in F6 field trials
H2F7=0.6 #Heritability in F7 field trials

# Number of crosses per year
nCrosses = 50

# Individuals produced per cross
nind = 16

# Degree of selfing when accession allowed to freely open pollinate in isolation from other accessions
selfing=0.5


###Number individuals selected per selection stage
nF1sel=200 #Random plants selected for seed multiplication
nF2sel=20 #Individual Plants selected per cross in greenhouse
nF3fam=300 #Number of families selected in F3
nF4sel=100 #Lines selected in F4
nF5sel=25 #Lines selected in F5
nF6sel=10 #Lines selected in F6
nF7sel=2 #Lines selected in F7



# --- Set initial parent population ---------------------------------------------
Parents <- Divpanel_D

# --- Output data frame for tracking progress -----------------------------------
output <- data.frame(
  year     = 1:nBurnin,
  meanG    = numeric(nBurnin),
  varG     = numeric(nBurnin),
  accSel   = numeric(nBurnin)
)



# Populate breeding programme
for (year in 1:7) {
  
  
  # F1
  F1 = vector("list",nCrosses) # Keep crosses seperate
  for(i in 1:nCrosses){ # Loop over crosses
    F1_i = randCross(Parents, nCrosses=1, nProgeny=nind)
    F1[[i]] = selectInd(F1_i, nInd = nind, use="rand")
  }
  #Create F2
  F2 = vector("list",nCrosses) # Keep crosses seperate while creating F2 seeds
  for(i in 1:nCrosses){ # Loop over crosses
    F2_i = selectOP(
      F1[[i]],
      nInd=nind,
      nSeeds=16,
      probSelf = 0.95,###
      pollenControl = FALSE,
      use = "rand")
    F2_i= setPheno(F2_i, H2=H2F2)
    F2[[i]] = selectInd(F2_i, nInd = nF1sel, use="rand") 
  }
  
  
  
  if (year < 7) {
    
    # Create F3
    F3 = vector("list",nCrosses) # Keep crosses separate
    for(i in 1:nCrosses){ # Loop over crosses
      F3_i = selectOP(
        F2[[i]],
        nInd=nF2sel,
        nSeeds=40,
        probSelf = 0.95,###
        pollenControl = FALSE,
        use = "pheno")
      F3_i = setPheno(F3_i, H2=H2F3)###Phenotype in Chile
      F3[[i]] = selectInd(F3_i, nInd = 800, use="rand")
    }
    
  }
  
  
  if (year < 6) {
    
    ###### Select among F3 to create the F4 to be trialled in Svalov. 
    F3_merged=mergePops(F3)
    
    F3_select = selectFam(
      F3_merged,
      nFam=nF3fam,
      trait = 1,
      use = "pheno",
      famType = "F")
    
    ##Split population by half-sib and collect within-family open-pollinated seeds for each "line" in Chile
    
    motherSval_F4= unique(F3_select@mother)
    F4_split=vector("list", length(motherSval_F4))
    
    for(i in 1:length(motherSval_F4)) {  #Loop over maternal parent
      motherpop= subset(F3_select, F3_select@mother==paste(motherSval_F4[i]))
      F4_split_i = selectOP(
        pop=motherpop,
        nInd=40,##40 seeds from each plant
        nSeeds=40,
        probSelf = 0.95,###
        pollenControl = FALSE,
        use = "rand") 
      F4_split_i = setPheno(F4_split_i, H2=H2F4) ###Corresponds to the phenotype of the F4 in Svalov
      F4_split[[i]] = selectInd(F4_split_i, nInd = 1600)
    } 
    
    
    ###Select in F4 for seed multiplication in Chile - assumed  seed  used for F4 Svalöv and F4 Chile is same
    
    ############Select top 50 lines
    F4_fampheno=data.frame(1:length(F4_split),NA)
    for(i in 1:length(F4_split)) {  
      
      F4_fampheno_i = mean(F4_split[[i]]@pheno) 
      F4_fampheno[i,2]= F4_fampheno_i
    }
    
    top100_F4 <-top_n(F4_fampheno, 100, NA. )
    select_F4=F4_split[top100_F4[,1]]
    
    
  }
  
  if (year < 5) {
    
    
    ####Multiply F4 seeds in Chile to create F5
    
    F4_Chile = vector("list",length(select_F4)) # Keep lines separate
    for(i in 1:length(select_F4)){ 
      F4Chile_i = selectOP(
        select_F4[[i]],
        nInd=200, ####Assumption this many seeds are held back in Chile from the F3 generation to be used for seed multiplication
        nSeeds=10,
        probSelf = selfing,###
        pollenControl = FALSE,
        use = "rand")
      F4Chile_i = setPheno(F4Chile_i, H2=H2F5) ##########Refers to phenotyping in Svalöv/Sweden of F5s
      F4_Chile[[i]] = selectInd(F4Chile_i, nInd = 2000, use="rand")
    }
    
    ######Select 25 lines from F5 based on trial in Svalov
    
    ############Select top 25 lines
    F5_fampheno=data.frame(1:length(F4_Chile),NA)
    for(i in 1:length(F4_Chile)) {  
      
      F5_fampheno_i = mean(F4_Chile[[i]]@pheno) 
      F5_fampheno[i,2]= F5_fampheno_i
    }
    
    top25_F5 <-top_n(F5_fampheno, 25, NA.)
    
    
    select_F5=F4_Chile[ top25_F5[,1]]
    
  }
  if (year < 4) {
    
    ###Multiply F5 in off-season loca
    F5_Chile = vector("list",length(select_F5)) # Keep lines separate
    for(i in 1:length(select_F5)){ 
      F5Chile_i = selectOP(
        select_F5[[i]],
        nInd=500,##assumed 500 seeds are used for multiplication
        nSeeds=10,
        probSelf = selfing,###
        pollenControl = FALSE,
        use = "rand")
      F5Chile_i = setPheno(F5Chile_i, H2=H2F6) ############Refers to phenotyping in Svalöv/Sweden at F6 stage
      F5_Chile[[i]] = selectInd(F5Chile_i, nInd = 5000, use="rand")
    }
    
    ######Select 10 lines from F6 based on trial in Svalov
    F6_fampheno=data.frame(1:length(F5_Chile),NA)
    for(i in 1:length(F5_Chile)) {  
      
      F6_fampheno_i = mean(F5_Chile[[i]]@pheno) 
      F6_fampheno[i,2]= F6_fampheno_i
    }
    
    top10_F6 <-top_n(F6_fampheno, 10, NA.)
    
    
    select_F6=F5_Chile[ top10_F6[,1]]
  }
  if (year < 3) {
    ###Multiply seeds for F7  and phenotype in multilocation trials
    F7 = vector("list",length(select_F6)) # Keep lines separate
    for(i in 1:length(select_F6)){ 
      F7_i = selectOP(
        select_F6[[i]],
        nInd=500,
        nSeeds=10,
        probSelf = selfing,###
        pollenControl = FALSE,
        use = "rand")
      F7_i = setPheno(F7_i, H2=H2F7) ############Refers to phenotyping in Svalöv/Sweden at F6 stage
      F7[[i]] = selectInd(F7_i, nInd = 5000, use="rand")
    }
    
    
    ############Select top 2 lines
    F7_fampheno=data.frame(1:length(F7),NA)
    for(i in 1:length(F7)) {  
      
      F7_fampheno_i = mean(F7[[i]]@pheno) 
      F7_fampheno[i,2]=  F7_fampheno_i 
    }
    
    top2_F7 <-top_n(F7_fampheno, nF7sel, NA.)
    
    
    select_F7=F7[top2_F7[,1]]
  }
  if (year < 2) {
    variety = vector("list",length(select_F7)) # Keep lines separate
    for(i in 1:length(select_F7)){ 
      variety_i = selectOP(
        select_F7[[i]],
        nInd=500,
        nSeeds=10,
        probSelf = selfing,###
        pollenControl = FALSE,
        use = "rand")
      variety_i = setPheno(variety_i, H2=H2F6) ############not done but I believe needed for selectInd function
      variety[[i]] = selectInd(variety_i, nInd = 5000, use="rand")
    }
  }
}




# Creating empty vectors to store genetic values

nYears = 20


output = data.frame(year = 1:nYears,
                    meanGF3 = numeric(nYears),
                    varGF3 = numeric(nYears),
                    meanGvariety = numeric(nYears),
                    varGvariety = numeric(nYears),
                    meanGParents =numeric(nYears),
                    varGParents=numeric(nYears),
                    varAParents=numeric(nYears),
                    varDParents=numeric(nYears),
                    accF3 = numeric(nYears),
                    accF4 = numeric(nYears),
                    accF5 = numeric(nYears),
                    accF6 = numeric(nYears))

for (year in 1:nYears) {
  
  #####################Year 7:  Create variety seed (F8) 
  ###########################################
  
  variety = vector("list",length(select_F7)) # Keep lines separate
  for(i in 1:length(select_F7)){ 
    variety_i = selectOP(
      select_F7[[i]],
      nInd=500,
      nSeeds=10,
      probSelf = selfing,###
      pollenControl = FALSE,
      use = "rand")
    variety[[i]] = selectInd(variety_i, nInd = 5000, use="rand")
  }
  
  outputvariety =mergePops(variety)
  
  ####################################################
  ###Year 6: Phenotype F7 in Svalöv/Sweden and select two lines that will be developed into varieties.
  ########Create parents for crossing
  ###################################
  
  
  ############Select top 2 lines
  F7_fampheno=data.frame(1:length(F7),NA)
  for(i in 1:length(F7)) {  
    
    F7_fampheno_i = mean(F7[[i]]@pheno) 
    F7_fampheno[i,2]=  F7_fampheno_i 
  }
  
  top2_F7 <-top_n(F7_fampheno, nF7sel, NA.)
  
  
  select_F7=F7[top2_F7[,1]]
  
  outputvariety =mergePops(F7)
  
  ############################
  #create parents##################
  Parents = c(select_F5, F7)### F7 are seeds produced by F6 in Chile
  
  ##Select two parental individuals from each F5 and F6 accession randomly
  Parentssel = vector("list",length(Parents)) # Keep crosses seperate
  for(i in 1:length(Parents)){ 
    Parentssel[[i]]<-selectInd(Parents[[i]], nInd = 2, use = "rand")
  }
  Parentssel=mergePops(Parentssel)
  
  ######################## 
  ######Year 5: Multiply seeds from lines selected in F5 and trial F6 in Svalöv 
  
  
  ###Multiply F5 in Chile 
  F5_Chile = vector("list",length(select_F5)) # Keep lines separate
  for(i in 1:length(select_F5)){ 
    F5Chile_i = selectOP(
      select_F5[[i]],
      nInd=500,
      nSeeds=10,
      probSelf = selfing,###
      pollenControl = FALSE,
      use = "rand")
    F5Chile_i = setPheno(F5Chile_i, H2=H2F6) ############Refers to phenotyping in Svalöv/Sweden at F6 stage
    F5_Chile[[i]] = selectInd(F5Chile_i, nInd = 5000, use="rand")
  }
  
  output_F6 =mergePops(F5_Chile)
  
  ######Select 10 lines from F6 based on trial in Svalov
  F6_fampheno=data.frame(1:length(F5_Chile),NA)
  for(i in 1:length(F5_Chile)) {  
    
    F6_fampheno_i = mean(F5_Chile[[i]]@pheno) 
    F6_fampheno[i,2]= F6_fampheno_i
  }
  
  top10_F6 <-top_n(F6_fampheno, 10, NA.)
  
  
  select_F6=F5_Chile[ top10_F6[,1]]
  
  ###Multiply seeds for F7  and phenotype in multilocation trials
  F7 = vector("list",length(select_F6)) # Keep lines separate
  for(i in 1:length(select_F6)){ 
    F7_i = selectOP(
      select_F6[[i]],
      nInd=500,
      nSeeds=10,
      probSelf = selfing,###
      pollenControl = FALSE,
      use = "rand")
    F7_i = setPheno(F7_i, H2=H2F7) ############Refers to phenotyping in Svalöv/Sweden at F7 stage
    F7[[i]] = selectInd(F7_i, nInd = 5000, use="rand")
  }
  
  
  #######################Year 4: Multiplication of seeds from F4 in Chile and F5 trial in Svalöv/Sweden       ##################
  ######################## 
  
  ####Multiply F4 seeds in Chile to create F5
  
  F4_Chile = vector("list",length(select_F4)) # Keep lines separate
  for(i in 1:length(select_F4)){ 
    F4Chile_i = selectOP(
      select_F4[[i]],
      nInd=200, ####Assumption this many seeds are held back in Chile from the F3 generation to be used for seed multiplication
      nSeeds=10,
      probSelf = selfing,###
      pollenControl = FALSE,
      use = "rand")
    F4Chile_i = setPheno(F4Chile_i, H2=H2F5) ##########Refers to phenotyping in Svalöv/Sweden of F5s
    F4_Chile[[i]] = selectInd(F4Chile_i, nInd = 2000, use="rand")
  }
  
  output_F5 =mergePops(F4_Chile)
  ######Select 25 lines from F5 based on trial in Svalov
  
  ############Select top 25 lines
  F5_fampheno=data.frame(1:length(F4_Chile),NA)
  for(i in 1:length(F4_Chile)) {  
    
    F5_fampheno_i = mean(F4_Chile[[i]]@pheno) 
    F5_fampheno[i,2]= F5_fampheno_i
  }
  
  top25_F5 <-top_n(F5_fampheno, 25, NA.)
  
  
  select_F5=F4_Chile[ top25_F5[,1]]
  
  
  ################################################################
  ###############################Year 3: F3 selection in Chile and trial of F4 in Svalöv
  
  
  ###### Select among F3 to create the F4 to be trialled in Svalov. 
  F3_merged=mergePops(F3)
  
  F3_select = selectFam(
    F3_merged,
    nFam=300,
    trait = 1,
    use = "pheno",
    famType = "F")
  
  ##Split population by half-sib and collect within-family open-pollinated seeds for each "line" in Chile
  
  motherSval_F4= unique(F3_select@mother)
  F4_split=vector("list", length(motherSval_F4))
  
  for(i in 1:length(motherSval_F4)) {  #Loop over maternal parent
    motherpop= subset(F3_select, F3_select@mother==paste(motherSval_F4[i]))
    F4_split_i = selectOP(
      pop=motherpop,
      nInd=40,
      nSeeds=40,
      probSelf = 0.95,###
      pollenControl = FALSE,
      use = "rand") 
    F4_split_i = setPheno(F4_split_i, H2=H2F4) ###Corresponds to the phenotype of the F4 in Svalov
    F4_split[[i]] = selectInd(F4_split_i, nInd = 1600)
  } 
  
  output_F4 =mergePops(F4_split)
  ###Select in F4 for seed multiplication in Chile - assumed  seed  used for F4 Svalöv and F4 Chile is same
  
  ############Select top 100 lines
  F4_fampheno=data.frame(1:length(F4_split),NA)
  for(i in 1:length(F4_split)) {  
    
    F4_fampheno_i = mean(F4_split[[i]]@pheno) 
    F4_fampheno[i,2]= F4_fampheno_i
  }
  
  top100_F4 <-top_n(F4_fampheno, 100, NA. )
  
  select_F4=F4_split[top100_F4[,1]]
  
  ########################################
  ###########Year 2: F2 planted to create seed for F3
  
  # Create F3
  F3 = vector("list",nCrosses) # Keep crosses separate
  for(i in 1:nCrosses){ # Loop over crosses
    F3_i = selectOP(
      F2[[i]],
      nInd=20,
      nSeeds=40,
      probSelf = 0.95,###
      pollenControl = FALSE,
      use = "pheno")
    F3_i = setPheno(F3_i, H2=H2F3)
    F3[[i]] = selectInd(F3_i, nInd = 800, use="rand")
  }
  ####################################################
  #######################Year 1: Crosses, F1 and create seed for F2
  
  
  F1 = vector("list",nCrosses) # Keep crosses seperate
  for(i in 1:nCrosses){ # Loop over crosses
    F1[[i]] = randCross(Parentssel, nCrosses=1, nProgeny=nind)
  }
  
  F2 = vector("list",nCrosses) # Keep crosses separate
  for(i in 1:nCrosses){ # Loop over crosses
    F2_i = selectOP(
      F1[[i]],
      nInd=nind,
      nSeeds=16,
      probSelf = 0.95,###
      pollenControl = FALSE,
      use = "rand")
    F2_i = setPheno(F2_i, H2=H2F2)
    F2[[i]] = selectInd(F2_i, nInd = 200, use="rand")
  }
  
  
  # Report results
  
  output$meanGF3[year] = meanG(F3_merged)
  output$varGF3[year] = varG(F3_merged)
  output$meanGvariety[year] = meanG(outputvariety)
  output$varGvariety[year] = varG(outputvariety)
  output$meanGParents[year]=meanG(Parentssel)
  output$varGParents[year]=varG(Parentssel)
  output$varAParents[year]=varA(Parentssel)
  output$varDParents[year]=varD(Parentssel)
  output$accF3[year] = cor(F3_merged@gv, F3_merged@pheno)
  output$accF4[year] = cor(output_F4@gv, output_F4@pheno)
  output$accF5[year] = cor(output_F5@gv, output_F5@pheno)
  output$accF6[year] = cor(output_F6@gv, output_F6@pheno)
  start_pop_1 =Parentssel[year==20]
  
}


##Phenotype starting pop rigorously
start_pop_1= setPheno(start_pop_1, H2=0.2, reps=6)
# Save burn-in to load later for subsequent scenarios
save(start_pop_1,file="start_pop_ADG_DD1.RData")

write.csv(output, "Burninperiod_output_DD1.csv")



