# Reproducibility
set.seed(123)

# Libraries
library(dplyr) # case_when and others
library(ggplot2)
library(tidyr) # pivot
library(grf)
library(glmnet)
library(ranger) # efficient forest
library(splines) # function bs() for splines
library(mvtnorm) # rmvnorm
library(tmle)
library(SuperLearner)
library(AIPW)

source("estimators.R")
source("generate_data_models.R")

results.linear <- data.frame("sample.size" = c(),
                             "estimate" = c(),
                             "estimator" = c(),
                             "subset" = c(),
                             "nuisance" = c(),
                             "term.A" = c(), 
                             "term.B" = c(), 
                             "term.C" = c(),
                             "term.D" = c(), 
                             "term.E" = c(), 
                             "term.F" = c())

different_subset_tested <- c("extended",
                             "smart",
                             "minimal")

for (sample.size in c(100, 300, 1000, 3000, 10000)){
  print(paste0("Starting sample size ", sample.size))
  for (i in 1:50){
    
    # generate a simulation
    a_simulation <- generate_simulation_wager_nie(n = sample.size, setup = "A", all_covariates_output = TRUE)
    
    # choose subset
    for (method in different_subset_tested){
      if (method == "extended"){
        X_treatment <- paste0("X.", 1:6)
        X_outcome <- paste0("X.", 1:6)
        
        custom_ipw <- ipw_forest(covariates_names = X_treatment, 
                                 dataframe = a_simulation,
                                 min.node.size.if.forest = 1,
                                 return.decomposition = TRUE)
        
        new.row <- data.frame("sample.size" = sample.size,
                              "estimate" = custom_ipw,
                              "estimator" = "ipw",
                              "subset" = method,
                              "nuisance" = "forest",
                              "term.A" = NA, 
                              "term.B" = NA, 
                              "term.C" = NA,
                              "term.D" = NA, 
                              "term.E" = NA, 
                              "term.F" = NA)
        
        results.linear <- rbind(results.linear, new.row)
      } else if (method == "smart"){
        X_treatment <- paste0("X.", 1:2)
        X_outcome <- paste0("X.", 1:6)
      } else if (method == "minimal"){
        X_treatment <- paste0("X.", 1:2)
        X_outcome <- paste0("X.", 1:2)
        
        custom_ipw <- ipw_forest(covariates_names = paste0("X.", 1:2), 
                                 dataframe = a_simulation,
                                 min.node.size.if.forest = 1,
                                 return.decomposition = TRUE)
        
        new.row <- data.frame("sample.size" = sample.size,
                              "estimate" = custom_ipw,
                              "estimator" = "ipw",
                              "subset" = method,
                              "nuisance" = "forest",
                              "term.A" = NA, 
                              "term.B" = NA, 
                              "term.C" = NA,
                              "term.D" = NA, 
                              "term.E" = NA, 
                              "term.F" = NA)
        
        results.linear <- rbind(results.linear, new.row)
        
      } else {
        stop("error in subset.")
      }
      
      custom_aipw <- aipw_forest(X_treatment, 
                                 X_outcome, 
                                 dataframe = a_simulation,
                                 min.node.size.if.forest = 1,
                                 n.folds = 5,
                                 return.decomposition = TRUE,
                                 with.weights = TRUE)
      
      new.row <- data.frame("sample.size" = rep(sample.size, 3),
                            "estimate" = c(custom_aipw["aipw"], custom_aipw["t.learner"], custom_aipw["ipw"]),
                            "estimator" = c("aipw", "t-learner", "ipw.cross.fit"),
                            "subset" = rep(method, 3),
                            "nuisance" = rep("forest", 3),
                            "term.A" = c(custom_aipw["term.A"], NA, NA), 
                            "term.B" = c(custom_aipw["term.B"], NA, NA), 
                            "term.C" = c(custom_aipw["term.C"], NA, NA),
                            "term.D" = c(custom_aipw["term.D"], NA, NA),
                            "term.E" = c(custom_aipw["term.E"], NA, NA), 
                            "term.F" = c(custom_aipw["term.F"], NA, NA))
      results.linear <- rbind(results.linear, new.row)
    }
  }
}

write.csv(x=results.linear, file="./data/A.with.weights.csv")
