# Reproducibility
set.seed(123)

# Libraries
library(dplyr) # case_when and others
library(ggplot2)
library(tidyr) # pivot
library(glmnet)
library(mvtnorm) # rmvnorm
library(tmle)
library(SuperLearner)
library(AIPW)
library(grf)

source("estimators.R")
source("generate_data_models.R")

results.linear <- data.frame("sample.size" = c(),
                             "estimate" = c(),
                             "estimator" = c(),
                             "subset" = c(),
                             "simulation" = c(),
                             "cross-fitting" = c(),
                             "independence" = c(),
                             "nuisance" = c())

different_subset_tested <- c("outcome.and.instruments",
                             "outcome.wo.instruments",
                             "smart",
                             "minimal.set")

for (sample.size in c(300, 1000, 3000, 9000, 30000)){
  print(paste0("Starting sample size ", sample.size))
  for (i in 1:30){
    for (independence in c(TRUE, FALSE)){
      
      # generate a simulation
      a_simulation <- generate_simulation_linear(n_obs = sample.size, independent_covariate = independence, constant_cate = FALSE)
      
      # choose subset
      for (method in different_subset_tested){
        if (method == "outcome.and.instruments"){
          X_treatment <- paste0("X.", 1:10)
          X_outcome <- paste0("X.", 1:10)
        } else if (method == "outcome.wo.instruments"){
          X_treatment <- paste0("X.", 4:10)
          X_outcome <- paste0("X.", 4:10)
        } else if (method == "smart"){
          X_treatment <- paste0("X.", 4:7)
          X_outcome <- paste0("X.", 4:10)
        } else if (method == "minimal.set"){
          X_treatment <- paste0("X.", 4:7)
          X_outcome <- paste0("X.", 4:7)
        } else {
          stop("error in subset.")
        }
        
        custom_aipw_1 <- aipw_splines(X_treatment, X_outcome, dataframe = a_simulation, n.folds = 1)
        custom_aipw_2 <- aipw_splines(X_treatment, X_outcome, dataframe = a_simulation, n.folds = 2)
        
        
        new.row <- data.frame("sample.size" = rep(sample.size, 6),
                              "estimate" = c(custom_aipw_1["ipw"],
                                             custom_aipw_1["t.learner"],
                                             custom_aipw_1["aipw"],
                                             custom_aipw_2["ipw"],
                                             custom_aipw_2["t.learner"],
                                             custom_aipw_2["aipw"]),
                              "estimator" = rep(c("ipw",
                                                  "t-learner",
                                                  "aipw"),2),
                              "subset" = rep(method, 6),
                              "simulation" = rep("linear.hte", 6),
                              "cross-fitting" = c(1,1,1,2,2,2),
                              "independence" = rep(independence, 6),
                              "nuisance" = rep("splines", 6))
        
        results.linear <- rbind(results.linear, new.row)
        
      }
    }
  }
}

results.linear$sample.size <- as.factor(results.linear$sample.size)

write.csv(x=results.linear, file="./data/2021-11-01-linear-hte-splines-nuisance.csv")