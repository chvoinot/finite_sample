# Reproducibility
set.seed(123)

# Libraries
library(dplyr) # case_when and others
library(tidyr) # pivot
library(grf)
source("estimators.R")
source("generate_data_models.R")

# load data
PATH <- "../data/cohort_2019_imputed_after_2_composite_covariate.RData"
load(PATH)

data_depp <- data_depp[data_depp$Taille_Classe > 8,]
data_depp <- data_depp[data_depp$Taille_Classe < 13 | data_depp$Taille_Classe > 18,]
data_depp$Treatment <- ifelse(data_depp$Taille_Classe < 13, 1, 0)
print(nrow(data_depp))

# Categ etab into one hot encoder
data_depp$REP <- ifelse(data_depp$Categ_Etab_CP == "REP", 1, 0)
data_depp$REPp <- ifelse(data_depp$Categ_Etab_CP == "REP+", 1, 0)
data_depp$Public <- ifelse(data_depp$Categ_Etab_CP == "Public", 1, 0)
data_depp$Private <- ifelse(data_depp$Categ_Etab_CP == "Private", 1, 0)

categ <- c("REP", "REPp", "Public", "Private", "IPS_Etab_CP")
minimal_set <- c(categ)
extended_set <- c(minimal_set, "Age_CP", "Sexe_Num") #"T1_Language" "T1_Math" are post-treatment covariates

results <- data_frame("estimator" = c(),
                      "estimate" = c(),
                      "sample.size" = c(),
                      "extended.set" = c(),
                      "nuisance" = c())

for (sample.size in c(300, 1000, 3000, 10000, 30000)){
  print(paste0("starting sample size ", str(sample.size)))
  for (i in 1:30){
    if(i == 5){
      print("starting 5")
    } else if (i == 15){
      print("starting 15")
    }
    workind_df <- data_depp[sample(nrow(data_depp), sample.size), ]
    estimate.with.minimal.set <- aipw_forest(covariates_names_vector_treatment = minimal_set,
                                             covariates_names_vector_outcome = minimal_set,
                                             dataframe = workind_df,
                                             outcome_name = "T3_Language",
                                             treatment_name = "Treatment",
                                             n.folds = 2,
                                             min.node.size.if.forest = 3)
    estimate.with.minimal.set.linear <- aipw_linear(covariates_names_vector_treatment = minimal_set,
                                             covariates_names_vector_outcome = minimal_set,
                                             dataframe = workind_df,
                                             outcome_name = "T3_Language",
                                             treatment_name = "Treatment",
                                             n.folds = 2)
    estimate.with.extended.set <- aipw_forest(covariates_names_vector_treatment = minimal_set,
                                              covariates_names_vector_outcome = extended_set,
                                              dataframe = workind_df,
                                              outcome_name = "T3_Language",
                                              treatment_name = "Treatment",
                                              n.folds = 2,
                                              min.node.size.if.forest = 3)
    estimate.with.extended.set.linear <- aipw_linear(covariates_names_vector_treatment = minimal_set,
                                              covariates_names_vector_outcome = extended_set,
                                              dataframe = workind_df,
                                              outcome_name = "T3_Language",
                                              treatment_name = "Treatment",
                                              n.folds = 2)

  
    new_row <- data.frame("estimator" = rep(c("ipw", "t-learner", "aipw"),4),
                        "estimate" = c(estimate.with.minimal.set["ipw"],
                                       estimate.with.minimal.set["t.learner"],
                                       estimate.with.minimal.set["aipw"],
                                       estimate.with.extended.set["ipw"],
                                       estimate.with.extended.set["t.learner"],
                                       estimate.with.extended.set["aipw"],
                                       estimate.with.minimal.set.linear["ipw"],
                                       estimate.with.minimal.set.linear["t.learner"],
                                       estimate.with.minimal.set.linear["aipw"],
                                       estimate.with.extended.set.linear["ipw"],
                                       estimate.with.extended.set.linear["t.learner"],
                                       estimate.with.extended.set.linear["aipw"]),
                        "sample.size" = rep(sample.size,12),
                        "extended.set" = rep(c("no", "no", "no", "yes", "yes", "yes"),2),
                        "nuisance" = c(rep("forest",6), rep("linear", 6)))
  
    results <- rbind(results, new_row)
  }
}



write.csv(x=results, file="./data/semi-synthetic-education.csv")