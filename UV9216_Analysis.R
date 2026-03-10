####
#### Packages and Data Prep ####
####
setwd() # IMPORTANT! SET WORKING DIRECTORY WHERE DATASET IS

# Run this to a) Download the packages if not done and b) Open them
required_packages <- c("metaforest", "caret", "readxl", "dplyr",
                       "ggplot2", "robumeta", "stringr")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}

# Data
dat_raw <- read_xls("Data.xls", sheet = "Data - EI only")

# Computing effect sizes and sampling variance
dat_raw <- dat_raw %>%
  mutate(
    yi  = EI_GPA_cor2,        # our effect size
    vi  = (1 - yi^2)^2 / (N - 1)  # sampling variance
  )


# Sanity check: any impossible values?
cat("Effect sizes out of [-1, 1]:", sum(abs(dat_raw$yi) > 1, na.rm = TRUE), "\n") # One correlation outside -1, +1
cat("Negative or zero N:", sum(dat_raw$N <= 1, na.rm = TRUE), "\n")
cat("Missing yi:", sum(is.na(dat_raw$yi)), "\n")
cat("Missing N:", sum(is.na(dat_raw$N)), "\n")

####
#### Outlier removal ####
####

# MacCann et al. removed 30 outlier effect sizes using the Hoaglin &
# Iglewicz (1987) labeling rule with the more conservative multiplier g = 2.2.
# We replicate this step on cor2. The rule flags values beyond
# Q1 - g*IQR and Q3 + g*IQR as outliers.

g <- 2.2
Q1 <- quantile(dat_raw$yi, 0.25, na.rm = TRUE)
Q3 <- quantile(dat_raw$yi, 0.75, na.rm = TRUE)
IQR_val <- Q3 - Q1
lower_bound <- Q1 - g * IQR_val
upper_bound <- Q3 + g * IQR_val

cat("Outlier bounds: [", round(lower_bound, 3), ",", round(upper_bound, 3), "]\n")
n_outliers <- sum(dat_raw$yi < lower_bound | dat_raw$yi > upper_bound,
                  na.rm = TRUE)
cat("Outliers identified:", n_outliers, "\n")

dat_clean <- dat_raw %>%
  filter(!is.na(yi), !is.na(N), N > 1,
         yi >= lower_bound, yi <= upper_bound,
         vi > 0, is.finite(vi))

cat("After outlier removal:", nrow(dat_clean), "effect sizes,",
    length(unique(dat_clean$PaperNo)), "papers\n")

####
#### Moderators Prep ####
####
dat <- dat_clean %>%
  mutate(
    Ability_EI   = as.integer(Stream == 1), # We code the three EI streams as dummy variables, 
    SelfRated_EI = as.integer(Stream == 2), # to allow metaforest to check how they differenlty predict
    Mixed_EI     = as.integer(Stream == 3), # academic performance, enabling comparison to MacCann analysis
    Ed_Level_Num = as.numeric(Ed_Level_Num), # Educational level (3 levels)
    GPA_type_num = factor(GPA_type_num), # GPA type (self reported vs standardized)
    GPA_gen_sci_art_num = factor(GPA_gen_sci_art_num), # GPA Area 
    Pub_type_Num = factor(Pub_type_Num), # Publication type (4 levels)
    PcFemale = ifelse(is.na(PcFemale), # Gender. Note - Imputing missing values with median.
                      median(PcFemale, na.rm = TRUE), # Continuous because it's the % of female in a study
                      PcFemale),
    Mean_age = as.numeric(Mean_age), # Age
    Year = as.numeric(str_extract(as.character(Year), "\\d{4}")), # Year of publication. Taking the first 4 string (year)
    GPA_source = factor(GPA_source), # Self reported vs Official record
    EI_reliability = as.numeric(EI_reliability), # Reliability of EI. In the study used as a correction, here as a predictor
  )

# Summary
# Stream dummies (binary)
cat("\nEI Streams (binary dummies):\n")
cat("  Ability EI:    n =", sum(dat$Ability_EI), 
    "(", round(mean(dat$Ability_EI)*100, 1), "%)\n")
cat("  Self-rated EI: n =", sum(dat$SelfRated_EI),
    "(", round(mean(dat$SelfRated_EI)*100, 1), "%)\n")
cat("  Mixed EI:      n =", sum(dat$Mixed_EI),
    "(", round(mean(dat$Mixed_EI)*100, 1), "%)\n")

# Education level
cat("\nEducation level (0=Primary, 1=Secondary, 2=Tertiary, 3=Mixed):\n")
print(table(dat$Ed_Level_Num, useNA = "ifany"))

# GPA type
cat("\nGPA type:\n")
print(table(dat$GPA_type_num, useNA = "ifany"))

# GPA subject
cat("\nGPA subject area:\n")
print(table(dat$GPA_gen_sci_art_num, useNA = "ifany"))

# Publication type
cat("\nPublication type:\n")
print(table(dat$Pub_type_Num, useNA = "ifany"))

# GPA source
cat("\nGPA source:\n")
print(table(dat$GPA_source, useNA = "ifany"))

# Continuous moderators
cat("\nContinuous moderators:\n")
cont_vars <- c("PcFemale", "Mean_age", "Year", "EI_reliability")
for (v in cont_vars) {
  x <- dat[[v]]
  cat(sprintf("  %-20s M = %.2f, SD = %.2f, range = [%.2f, %.2f], missing = %d\n",
              v,
              mean(x, na.rm=TRUE),
              sd(x, na.rm=TRUE),
              min(x, na.rm=TRUE),
              max(x, na.rm=TRUE),
              sum(is.na(x))))
}


####
#### Dataset for analysis ####
####
dat_mf <- dat %>%
  select(
    PaperNo,               # study cluster identifier for clustered bootstrap
    yi,                    # effect size (cor2)
    vi,                    # sampling variance
    Ability_EI,            # Stream 1 dummy
    SelfRated_EI,          # Stream 2 dummy
    Mixed_EI,              # Stream 3 dummy
    Ed_Level_Num,          # education level (0=Primary, 1=Secondary, 2=Tertiary, 3=Mixed)
    GPA_type_num,          # achievement type (0=Course grade, 1=Standardised test)
    GPA_gen_sci_art_num,   # subject area collapsed (0=General, 1=Maths/Sciences, 2=Humanities)
    Pub_type_Num,          # publication type (0=Journal, 1=Dissertation, 3=Conference, 4=Unpublished)
    PcFemale,              # % female in sample (continuous, 0-100)
    Mean_age,              # mean age of sample (continuous)
    Year,                  # publication year (continuous, 1997-2019)
    GPA_source,            # GPA source (records vs. self-report)
    EI_reliability         # reliability of EI instrument (continuous, Cronbach's alpha)
  ) %>%
  as.data.frame()

# Handle missing values. Missing proportions are negligible, median/mode used

# GPA_gen_sci_art_num: 16 missing — impute with mode (most frequent category)
mode_subj <- names(which.max(table(dat_mf$GPA_gen_sci_art_num)))
dat_mf$GPA_gen_sci_art_num[is.na(dat_mf$GPA_gen_sci_art_num)] <- mode_subj

# Mean_age: 1 missing — impute with median
# Note: MacCann already imputed most missing ages at source (using education-
# level means), so this residual missing is a genuine outlier case.
dat_mf$Mean_age[is.na(dat_mf$Mean_age)] <- median(dat_mf$Mean_age, na.rm = TRUE)

# GPA_source: 1 missing — impute with mode
mode_source <- names(which.max(table(dat_mf$GPA_source)))
dat_mf$GPA_source[is.na(dat_mf$GPA_source)] <- mode_source

cat("\nFinal analysis dataset:", nrow(dat_mf), "effect sizes,",
    length(unique(dat_mf$PaperNo)), "papers\n")
cat("Missing values after imputation:\n")
print(colSums(is.na(dat_mf)))

####    
#### Convergence of Metaforest ####
####
set.seed(2025)
check_conv <- MetaForest(
  yi ~.,
  data         = dat_mf,
  vi           = "vi",
  study        = "PaperNo",
  whichweights = "random",
  num.trees    = 20000
)

# Plot convergence: look for where the cumulative OOB MSE curve flattens.
# Record this number as conv_trees for use in subsequent steps.
plot(check_conv) # 5000 trees


####
#### Pre - selection ####
####

# Pre - selection of variables done as per Van Lissa (2023) reccomended workflow
set.seed(2025)
mf_rep <- MetaForest(
  yi ~ .,
  data         = dat_mf,
  study        = "PaperNo",
  whichweights = "random",
  num.trees    = 5000
)

preselected <- preselect(mf_rep, 
                         replications = 100, 
                         algorithm = "recursive")


plot(preselected) # Distribution of variable importance across 100 replications

# Rename for clarity and replot
var_labels <- c(
  "Mean_age"            = "Mean Age",
  "Ability_EI"          = "Ability EI",
  "Ed_Level_Num"        = "Education Level",
  "PcFemale"            = "Gender Composition (% Female)",
  "GPA_gen_sci_art_num" = "Subject Area",
  "EI_reliability"      = "EI Instrument Reliability",
  "Mixed_EI"            = "Mixed EI",
  "SelfRated_EI"        = "Self-Rated EI",
  "Pub_type_Num"        = "Publication Type",
  "GPA_source"          = "GPA Source",
  "GPA_type_num"        = "Achievement Type",
  "Year"                = "Publication Year"
)
preselected_renamed <- preselected
colnames(preselected_renamed$selected) <- var_labels[colnames(preselected_renamed$selected)]
plot(preselected_renamed)



retained_vars <- preselect_vars(preselected, cutoff = 0.5) # Retain variables with positive importance in more than 50% of the replications
print(retained_vars)


####
#### Cross - validation ###
####
# Define cross-validation strategy (clustered by paper)
grouped_cv <- trainControl(method = "cv",
                           index  = groupKFold(dat_mf$PaperNo, k = 10)
)

# Define hyperparameter grid
# whichweights: how effect sizes are weighted during training
# mtry: number of moderators available at each split
# min.node.size: minimum effect sizes per terminal node
tuning_grid <- expand.grid(
  whichweights  = c("random", "fixed", "unif"),
  mtry          = 2:6,
  min.node.size = 2:6
)

# Prepare data: retained moderators + clustering variable + vi only
X <- dat_mf[, c("PaperNo", "vi", retained_vars)]

# Train with cross-validation
mf_cv <- train(
  y             = dat_mf$yi,
  x             = X,
  method        = ModelInfo_mf(),
  trControl     = grouped_cv,
  tuneGrid      = tuning_grid,
  num.trees     = 5000,
  study         = "PaperNo"
)

# So, which hyperparameters?
print(mf_cv$bestTune)

r2_cv <- mf_cv$results$Rsquared[which.min(mf_cv$results$RMSE)] # Cross-validated R²(the lowest across all models ran)
r2_cv
final_mf <- mf_cv$finalModel # Extract final model 
r2_oob <- final_mf$forest$r.squared # R2 oob for final model
r2_oob

plot(final_mf) # Has the final model converged?
VarImpPlot(final_mf) # Final variable selection and importance in predicting the outcome


