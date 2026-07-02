# Project 3: Sex-specific differences in HFpEF
#
# Multi-omic (proteomics + phenomics) patient stratification in heart failure
# with preserved ejection fraction (HFpEF), analyzed separately by sex.
# Pipeline: similarity network fusion (SNF) + spectral clustering to find
# subgroups -> Cox/Kaplan-Meier survival validation -> XGBoost to project
# cluster labels onto an independent cohort -> RNA differential expression
# and pathway enrichment on the resulting clusters.
#
# NOTE ON DATA ACCESS
# This script operates on patient-level clinical, proteomic, and RNA-seq
# data from the MYOVASC and Gutenberg Health Study (GHS) cohorts. Because
# this is identifiable/sensitive health data, none of it is included in
# this repository and the raw data loading calls below (read_SQL_MyoVasc_BL,
# read_SQL_A6_data, etc.) depend on internal database access that is not
# publicly available. This script is shared to document the analysis
# methodology, not to be run end-to-end without institutional data access.
#
# ---- Configuration ------------------------------------------------------
# Set these paths/credentials for your own environment (e.g. via a local
# ---------------------------------------------------------------------------

setwd(project_dir)
source("Imputa.R")
#Libraries
library(SNFtool)
library(survival)
library(survminer)
library(survcomp)
library(flexclust)
library(ggplot2)
library(Rtsne)
library(ggplot2)
library(rstatix)
library(stringr)
library(rstatix)
library(janitor)
library(ggsurvfit)
library(readxl)
library(caret)
library(Hmisc)
# load baseline dataset MYOVASC

dall1 = read_SQL_MyoVasc_BL(proteins=T, lipids = T, user = db_user, ATC_7_digits = 'B01AC06')
baseline = read.csv2(file.path(data_dir, "file.csv"))
setwd(project_dir)
id_lid = data.frame(dall1$v11_sid01,dall1$lid)   
colnames(id_lid) = c("ID", "LID")

# Data loading GHS
Sys.setlocale("LC_ALL", "GER") # switching to German locale

dattr_ghs = attributes(ghs_bas)
whf_ghs <- read.csv(file.path(ghs_sql_dir, "whf_ghs.csv"), sep = ";")
ghs_bas = merge(ghs_bas, whf_ghs, by.x = "v11_sid01", by.y = "v11_sid01")

# Proteins
dattr = attributes(dall1)
proteins = baseline[,colnames(baseline) %in% dattr_ghs$vars_prot]
proteins = data.frame(v11_sid01 = baseline$v11_sid01,  proteins)
rownames(proteins) = proteins$v11_sid01
proteins = proteins[,-1]

#Phenomics
cardiac = c("sbp", "dbp", "ee", "ef", "lvm", "rwt")
blood = c("tni",  "natrium", "kalium", "chlorid", "calcium",  "v34_lbl07", "v34_lbl08", 
                             "got",  "ggt", "v34_lbl13", "v34_lbl14", "v34_lbl15", "v34_lbl40", "v34_lbl41", "crp_ln", "v34_lbl50", "glucose",
                             "v34_lbl54", "v34_lbl55", "v34_lbl56", "v34_lbl57",   "v34_lbl67", "egfr", 
                             "v34_lbl18",  "v34_lbl19", "v34_lbl20", "v34_lbl21", "v34_lbl22", "v34_lbl23", "v34_lbl24", "v34_lbl25",
                             "v34_lbl26", "v34_lbl81", "v34_lbl27", "v34_lbl28","v34_lbl29","v34_lbl30","v34_lbl31")
blood_variables_names  <- c(  "Troponin I",   "Na", "K", "Cl", "Ca",  "Creatinine", "Urea", "GOT",  "GGT", "Alkalinephosphatase", "Bilirubin", "Albumin",
                              "LDH", "Lipase","CRP", "HbA1c", "Glucose", "Cholesterine", "Triglycerides", "HDL", "LDL","TSH", "eGFR",
                              "Erythrocytes",   "Leukocytes", "Hemoglobin", "Hematocrit", "MCV", "MCH", "MCHC", "EVB", "Thrombocytes", "MTV",
                              "Neutrophils", "Lymphocytes", "Monocytes", "Eosinophils", "Basophils")

anthro = c("age", "weight", "height", "waist","hip","wthr")

phenomics = dall1[, colnames(dall1) %in% c(cardiac, blood, anthro)]
rownames(phenomics) = dall1$v11_sid01

var = c(dattr$vars_prot, cardiac, blood, anthro)
var_names = c(dattr$nams_prot, cardiac, blood_variables_names, anthro)
var_data_match = data.frame(var = var, names = var_names)


# GHS Females
colnames(ghs_bas)[colnames(ghs_bas) == "whtr"] = "wthr"
colnames(ghs_bas)[colnames(ghs_bas) == "v34_lbl01"] = "natrium"
colnames(ghs_bas)[colnames(ghs_bas) == "v34_lbl03"] = "chlorid"
colnames(ghs_bas)[colnames(ghs_bas) == "v34_lbl04"] = "calcium"

proteins_ghs = ghs_bas[,colnames(ghs_bas) %in% dattr_ghs$vars_prot]
rownames(proteins_ghs) = ghs_bas$v11_sid01
phenomics_ghs = ghs_bas[, colnames(ghs_bas) %in% c(cardiac, blood, anthro)]

X_test_df = cbind(proteins_ghs, phenomics_ghs)
X_test_df = X_test_df[ghs_bas$sex == "Women" & ghs_bas$hf_a3 == "HFPEF",]

#Intersection of omics (proteomics and phenomics)
panels_list = list(proteins, phenomics)
intersection <- intersect(rownames(panels_list[[1]]), rownames(panels_list[[2]]))
d = dall1[dall1$v11_sid01 %in% intersection,]
d = d[d$hf_a3 == "HFPEF",]
hfpef_id_female = d$v11_sid01[which(d$hf_a3 =="HFPEF" & d$p002 =="2")]
d = d[d$v11_sid01 %in% hfpef_id_female,]
intersection = intersection[intersection %in% hfpef_id_female]
panels_list_intersected = list()
for (i in 1:length(panels_list)){
  panels_list_intersected[[i]] =  panels_list[[i]][(rownames(panels_list[[i]]) %in% intersection),]
}
panels_list_intersected <- lapply(panels_list_intersected, function(x) x[order(rownames(x)), , drop = FALSE])

panels_list_intersected[[2]] <- apply(panels_list_intersected[[2]], 2, as.numeric)
rownames(panels_list_intersected[[2]]) = rownames(panels_list_intersected[[1]])

W = list()
Normalised_list = list()
train_means = list()
train_sds = list()
for (i in 1:length(panels_list_intersected)){  
  
  #if (data_integration == 'SNF'){
  # as suggested by authors apply data normalization
  
  data_imputed = imputa(panels_list_intersected[[i]])
  nzv = nearZeroVar(data_imputed, names = TRUE)
  print(length(nzv))
  
  if  (length(nzv) > 0) {
    data_imputed= data_imputed[,-which(colnames(data_imputed) %in% nzv)]
  }
  
  train_means[[i]] <- apply(data_imputed, 2, mean, na.rm = TRUE)
  train_sds[[i]]   <- apply(data_imputed, 2, sd, na.rm = TRUE)
  
  data_norm = scale(data_imputed, center = train_means[[i]], scale = train_sds[[i]])
  
  cc = cor(data_norm, use = "pairwise.complete.obs", method = "pearson")
  select_corr = caret::findCorrelation(cc, cutoff = 0.9, exact = FALSE)
  print(length(select_corr))
  
  if (length(select_corr)>0){
    data_norm_corr= data_norm[, -select_corr]
  }
 
  ## Calculate the pair-wise distance;
  dist_mat = (dist2(as.matrix(data_norm),as.matrix(data_norm)))^(1/2)
  W_mat = affinityMatrix(dist_mat, K = 20, sigma = 0.5) 
  W[[i]] = W_mat
  Normalised_list[[i]] = data_norm_corr
  print(dim(W[[i]]))
}

setwd(project_dir)

### SNF & Spectral Clustering
set.seed(256)
snf_females = SNF(W, 20, 20)
memberships_female = data.frame()
for (i in 2:10){
  memberships_female = rbind(memberships_female, SNFtool::spectralClustering(snf_females, K = i))
}
memberships_female = t(memberships_female)
colnames(memberships_female) = c("cl2", "cl3", "cl4", "cl5", "cl6", "cl7", "cl8", "cl9", "cl10")
rownames(memberships_female) = rownames(snf_females)
memberships_female = as.data.frame(memberships_female)

saveRDS(list(Normalised_list = Normalised_list, memberships_female = memberships_female), file.path(output_dir, "FEMALES.RDS"))

#### Clinical check 
d$memberships = as.factor(memberships_female$cl2)
cox_age <- coxph(Surv(whf_all_time,whf_all_event) ~  memberships + age , data =  d)
summary(cox_age)

# Unadjusted Kaplan-Meier
library(survival)
library(survminer)
sf_km <- survfit(Surv(whf_all_time,whf_all_event) ~ memberships, data = d)
#d$age_event = round(d$age+d$whf_all_time)
#sf_km = survfit(Surv(age_event,whf_all_event) ~ memberships, data = d)
female_cumulative = ggsurvplot(
  sf_km,
  data = d,
  fun = "event",
  conf.int = TRUE,
  pval = TRUE,
  risk.table = TRUE,
  cumevents = TRUE,
  xlim = c(0, 4),
  
  # Rename legend labels
  legend.labs = c("High risk cluster", "Low risk cluster"),
  
  legend = "right",
  title = "Females",
  xlab = "Time (years)",
  
  # Increase general font sizes
  fontsize = 6,                 # overall scaling (affects curves & tables)
  risk.table.fontsize = 5,      # risk table text
  cumevents.fontsize = 5,       # cumulative events table
  
  ggtheme = theme_bw() +
    theme(
      legend.box = "vertical",
      legend.direction = "vertical",
      
      # Title
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
      
      # Axis text
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      
      # Legend text
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 12)
    )
)

# Clinical characteristics
source(file.path(shared_functions_dir, "baseline13.r"))
setwd(shared_functions_dir)
suppressWarnings(eval(parse(file = "ade.fun.run.r")))

vars_1 <- c('age', 'sex', 'maggic')
nams_1 <- c('Age [y]', 'Sex (Women)', 'MAGGIC Score')
vars_2 <- c( vars_1, 'bmi',		     'hba1c',		  'fli_cat'	,	'nic',     'hyper',  "lvh_ord" , 'diab',     'dyslip',   'adipos', "ckd", "cvd", "chf", "cancer",  "plaq")
nams_2 <- c( nams_1, 'BMI [kg/m2]', 'HbA1c [%]',	'FLI'	,  'Smoking', 'Hypertension', "Hypertrophy type",
             'Diabetes', 'Dyslipidemia',  'Obesity', "Chronic Kidney Disease", "CVD", "Congestive heart failure", "Cancer", "Plaques (yes)")
vars_3 <- c( vars_2, "nyha_ord2" ,'afib',  "mi", 'stroke', 'cad', 'pad')
nams_3 <- c(nams_2,   "NYHA", 'AF' ,   'MI', 'Stroke', 'CAD', 'PAD')

vars_4 <- c( vars_3, "sbp", "dbp", "heart_rate_pre", "peak_vo2", "ee", "ef", "lvm", "rwt", "bnp" )
nams_4 <- c(nams_3, "SBP", "DBP", "Heart rate", "peak_vO2","E/E'","EF [%]",  "LVM", "RWT", 'NT-proBNP [pg/ml]')

vars_5 <- c(vars_4, "tni",  "natrium", "kalium", "chlorid", "calcium", "v34_lbl05", "v34_lbl06", "v34_lbl07", "v34_lbl08", 
            "got", "gpt", "ggt", "v34_lbl13", "v34_lbl14", "v34_lbl15", "v34_lbl40", "v34_lbl41", "crp_ln", "v34_lbl50", "glucose",
            "v34_lbl54", "v34_lbl55", "v34_lbl56", "v34_lbl57",   "v34_lbl67", "egfr")
nams_5 <- c( nams_4, "Troponin I",   "Na", "K", "Cl", "Ca", "Mg", "P", "Creatinine", "Urea", "GOT", "GPT", "GGT", "Alkalinephosphatase", "Bilirubin", "Albumin",
             "LDH", "Lipase","CRP", "HbA1c", "Glucose", "Cholesterine", "Triglycerides", "HDL", "LDL","TSH", "GFR")

vars_6 <- c(vars_5, "v34_lbl18",  "v34_lbl19", "v34_lbl20", "v34_lbl21", "v34_lbl22", "v34_lbl23", "v34_lbl24", "v34_lbl25",
            "v34_lbl26", "v34_lbl81", "v34_lbl27", "v34_lbl28","v34_lbl29","v34_lbl30","v34_lbl31")
nams_6 <- c( nams_5, "Erythrocytes",   "Leukocytes", "Hemoglobin", "Hematocrit", "MCV", "MCH", "MCHC", "EVB", "Thrombocytes", "MTV",
             "Neutrophils", "Lymphocytes", "Monocytes", "Eosinophils", "Basophils")

d$v34_lbl34 = as.numeric(d$v34_lbl34)
d$v34_lbl50 = as.numeric(d$v34_lbl50)
d$v34_lbl56 = as.numeric(d$v34_lbl56)
d$v34_lbl57 = as.numeric(d$v34_lbl57)
d$v34_lbl14 = as.numeric(d$v34_lbl14)
d$v34_lbl41 = as.numeric(d$v34_lbl41)
d$v34_lbl67 = as.numeric(d$v34_lbl67)
d$nyha_ord2 = as.factor(d$nyha_ord2)
d$egfr = as.numeric(d$egfr)
d$v34_lbl50 = as.numeric(d$v34_lbl50)
d$v34_lbl50 = as.numeric(d$v34_lbl50)
d$v34_lbl56 = as.numeric(d$v34_lbl56)
d$v34_lbl57 = as.numeric(d$v34_lbl57)
d$nyha_ord2 = as.factor(d$nyha_ord2)
d$afib_type <- d$myo029
d$afib_type[d$myo029 == 88 | d$myo029 == 99] <- NA
d$afib_type <- factor(d$afib_type, labels = c('paroxysmal', 'permanent', 'persistent'))
d$egfr = as.numeric(d$egfr)
d$fli_cat = ifelse(d$fli>=60, "yes", "no")
setwd(project_dir)
BT1_females <- baseline13(vars_6, nams_6, data = d, skeew=1, group = 'memberships', pv=T, style=2, minuniq = 3)

#XGBoost
X_train_df  = cbind(Normalised_list[[1]], Normalised_list[[2]])
X_test_df = X_test_df[,colnames(X_test_df) %in% colnames(X_train_df)]
cluster_train = memberships_female$cl2
library(xgboost)
library(caret)
library(pROC)
library(PRROC)

# 1) Align columns train/test (critical)
common_features <- intersect(colnames(X_train_df), colnames(X_test_df))
X_train_df <- X_train_df[, common_features, drop = FALSE]
X_test_df  <- X_test_df[,  common_features, drop = FALSE]

X_train <- data.matrix(X_train_df)
X_test  <- data.matrix(X_test_df)

# 2) Binary labels: cluster 1/2 -> y 0/1 (1 = cluster 2)
y_train <- ifelse(cluster_train == 2, 1, 0)
stopifnot(all(y_train %in% c(0,1)))

dtrain <- xgb.DMatrix(X_train, label = y_train)

# 3) Parameters (binary)
params <- list(
  booster = "gbtree",
  objective = "binary:logistic",
  eval_metric = c("auc", "aucpr"),
  eta = 0.05,
  max_depth = 4,
  min_child_weight = 1,
  subsample = 0.8,
  colsample_bytree = 0.8,
  lambda = 1,
  alpha = 0
)

# 4) CV to pick best nrounds
set.seed(256)
cv <- xgb.cv(
  params = params,
  data = dtrain,
  nrounds = 5000,
  nfold = 5,
  stratified = TRUE,
  early_stopping_rounds = 50,
  maximize = TRUE,
  verbose = 0
)

best_nrounds <- cv$best_iteration
if (is.null(best_nrounds) || !is.finite(best_nrounds) || best_nrounds < 1) {
  elog <- cv$evaluation_log
  best_nrounds <- if ("test_aucpr_mean" %in% names(elog)) which.max(elog$test_aucpr_mean) else which.max(elog$test_auc_mean)
}

# 5) Fit final model on full training set
model <- xgb.train(
  params  = params,
  data    = dtrain,
  nrounds = best_nrounds,
  verbose = 0
)

dtest <- xgb.DMatrix(X_test)
p_test <- predict(model, dtest)        # probability of y=1 => cluster 2

# Hard assignment (0.5 threshold)
pred_cluster_test <- ifelse(p_test >= 0.5, 2, 1)

# Confidence = distance from 0.5, or just max prob
conf <- pmax(p_test, 1 - p_test)

table(pred_cluster_test)
summary(conf)

library(pROC)
pred_cluster_test_01 = ifelse(pred_cluster_test==2,0,1)
# roc_obj <- roc(response = test_surv_df$event, predictor = pred_cluster_test_01)
# auc(roc_obj)
# plot(roc_obj)

imp_fem <- xgb.importance(feature_names = colnames(X_train), model = model)
head(imp_fem, 30)
imp_fem$Feature = var_data_match$names[match(imp_fem$Feature,var_data_match$var)]
imp_fem$Feature_clean <- ifelse(
  grepl("\\(", imp_fem$Feature),
  sub(".*\\((.*)\\)", "\\1", imp_fem$Feature),
  imp_fem$Feature
)
library(dplyr)

imp_fem <- imp_fem %>%
  mutate(
    Feature_clean = ifelse(
      grepl("\\(", Feature),
      sub(".*\\((.*)\\)", "\\1", Feature),
      Feature
    ),
    Type = ifelse(
      grepl("\\(", Feature),
      "Proteomics",
      "Phenomics"
    )
  )

library(ggplot2)

imp_fem %>%
  slice_max(order_by = Gain, n = 20) %>%
  mutate(Feature_clean = factor(Feature_clean, levels = rev(Feature_clean))) %>%
  ggplot(aes(x = Gain, y = Feature_clean, fill = Type)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = c(
    "Proteomics" = "#1f78b4",
    "Phenomics" = "#33a02c"
  )) +
  labs(
    x = "Importance (Gain)",
    y = NULL,
    fill = NULL,
    title = "Females"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "left",        # move legend to the left
    legend.direction = "vertical",   # arrange legend items vertically
    panel.grid.major.y = element_blank()
  )



top20_fem = imp_fem$Feature[1:20]
top20_fem_clean <- ifelse(
  grepl("^[OPQ][0-9][A-Z0-9]{3}[0-9]", top20_fem),
  tolower(sub(" .*", "", top20_fem)),
  top20_fem
)
a = Normalised_list[[1]][,colnames(Normalised_list[[1]]) %in% top20_fem_clean]
b = Normalised_list[[2]][,colnames(Normalised_list[[2]]) %in% top20_fem_clean]


#GHS FEMALES

ghs_females = ghs_bas[ghs_bas$sex == "Women" & ghs_bas$hf_a3 == "HFPEF",]
ghs_females$whf_strict_time =  gsub(",", ".", ghs_females$whf_strict_time, fixed = TRUE)
ghs_females$whf_strict_time = as.numeric(ghs_females$whf_strict_time)
test_surv_df <- data.frame(
 # time  = ghs_females$whf_strict_time,
  time  = ghs_females$whf_strict_time,
  event = ghs_females$whf_strict_event,
  pred_cluster = factor(pred_cluster_test)
)
cox <- coxph(Surv(time, event) ~ pred_cluster, data = test_surv_df)
summary(cox)
# Unadjusted Kaplan-Meier
library(survival)
library(survminer)


sf_km <- survfit(Surv(time, event) ~ pred_cluster, data = test_surv_df)
female_cumulative_ghs = ggsurvplot(
  sf_km,
  data = test_surv_df,
  fun = "event",
  conf.int = TRUE,
  pval = TRUE,
  risk.table = TRUE,
  cumevents = TRUE,
  xlim = c(0, 4),
  
  # Rename legend labels
  legend.labs = c("Predicted high risk cluster", "Predicted low risk cluster"),
  
  legend = "left",
  title = "Females",
  xlab = "Time (years)",
  
  # Increase general font sizes
  fontsize = 6,                 # overall scaling (affects curves & tables)
  risk.table.fontsize = 5,      # risk table text
  cumevents.fontsize = 5,       # cumulative events table
  
  ggtheme = theme_bw() +
    theme(
      legend.box = "vertical",
      legend.direction = "vertical",
      
      # Title
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
      
      # Axis text
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      
      # Legend text
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 12)
    )
)

### Baseline characteristics

ghs_females$cluster = pred_cluster_test
# Clinical characteristics
source(file.path(shared_functions_dir, "baseline13.r"))
setwd(shared_functions_dir)
suppressWarnings(eval(parse(file = "ade.fun.run.r")))

vars_1 <- c('age', 'sex')
nams_1 <- c('Age [y]', 'Sex (Women)')
vars_2 <- c( vars_1, 'bmi',		     'hba1c',	'nic',     'hyper', 'diab',     'dyslip',   'adipos', "ckd", "cvd", "chf", "cancer",  "plaq")
nams_2 <- c( nams_1, 'BMI [kg/m2]', 'HbA1c [%]',	 'Smoking', 'Hypertension', 
             'Diabetes', 'Dyslipidemia',  'Obesity', "Chronic Kidney Disease", "CVD", "Congestive heart failure", "Cancer", "Plaques (yes)")
vars_3 <- c( vars_2, 'afib',  "mi", 'stroke', 'cad', 'pad')
nams_3 <- c(nams_2,    'AF' ,   'MI', 'Stroke', 'CAD', 'PAD')

vars_4 <- c( vars_3, "sbp", "dbp", "heart_rate_pre", "peak_vo2", "ee", "ef", "lvm", "rwt", "bnp" )
nams_4 <- c(nams_3, "SBP", "DBP", "Heart rate", "peak_vO2","E/E'","EF [%]",  "LVM", "RWT", 'NT-proBNP [pg/ml]')

vars_5 <- c(vars_4, "tni",  "natrium", "kalium", "chlorid", "calcium", "v34_lbl05", "v34_lbl06", "v34_lbl07", "v34_lbl08", 
            "got", "gpt", "ggt", "v34_lbl13", "v34_lbl14", "v34_lbl15", "v34_lbl40", "v34_lbl41", "crp_ln", "v34_lbl50", "glucose",
            "v34_lbl54", "v34_lbl55", "v34_lbl56", "v34_lbl57",   "v34_lbl67", "egfr")
nams_5 <- c( nams_4, "Troponin I",   "Na", "K", "Cl", "Ca", "Mg", "P", "Creatinine", "Urea", "GOT", "GPT", "GGT", "Alkalinephosphatase", "Bilirubin", "Albumin",
             "LDH", "Lipase","CRP", "HbA1c", "Glucose", "Cholesterine", "Triglycerides", "HDL", "LDL","TSH", "GFR")

vars_6 <- c(vars_5, "v34_lbl18",  "v34_lbl19", "v34_lbl20", "v34_lbl21", "v34_lbl22", "v34_lbl23", "v34_lbl24", "v34_lbl25",
            "v34_lbl26", "v34_lbl81", "v34_lbl27", "v34_lbl28","v34_lbl29","v34_lbl30","v34_lbl31")
nams_6 <- c( nams_5, "Erythrocytes",   "Leukocytes", "Hemoglobin", "Hematocrit", "MCV", "MCH", "MCHC", "EVB", "Thrombocytes", "MTV",
             "Neutrophils", "Lymphocytes", "Monocytes", "Eosinophils", "Basophils")

setwd(project_dir)
BT1_females_ghs <- baseline13(vars_3, nams_3, data = ghs_females, skeew=1, group = 'cluster', pv=T, style=2, minuniq = 3)


# RNA
counts_vsd_tr <- read.csv(file.path(data_dir, "counts_vsd_tr.csv"))
rna_lids <- readRDS(file.path(data_dir, "rna_lids.rds"))
dim(counts_vsd_tr)
rna_myovasc = merge(rna_lids, counts_vsd_tr, by.x = "Data", by.y = "X")
rna_myovasc = merge(rna_myovasc, id_lid,  by ="LID")
rownames(rna_myovasc) = rna_myovasc$ID
rna_myovasc = rna_myovasc[,-c(1:3,13903)]


all_rna_females = rna_myovasc[rownames(rna_myovasc) %in% hfpef_id_female,]
all_rna_females$cluster = memberships_female$cl2[na.omit(match(rownames(memberships_female),rownames(all_rna_females)))]
df <- all_rna_females
df$cluster <- factor(df$cluster)

df <- all_rna_females
df$cluster <- factor(df$cluster)


# Assuming df has cluster in column "cluster"
# and all other columns are genes
genes <- setdiff(colnames(df), "cluster")
df$cluster <- factor(df$cluster, levels = c(1,2), labels = c("High", "Low"))
table(df$cluster)
results <- data.frame(
  Gene = genes,
  p.value = NA,
  log2FC = NA
)

for (gene in genes) {
  high <- df[[gene]][df$cluster == "High"]
  low  <- df[[gene]][df$cluster == "Low"]

  t_res <- t.test(high, low)

  results$p.value[results$Gene == gene] <- t_res$p.value
  results$log2FC[results$Gene == gene] <- mean(high) - mean(low)
}

results$adj.p <- p.adjust(results$p.value, method = "fdr")
head(results)


library(clusterProfiler)
library(org.Hs.eg.db)  # human gene annotations
gene_df <- results %>%
  dplyr::select(Gene, log2FC)   # assuming you have a log2FC column

# Map symbols to Entrez IDs
gene_df$EntrezID <- mapIds(org.Hs.eg.db,
                           keys = gene_df$Gene,
                           column = "ENTREZID",
                           keytype = "SYMBOL",
                           multiVals = "first")

# Remove genes that could not be mapped
gene_df <- gene_df[!is.na(gene_df$EntrezID), ]

gene_list <- gene_df$log2FC
names(gene_list) <- gene_df$EntrezID
gene_list <- sort(gene_list, decreasing = TRUE)  # required by GSEA

gsea_res <- gseKEGG(
  geneList = gene_list,
  organism = "hsa",
  minGSSize = 10,
  pvalueCutoff = 0.5      # keep soft threshold
  # nPerm removed!
)

library(ReactomePA)
gene_list <- gene_df$log2FC
names(gene_list) <- gene_df$EntrezID

# clean
gene_list <- gene_list[!is.na(names(gene_list))]
gene_list <- gene_list[!duplicated(names(gene_list))]

# IMPORTANT: names must be character
names(gene_list) <- as.character(names(gene_list))

# sort (required for GSEA)
gene_list <- sort(gene_list, decreasing = TRUE)

gsea_res <- gsePathway(
  geneList = gene_list,
  organism = "human",   # or "mouse"
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  minGSSize = 10,
  maxGSSize = 500,
  verbose = FALSE
)
library(enrichplot)
library(ggplot2)

dotplot(gsea_res, showCategory = 20) +
  scale_x_continuous(limits = c(0, NA)) +
  ggtitle("Significantly Enriched RNA Pathways in Females")


################################ Males
intersection <- intersect(rownames(panels_list[[1]]), rownames(panels_list[[2]]))
d = dall1[dall1$v11_sid01 %in% intersection,]
d = d[d$hf_a3 == "HFPEF",]
hfpef_id_male = d$v11_sid01[which(d$hf_a3 =="HFPEF" & d$p002 =="1")]
d = d[d$v11_sid01 %in% hfpef_id_male,]
intersection = intersection[intersection %in% hfpef_id_male]
panels_list_intersected = list()
for (i in 1:length(panels_list)){
  panels_list_intersected[[i]] =  panels_list[[i]][(rownames(panels_list[[i]]) %in% intersection),]
}
panels_list_intersected <- lapply(panels_list_intersected, function(x) x[order(rownames(x)), , drop = FALSE])

saveRDS(panels_list_intersected, file.path(output_dir, "panels_list_intersected_male.RDS"))
panels_list_intersected[[2]] <- apply(panels_list_intersected[[2]], 2, as.numeric)
rownames(panels_list_intersected[[2]]) = rownames(panels_list_intersected[[1]])

W = list()
Normalised_list = list()
train_means = list()
train_sds = list()
for (i in 1:length(panels_list_intersected)){  
  
  #if (data_integration == 'SNF'){
  # as suggested by authors apply data normalization
  
  data_imputed = imputa(panels_list_intersected[[i]])
  nzv = nearZeroVar(data_imputed, names = TRUE)
  print(length(nzv))
  
  if  (length(nzv) > 0) {
    data_imputed= data_imputed[,-which(colnames(data_imputed) %in% nzv)]
  }
  
  train_means[[i]] <- apply(data_imputed, 2, mean, na.rm = TRUE)
  train_sds[[i]]   <- apply(data_imputed, 2, sd, na.rm = TRUE)
  
  data_norm = scale(data_imputed, center = train_means[[i]], scale = train_sds[[i]])
  
  cc = cor(data_norm, use = "pairwise.complete.obs", method = "pearson")
  select_corr = caret::findCorrelation(cc, cutoff = 0.9, exact = FALSE)
  print(length(select_corr))
  
  if (length(select_corr)>0){
    data_norm_corr= data_norm[, -select_corr]
  }
  
  ## Calculate the pair-wise distance;
  dist_mat = (dist2(as.matrix(data_norm),as.matrix(data_norm)))^(1/2)
  W_mat = affinityMatrix(dist_mat, K = 20, sigma = 0.5) 
  W[[i]] = W_mat
  Normalised_list[[i]] = data_norm_corr
  print(dim(W[[i]]))
}

setwd(project_dir)

### SNF & Spectral Clustering
snf_males = SNF(W, 20, 20)
memberships_male = data.frame()
for (i in 2:3){
  memberships_male = rbind(memberships_male, SNFtool::spectralClustering(snf_males, K = i))
}
memberships_male = t(memberships_male)
colnames(memberships_male) = c("cl2", "cl3")
rownames(memberships_male) = rownames(snf_males)
memberships_male = as.data.frame(memberships_male)


#### Clinical check 
d$memberships = as.factor(memberships_male$cl2)
cox_age <- coxph(Surv(whf_all_time,whf_all_event) ~  memberships + age , data =  d)
summary(cox_age)

# Unadjusted Kaplan-Meier
library(survival)
library(survminer)

sf_km <- survfit(Surv(whf_all_time,whf_all_event) ~ memberships, data = d)
male_cumulative = ggsurvplot(
  sf_km,
  data = d,
  fun = "event",
  conf.int = TRUE,
  pval = TRUE,
  risk.table = TRUE,
  cumevents = TRUE,
  xlim = c(0, 4),
  
  # Rename legend labels
  legend.labs = c("High risk cluster", "Low risk cluster"),
  
  legend = "left",
  title = "Males",
  xlab = "Time (years)",
  
  # Increase general font sizes
  fontsize = 6,                 # overall scaling (affects curves & tables)
  risk.table.fontsize = 5,      # risk table text
  cumevents.fontsize = 5,       # cumulative events table
  
  ggtheme = theme_bw() +
    theme(
      legend.box = "vertical",
      legend.direction = "vertical",
      
      # Title
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
      
      # Axis text
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      
      # Legend text
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 12)
    )
)


# Clinical characteristics
source(file.path(shared_functions_dir, "baseline13.r"))
setwd(shared_functions_dir)
suppressWarnings(eval(parse(file = "ade.fun.run.r")))

vars_1 <- c('age', 'sex', 'maggic')
nams_1 <- c('Age [y]', 'Sex (Women)')
vars_2 <- c( vars_1,   'nic',     'hyper',  'diab',     'dyslip',   'adipos', "ckd", "cvd", "chf", "cancer")
nams_2 <- c( nams_1,   'Smoking', 'Hypertension',
             'Diabetes', 'Dyslipidemia',  'Obesity', "Chronic Kidney Disease", "CVD", "Congestive heart failure", "Cancer")
vars_3 <- c( vars_2, "nyha_ord2" ,'afib',  "mi", 'stroke', 'cad', 'pad')
nams_3 <- c(nams_2,   "NYHA", 'AF' ,   'MI', 'Stroke', 'CAD', 'PAD')

vars_4 <- c( vars_3, "sbp", "dbp", "heart_rate_pre", "peak_vo2", "ee", "ef", "lvm", "rwt", "bnp" )
nams_4 <- c(nams_3, "SBP", "DBP", "Heart rate", "peak_vO2","E/E'","EF [%]",  "LVM", "RWT", 'NT-proBNP [pg/ml]')

vars_5 <- c(vars_4, "tni",  "natrium", "kalium", "chlorid", "calcium", "v34_lbl05", "v34_lbl06", "v34_lbl07", "v34_lbl08", 
            "got", "gpt", "ggt", "v34_lbl13", "v34_lbl14", "v34_lbl15", "v34_lbl40", "v34_lbl41", "crp_ln", "v34_lbl50", "glucose",
            "v34_lbl54", "v34_lbl55", "v34_lbl56", "v34_lbl57",   "v34_lbl67", "egfr")
nams_5 <- c( nams_4, "Troponin I",   "Na", "K", "Cl", "Ca", "Mg", "P", "Creatinine", "Urea", "GOT", "GPT", "GGT", "Alkalinephosphatase", "Bilirubin", "Albumin",
             "LDH", "Lipase","CRP", "HbA1c", "Glucose", "Cholesterine", "Triglycerides", "HDL", "LDL","TSH", "GFR")

vars_6 <- c(vars_5, "v34_lbl18",  "v34_lbl19", "v34_lbl20", "v34_lbl21", "v34_lbl22", "v34_lbl23", "v34_lbl24", "v34_lbl25",
            "v34_lbl26", "v34_lbl81", "v34_lbl27", "v34_lbl28","v34_lbl29","v34_lbl30","v34_lbl31")
nams_6 <- c( nams_5, "Erythrocytes",   "Leukocytes", "Hemoglobin", "Hematocrit", "MCV", "MCH", "MCHC", "EVB", "Thrombocytes", "MTV",
             "Neutrophils", "Lymphocytes", "Monocytes", "Eosinophils", "Basophils")

d$v34_lbl34 = as.numeric(d$v34_lbl34)
d$v34_lbl50 = as.numeric(d$v34_lbl50)
d$v34_lbl56 = as.numeric(d$v34_lbl56)
d$v34_lbl57 = as.numeric(d$v34_lbl57)
d$v34_lbl14 = as.numeric(d$v34_lbl14)
d$v34_lbl41 = as.numeric(d$v34_lbl41)
d$v34_lbl67 = as.numeric(d$v34_lbl67)
d$nyha_ord2 = as.factor(d$nyha_ord2)
d$egfr = as.numeric(d$egfr)
d$v34_lbl50 = as.numeric(d$v34_lbl50)
d$v34_lbl50 = as.numeric(d$v34_lbl50)
d$v34_lbl56 = as.numeric(d$v34_lbl56)
d$v34_lbl57 = as.numeric(d$v34_lbl57)
d$nyha_ord2 = as.factor(d$nyha_ord2)
d$afib_type <- d$myo029
d$afib_type[d$myo029 == 88 | d$myo029 == 99] <- NA
d$afib_type <- factor(d$afib_type, labels = c('paroxysmal', 'permanent', 'persistent'))
d$egfr = as.numeric(d$egfr)
d$fli_cat = ifelse(d$fli>=60, "yes", "no")
setwd(project_dir)
#d$memberships = memberships_male$cl2
#BT1_males <- baseline13(vars_6, nams_6, data = d, skeew=1, group = 'memberships', pv=T, style=2, minuniq = 3)
ghs_males$memberships = pred_cluster_test
BT1_males <- baseline13(vars_6, nams_6, data = ghs_males, skeew=1, group = 'memberships', pv=T, style=2, minuniq = 3)

#XGBoost
# GHS male
colnames(ghs_bas)[colnames(ghs_bas) == "whtr"] = "wthr"
colnames(ghs_bas)[colnames(ghs_bas) == "v34_lbl01"] = "natrium"
colnames(ghs_bas)[colnames(ghs_bas) == "v34_lbl03"] = "chlorid"
colnames(ghs_bas)[colnames(ghs_bas) == "v34_lbl04"] = "calcium"

proteins_ghs = ghs_bas[,colnames(ghs_bas) %in% dattr_ghs$vars_prot]
rownames(proteins_ghs) = ghs_bas$v11_sid01
phenomics_ghs = ghs_bas[, colnames(ghs_bas) %in% c(cardiac, blood, anthro)]

X_test_df = cbind(proteins_ghs, phenomics_ghs)
X_test_df = X_test_df[ghs_bas$sex == "Men" & ghs_bas$hf_a3 == "HFPEF",]

X_train_df  = cbind(Normalised_list[[1]], Normalised_list[[2]])
X_test_df = X_test_df[,colnames(X_test_df) %in% colnames(X_train_df)]
cluster_train = memberships_male$cl2
library(xgboost)
library(caret)
library(pROC)
library(PRROC)

# 1) Align columns train/test (critical)
common_features <- intersect(colnames(X_train_df), colnames(X_test_df))
X_train_df <- X_train_df[, common_features, drop = FALSE]
X_test_df  <- X_test_df[,  common_features, drop = FALSE]

X_train <- data.matrix(X_train_df)
X_test  <- data.matrix(X_test_df)

# 2) Binary labels: cluster 1/2 -> y 0/1 (1 = cluster 2)
y_train <- ifelse(cluster_train == 2, 1, 0)
stopifnot(all(y_train %in% c(0,1)))

dtrain <- xgb.DMatrix(X_train, label = y_train)

# 3) Parameters (binary)
params <- list(
  booster = "gbtree",
  objective = "binary:logistic",
  eval_metric = c("auc", "aucpr"),
  eta = 0.05,
  max_depth = 4,
  min_child_weight = 1,
  subsample = 0.8,
  colsample_bytree = 0.8,
  lambda = 1,
  alpha = 0
)

# 4) CV to pick best nrounds
set.seed(256)
cv <- xgb.cv(
  params = params,
  data = dtrain,
  nrounds = 5000,
  nfold = 5,
  stratified = TRUE,
  early_stopping_rounds = 50,
  maximize = TRUE,
  verbose = 0
)

best_nrounds <- cv$best_iteration
if (is.null(best_nrounds) || !is.finite(best_nrounds) || best_nrounds < 1) {
  elog <- cv$evaluation_log
  best_nrounds <- if ("test_aucpr_mean" %in% names(elog)) which.max(elog$test_aucpr_mean) else which.max(elog$test_auc_mean)
}

# 5) Fit final model on full training set
model <- xgb.train(
  params  = params,
  data    = dtrain,
  nrounds = best_nrounds,
  verbose = 0
)

dtest <- xgb.DMatrix(X_test)
p_test <- predict(model, dtest)        # probability of y=1 => cluster 2

# Hard assignment (0.5 threshold)
pred_cluster_test <- ifelse(p_test >= 0.5, 2, 1)

# Confidence = distance from 0.5, or just max prob
conf <- pmax(p_test, 1 - p_test)

table(pred_cluster_test)
summary(conf)

imp_male <- xgb.importance(feature_names = colnames(X_train), model = model)
head(imp_male, 30)
imp_male$Feature = var_data_match$names[match(imp_male$Feature,var_data_match$var)]
imp_male$Feature_clean <- ifelse(
  grepl("\\(", imp_male$Feature),
  sub(".*\\((.*)\\)", "\\1", imp_male$Feature),
  imp_male$Feature
)
library(dplyr)

imp_male <- imp_male %>%
  mutate(
    Feature_clean = ifelse(
      grepl("\\(", Feature),
      sub(".*\\((.*)\\)", "\\1", Feature),
      Feature
    ),
    Type = ifelse(
      grepl("\\(", Feature),
      "Proteomics",
      "Phenomics"
    )
  )

library(ggplot2)

imp_male %>%
  slice_max(order_by = Gain, n = 20) %>%
  mutate(Feature_clean = factor(Feature_clean, levels = rev(Feature_clean))) %>%
  ggplot(aes(x = Gain, y = Feature_clean, fill = Type)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = c(
    "Proteomics" = "#1f78b4",
    "Phenomics" = "#33a02c"
  )) +
  labs(
    x = "Importance (Gain)",
    y = NULL,
    fill = NULL,
    title = "Males"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "left",        # move legend to the left
    legend.direction = "vertical",   # arrange legend items vertically
    panel.grid.major.y = element_blank()
  )




library(survival)
library(survminer)

ghs_males = ghs_bas[ghs_bas$sex == "Men" & ghs_bas$hf_a3 == "HFPEF",]
ghs_males$whf_strict_time =  gsub(",", ".", ghs_males$whf_strict_time, fixed = TRUE)
ghs_males$whf_strict_time = as.numeric(ghs_males$whf_strict_time)
test_surv_df <- data.frame(
  time  = ghs_males$whf_strict_time,
  event = ghs_males$whf_strict_event,
  pred_cluster = factor(pred_cluster_test)
)
cox <- coxph(Surv(time, event) ~ pred_cluster, data = test_surv_df)
summary(cox)


# Unadjusted Kaplan-Meier
library(survival)
library(survminer)

sf_km <- survfit(Surv(time, event) ~ pred_cluster, data = test_surv_df)
ggsurvplot(
  sf_km,
  data = test_surv_df,
  fun = "event",        # cumulative incidence = 1 - S(t)
  conf.int = TRUE,
  pval = TRUE,
  
  risk.table = TRUE,    # number at risk
  cumevents = TRUE,     # cumulative number of events ("cases")
  
  risk.table.title = "Number at risk",
  cumevents.title   = "Cumulative cases",
  
  legend.title = "Memberships",
  xlab = "Follow-up time",
  ylab = "Cumulative incidence",
  ggtheme = theme_bw()
)



# Important variables

renamed_males = merge(var_data_match, imp, by.x = "var" , by.y = "Feature")
renamed_females = merge(var_data_match, imp_fem, by.x = "var" , by.y = "Feature")

library(dplyr)
library(ggplot2)

top_n <- 20

imp_top <- renamed_males %>%
  arrange(desc(Gain)) %>%
  slice_head(n = top_n) %>%
  mutate(feature = factor(.data$names, levels = rev(.data$names)))

ggplot(imp_top, aes(x = feature, y = Gain)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(
    x = NULL,
    y = "Gain",
    title = "Defining Features in Males"
  ) +
  theme_minimal(base_size = 13)
library(dplyr)
library(ggplot2)

top_n <- 20

imp_top <- renamed_females %>%
  arrange(desc(Gain)) %>%
  slice_head(n = top_n) %>%
  mutate(feature = factor(.data$names, levels = rev(.data$names)))

ggplot(imp_top, aes(x = feature, y = Gain)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(
    x = NULL,
    y = "Gain",
    title = "Defining Features in Females"
  ) +
  theme_minimal(base_size = 13)


boxplot(
  CD40 ~ cluster,
  data = df,
  main = "CD40 expression by cluster",
  xlab = "Cluster",
  ylab = "Expression",
  col = "lightgray"
)


# Figure for publication
male_cumulative <- ggsurvplot(
  sf_km,
  data = d,
  fun = "event",
  conf.int = TRUE,
  pval = TRUE,
  risk.table = TRUE,
  cumevents = TRUE,
  legend = "right",   # place legend on right
  title = "Males",
  ggtheme = theme_bw() +
    theme(
      legend.box = "vertical",
      legend.direction = "vertical",
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold")  # center title
    )
)

test_cumulative <- ggsurvplot(
  sf_km,
  data = test_surv_df,
  fun = "event",
  conf.int = TRUE,
  pval = TRUE,
  risk.table = TRUE,
  cumevents = TRUE,
  legend = "none",   # remove duplicate legend
  ggtheme = theme_bw()
)

arrange_ggsurvplots(
  list(male_cumulative, test_cumulative),
  ncol = 2,
  nrow = 1
)


# RNA
counts_vsd_tr <- read.csv(file.path(data_dir, "counts_vsd_tr.csv"))
rna_lids <- readRDS(file.path(data_dir, "rna_lids.rds"))
dim(counts_vsd_tr)
rna_myovasc = merge(rna_lids, counts_vsd_tr, by.x = "Data", by.y = "X")
rna_myovasc = merge(rna_myovasc, id_lid,  by ="LID")
rownames(rna_myovasc) = rna_myovasc$ID
rna_myovasc = rna_myovasc[,-c(1:3,13903)]


all_rna_males = rna_myovasc[rownames(rna_myovasc) %in% hfpef_id_male,]
all_rna_males$cluster = memberships_male$cl2[na.omit(match(rownames(memberships_male),rownames(all_rna_males)))]
df <- all_rna_males
df$cluster <- factor(df$cluster)

df <- all_rna_males
df$cluster <- factor(df$cluster)


# Assuming df has cluster in column "cluster"
# and all other columns are genes
genes <- setdiff(colnames(df), "cluster")
df$cluster <- factor(df$cluster, levels = c(1,2), labels = c("High", "Low"))
table(df$cluster)
results <- data.frame(
  Gene = genes,
  p.value = NA,
  log2FC = NA
)

for (gene in genes) {
  high <- df[[gene]][df$cluster == "High"]
  low  <- df[[gene]][df$cluster == "Low"]
  
  t_res <- t.test(high, low)
  
  results$p.value[results$Gene == gene] <- t_res$p.value
  results$log2FC[results$Gene == gene] <- mean(high) - mean(low)
}

results$adj.p <- p.adjust(results$p.value, method = "fdr")
head(results)


library(clusterProfiler)
library(org.Hs.eg.db)  # human gene annotations
gene_df <- results %>%
  dplyr::select(Gene, log2FC)   # assuming you have a log2FC column

# Map symbols to Entrez IDs
gene_df$EntrezID <- mapIds(org.Hs.eg.db,
                           keys = gene_df$Gene,
                           column = "ENTREZID",
                           keytype = "SYMBOL",
                           multiVals = "first")

# Remove genes that could not be mapped
gene_df <- gene_df[!is.na(gene_df$EntrezID), ]

gene_list <- gene_df$log2FC
names(gene_list) <- gene_df$EntrezID
gene_list <- sort(gene_list, decreasing = TRUE)  # required by GSEA

gsea_res <- gseKEGG(
  geneList = gene_list,
  organism = "hsa",
  minGSSize = 10,
  pvalueCutoff = 0.5      # keep soft threshold
  # nPerm removed!
)

library(ReactomePA)
gene_list <- gene_df$log2FC
names(gene_list) <- gene_df$EntrezID

# clean
gene_list <- gene_list[!is.na(names(gene_list))]
gene_list <- gene_list[!duplicated(names(gene_list))]

# IMPORTANT: names must be character
names(gene_list) <- as.character(names(gene_list))

# sort (required for GSEA)
gene_list <- sort(gene_list, decreasing = TRUE)

gsea_res <- gsePathway(
  geneList = gene_list,
  organism = "human",   
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  minGSSize = 10,
  maxGSSize = 500,
  verbose = FALSE
)
library(enrichplot)
library(ggplot2)

dotplot(gsea_res, showCategory = 20) +
  scale_x_continuous(limits = c(0, NA)) +
  ggtitle("Significantly Enriched RNA Pathways in Males")




save(memberships_male, Normalised_list, imp_male, file = file.path(output_dir, "males.RData"))
save(memberships_female, Normalised_list, imp_fem, file = file.path(output_dir, "females.RData"))
