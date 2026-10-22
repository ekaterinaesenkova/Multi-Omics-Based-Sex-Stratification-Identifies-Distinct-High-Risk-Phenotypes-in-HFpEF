# Main project for Cardiovascular Diabetology
# ===== BLOCK 0: setup =====
library(SNFtool); library(survival); library(caret); library(missRanger)
library(stringr); library(janitor); library(readxl); library(Hmisc); library(dplyr)
project_dir <- "Q:\\DIASyM_MSCoreSys\\Workfolder\\E.Esenkova\\Esenkova\\R\\MYOVASC_PROJECT_3"

setwd("T:\\Myovasc\\R_SQL"); source("Read_SQL_MyoVasc_BL_v11.R")
dall1    <- read_SQL_MyoVasc_BL(proteins = TRUE, lipids = TRUE, user = "ese0e", ATC_7_digits = "B01AC06")
olink_mv <- read.csv2("N:/Transfer/gie2a/MyoVasc/Olink transformed/olink_6_BL_FU2_FU4_co-transformed_imp_BL_2024-02-29.csv")
Sys.setlocale("LC_ALL", "GER")
suppressWarnings(eval(parse(file = "T:/GHS/R_SQL/Read_SQL_Data_1_S_jch.r")))
ghs_bas <- read_SQL_A6_data(reload = TRUE, read_proteins = TRUE, u_id = "ese0e")
ghs_bas <- merge(ghs_bas, read.csv("T:/GHS/R_SQL/whf_ghs.csv", sep = ";"), by = "v11_sid01")
setwd(project_dir)

num <- function(x) suppressWarnings(as.numeric(gsub(",", ".", as.character(x))))
yes <- function(x) x %in% c(1, "1", "yes", "Yes", TRUE)

prot_vars <- intersect(intersect(attributes(dall1)$vars_prot, colnames(olink_mv)), colnames(ghs_bas))
for (o in c(whtr = "wthr", v34_lbl01 = "natrium", v34_lbl03 = "chlorid", v34_lbl04 = "calcium") |> names())
  names(ghs_bas)[names(ghs_bas) == o] <- c(whtr = "wthr", v34_lbl01 = "natrium", v34_lbl03 = "chlorid", v34_lbl04 = "calcium")[o]

add_vars <- function(d, female) {
  d$female <- as.integer(female)
  h <- num(d$height) * 100
  d$lvmi <- num(d$lvm) / sqrt(h * num(d$weight) / 3600)
  d$tyg  <- log(num(d$v34_lbl55) * num(d$glucose) / 2)
  d
}
mv  <- merge(dall1[, setdiff(colnames(dall1), prot_vars)], olink_mv[, c("v11_sid01", prot_vars)], by = "v11_sid01")
mv  <- add_vars(mv, mv$p002 == "2")
ghs_ol <- rowMeans(!is.na(ghs_bas[, prot_vars])) > 0.5
ghs <- add_vars(ghs_bas[ghs_ol, ], ghs_bas$sex[ghs_ol] == "Women")

# preprocessing + SNF (ID-safe, correlation filter applied BEFORE SNF)
impute_rf <- function(X) if (!anyNA(X)) X else missRanger(X, num.trees = 100, pmm.k = 0, verbose = 0, seed = 256)
prep_view <- function(d, vars, view_name = "") {
  X <- as.data.frame(lapply(d[, vars, drop = FALSE], num)); rownames(X) <- as.character(d$v11_sid01)
  X <- X[, colMeans(is.na(X)) <= 0.3, drop = FALSE]
  keep <- rowMeans(is.na(X)) <= 0.5
  if (any(!keep)) cat("  ", view_name, ": removed", sum(!keep), "participants (>50% missing)\n")
  X <- X[keep, , drop = FALSE]; Xi <- impute_rf(X); rownames(Xi) <- rownames(X)
  nzv <- nearZeroVar(Xi, names = TRUE); if (length(nzv)) Xi <- Xi[, setdiff(colnames(Xi), nzv), drop = FALSE]
  Z <- scale(as.matrix(Xi))
  Z[, setdiff(colnames(Z), findCorrelation(cor(Z), cutoff = 0.9, names = TRUE)), drop = FALSE]
}
run_snf <- function(d) {
  Z <- lapply(names(views), function(v) prep_view(d, intersect(views[[v]], colnames(d)), v)); names(Z) <- names(views)
  ids <- Reduce(intersect, lapply(Z, rownames)); Z <- lapply(Z, function(z) z[ids, , drop = FALSE])
  W <- lapply(Z, function(z) affinityMatrix(dist2(z, z)^(1/2), K = 20, sigma = 0.5))
  fused <- SNF(W, 20, 20); dimnames(fused) <- list(ids, ids)
  list(Z = Z, W = W, fused = fused)
}

# original 2 views (as in the manuscript)
clin_vars <- c("sbp", "dbp", "ee", "ef", "lvm", "rwt",
               "tni", "natrium", "kalium", "chlorid", "calcium", "v34_lbl07", "v34_lbl08", "got", "ggt",
               "v34_lbl13", "v34_lbl14", "v34_lbl15", "v34_lbl40", "v34_lbl41", "crp_ln", "v34_lbl50", "glucose",
               "v34_lbl54", "v34_lbl55", "v34_lbl56", "v34_lbl57", "v34_lbl67", "egfr",
               "v34_lbl18", "v34_lbl19", "v34_lbl20", "v34_lbl21", "v34_lbl22", "v34_lbl23", "v34_lbl24", "v34_lbl25",
               "v34_lbl26", "v34_lbl81", "v34_lbl27", "v34_lbl28", "v34_lbl29", "v34_lbl30", "v34_lbl31",
               "age", "weight", "height", "waist", "hip", "wthr")
views <- list(proteins = prot_vars, clinical = clin_vars)

cat("Proteins:", length(prot_vars), "| clinical vars found in MyoVasc:", sum(clin_vars %in% names(mv)),
    "of", length(clin_vars), "| in GHS:", sum(clin_vars %in% names(ghs)), "\n")
cat("HFpEF with Olink - MyoVasc women/men:", sum(mv$hf_a3 == "HFPEF" & mv$female == 1, na.rm = TRUE), "/",
    sum(mv$hf_a3 == "HFPEF" & mv$female == 0, na.rm = TRUE),
    "| GHS women/men:", sum(ghs$hf_a3 == "HFPEF" & ghs$female == 1, na.rm = TRUE), "/",
    sum(ghs$hf_a3 == "HFPEF" & ghs$female == 0, na.rm = TRUE), "\n")
cat("Candidate outcome/death variables:\n")
print(grep("whf|death|tod|cardiac", names(mv), value = TRUE))

# ===== BLOCK 1: discovery (HFpEF, bug-fixed) =====
hf_mv <- mv[which(mv$hf_a3 == "HFPEF"), ]
res <- list()
for (s in c("women", "men")) {
  set.seed(256)
  d0 <- hf_mv[hf_mv$female == as.integer(s == "women"), ]
  r  <- run_snf(d0)
  lab <- setNames(spectralClustering(r$fused, K = 2), rownames(r$fused))
  d  <- d0[match(names(lab), as.character(d0$v11_sid01)), ]
  stopifnot(identical(as.character(d$v11_sid01), names(lab)))
  d$time <- num(d$whf_all_time); d$event <- num(d$whf_all_event)
  rate <- tapply(d$event == 1, lab, mean, na.rm = TRUE)
  d$cluster <- factor(ifelse(lab == names(which.max(rate)), "High", "Low"), levels = c("Low", "High"))
  r$data <- d; r$k_est <- estimateNumberOfClustersGivenGraph(r$fused, NUMC = 2:6); res[[s]] <- r
  
  m1 <- coxph(Surv(time, event) ~ cluster + age, data = d)
  m2 <- coxph(Surv(time, event) ~ cluster + age + nyha, data = d)
  ci <- function(m) { x <- summary(m)$conf.int["clusterHigh", c(1, 3, 4)]; p <- summary(m)$coef["clusterHigh", 5]
  sprintf("HR %.2f (%.2f-%.2f), p = %.2g", x[1], x[2], x[3], p) }
  cat("\n=====", s, "=====  n =", nrow(d), "\n")
  print(table(cluster = d$cluster, WHF_event = d$event))
  cat("K suggested: eigengap", r$k_est[[1]], "| rotation cost", r$k_est[[3]], "\n")
  cat("Age-adjusted:        ", ci(m1), "\n")
  cat("Age + NYHA-adjusted: ", ci(m2), "\n")
  cat("Log-rank p:          ", signif(1 - pchisq(survdiff(Surv(time, event) ~ cluster, data = d)$chisq, 1), 2), "\n")
  cat("Diabetes % (High/Low):", round(100 * tapply(yes(d$diab), d$cluster, mean)[c("High", "Low")], 1), "\n")
}
saveRDS(res, file.path(project_dir, "WP1_discovery.rds"))

# Did the old row-position bug hit? (old script assigned clusters by row position)
old_f <- tryCatch(readRDS(file.path(project_dir, "FEMALES.RDS"))$memberships_female, error = function(e) NULL)
if (!is.null(old_f)) {
  old_order <- as.character(dall1$v11_sid01[dall1$v11_sid01 %in% rownames(old_f)])
  cat("\nOld script - row order of d identical to cluster order?", identical(old_order, rownames(old_f)),
      " (FALSE = labels were attached to the wrong patients)\n")
}
# ===== BLOCK 2: HFpEF documentation, stability, hyperparameters =====
med_iqr <- function(x) { x <- num(x); if (all(is.na(x))) return(NA)
sprintf("%.1f (%.1f-%.1f)", median(x, na.rm = TRUE), quantile(x, .25, na.rm = TRUE), quantile(x, .75, na.rm = TRUE)) }
cat("Possible LA volume variables:", grep("lavi|la_vol|lav|lae", names(mv), value = TRUE), "\n")

echo_vars <- c(EF = "ef", Ee = "ee", LVMI = "lvmi", RWT = "rwt", NTproBNP = "bnp", TnI = "tni")
doc <- list()
for (s in names(res)) {
  d <- res[[s]]$data
  rows <- lapply(names(echo_vars), function(v) {
    x <- num(d[[echo_vars[v]]])
    c(variable = v, High = med_iqr(x[d$cluster == "High"]), Low = med_iqr(x[d$cluster == "Low"]),
      p = signif(wilcox.test(x ~ d$cluster)$p.value, 2))
  })
  extra <- rbind(
    c("EF >= 50 %", sprintf("%.1f%%", 100 * tapply(num(d$ef) >= 50, d$cluster, mean, na.rm = TRUE)[c("High", "Low")]), NA),
    c("chf (history) %", sprintf("%.1f%%", 100 * tapply(yes(d$chf), d$cluster, mean)[c("High", "Low")]), NA))
  doc[[s]] <- rbind(do.call(rbind, rows), extra)
  cat("\n=====", s, ": cardiac characterisation by cluster =====\n"); print(doc[[s]], quote = FALSE)
  cat("NYHA by cluster:\n"); print(table(NYHA = d$nyha, d$cluster))
}
write.csv(do.call(rbind, lapply(names(doc), function(s) data.frame(sex = s, doc[[s]]))),
          file.path(project_dir, "WP2_cardiac_by_cluster_MyoVasc.csv"), row.names = FALSE)

# Stability (200 subsamples of 80%) and SNF hyperparameter grid
agree2 <- function(a, b) max(mean(a == b), mean(a != b))
stab <- list(); grid <- list()
for (s in names(res)) {
  Z <- res[[s]]$Z
  ref <- setNames(as.integer(res[[s]]$data$cluster), as.character(res[[s]]$data$v11_sid01))
  ids <- rownames(Z[[1]]); set.seed(256)
  stab[[s]] <- t(replicate(200, {
    sub <- sample(ids, floor(0.8 * length(ids)))
    W <- lapply(Z, function(z) affinityMatrix(dist2(z[sub, ], z[sub, ])^(1/2), K = 20, sigma = 0.5))
    lb <- spectralClustering(SNF(W, 20, 20), K = 2)
    c(agreement = agree2(lb, ref[sub]), NMI = calNMI(lb, ref[sub]))
  }))
  grid[[s]] <- expand.grid(K = c(10, 15, 20, 25, 30), sigma = c(0.3, 0.5, 0.8))
  grid[[s]]$agreement <- apply(grid[[s]], 1, function(g) {
    W <- lapply(Z, function(z) affinityMatrix(dist2(z, z)^(1/2), K = g["K"], sigma = g["sigma"]))
    agree2(spectralClustering(SNF(W, g["K"], 20), K = 2), ref[ids])
  })
  cat("\n=====", s, ": stability (200 x 80%) =====\n")
  print(round(apply(stab[[s]], 2, quantile, c(0.025, 0.5, 0.975)), 3))
  cat("Hyperparameter grid - agreement with main solution: median",
      round(median(grid[[s]]$agreement), 3), "| min", round(min(grid[[s]]$agreement), 3), "\n")
}
saveRDS(list(doc = doc, stab = stab, grid = grid), file.path(project_dir, "WP2_doc_stability.rds"))

# ===== BLOCK 3: GHS validation =====
library(xgboost)
auc_fn <- function(p, y) { r <- rank(p); n1 <- sum(y == 1); n0 <- sum(y == 0); (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0) }
metrics <- function(p, y, thr = 0.5) {
  pr <- as.integer(p >= thr); tp <- sum(pr == 1 & y == 1); tn <- sum(pr == 0 & y == 0)
  fp <- sum(pr == 1 & y == 0); fn <- sum(pr == 0 & y == 1)
  mcc <- (tp * tn - fp * fn) / sqrt(as.numeric(tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
  c(AUC = auc_fn(p, y), sens = tp / (tp + fn), spec = tn / (tn + fp), PPV = tp / (tp + fp), NPV = tn / (tn + fn),
    acc = (tp + tn) / length(y), bal_acc = (tp / (tp + fn) + tn / (tn + fp)) / 2, MCC = mcc)
}
scale_ghs <- function(gd, Zref) {   # GHS features: same columns, imputed + z-scored WITHIN GHS
  lapply(names(Zref), function(v) {
    f <- intersect(colnames(Zref[[v]]), colnames(gd))
    X <- as.data.frame(lapply(gd[, f, drop = FALSE], num)); rownames(X) <- as.character(gd$v11_sid01)
    scale(as.matrix(impute_rf(X)))
  })
}
hf_ghs <- ghs[which(ghs$hf_a3 == "HFPEF"), ]
hf_ghs$time <- num(hf_ghs$whf_strict_time); hf_ghs$event <- num(hf_ghs$whf_strict_event)

# 1) Olink harmonisation: raw protein medians MyoVasc vs GHS (HFpEF)
med_mv  <- sapply(prot_vars, function(p) median(num(hf_mv[[p]]),  na.rm = TRUE))
med_ghs <- sapply(prot_vars, function(p) median(num(hf_ghs[[p]]), na.rm = TRUE))
cat("Olink harmonisation: correlation of protein medians MyoVasc vs GHS =", round(cor(med_mv, med_ghs), 3),
    "| median |difference| =", round(median(abs(med_mv - med_ghs)), 2), "NPX\n")

val <- list()
for (s in names(res)) {
  Zmv <- res[[s]]$Z
  y   <- as.integer(res[[s]]$data$cluster[match(rownames(Zmv[[1]]), as.character(res[[s]]$data$v11_sid01))] == "High")  # High = 1
  gd  <- hf_ghs[hf_ghs$female == as.integer(s == "women"), ]
  Zg  <- scale_ghs(gd, Zmv); names(Zg) <- names(Zmv)
  feat <- unlist(lapply(names(Zmv), function(v) intersect(colnames(Zmv[[v]]), colnames(Zg[[v]]))))
  Xtr <- do.call(cbind, Zmv)[, feat]; Xte <- do.call(cbind, Zg)[, feat]
  cat("\n=====", s, "=====  features used:", length(feat), "of", sum(sapply(Zmv, ncol)), "\n")
  
  # 2) internal performance: 20 x stratified 10-fold CV, out-of-fold predictions pooled per repetition
  params <- list(objective = "binary:logistic", eval_metric = "auc", eta = 0.05, max_depth = 4,
                 subsample = 0.8, colsample_bytree = 0.8, min_child_weight = 1)
  set.seed(256)
  cvb <- xgb.cv(params, xgb.DMatrix(Xtr, label = y), nrounds = 2000, nfold = 5, stratified = TRUE,
                early_stopping_rounds = 50, verbose = 0)
  nr <- cvb$best_iteration
  cvm <- t(sapply(1:20, function(rep) {
    set.seed(rep); folds <- caret::createFolds(factor(y), k = 10); oof <- numeric(length(y))
    for (f in folds) oof[f] <- predict(xgb.train(params, xgb.DMatrix(Xtr[-f, ], label = y[-f]), nrounds = nr, verbose = 0),
                                       xgb.DMatrix(Xtr[f, , drop = FALSE]))
    if (rep == 1) { cat("Confusion matrix (repetition 1, threshold 0.5, High = positive):\n")
      print(table(predicted = ifelse(oof >= 0.5, "High", "Low"), true = ifelse(y == 1, "High", "Low"))) }
    metrics(oof, y)
  }))
  cat("CV performance, mean (2.5-97.5 percentile over 20 repetitions):\n")
  print(round(rbind(mean = colMeans(cvm), lo = apply(cvm, 2, quantile, .025), hi = apply(cvm, 2, quantile, .975)), 3))
  
  # 3) transfer to GHS
  mdl <- xgb.train(params, xgb.DMatrix(Xtr, label = y), nrounds = nr, verbose = 0)
  g <- gd[match(rownames(Xte), as.character(gd$v11_sid01)), ]
  g$p_high <- predict(mdl, xgb.DMatrix(Xte))
  g$cluster_tr <- factor(ifelse(g$p_high >= 0.5, "High", "Low"), levels = c("Low", "High"))
  
  # 4) de novo SNF in GHS
  set.seed(256)
  rg <- run_snf(gd); lab <- setNames(spectralClustering(rg$fused, K = 2), rownames(rg$fused))
  rate <- tapply(g$event[match(names(lab), as.character(g$v11_sid01))] == 1, lab, mean, na.rm = TRUE)
  g$cluster_dn <- factor(ifelse(lab[as.character(g$v11_sid01)] == names(which.max(rate)), "High", "Low"), levels = c("Low", "High"))
  smd <- function(Z, cl) colMeans(Z[cl == "High", , drop = FALSE]) - colMeans(Z[cl == "Low", , drop = FALSE])
  cl_mv <- ifelse(y == 1, "High", "Low")
  Zg_dn <- do.call(cbind, rg$Z); fdn <- intersect(feat, colnames(Zg_dn))
  cat("\nHigh-risk share: MyoVasc", round(100 * mean(y)), "% | GHS transferred", round(100 * mean(g$cluster_tr == "High")),
      "% | GHS de novo", round(100 * mean(g$cluster_dn == "High", na.rm = TRUE)), "%\n")
  cat("Profile correlation (High-minus-Low differences across features):\n",
      "  MyoVasc vs GHS transferred:", round(cor(smd(Xtr, cl_mv), smd(Xte, as.character(g$cluster_tr))), 2), "\n",
      "  MyoVasc vs GHS de novo:    ", round(cor(smd(Xtr[, fdn], cl_mv),
                                                 smd(Zg_dn[, fdn], as.character(g$cluster_dn[match(rownames(Zg_dn), as.character(g$v11_sid01))]))), 2), "\n")
  cat("Agreement transferred vs de novo labels:\n"); print(table(transferred = g$cluster_tr, de_novo = g$cluster_dn))
  
  # 5) outcomes in GHS
  for (cl in c("cluster_tr", "cluster_dn")) {
    m <- coxph(as.formula(paste("Surv(time, event) ~", cl, "+ age")), data = g)
    x <- summary(m)$conf.int[1, c(1, 3, 4)]
    lr <- 1 - pchisq(survdiff(as.formula(paste("Surv(time, event) ~", cl)), data = g)$chisq, 1)
    cat(sprintf("GHS %-11s n = %d, events High/Low = %d/%d, HR %.2f (%.2f-%.2f), log-rank p = %.2g\n", cl, sum(!is.na(g[[cl]])),
                sum(g$event[g[[cl]] == "High"], na.rm = TRUE), sum(g$event[g[[cl]] == "Low"], na.rm = TRUE), x[1], x[2], x[3], lr))
  }
  cat("Diabetes % High/Low - transferred:", round(100 * tapply(yes(g$diab), g$cluster_tr, mean)[c("High", "Low")], 1),
      "| de novo:", round(100 * tapply(yes(g$diab), g$cluster_dn, mean)[c("High", "Low")], 1), "\n")
  val[[s]] <- list(cv = cvm, model = mdl, features = feat, ghs = g, ghs_snf = rg)
}
saveRDS(val, file.path(project_dir, "WP4_GHS_validation.rds"))

# ===== BLOCK 3b =====
smd <- function(Z, cl) colMeans(Z[cl == "High", , drop = FALSE]) - colMeans(Z[cl == "Low", , drop = FALSE])
for (s in names(val)) {
  feat <- val[[s]]$features; g <- val[[s]]$ghs; rg <- val[[s]]$ghs_snf
  Xtr  <- do.call(cbind, res[[s]]$Z)[, feat]
  cl_mv <- as.character(res[[s]]$data$cluster[match(rownames(Xtr), as.character(res[[s]]$data$v11_sid01))])
  Zdn  <- do.call(cbind, rg$Z); f2 <- intersect(feat, colnames(Zdn))
  lab  <- as.character(g$cluster_dn[match(rownames(Zdn), as.character(g$v11_sid01))])
  r0   <- cor(smd(Xtr[, f2], cl_mv), smd(Zdn[, f2], lab))
  if (r0 < 0) g$cluster_dn <- factor(ifelse(g$cluster_dn == "High", "Low", "High"), levels = c("Low", "High"))   # outcome-blind orientation
  lab  <- as.character(g$cluster_dn[match(rownames(Zdn), as.character(g$v11_sid01))])
  cat("\n=====", s, "===== de novo profile correlation after orientation:", round(cor(smd(Xtr[, f2], cl_mv), smd(Zdn[, f2], lab)), 2),
      "| agreement transferred vs de novo:", round(100 * mean(g$cluster_tr == g$cluster_dn, na.rm = TRUE)), "%\n")
  m <- coxph(Surv(time, event) ~ cluster_dn + age, data = g); x <- summary(m)$conf.int[1, c(1, 3, 4)]
  cat(sprintf("De novo (profile-oriented): events High/Low = %d/%d, HR %.2f (%.2f-%.2f), log-rank p = %.2g\n",
              sum(g$event[g$cluster_dn == "High"], na.rm = TRUE), sum(g$event[g$cluster_dn == "Low"], na.rm = TRUE), x[1], x[2], x[3],
              1 - pchisq(survdiff(Surv(time, event) ~ cluster_dn, data = g)$chisq, 1)))
  cat("Diabetes % High/Low de novo:", round(100 * tapply(yes(g$diab), g$cluster_dn, mean)[c("High", "Low")], 1), "\n")
  cat("GHS cardiac characterisation by transferred cluster (High | Low):\n")
  for (v in c("ef", "ee", "lvmi", "bnp")) cat(sprintf("  %-5s %s | %s\n", v, med_iqr(g[[v]][g$cluster_tr == "High"]), med_iqr(g[[v]][g$cluster_tr == "Low"])))
  cat(sprintf("  EF>=50%%: %.0f%% | %.0f%%   chf history: %.0f%% | %.0f%%\n",
              100 * mean(num(g$ef[g$cluster_tr == "High"]) >= 50, na.rm = TRUE), 100 * mean(num(g$ef[g$cluster_tr == "Low"]) >= 50, na.rm = TRUE),
              100 * mean(yes(g$chf[g$cluster_tr == "High"])), 100 * mean(yes(g$chf[g$cluster_tr == "Low"]))))
  val[[s]]$ghs <- g
}
# Sex-pooled GHS estimate (stratified by sex)
gg <- rbind(val$women$ghs[, c("time", "event", "age", "female", "cluster_tr", "cluster_dn")],
            val$men$ghs[,   c("time", "event", "age", "female", "cluster_tr", "cluster_dn")])
for (cl in c("cluster_tr", "cluster_dn")) {
  x <- summary(coxph(as.formula(paste("Surv(time, event) ~", cl, "+ age + strata(female)")), data = gg))
  cat(sprintf("\nGHS pooled (strata = sex) %s: HR %.2f (%.2f-%.2f), p = %.2g, events = %d\n", cl,
              x$conf.int[1, 1], x$conf.int[1, 3], x$conf.int[1, 4], x$coef[1, 5], sum(gg$event, na.rm = TRUE)))
}
cat("\nGHS HFpEF without Olink (excluded now, included in old analysis): women",
    sum(ghs_bas$hf_a3 == "HFPEF" & ghs_bas$sex == "Women" & !ghs_ol, na.rm = TRUE), "| men",
    sum(ghs_bas$hf_a3 == "HFPEF" & ghs_bas$sex == "Men" & !ghs_ol, na.rm = TRUE), "\n")
saveRDS(val, file.path(project_dir, "WP4_GHS_validation.rds"))

# ===== BLOCK 4: diabetes + parsimonious comorbidity sensitivity =====
cat("Medication-related variables:\n"); print(grep("sglt|metf|insul|a10|antidiab|glp|diuret|atc|ace|arb|mra|statin", names(mv), value = TRUE, ignore.case = TRUE))

hr_txt <- function(m, term) { x <- summary(m)$conf.int[term, c(1, 3, 4)]
sprintf("HR %.2f (%.2f-%.2f), p = %.2g, n = %d", x[1], x[2], x[3], summary(m)$coef[term, 5], m$n) }
cx <- function(rhs, d) coxph(as.formula(paste("Surv(time, event) ~", rhs)), data = d)
dia <- list()
for (s in names(res)) {
  d <- res[[s]]$data
  d$dm <- yes(d$diab); d$egfr_n <- num(d$egfr); d$bmi_n <- num(d$bmi); d$af <- yes(d$afib)
  d$hba1c <- num(d$v34_lbl50)
  cat("\n=====", s, "=====\n")
  cat("Cluster x diabetes, n (events):\n")
  print(matrix(paste0(table(d$cluster, d$dm), " (", tapply(d$event, list(d$cluster, d$dm), sum), ")"), 2,
               dimnames = list(cluster = levels(d$cluster), diabetes = c("no", "yes"))), quote = FALSE)
  cat("HbA1c median High/Low:", round(tapply(d$hba1c, d$cluster, median, na.rm = TRUE)[c("High", "Low")], 1),
      "| TyG:", round(tapply(num(d$tyg), d$cluster, median, na.rm = TRUE)[c("High", "Low")], 2), "\n")
  
  cat("\n-- One-at-a-time adjustment (base: age + NYHA) --\n")
  for (adj in c("", "egfr_n", "dm", "bmi_n", "af")) {
    rhs <- paste("cluster + age + nyha", if (adj != "") paste("+", adj) else "")
    cat(sprintf("  + %-7s %s\n", ifelse(adj == "", "(none)", adj), hr_txt(cx(rhs, d), "clusterHigh")))
  }
  cat("\n-- Diabetes vs cluster --\n")
  cat("  Diabetes alone (age):        ", hr_txt(cx("dm + age", d), "dmTRUE"), "\n")
  mb <- cx("cluster + dm + age", d)
  cat("  Both - cluster:              ", hr_txt(mb, "clusterHigh"), "\n")
  cat("  Both - diabetes:             ", hr_txt(mb, "dmTRUE"), "\n")
  cat("  Interaction p (cluster x dm):", signif(summary(cx("cluster * dm + age", d))$coef["clusterHigh:dmTRUE", 5], 2), "\n")
  cat("  Within NON-diabetics:        ", hr_txt(cx("cluster + age", d[!d$dm, ]), "clusterHigh"), "\n")
  cat("  Within diabetics:            ", hr_txt(cx("cluster + age", d[d$dm, ]), "clusterHigh"), "\n")
  cat("  C-index age+dm / age+cluster / age+dm+cluster:",
      round(concordance(cx("age + dm", d))$concordance, 3), "/", round(concordance(cx("age + cluster", d))$concordance, 3), "/",
      round(concordance(cx("age + dm + cluster", d))$concordance, 3), "\n")
  
  # Glycaemia-free clustering: remove HbA1c and glucose from the inputs
  views_ng <- list(proteins = prot_vars, clinical = setdiff(clin_vars, c("v34_lbl50", "glucose")))
  views_bak <- views; views <- views_ng; set.seed(256)
  d0 <- hf_mv[hf_mv$female == as.integer(s == "women"), ]
  rng <- run_snf(d0); views <- views_bak
  lab <- setNames(spectralClustering(rng$fused, K = 2), rownames(rng$fused))
  main <- setNames(as.character(d$cluster), as.character(d$v11_sid01))[names(lab)]
  tb <- table(lab, main); hi <- rownames(tb)[which.max(tb[, "High"] / rowSums(tb))]
  d$cluster_ng <- factor(ifelse(lab[as.character(d$v11_sid01)] == hi, "High", "Low"), levels = c("Low", "High"))
  cat("\n-- Glycaemia-free clustering (no HbA1c, no glucose) --\n")
  cat("  Agreement with main clusters:", round(100 * mean(d$cluster_ng == d$cluster)), "%\n")
  cat("  Diabetes % High/Low:", round(100 * tapply(d$dm, d$cluster_ng, mean)[c("High", "Low")], 1), "\n")
  cat("  WHF (age + NYHA):   ", hr_txt(cx("cluster_ng + age + nyha", d), "cluster_ngHigh"), "\n")
  dia[[s]] <- d
}

# Cohort meta-analysis per sex (MyoVasc + GHS transferred, age-adjusted, fixed effect)
for (s in names(res)) {
  b <- sapply(list(MyoVasc = cx("cluster + age", res[[s]]$data), GHS = cx("cluster_tr + age", val[[s]]$ghs)),
              function(m) c(coef(m)[1], sqrt(vcov(m)[1, 1])))
  w <- 1 / b[2, ]^2; est <- sum(w * b[1, ]) / sum(w); se <- sqrt(1 / sum(w))
  Q <- sum(w * (b[1, ] - est)^2)
  cat(sprintf("\nMeta-analysis %s: HR %.2f (%.2f-%.2f), p = %.2g | heterogeneity p = %.2g\n", s,
              exp(est), exp(est - 1.96 * se), exp(est + 1.96 * se), 2 * pnorm(-abs(est / se)), 1 - pchisq(Q, 1)))
}
saveRDS(dia, file.path(project_dir, "WP5_diabetes.rds"))

# ===== BLOCK 4b: treatment adjustment + HOMA-IR =====
cat("Other cardiovascular drug variables (C03 diuretics, C07 beta-blockers, C09 RAAS, C10 lipid):\n")
print(grep("^c0[3-9]|^c10|loop|mra|betabl|diur", names(mv), value = TRUE, ignore.case = TRUE)[1:40])

meds <- c(any_glucose_lowering = "a10", insulin_tx = "a10a", metformin = "a10ba", DPP4i = "a10bh",
          GLP1RA = "a10bj", SGLT2i = "a10bk", ACEi = "acei", ARB = "arbs")
for (s in names(dia)) {
  d <- dia[[s]]
  d$homa_ir <- num(d$insulin) * num(d$glucose) / 405            # insulin [uU/mL] x glucose [mg/dL] / 405
  d$glc_tx  <- yes(d$a10); d$raas <- yes(d$acei) | yes(d$arbs)
  cat("\n=====", s, "=====\n")
  cat("Treatment prevalence % (High | Low):\n")
  for (m in names(meds)) if (meds[m] %in% names(d))
    cat(sprintf("  %-22s %5.1f | %5.1f\n", m, 100 * mean(yes(d[[meds[m]]][d$cluster == "High"])), 100 * mean(yes(d[[meds[m]]][d$cluster == "Low"]))))
  cat("HOMA-IR median (High | Low):", med_iqr(d$homa_ir[d$cluster == "High"]), "|", med_iqr(d$homa_ir[d$cluster == "Low"]),
      " p =", signif(wilcox.test(homa_ir ~ cluster, data = d)$p.value, 2), "\n")
  cat("Among NON-diabetics, HOMA-IR (High | Low):", med_iqr(d$homa_ir[d$cluster == "High" & !d$dm]), "|",
      med_iqr(d$homa_ir[d$cluster == "Low" & !d$dm]), "\n")
  cat("Cluster HR + glucose-lowering treatment:", hr_txt(cx("cluster + age + nyha + glc_tx", d), "clusterHigh"), "\n")
  cat("Cluster HR + RAAS inhibitor:            ", hr_txt(cx("cluster + age + nyha + raas", d), "clusterHigh"), "\n")
  cat("Cluster HR + log HOMA-IR:               ", hr_txt(cx("cluster + age + nyha + log(homa_ir)", d[!is.na(d$homa_ir) & d$homa_ir > 0, ]), "clusterHigh"), "\n")
  dia[[s]] <- d
}
saveRDS(dia, file.path(project_dir, "WP5_diabetes.rds"))

# ===== BLOCK 5: sex interactions, pooled comparator, protein x sex =====
pool <- rbind(dia$women, dia$men[, names(dia$women)])
pool$hi <- as.integer(pool$cluster == "High")
rhr <- function(m, term) { b <- coef(m)[term]; se <- sqrt(vcov(m)[term, term])
sprintf("ratio women:men %.2f (%.2f-%.2f), p-interaction = %.2g", exp(b), exp(b - 1.96 * se), exp(b + 1.96 * se), summary(m)$coef[term, 5]) }

cat("== Outcome interactions (MyoVasc HFpEF, pooled) ==\n")
cat("Cluster x sex (age + NYHA):  ", rhr(cx("hi * female + age + nyha", pool), "hi:female"), "\n")
cat("Diabetes x sex (age + NYHA): ", rhr(cx("dm * female + age + nyha", pool), "dmTRUE:female"), "\n")
cat("Diabetes x sex, + eGFR:      ", rhr(cx("dm * female + age + nyha + egfr_n", pool), "dmTRUE:female"), "\n")
gpool <- rbind(val$women$ghs, val$men$ghs[, names(val$women$ghs)]); gpool$dm <- yes(gpool$diab)
cat("GHS diabetes x sex (age):    ", rhr(cx("dm * female + age", gpool), "dmTRUE:female"), "\n")
cat("GHS diabetes HR women / men: ",
    round(exp(coef(cx("dm + age", gpool[gpool$female == 1, ]))[1]), 2), "/",
    round(exp(coef(cx("dm + age", gpool[gpool$female == 0, ]))[1]), 2), "\n")

# Pooled SNF comparator (sexes together, variables z-scored within sex)
set.seed(256)
Zp <- lapply(names(views), function(v) {
  zz <- lapply(c(1, 0), function(f) prep_view(hf_mv[hf_mv$female == f, ], intersect(views[[v]], names(hf_mv)), v))
  fc <- Reduce(intersect, lapply(zz, colnames)); do.call(rbind, lapply(zz, function(z) z[, fc, drop = FALSE]))
})
ids <- Reduce(intersect, lapply(Zp, rownames)); Zp <- lapply(Zp, function(z) z[ids, , drop = FALSE])
Wp  <- SNF(lapply(Zp, function(z) affinityMatrix(dist2(z, z)^(1/2), K = 20, sigma = 0.5)), 20, 20)
lab_p <- setNames(spectralClustering(Wp, K = 2), ids)
for (s in names(dia)) {
  d <- dia[[s]]; lp <- lab_p[as.character(d$v11_sid01)]
  cat(sprintf("\nPooled vs sex-specific clustering, %s: agreement %.0f%%, NMI %.2f\n", s,
              100 * max(mean(lp == as.integer(d$cluster)), mean(lp != as.integer(d$cluster))), calNMI(lp, as.integer(d$cluster))))
}

# Protein x sex interactions (FDR across 92 proteins)
pz <- pool; for (p in prot_vars) pz[[p]] <- as.numeric(scale(num(pz[[p]])))
int_tab <- do.call(rbind, lapply(prot_vars, function(p) {
  m1 <- glm(as.formula(paste("hi ~", p, "* female + age")), data = pz, family = binomial)
  m2 <- lm(as.formula(paste(p, "~ dm * female + age + egfr_n")), data = pz)
  data.frame(protein = p,
             OR_women = exp(sum(coef(m1)[c(p, paste0(p, ":female"))])), OR_men = exp(coef(m1)[p]),
             p_int_cluster = summary(m1)$coef[paste0(p, ":female"), 4],
             dm_effect_women = sum(coef(m2)[c("dmTRUE", "dmTRUE:female")]), dm_effect_men = coef(m2)["dmTRUE"],
             p_int_diabetes = summary(m2)$coef["dmTRUE:female", 4])
}))
int_tab$fdr_cluster  <- p.adjust(int_tab$p_int_cluster,  "BH")
int_tab$fdr_diabetes <- p.adjust(int_tab$p_int_diabetes, "BH")
nm <- attributes(dall1)$nams_prot; names(nm) <- attributes(dall1)$vars_prot
int_tab$name <- nm[int_tab$protein]
cat("\nProteins with sex-different association to HIGH-RISK cluster (FDR < 0.05):", sum(int_tab$fdr_cluster < 0.05), "\n")
print(head(int_tab[order(int_tab$p_int_cluster), c("name", "OR_women", "OR_men", "p_int_cluster", "fdr_cluster")], 10), digits = 2)
cat("\nProteins with sex-different association to DIABETES (FDR < 0.05):", sum(int_tab$fdr_diabetes < 0.05), "\n")
print(head(int_tab[order(int_tab$p_int_diabetes), c("name", "dm_effect_women", "dm_effect_men", "p_int_diabetes", "fdr_diabetes")], 10), digits = 2)
write.csv(int_tab, file.path(project_dir, "WP6_protein_sex_interactions.csv"), row.names = FALSE)

# ---- A) which reproductive variables exist? ----
cand <- grep("menop|menar|mens|period|regel|letzt|hrt|hormon|oestr|estr|schwang|geburt|parit|gravid|ovar|hyster|gestat|praekl|preecl|wechsel|klimak",
             names(mv), value = TRUE, ignore.case = TRUE)
print(cand)
for (v in head(cand, 25)) {                     # quick look: values in HFpEF women
  x <- res$women$data[[v]]
  cat("\n", v, ":"); if (is.numeric(num(x)) && length(unique(na.omit(x))) > 8) print(summary(num(x))) else print(table(x, useNA = "ifany"))
}

# ---- B) protein concordance women vs men ----
conc <- int_tab[, c("name", "OR_women", "OR_men", "p_int_cluster", "fdr_cluster")]
conc$logOR_w <- log(conc$OR_women); conc$logOR_m <- log(conc$OR_men)
cat("\nCorrelation of protein-cluster associations women vs men (log OR):", round(cor(conc$logOR_w, conc$logOR_m), 2), "\n")
top_w <- head(shap$women$label, 20); top_m <- head(shap$men$label, 20)
cat("Top-20 SHAP shared:", paste(intersect(top_w, top_m), collapse = ", "), "\n")
cat("Top-20 only women: ", paste(setdiff(top_w, top_m), collapse = ", "), "\n")
cat("Top-20 only men:   ", paste(setdiff(top_m, top_w), collapse = ", "), "\n")

library(ggplot2)
conc$label <- ifelse(conc$p_int_cluster < 0.05, trimws(sub(".*\\((.*)\\).*", "\\1", conc$name)), NA)
p_conc <- ggplot(conc, aes(logOR_m, logOR_w)) +
  geom_abline(linetype = 2, colour = "grey50") + geom_hline(yintercept = 0, colour = "grey85") + geom_vline(xintercept = 0, colour = "grey85") +
  geom_point(aes(colour = p_int_cluster < 0.05), size = 2) +
  geom_text(aes(label = label), vjust = -0.7, size = 3, na.rm = TRUE) +
  scale_colour_manual(values = c(`FALSE` = "grey60", `TRUE` = "#D85A30"), labels = c("concordant", "nominal sex interaction (p<0.05)"), name = NULL) +
  labs(x = "Men: log OR for high-risk cluster (per SD)", y = "Women: log OR for high-risk cluster (per SD)") +
  theme_bw(base_size = 12) + theme(legend.position = "bottom")
ggsave(file.path(project_dir, "Fig_protein_concordance.png"), p_conc, width = 6, height = 6, dpi = 300)
print(p_conc)

# ===== BLOCK 7: reproductive factors in HFpEF women =====
d <- dia$women
cat("Raw coding check:\n")
for (v in c("f002", "f005", "f007", "f011", "f015", "f022")) { cat("\n", v, ":"); print(table(d[[v]], useNA = "ifany")) }

rng   <- function(x, lo, hi) { x <- num(x); x[x < lo | x > hi] <- NA; x }      # removes codes like 88/99
yes2  <- function(x) ifelse(is.na(x), NA, x %in% c(1, "1", "yes", "Yes", "ja", "Ja", TRUE))
d$menarche   <- rng(d$f001a, 8, 20)
d$meno_age   <- rng(d$f013_1, 30, 65)
d$early_meno <- d$meno_age < 45
d$repro_span <- d$meno_age - d$menarche
d$yrs_since_meno <- ifelse(num(d$age) >= d$meno_age, num(d$age) - d$meno_age, NA)
d$ever_preg  <- yes2(d$f002)
d$gravidity  <- rng(d$f003, 0, 20)
d$parity     <- rng(d$f004, 0, 15)
d$parity3    <- d$parity >= 3
d$abortion   <- yes2(d$f005)
d$breastfeed <- rng(d$f006, 0, 240)
d$oc_ever    <- yes2(d$f007)
d$hysterect  <- yes2(d$f015)
d$homa_ok    <- ifelse(yes(d$a10a), NA, d$homa_ir)          # HOMA-IR only without insulin therapy

cont <- c(menarche = "Age at menarche", meno_age = "Age at menopause", repro_span = "Reproductive lifespan (y)",
          yrs_since_meno = "Years since menopause", gravidity = "Pregnancies (n)", parity = "Births (n)", breastfeed = "Breastfeeding (months)")
bin  <- c(early_meno = "Early menopause (<45 y)", ever_preg = "Ever pregnant", parity3 = "Births >= 3",
          abortion = "Any abortion", oc_ever = "Ever oral contraceptives", hysterect = "Hysterectomy")

cat("\n===== Reproductive factors by cluster (HFpEF women) =====\n")
for (v in names(cont)) cat(sprintf("%-28s High %-18s Low %-18s p = %.2g  (n = %d)\n", cont[v], med_iqr(d[[v]][d$cluster == "High"]),
                                   med_iqr(d[[v]][d$cluster == "Low"]), wilcox.test(d[[v]] ~ d$cluster)$p.value, sum(!is.na(d[[v]]))))
for (v in names(bin)) { x <- d[[v]]
cat(sprintf("%-28s High %5.1f%%            Low %5.1f%%            p = %.2g  (n = %d)\n", bin[v],
            100 * mean(x[d$cluster == "High"], na.rm = TRUE), 100 * mean(x[d$cluster == "Low"], na.rm = TRUE),
            fisher.test(table(x, d$cluster))$p.value, sum(!is.na(x)))) }

# Associations (each factor separately, age-adjusted)
ef <- function(m, term) { b <- coef(m)[term]; se <- sqrt(vcov(m)[term, term]); sprintf("%.2f (%.2f-%.2f), p = %.2g", exp(b), exp(b - 1.96 * se), exp(b + 1.96 * se), summary(m)$coef[term, ncol(summary(m)$coef)]) }
cat("\n===== Associations (age-adjusted; continuous per 5 years / per unit) =====\n")
cat(sprintf("%-28s %-30s %-30s %-30s %-30s\n", "Factor", "OR High cluster", "OR diabetes", "ratio HOMA-IR (+BMI)", "HR WHF"))
d$meno_age5 <- d$meno_age / 5; d$repro_span5 <- d$repro_span / 5
for (v in c("early_meno", "meno_age5", "repro_span5", "parity", "parity3", "hysterect", "oc_ever")) {
  dd <- d[!is.na(d[[v]]), ]; f <- function(lhs, extra = "") as.formula(paste(lhs, "~", v, "+ age", extra))
  term <- if (is.logical(dd[[v]])) paste0(v, "TRUE") else v
  o1 <- ef(glm(f("I(cluster == 'High')"), data = dd, family = binomial), term)
  o2 <- ef(glm(f("dm"), data = dd, family = binomial), term)
  o3 <- ef(lm(f("log(homa_ok)", "+ num(bmi)"), data = dd[!is.na(dd$homa_ok) & dd$homa_ok > 0, ]), term)
  o4 <- ef(coxph(f("Surv(time, event)"), data = dd), term)
  cat(sprintf("%-28s %-30s %-30s %-30s %-30s\n", v, o1, o2, o3, o4))
}
cat("\nCluster HR within early vs not-early menopause (exploratory):\n")
for (g in c(TRUE, FALSE)) { dd <- d[which(d$early_meno == g), ]
cat(sprintf("  early_meno = %-5s n = %d, events = %d: %s\n", g, nrow(dd), sum(dd$event), hr_txt(cx("cluster + age", dd), "clusterHigh"))) }
dia$women <- d
saveRDS(dia, file.path(project_dir, "WP5_diabetes.rds"))


##### Tables and Figures
# ===== STEP 1: Supplementary Tables 1 and 2 =====
source("T:\\R\\Gieswinkel\\R_my_functions\\baseline13.r")
setwd("T:\\R\\Gieswinkel\\R_my_functions"); suppressWarnings(eval(parse(file = "ade.fun.run.r"))); setwd(project_dir)

prep_tab <- function(d) {
  d$Cluster    <- factor(as.character(d$cluster), levels = c("High", "Low"))
  d$homa_ir    <- num(d$insulin) * num(d$glucose) / 405
  d$glc_tx     <- if ("a10"  %in% names(d)) ifelse(yes(d$a10),  "yes", "no") else NA
  d$insulin_tx <- if ("a10a" %in% names(d)) ifelse(yes(d$a10a), "yes", "no") else NA
  for (v in c("age", "maggic", "bmi", "v34_lbl50", "tyg", "homa_ir", "sbp", "dbp", "ef", "ee", "lvmi", "rwt", "bnp", "tni", "egfr", "crp_ln"))
    if (v %in% names(d)) d[[v]] <- num(d[[v]])
  if ("nyha" %in% names(d)) d$nyha_f <- factor(d$nyha)
  d
}
vars <- c("age", "maggic", "bmi", "v34_lbl50", "tyg", "homa_ir", "hyper", "diab", "glc_tx", "insulin_tx", "dyslip", "adipos", "nic",
          "ckd", "afib", "cad", "mi", "stroke", "pad", "cancer", "nyha_f", "sbp", "dbp", "ef", "ee", "lvmi", "rwt", "bnp", "tni", "egfr", "crp_ln")
nams <- c("Age [y]", "MAGGIC score", "BMI [kg/m2]", "HbA1c [%]", "TyG index", "HOMA-IR", "Arterial hypertension", "Diabetes mellitus",
          "Glucose-lowering therapy", "Insulin therapy", "Dyslipidaemia", "Obesity", "Current smoking", "Chronic kidney disease",
          "Atrial fibrillation", "Coronary artery disease", "History of MI", "History of stroke", "Peripheral artery disease",
          "History of cancer", "NYHA class", "SBP [mmHg]", "DBP [mmHg]", "LVEF [%]", "E/e'", "LV mass index [g/m2]", "RWT",
          "NT-proBNP [pg/mL]", "Troponin I [ng/L]", "eGFR [mL/min/1.73m2]", "CRP (log)")
names(nams) <- vars

tabs <- list()
for (s in c("women", "men")) {
  for (co in c("MyoVasc", "GHS")) {
    d <- if (co == "MyoVasc") dia[[s]] else { g <- val[[s]]$ghs; g$cluster <- g$cluster_tr; g }
    d <- prep_tab(d)
    v <- vars[vars %in% names(d) & sapply(vars, function(x) x %in% names(d) && !all(is.na(d[[x]])))]
    cat("\n=====", s, co, "| n High/Low:", paste(table(d$Cluster), collapse = "/"), "| variables:", length(v), "of", length(vars), "\n")
    if (length(setdiff(vars, v))) cat("  not available:", paste(setdiff(vars, v), collapse = ", "), "\n")
    bt <- baseline13(v, unname(nams[v]), data = d, skeew = 1, group = "Cluster", pv = TRUE, style = 2, minuniq = 3)
    tabs[[paste(s, co)]] <- bt
    write.csv(as.data.frame(bt), file.path(project_dir, sprintf("SuppTable_baseline_%s_%s.csv", s, co)), row.names = FALSE)
  }
}
print(head(as.data.frame(tabs[["women MyoVasc"]]), 12))

# Missing value for the text: concordance of protein-cluster associations women vs men
cat("\nConcordance r (log OR women vs men):", round(cor(log(int_tab$OR_women), log(int_tab$OR_men)), 2), "\n")


# ===== STEP 2 (fix): p-value formatting + adjusted HR in each panel =====
fmt_p <- function(p) if (p < 0.0001) "p < 0.0001" else if (p < 0.001) sprintf("p = %.4f", p) else if (p < 0.01) sprintf("p = %.3f", p) else sprintf("p = %.2f", p)
cif_plot <- function(d, title, cl = "cluster", adj = "age") {
  d$grp <- factor(ifelse(d[[cl]] == "High", "High-risk", "Low-risk"), levels = c("Low-risk", "High-risk"))
  fit <- survfit(Surv(time, event) ~ grp, data = d)
  lr  <- 1 - pchisq(survdiff(Surv(time, event) ~ grp, data = d)$chisq, 1)
  hr  <- summary(coxph(as.formula(paste("Surv(time, event) ~ grp +", adj)), data = d))$conf.int["grpHigh-risk", c(1, 3, 4)]
  lab <- sprintf("Log-rank %s\nHR %.2f (95%% CI %.2f–%.2f)\nadjusted for %s", fmt_p(lr), hr[1], hr[2], hr[3], gsub("\\+", "and", gsub("nyha", "NYHA", adj)))
  ggsurvplot(fit, data = d, fun = "event", conf.int = TRUE, palette = c("#2E86C1", "#C0392B"), legend.labs = c("Low-risk", "High-risk"),
             legend.title = "", xlim = c(0, 4), break.time.by = 1, ylim = c(0, 0.5),
             xlab = "Time (years)", ylab = "Cumulative incidence of WHF", title = title,
             risk.table = TRUE, risk.table.height = 0.25, risk.table.y.text = FALSE, tables.y.text = FALSE,
             pval = lab, pval.size = 3.6, pval.coord = c(0.05, 0.43),
             ggtheme = theme_classic(base_size = 12) + theme(plot.title = element_text(hjust = 0.5, face = "bold")),
             tables.theme = theme_cleantable(), fontsize = 3.5)
}
save_pair(cif_plot(res$women$data, sprintf("Women (n = %d)", nrow(res$women$data)), adj = "age + nyha"),
          cif_plot(res$men$data,   sprintf("Men (n = %d)",   nrow(res$men$data)),   adj = "age + nyha"), "Figure1_MyoVasc_CIF")
save_pair(cif_plot(val$women$ghs, sprintf("Women (n = %d)", nrow(val$women$ghs)), "cluster_tr"),
          cif_plot(val$men$ghs,   sprintf("Men (n = %d)",   nrow(val$men$ghs)),   "cluster_tr"), "Figure2_GHS_CIF")

# ===== STEP 3: Figure 3 - SHAP =====
for (s in names(shap)) {
  shap[[s]]$label <- gsub("TGF-alpha", "TGF-α", shap[[s]]$label)
  shap[[s]]$label[shap[[s]]$label == "K"] <- "Potassium"
}
library(patchwork)
shared <- intersect(head(shap$women$label, 20), head(shap$men$label, 20))
shap_plot <- function(s, title) {
  d <- head(shap[[s]], 20)
  d$type <- ifelse(d$feature %in% prot_vars, "Protein", "Clinical")
  d$label <- ifelse(d$label %in% shared, paste0(d$label, " *"), d$label)
  d$label <- factor(d$label, levels = rev(d$label))
  d$dir <- ifelse(d$direction == "higher in High", "↑", "↓")
  ggplot(d, aes(mean_abs_SHAP, label, fill = type)) + geom_col(width = 0.7) +
    geom_text(aes(label = dir), hjust = -0.4, size = 3.5) +
    scale_fill_manual(values = c(Protein = "#2E86C1", Clinical = "#7DCEA0"), name = NULL) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
    labs(x = "Mean |SHAP value|", y = NULL, title = title) +
    theme_classic(base_size = 11) + theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.position = "bottom")
}
fig3 <- shap_plot("women", "Women") + shap_plot("men", "Men") + plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
fig3 <- fig3 + plot_annotation(caption = "* among the top 20 in both sexes;  ↑ higher / ↓ lower in the high-risk cluster")
ggsave(file.path(project_dir, "Figure3_SHAP.png"),  fig3, width = 10, height = 6.5, dpi = 300)
ggsave(file.path(project_dir, "Figure3_SHAP.tiff"), fig3, width = 10, height = 6.5, dpi = 300, compression = "lzw")

# ===== STEP 4: Figures 4 and 5 - protein boxplots =====
lab_prot <- gsub("TGF-alpha", "TGF-α", trimws(sub(".*\\((.*)\\).*", "\\1", nm)))
prot_box <- function(s, fig_name) {
  d <- res[[s]]$data
  d$grp <- factor(ifelse(d$cluster == "High", "High-risk", "Low-risk"), levels = c("Low-risk", "High-risk"))
  # all 92 proteins: Wilcoxon + BH
  tt <- do.call(rbind, lapply(prot_vars, function(p) { x <- num(d[[p]])
  data.frame(protein = p, label = lab_prot[p], median_high = median(x[d$grp == "High-risk"], na.rm = TRUE),
             median_low = median(x[d$grp == "Low-risk"], na.rm = TRUE), p = wilcox.test(x ~ d$grp)$p.value) }))
  tt$p_adj <- p.adjust(tt$p, "BH"); tt$diff <- tt$median_high - tt$median_low
  write.csv(tt[order(tt$p_adj), ], file.path(project_dir, sprintf("SuppTable_protein_differences_%s.csv", s)), row.names = FALSE)
  # top 15 proteins by SHAP
  top <- head(shap[[s]]$feature[shap[[s]]$feature %in% prot_vars], 15)
  long <- do.call(rbind, lapply(top, function(p) data.frame(protein = lab_prot[p], grp = d$grp, npx = num(d[[p]]))))
  long$protein <- factor(long$protein, levels = lab_prot[top])
  pl <- tt[match(top, tt$protein), ]
  pl$txt <- ifelse(pl$p_adj < 0.0001, "p.adj < 0.0001", sprintf("p.adj = %.2g", pl$p_adj))
  pl$protein <- factor(lab_prot[top], levels = lab_prot[top])
  pl$y <- sapply(top, function(p) { x <- num(d[[p]]); max(x, na.rm = TRUE) + 0.06 * diff(range(x, na.rm = TRUE)) })
  g <- ggplot(long, aes(grp, npx, fill = grp)) +
    geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.85) +
    geom_jitter(width = 0.15, size = 0.5, alpha = 0.35) +
    geom_text(data = pl, aes(x = 1.5, y = y, label = txt), inherit.aes = FALSE, size = 2.8, vjust = 0) +
    facet_wrap(~ protein, scales = "free_y", ncol = 5) +
    scale_fill_manual(values = c("Low-risk" = "#2E86C1", "High-risk" = "#C0392B"), name = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    labs(x = NULL, y = "Protein level (NPX)") +
    theme_bw(base_size = 10) + theme(legend.position = "bottom", strip.background = element_rect(fill = "grey92"),
                                     strip.text = element_text(face = "bold"), axis.text.x = element_blank(), axis.ticks.x = element_blank())
  ggsave(file.path(project_dir, paste0(fig_name, ".png")),  g, width = 10, height = 7, dpi = 300)
  ggsave(file.path(project_dir, paste0(fig_name, ".tiff")), g, width = 10, height = 7, dpi = 300, compression = "lzw")
  cat(s, ": proteins different at FDR < 0.05:", sum(tt$p_adj < 0.05), "of 92 | higher in High:", sum(tt$p_adj < 0.05 & tt$diff > 0),
      "| lower in High:", sum(tt$p_adj < 0.05 & tt$diff < 0), "\n")
}
prot_box("women", "Figure4_proteins_women")
prot_box("men",   "Figure5_proteins_men")

# ===== STEP 5 (fix): Figure 6 with direction colours and shared scales =====
all_g <- do.call(rbind, lapply(names(rna_out), function(s) rna_out[[s]]$gsea))
fdr_lim <- range(-log10(all_g$p.adjust))
rna_plot <- function(s, title, n_each = 8) {
  g <- rna_out[[s]]$gsea
  up <- head(g[g$NES > 0, ][order(g$p.adjust[g$NES > 0]), ], n_each)
  dn <- head(g[g$NES < 0, ][order(g$p.adjust[g$NES < 0]), ], n_each)
  d  <- rbind(up, dn)
  d$Description <- ifelse(nchar(d$Description) > 55, paste0(substr(d$Description, 1, 52), "…"), d$Description)
  d$Description <- factor(d$Description, levels = d$Description[order(d$NES)])
  d$dir <- factor(ifelse(d$NES > 0, "Higher in high-risk", "Lower in high-risk"), levels = c("Lower in high-risk", "Higher in high-risk"))
  ggplot(d, aes(NES, Description, size = -log10(p.adjust), colour = dir)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey60") + geom_point() +
    scale_colour_manual(values = c("Lower in high-risk" = "#2E86C1", "Higher in high-risk" = "#C0392B"), name = NULL, drop = FALSE) +
    scale_size(range = c(2, 6), limits = fdr_lim, name = expression(-log[10]~FDR)) +
    scale_x_continuous(limits = c(-3.3, 3.3)) +
    labs(x = "Normalized enrichment score (NES)", y = NULL, title = title) +
    theme_bw(base_size = 10) + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
}
fig6 <- (rna_plot("women", "Women") + rna_plot("men", "Men")) + plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
fig6 <- fig6 + plot_annotation(caption = "Preranked GSEA (Reactome) on leukocyte-adjusted t-statistics; top 8 pathways per direction (FDR < 0.05).")
ggsave(file.path(project_dir, "Figure6_RNA_pathways.png"),  fig6, width = 13, height = 6.5, dpi = 300)
ggsave(file.path(project_dir, "Figure6_RNA_pathways.tiff"), fig6, width = 13, height = 6.5, dpi = 300, compression = "lzw")

# ===== STEP 6: Supplementary Table 4 (classifier), Table 6 (reproductive), Figure 7 =====
# Supp Table 4: cross-validated classifier performance + confusion matrix (repetition-averaged)
st4 <- do.call(rbind, lapply(names(val), function(s) { m <- val[[s]]$cv
data.frame(sex = s, metric = colnames(m), mean = round(colMeans(m), 3),
           lower_2.5 = round(apply(m, 2, quantile, .025), 3), upper_97.5 = round(apply(m, 2, quantile, .975), 3)) }))
write.csv(st4, file.path(project_dir, "SuppTable4_classifier_performance.csv"), row.names = FALSE)

# Supp Table 6: reproductive factors in women by cluster (descriptive + age-adjusted OR for high-risk cluster)
Sys.setlocale("LC_ALL", "GER")
setwd("T:\\Myovasc\\R_SQL"); source("Read_SQL_MyoVasc_BL_v11.R")
dall1 <- read_SQL_MyoVasc_BL(proteins = TRUE, lipids = TRUE, user = "ese0e", ATC_7_digits = "B01AC06")
setwd(project_dir)
res <- readRDS(file.path(project_dir, "WP1_discovery.rds"))          # cluster labels
cl  <- setNames(as.character(res$women$data$cluster), as.character(res$women$data$v11_sid01))
d   <- dall1[as.character(dall1$v11_sid01) %in% names(cl), ]
d$cluster <- factor(cl[as.character(d$v11_sid01)], levels = c("Low", "High"))
d$age <- num(d$age)
cat("Women:", nrow(d), "| High/Low:", paste(table(d$cluster), collapse = "/"), "\n")
#reproductive + hormonal factors with completeness =====
rng <- function(x, lo, hi) { x <- num(x); x[x < lo | x > hi] <- NA; x }
yn  <- function(x) { x <- num(x); ifelse(x %in% c(0, 1), x == 1, NA) }
never_preg <- yn(d$f002) %in% FALSE
premeno    <- yn(d$f022) %in% FALSE
d$menarche    <- rng(d$f001a, 8, 20)
d$ever_preg   <- yn(d$f002)
d$gravidity   <- ifelse(never_preg, 0, rng(d$f003, 0, 20))
d$parity      <- ifelse(never_preg, 0, rng(d$f004, 0, 15))
d$parity3     <- d$parity >= 3
d$miscarriage <- ifelse(never_preg, FALSE, yn(d$f005))
d$breastfeed  <- ifelse(never_preg | d$parity %in% 0, 0, rng(d$f006, 0, 240))
d$oc_ever     <- yn(d$f009)
d$oc_years    <- ifelse(d$oc_ever %in% FALSE, 0, rng(d$f010, 0, 50))
d$postmeno    <- yn(d$f022)
d$meno_age    <- rng(d$f013_1, 30, 65)
d$early_meno  <- d$meno_age < 45
d$repro_span  <- d$meno_age - d$menarche
d$yrs_since_meno <- ifelse(num(d$age) >= d$meno_age, num(d$age) - d$meno_age, NA)
d$hysterect   <- yn(d$f014)
d$hrt_ever    <- yn(d$f016)
d$hrt_years   <- ifelse(d$hrt_ever %in% FALSE, 0, rng(d$f017, 0, 60))
d$hrt_current <- ifelse(d$hrt_ever %in% FALSE, FALSE, yn(d$f019))

# variable, label, raw source, reason-for-missing variable, not-applicable mask
specs <- list(
  list("menarche", "Age at menarche, y", "f001a"),
  list("ever_preg", "Ever pregnant", "f002"),
  list("gravidity", "Pregnancies, n", "f003"),
  list("parity", "Live births, n", "f004"),
  list("parity3", "Live births >= 3", "f004"),
  list("miscarriage", "Any miscarriage", "f005"),
  list("breastfeed", "Breastfeeding, months", "f006"),
  list("oc_ever", "Ever oral contraceptive use", "f009"),
  list("oc_years", "Oral contraceptive use, years", "f010"),
  list("postmeno", "Postmenopausal", "f022"),
  list("meno_age", "Age at menopause, y", "f013_1", "f013_2", premeno),
  list("early_meno", "Early menopause (<45 y)", "f013_1", "f013_2", premeno),
  list("repro_span", "Reproductive lifespan, y", "f013_1", "f013_2", premeno),
  list("yrs_since_meno", "Years since menopause", "f013_1", "f013_2", premeno),
  list("hysterect", "Hysterectomy", "f014"),
  list("hrt_ever", "Ever hormone replacement therapy", "f016"),
  list("hrt_years", "Hormone replacement therapy, years", "f017"),
  list("hrt_current", "Current hormone replacement therapy", "f019"))

miss_tab <- function(x, raw, reason = NULL, na_app = NULL) {
  r  <- num(d[[raw]]); rr <- if (!is.null(reason) && reason %in% names(d)) num(d[[reason]]) else rep(NA, nrow(d))
  code <- ifelse(r %in% c(88, 8888) | rr %in% 88, "no_answer", ifelse(r %in% c(99, 9999) | rr %in% 99, "dont_know", "blank"))
  if (!is.null(na_app)) code[na_app] <- "not_applicable"
  m <- is.na(x)
  c(available = sum(!m), not_applicable = sum(m & code == "not_applicable"), no_answer_88 = sum(m & code == "no_answer"),
    dont_know_99 = sum(m & code == "dont_know"), blank_or_implausible = sum(m & code == "blank"))
}

st6 <- do.call(rbind, lapply(specs, function(sp) {
  v <- sp[[1]]; x <- d[[v]]; bin <- is.logical(x)
  smry <- function(z) if (bin) sprintf("%.1f%% (%d)", 100 * mean(z, na.rm = TRUE), sum(z, na.rm = TRUE)) else med_iqr(z)
  p <- tryCatch(if (bin) fisher.test(table(x, d$cluster))$p.value else wilcox.test(x ~ d$cluster)$p.value, error = function(e) NA)
  or <- tryCatch({ dd <- data.frame(y = d$cluster == "High", x = as.numeric(x), age = num(d$age))
  m <- glm(y ~ x + age, data = dd, family = binomial); b <- coef(m)["x"]; se <- sqrt(vcov(m)["x", "x"])
  sprintf("%.2f (%.2f-%.2f)", exp(b), exp(b - 1.96 * se), exp(b + 1.96 * se)) }, error = function(e) "n/a")
  mt <- miss_tab(x, sp[[3]], if (length(sp) >= 4) sp[[4]], if (length(sp) >= 5) sp[[5]])
  data.frame(factor = sp[[2]], t(mt), high = smry(x[d$cluster == "High"]), low = smry(x[d$cluster == "Low"]),
             p = signif(p, 2), OR_high_age_adj = or, check.names = FALSE)
}))
print(st6, row.names = FALSE)
write.csv(st6, file.path(project_dir, "SuppTable5_reproductive_hormonal_factors.csv"), row.names = FALSE)

# HRT preparations (free text) - descriptive overview
cat("\nHRT preparations ever (f018), top 15:\n");  print(head(sort(table(trimws(d$f018[d$f018 != ""])), decreasing = TRUE), 15))
cat("\nHRT preparations current (f021), top 15:\n"); print(head(sort(table(trimws(d$f021[d$f021 != ""])), decreasing = TRUE), 15))

# Figure 7: protein concordance (final; uses p_conc from earlier, re-created here if needed)
library(ggrepel)
conc <- data.frame(name = int_tab$name, logOR_w = log(int_tab$OR_women), logOR_m = log(int_tab$OR_men), p_int = int_tab$p_int_cluster)
conc$label <- ifelse(conc$p_int < 0.05, gsub("TGF-alpha", "TGF-α", trimws(sub(".*\\((.*)\\).*", "\\1", conc$name))), NA)
fig7 <- ggplot(conc, aes(logOR_m, logOR_w)) +
  geom_abline(linetype = 2, colour = "grey50") + geom_hline(yintercept = 0, colour = "grey85") + geom_vline(xintercept = 0, colour = "grey85") +
  geom_point(aes(colour = p_int < 0.05), size = 2) +
  geom_text_repel(aes(label = label), size = 3.2, na.rm = TRUE, min.segment.length = 0, box.padding = 0.4) +
  annotate("text", x = min(conc$logOR_m), y = max(conc$logOR_w), hjust = 0, vjust = 1, size = 3.5,
           label = sprintf("r = %.2f\nno protein x sex interaction at FDR < 0.05", cor(conc$logOR_w, conc$logOR_m))) +
  scale_colour_manual(values = c(`FALSE` = "grey60", `TRUE` = "#C0392B"), labels = c("concordant", "nominal sex interaction (p < 0.05)"), name = NULL) +
  labs(x = "Men: log OR for high-risk cluster (per SD)", y = "Women: log OR for high-risk cluster (per SD)") +
  theme_bw(base_size = 12) + theme(legend.position = "bottom")
ggsave(file.path(project_dir, "Figure7_protein_concordance.png"),  fig7, width = 6.5, height = 6.5, dpi = 300)
ggsave(file.path(project_dir, "Figure7_protein_concordance.tiff"), fig7, width = 6.5, height = 6.5, dpi = 300, compression = "lzw")
cat("Saved SuppTable4, SuppTable6, Figure7\n")


