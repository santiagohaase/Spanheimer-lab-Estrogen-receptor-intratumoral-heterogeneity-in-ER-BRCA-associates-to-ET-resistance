project_dir <- "/path/to/project_root"
# setwd(project_dir)  # Optional: uncomment for interactive use.
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
  library(readr)
  library(tibble)
  library(xgboost)
  library(survival)
  library(ggplot2)
  library(survminer)
})

# ============================
# USER INPUTS
# ============================
rfs_path      <- "/path/to/SC_RNA_ERHeterogeneity_QC/metabric_rfs_table.csv"
expr_path     <- "/path/to/SC_RNA_ERHeterogeneity_QC/data_mrna_illumina_microarray_zscores_ref_diploid_samples.txt"
clinical_path <- "/path/to/SC_RNA_ERHeterogeneity_QC/data_clinical_patient.txt"

out_dir <- "/path/to/SC_RNA_ERHeterogeneity_QC/xgb_survival_gene_risk_ERpos"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_qc                  <- file.path(out_dir, "QC.txt")
out_perf                <- file.path(out_dir, "performance.txt")
out_importance          <- file.path(out_dir, "gene_importance_robust.csv")
out_signed_rank         <- file.path(out_dir, "gene_signed_ranking.csv")
out_pred_test           <- file.path(out_dir, "test_predictions.csv")
out_model               <- file.path(out_dir, "xgb_aft_model.json")
out_signature_genes     <- file.path(out_dir, "derived_signature_genes.csv")
out_signature_scores    <- file.path(out_dir, "signature_scores_all_patients.csv")
out_signature_summary   <- file.path(out_dir, "signature_survival_summary.txt")
out_km_signature_pdf    <- file.path(out_dir, "KM_signature_top10_vs_bottom10.pdf")
out_km_signature_png    <- file.path(out_dir, "KM_signature_top10_vs_bottom10.png")
out_km_modelrisk_pdf    <- file.path(out_dir, "KM_modelrisk_top10_vs_bottom10.pdf")
out_km_modelrisk_png    <- file.path(out_dir, "KM_modelrisk_top10_vs_bottom10.png")

# ============================
# SETTINGS
# ============================
set.seed(1)

FILTER_TO_ER_POS <- TRUE

USE_ALL_GENES <- FALSE
TOP_VAR_GENES <- 8000

TEST_FRAC <- 0.20
N_BOOT    <- 100
BOOT_FRAC <- 0.85

NTHREAD <- 30

TOP_SIG_GENES <- 25
MIN_SELECTION_FREQ <- 0.10
EXTREME_FRAC <- 0.10

xgb_params <- list(
  booster = "gbtree",
  objective = "survival:aft",
  eval_metric = "aft-nloglik",
  aft_loss_distribution = "normal",
  aft_loss_distribution_scale = 1.0,
  eta = 0.05,
  max_depth = 3,
  min_child_weight = 10,
  subsample = 0.8,
  colsample_bytree = 0.6,
  lambda = 1.0,
  alpha  = 0.0,
  nthread = NTHREAD
)

NROUND_MAX <- 4000
EARLY_STOP <- 50

# ============================
# HELPERS
# ============================
`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}

norm_id <- function(x) {
  x %>%
    as.character() %>%
    str_trim() %>%
    str_replace_all("\\s+", "") %>%
    toupper()
}

log_line <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "|", paste(..., collapse = " "), "\n")
}

c_index <- function(time, event, risk_score) {
  cc <- survival::concordance(Surv(time, event) ~ risk_score)
  as.numeric(cc$concordance)
}

parse_er_status <- function(x) {
  er_norm <- x %>%
    as.character() %>%
    str_to_lower() %>%
    str_trim() %>%
    str_replace(regex("^\\s*[01]\\s*[:=-]\\s*"), "") %>%
    str_replace_all("[^a-z0-9\\+\\-]+", " ") %>%
    str_squish() %>%
    str_replace_all("\\bpositve\\b", "positive") %>%
    str_replace_all("\\bpositiv\\b", "positive") %>%
    str_replace_all("\\bnegtive\\b", "negative") %>%
    str_replace_all("\\bnegativ\\b", "negative")
  
  case_when(
    is.na(er_norm) | er_norm == "" ~ NA_character_,
    str_detect(er_norm, "\\bpositive\\b|\\bpos\\b|\\bplus\\b|\\ber\\s*\\+\\b|^1$") ~ "Positive",
    str_detect(er_norm, "\\bnegative\\b|\\bneg\\b|\\bminus\\b|\\ber\\s*\\-\\b|^0$") ~ "Negative",
    TRUE ~ "Unknown"
  )
}

pick_er_col <- function(df, id_col = "#Patient Identifier") {
  nm <- names(df)
  nm_clean <- nm %>%
    str_to_lower() %>%
    str_replace_all("[^a-z0-9]+", "_") %>%
    str_replace_all("_+", "_") %>%
    str_replace_all("^_|_$", "")
  
  keep <- nm != id_col
  nm <- nm[keep]
  nm_clean <- nm_clean[keep]
  
  preferred_exact <- c(
    "ER status measured by IHC",
    "ER Status",
    "ER_STATUS",
    "ER IHC Status",
    "Estrogen Receptor Status",
    "Estrogen Receptor"
  )
  
  hit_exact <- intersect(preferred_exact, nm)
  if (length(hit_exact) > 0) return(hit_exact[[1]])
  
  score <- rep(0L, length(nm))
  score <- score + 100L * str_detect(nm_clean, regex("er.*ihc", ignore_case = TRUE))
  score <- score +  90L * str_detect(nm_clean, regex("estrogen.*receptor", ignore_case = TRUE))
  score <- score +  80L * str_detect(nm_clean, regex("er(_|$).*status|status.*er(_|$)", ignore_case = TRUE))
  score <- score +  60L * str_detect(nm_clean, regex("er.*receptor|receptor.*er", ignore_case = TRUE))
  score <- score - 200L * str_detect(nm_clean, regex("identif|patient|case|sample|id$", ignore_case = TRUE))
  
  best <- which.max(score)
  if (length(best) == 0 || score[[best]] <= 0) {
    cand <- nm[str_detect(nm_clean, regex("er|estrogen", ignore_case = TRUE))]
    stop(
      "Could not confidently detect ER column.\nCandidates:\n  ",
      paste(cand, collapse = "\n  ")
    )
  }
  nm[[best]]
}

score_signature_from_matrix <- function(expr_mat, genes) {
  genes_use <- intersect(genes, rownames(expr_mat))
  if (length(genes_use) == 0) {
    return(rep(NA_real_, ncol(expr_mat)))
  }
  colMeans(expr_mat[genes_use, , drop = FALSE], na.rm = TRUE)
}

make_extreme_group <- function(score, frac = 0.10, low_label = "Low10", high_label = "High10") {
  q_low  <- quantile(score, probs = frac, na.rm = TRUE)
  q_high <- quantile(score, probs = 1 - frac, na.rm = TRUE)
  
  case_when(
    score <= q_low  ~ low_label,
    score >= q_high ~ high_label,
    TRUE ~ NA_character_
  )
}

plot_extreme_km <- function(df, score_col, time_col, event_col,
                            frac = 0.10,
                            legend_title = "Group",
                            low_label = "Low10",
                            high_label = "High10",
                            pdf_file, png_file, summary_file = NULL) {
  
  group_vec <- make_extreme_group(df[[score_col]], frac = frac,
                                  low_label = low_label, high_label = high_label)
  
  km_df <- df %>%
    mutate(.group = group_vec) %>%
    filter(!is.na(.group)) %>%
    mutate(.group = factor(.group, levels = c(low_label, high_label)))
  
  if (nrow(km_df) < 20 || length(unique(km_df$.group)) < 2) {
    log_line("Not enough patients to plot KM for", score_col)
    return(invisible(NULL))
  }
  
  surv_formula <- as.formula(paste0("Surv(", time_col, ", ", event_col, ") ~ .group"))
  fit_km <- survfit(surv_formula, data = km_df)
  lr <- survdiff(surv_formula, data = km_df)
  lr_p <- 1 - pchisq(lr$chisq, df = length(lr$n) - 1)
  
  cox_fit <- coxph(surv_formula, data = km_df)
  cox_sum <- summary(cox_fit)
  
  p <- ggsurvplot(
    fit_km,
    data = km_df,
    risk.table = TRUE,
    pval = TRUE,
    conf.int = TRUE,
    legend.title = legend_title,
    legend.labs = c(low_label, high_label),
    xlab = "Time",
    ylab = "Recurrence-free survival probability",
    ggtheme = theme_classic()
  )
  
  ggsave(pdf_file, p$plot, width = 7, height = 6)
  ggsave(png_file, p$plot, width = 7, height = 6, dpi = 300)
  
  if (!is.null(summary_file)) {
    writeLines(c(
      paste0("Score column: ", score_col),
      paste0("Extreme fraction: ", frac),
      paste0("n_extreme: ", nrow(km_df)),
      paste0(low_label, " n: ", sum(km_df$.group == low_label)),
      paste0(high_label, " n: ", sum(km_df$.group == high_label)),
      paste0("Log-rank chisq: ", signif(lr$chisq, 4)),
      paste0("Log-rank p: ", signif(lr_p, 4)),
      paste0("Cox HR ", high_label, " vs ", low_label, ": ", signif(cox_sum$conf.int[1, "exp(coef)"], 4)),
      paste0("Cox HR lower95: ", signif(cox_sum$conf.int[1, "lower .95"], 4)),
      paste0("Cox HR upper95: ", signif(cox_sum$conf.int[1, "upper .95"], 4)),
      paste0("Cox p: ", signif(cox_sum$coefficients[1, "Pr(>|z|)"], 4))
    ), summary_file)
  }
  
  invisible(list(
    km_df = km_df,
    fit_km = fit_km,
    cox_fit = cox_fit,
    logrank = lr
  ))
}

# ============================
# 1) READ RFS
# ============================
log_line("Reading RFS:", rfs_path)
rfs <- read_csv(rfs_path, show_col_types = FALSE) %>%
  mutate(patient_id = norm_id(patient_id)) %>%
  filter(
    str_detect(patient_id, "^(MB-|MTS-)"),
    is.finite(rfs_time), rfs_time > 0,
    rfs_event %in% c(0, 1)
  ) %>%
  distinct(patient_id, .keep_all = TRUE)

log_line("RFS kept:", nrow(rfs), "| events:", sum(rfs$rfs_event))

# ============================
# 2) READ CLINICAL and FILTER ER+
# ============================
log_line("Reading clinical:", clinical_path)
clin0 <- read_delim(clinical_path, delim = "\t", show_col_types = FALSE, guess_max = 50000)

id_col <- "#Patient Identifier"
if (!(id_col %in% names(clin0))) {
  stop("Clinical file missing column: ", id_col)
}

er_col <- pick_er_col(clin0, id_col = id_col)
log_line("Using ER column:", er_col)

nm0 <- names(clin0)
nm0_clean <- nm0 %>% str_to_lower() %>% str_replace_all("[^a-z0-9]+", "_")
er_like <- nm0[str_detect(nm0_clean, regex("er|estrogen", ignore_case = TRUE))]
log_line("Clinical ER-like columns:", paste(er_like, collapse = " | "))

clin <- clin0 %>%
  mutate(
    pid = as.character(.data[[id_col]]),
    er_raw = as.character(.data[[er_col]])
  ) %>%
  filter(
    !is.na(pid), pid != "",
    str_detect(pid, "^(MB-|MTS-)")
  ) %>%
  transmute(
    patient_id = norm_id(pid),
    er_status  = parse_er_status(er_raw)
  ) %>%
  distinct(patient_id, .keep_all = TRUE)

log_line("Clinical patients:", nrow(clin))
log_line("ER status counts:", paste(capture.output(print(table(clin$er_status, useNA = "ifany"))), collapse = " "))

if (FILTER_TO_ER_POS) {
  er_keep <- clin %>% filter(er_status == "Positive") %>% pull(patient_id)
  log_line("ER+ filter enabled. ER+ patients:", length(er_keep))
} else {
  er_keep <- clin %>% pull(patient_id)
  log_line("ER filter disabled. Using all clinical patients:", length(er_keep))
}

rfs <- rfs %>% filter(patient_id %in% er_keep)
log_line("RFS after ER filtering:", nrow(rfs), "| events:", sum(rfs$rfs_event))

if (nrow(rfs) < 200) {
  stop("Too few ER-filtered RFS patients. Check ER parsing / file.")
}

# ============================
# 3) READ EXPRESSION MATRIX
# ============================
log_line("Reading expression matrix:", expr_path)
expr_dt <- fread(expr_path, sep = "\t", header = TRUE, data.table = FALSE)

stopifnot(all(c("Hugo_Symbol", "Entrez_Gene_Id") %in% colnames(expr_dt)))

sample_cols <- setdiff(colnames(expr_dt), c("Hugo_Symbol", "Entrez_Gene_Id"))
sample_ids  <- norm_id(sample_cols)

keep_samples <- str_detect(sample_ids, "^(MB-|MTS-)")
sample_cols  <- sample_cols[keep_samples]
sample_ids   <- sample_ids[keep_samples]

log_line("Expression samples:", length(sample_ids), "| genes:", nrow(expr_dt))

common_ids <- intersect(rfs$patient_id, sample_ids)
log_line("Intersect patients (ER-filtered RFS ∩ expr):", length(common_ids))
if (length(common_ids) < 200) {
  stop("Too few intersected patients: ", length(common_ids))
}

rfs2 <- rfs %>% filter(patient_id %in% common_ids)

id_to_col <- setNames(sample_cols, sample_ids)
expr_cols_use <- id_to_col[rfs2$patient_id]
stopifnot(!any(is.na(expr_cols_use)))

X_gp <- as.matrix(expr_dt[, expr_cols_use, drop = FALSE])
mode(X_gp) <- "numeric"

gene_names <- expr_dt$Hugo_Symbol %>% as.character()
gene_names[is.na(gene_names) | gene_names == ""] <- paste0("ENTREZ_", expr_dt$Entrez_Gene_Id)
rownames(X_gp) <- make.unique(gene_names)

# ============================
# 4) GENE SELECTION
# ============================
if (!USE_ALL_GENES) {
  log_line("Filtering genes by variance; keeping top:", TOP_VAR_GENES)
  vars <- apply(X_gp, 1, var, na.rm = TRUE)
  vars[!is.finite(vars)] <- 0
  keep <- order(vars, decreasing = TRUE)[seq_len(min(TOP_VAR_GENES, length(vars)))]
  X_gp <- X_gp[keep, , drop = FALSE]
} else {
  log_line("USE_ALL_GENES=TRUE: using all genes:", nrow(X_gp))
}

genes_used <- rownames(X_gp)

X <- t(X_gp)
X[!is.finite(X)] <- 0
X <- matrix(as.single(X), nrow = nrow(X), ncol = ncol(X), dimnames = dimnames(X))

# ============================
# 5) AFT LABELS
# ============================
time  <- as.numeric(rfs2$rfs_time)
event <- as.numeric(rfs2$rfs_event)

label_lower <- time
label_upper <- ifelse(event == 1, time, Inf)

# ============================
# 6) TRAIN/TEST SPLIT
# ============================
n <- nrow(X)
test_n <- max(1, floor(TEST_FRAC * n))

test_idx <- sample.int(n, size = test_n)
train_idx <- setdiff(seq_len(n), test_idx)

X_train <- X[train_idx, , drop = FALSE]
X_test  <- X[test_idx,  , drop = FALSE]

time_train  <- time[train_idx]
event_train <- event[train_idx]
time_test   <- time[test_idx]
event_test  <- event[test_idx]

lb_train <- label_lower[train_idx]
ub_train <- label_upper[train_idx]
lb_test  <- label_lower[test_idx]
ub_test  <- label_upper[test_idx]

dtrain <- xgb.DMatrix(data = X_train, feature_names = colnames(X_train))
dtest  <- xgb.DMatrix(data = X_test,  feature_names = colnames(X_test))

setinfo(dtrain, "label_lower_bound", lb_train)
setinfo(dtrain, "label_upper_bound", ub_train)
setinfo(dtest,  "label_lower_bound", lb_test)
setinfo(dtest,  "label_upper_bound", ub_test)

watchlist <- list(train = dtrain, test = dtest)

# ============================
# 7) FIT MAIN MODEL
# ============================
log_line("Training XGBoost AFT model...")
fit <- xgb.train(
  params = xgb_params,
  data = dtrain,
  nrounds = NROUND_MAX,
  watchlist = watchlist,
  early_stopping_rounds = EARLY_STOP,
  verbose = 1
)

xgb.save(fit, out_model)
log_line("Saved model:", out_model)

pred_test <- predict(fit, dtest)
risk_test <- -pred_test

test_tbl <- tibble(
  patient_id = rfs2$patient_id[test_idx],
  rfs_time   = time_test,
  rfs_event  = event_test,
  pred_time  = pred_test,
  risk_score = risk_test
)

write_csv(test_tbl, out_pred_test)
log_line("Wrote test predictions:", out_pred_test)

cidx_test <- c_index(time_test, event_test, risk_test)
perf_lines <- c(
  paste0("best_iteration: ", fit$best_iteration),
  paste0("n_train: ", length(train_idx)),
  paste0("n_test: ", length(test_idx)),
  paste0("test_cindex: ", signif(cidx_test, 4))
)
writeLines(perf_lines, out_perf)
log_line("Wrote performance:", out_perf)
log_line("Test C-index:", signif(cidx_test, 4))

# ============================
# 8) BOOTSTRAP ROBUST IMPORTANCE
# ============================
log_line("Running bootstrap importance with N_BOOT =", N_BOOT)

boot_list <- vector("list", N_BOOT)

for (b in seq_len(N_BOOT)) {
  if (b %% 5 == 0 || b == 1 || b == N_BOOT) {
    log_line("Bootstrap", b, "of", N_BOOT)
  }
  
  boot_n <- max(20, floor(BOOT_FRAC * length(train_idx)))
  boot_rows <- sample(seq_len(nrow(X_train)), size = boot_n, replace = TRUE)
  
  Xb  <- X_train[boot_rows, , drop = FALSE]
  lbb <- lb_train[boot_rows]
  ubb <- ub_train[boot_rows]
  
  db <- xgb.DMatrix(data = Xb, feature_names = colnames(Xb))
  setinfo(db, "label_lower_bound", lbb)
  setinfo(db, "label_upper_bound", ubb)
  
  fit_b <- tryCatch(
    xgb.train(
      params = xgb_params,
      data = db,
      nrounds = max(100, fit$best_iteration %||% 300),
      verbose = 0
    ),
    error = function(e) NULL
  )
  
  if (is.null(fit_b)) next
  
  imp_b <- tryCatch(
    xgb.importance(model = fit_b),
    error = function(e) NULL
  )
  
  if (is.null(imp_b) || nrow(imp_b) == 0) next
  
  imp_b <- as_tibble(imp_b) %>%
    transmute(
      gene = Feature,
      gain = Gain,
      cover = Cover,
      frequency = Frequency,
      boot = b
    )
  
  boot_list[[b]] <- imp_b
}

boot_imp <- bind_rows(boot_list)

if (nrow(boot_imp) == 0) {
  stop("Bootstrap importance failed: no importance tables collected.")
}

imp_robust <- boot_imp %>%
  group_by(gene) %>%
  summarise(
    n_boot_present = n_distinct(boot),
    selection_freq = n_boot_present / N_BOOT,
    mean_gain = mean(gain, na.rm = TRUE),
    median_gain = median(gain, na.rm = TRUE),
    mean_cover = mean(cover, na.rm = TRUE),
    mean_frequency = mean(frequency, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_gain), desc(selection_freq))

write_csv(imp_robust, out_importance)
log_line("Wrote robust importance:", out_importance)

# ============================
# 9) SIGNED DIRECTION (Spearman gene vs predicted risk)
# ============================
log_line("Computing signed direction (Spearman gene vs predicted risk)...")

dfull <- xgb.DMatrix(data = X, feature_names = colnames(X))
setinfo(dfull, "label_lower_bound", label_lower)
setinfo(dfull, "label_upper_bound", label_upper)

pred_full <- predict(fit, dfull)
risk_full <- -pred_full

cors <- numeric(ncol(X))
names(cors) <- colnames(X)

pb <- txtProgressBar(min = 0, max = ncol(X), style = 3)
for (j in seq_len(ncol(X))) {
  cors[j] <- suppressWarnings(
    cor(X[, j], risk_full, method = "spearman", use = "pairwise.complete.obs")
  )
  if (j %% 50 == 0) {
    setTxtProgressBar(pb, j)
  }
}
setTxtProgressBar(pb, ncol(X))
close(pb)

signed_tbl <- imp_robust %>%
  mutate(
    spearman_risk = unname(cors[match(gene, names(cors))]),
    signed_score  = mean_gain * sign(spearman_risk)
  ) %>%
  arrange(desc(signed_score))

write_csv(signed_tbl, out_signed_rank)
log_line("Wrote signed ranking:", out_signed_rank)

# ============================
# 10) BUILD DERIVED SIGNATURE
# ============================
log_line("Building derived signature...")

sig_tbl_filt <- signed_tbl %>%
  filter(
    is.finite(signed_score),
    is.finite(selection_freq),
    selection_freq >= MIN_SELECTION_FREQ
  )

risk_high_genes <- sig_tbl_filt %>%
  arrange(desc(signed_score), desc(selection_freq), desc(mean_gain)) %>%
  slice_head(n = TOP_SIG_GENES) %>%
  pull(gene) %>%
  unique()

risk_low_genes <- sig_tbl_filt %>%
  arrange(signed_score, desc(selection_freq), desc(mean_gain)) %>%
  slice_head(n = TOP_SIG_GENES) %>%
  pull(gene) %>%
  unique()

if (length(risk_high_genes) < 5 || length(risk_low_genes) < 5) {
  stop("Too few genes passed filters for derived signature.")
}

sig_gene_tbl <- bind_rows(
  tibble(signature = "risk_high", gene = risk_high_genes),
  tibble(signature = "risk_low", gene = risk_low_genes)
)

write_csv(sig_gene_tbl, out_signature_genes)
log_line("Wrote signature genes:", out_signature_genes)

# ============================
# 11) SCORE ALL PATIENTS
# ============================
log_line("Scoring all patients...")

sig_high_score <- score_signature_from_matrix(X_gp, risk_high_genes)
sig_low_score  <- score_signature_from_matrix(X_gp, risk_low_genes)
sig_score      <- sig_high_score - sig_low_score

sig_scores_tbl <- tibble(
  patient_id = rfs2$patient_id,
  rfs_time = time,
  rfs_event = event,
  model_pred_time = pred_full,
  model_risk = risk_full,
  sig_high_score = sig_high_score,
  sig_low_score = sig_low_score,
  sig_score = sig_score
)

write_csv(sig_scores_tbl, out_signature_scores)
log_line("Wrote signature scores:", out_signature_scores)

# ============================
# 12) KM PLOT FOR SIGNATURE SCORE
# ============================
log_line("Plotting KM for derived signature extremes...")
plot_extreme_km(
  df = sig_scores_tbl,
  score_col = "sig_score",
  time_col = "rfs_time",
  event_col = "rfs_event",
  frac = EXTREME_FRAC,
  legend_title = "Derived risk signature",
  low_label = "Bottom 10%",
  high_label = "Top 10%",
  pdf_file = out_km_signature_pdf,
  png_file = out_km_signature_png,
  summary_file = out_signature_summary
)

# ============================
# 13) KM PLOT FOR MODEL RISK
# ============================
log_line("Plotting KM for model risk extremes...")
plot_extreme_km(
  df = sig_scores_tbl,
  score_col = "model_risk",
  time_col = "rfs_time",
  event_col = "rfs_event",
  frac = EXTREME_FRAC,
  legend_title = "Model risk",
  low_label = "Bottom 10%",
  high_label = "Top 10%",
  pdf_file = out_km_modelrisk_pdf,
  png_file = out_km_modelrisk_png,
  summary_file = NULL
)

# ============================
# 14) QC SUMMARY
# ============================
qc_lines <- c(
  paste0("XGB survival AFT | ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "INPUTS:",
  paste0("  rfs_path: ", rfs_path),
  paste0("  expr_path: ", expr_path),
  paste0("  clinical_path: ", clinical_path),
  "",
  "ER FILTER:",
  paste0("  FILTER_TO_ER_POS: ", FILTER_TO_ER_POS),
  paste0("  ER column used: ", er_col),
  paste0("  ER status counts: ", paste(capture.output(print(table(clin$er_status, useNA = 'ifany'))), collapse = " ")),
  "",
  "COUNTS:",
  paste0("  RFS kept (post-ER-filter): ", nrow(rfs2)),
  paste0("  Events: ", sum(event)),
  paste0("  Genes used: ", ncol(X)),
  "",
  "MODEL:",
  paste0("  best_iteration: ", fit$best_iteration),
  paste0("  test_cindex: ", signif(cidx_test, 4)),
  paste0("  N_BOOT: ", N_BOOT),
  paste0("  TOP_SIG_GENES: ", TOP_SIG_GENES),
  paste0("  MIN_SELECTION_FREQ: ", MIN_SELECTION_FREQ),
  paste0("  EXTREME_FRAC: ", EXTREME_FRAC),
  "",
  "DERIVED SIGNATURE:",
  paste0("  risk_high_genes_n: ", length(risk_high_genes)),
  paste0("  risk_low_genes_n: ", length(risk_low_genes))
)

writeLines(qc_lines, out_qc)
log_line("Wrote QC:", out_qc)
log_line("Done.")