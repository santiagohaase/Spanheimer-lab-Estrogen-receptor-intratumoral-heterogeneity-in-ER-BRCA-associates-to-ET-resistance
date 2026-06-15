


#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(survival)
  library(survminer)
  library(ggplot2)
  library(cmprsk)
})

# =========================================================
# MASTER SETTINGS
# =========================================================

# Choose panel: "nonprolif" or "prolif"
PANEL_NAME <- "nonprolif"

# Main analysis toggles
RUN_MAIN_RFS_ERPOS                <- TRUE
RUN_MAIN_RFS_ERPOS_EXTREMES       <- TRUE
RUN_COMPETING_RISKS_ERPOS         <- TRUE
RUN_SEPARATE_UP_DOWN_ERPOS        <- TRUE

RUN_ET_120MO                      <- TRUE
RUN_ET_120MO_GRID_SEARCH          <- TRUE
RUN_ET_120MO_REPLOT_BEST          <- TRUE
RUN_ET_120MO_WITH_ESR1_COVARIATE  <- TRUE

# General parameters
EXTREME_FRAC_MAIN <- 0.10   # top 10% vs bottom 10% for main ER+ script
EXTREME_FRAC_CR   <- 0.20   # top 20% vs bottom 20% for competing risks to match original intent

ET_TMAX_MONTHS <- 120
ET_CUT_TOP     <- 0.05
ET_CUT_BOTTOM  <- 0.50
CUT_GRID       <- c(0.05, 0.10, 0.20, 0.30, 0.40, 0.50)

COL_LOW  <- "#4C78A8"
COL_HIGH <- "#E45756"

# =========================================================
# PATHS
# =========================================================

base_dir <- "/path/to/SC_RNA_ERHeterogeneity_QC"

panel_dir <- switch(
  PANEL_NAME,
  nonprolif = file.path(base_dir, "NonProlif_DE_signature_panel"),
  prolif    = file.path(base_dir, "Prolif_DE_signature_panel"),
  stop("Unknown PANEL_NAME: ", PANEL_NAME)
)

up_file   <- file.path(panel_dir, "CONSENSUS_UP_supportGE3.txt")
down_file <- file.path(panel_dir, "CONSENSUS_DOWN_supportGE3.txt")

expr_file <- file.path(base_dir, "data_mrna_illumina_microarray_zscores_ref_diploid_samples.txt")
map_file  <- file.path(base_dir, "METABRIC_ERpos_heterogeneity_by_sample.csv")
clin_patient_file <- file.path(base_dir, "data_clinical_patient.txt")

out_root <- file.path(base_dir, paste0("metabric_master_signature_", PANEL_NAME))
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(out_root, "master_analysis.log")

# =========================================================
# LOGGING
# =========================================================
log_msg <- function(...) {
  line <- paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "|", paste(..., collapse = " "))
  cat(line, "\n", file = log_file, append = TRUE)
  cat(line, "\n")
}

log_msg("=== MASTER METABRIC signature analysis started ===")
log_msg("PANEL_NAME:", PANEL_NAME)
log_msg("Output root:", out_root)

# =========================================================
# HELPERS
# =========================================================

read_gene_list <- function(path) {
  stopifnot(file.exists(path))
  x <- readLines(path, warn = FALSE)
  x <- sub("^\ufeff", "", x)
  x <- gsub("\r$", "", x)
  x <- str_trim(x)
  x <- x[x != ""]
  x <- x[!str_starts(x, "#")]
  genes <- str_split(x, "\\s+|\t|,", simplify = TRUE)[, 1]
  genes <- unique(genes[genes != "" & !is.na(genes)])
  toupper(genes)
}

read_with_header_line <- function(path, header_startswith = "PATIENT_ID", sep = "\t") {
  lines <- readLines(path, warn = FALSE)
  lines[1] <- sub("^\ufeff", "", lines[1])
  lines <- gsub("\r$", "", lines)
  
  hdr_i <- which(grepl(paste0("^\\s*", header_startswith, "\\b"), lines))[1]
  if (is.na(hdr_i)) {
    writeLines(lines[1:min(250, length(lines))],
               file.path(out_root, paste0("HEAD_PREVIEW_", basename(path), ".txt")))
    stop("Could not find header starting with ", header_startswith, " in ", path)
  }
  
  txt <- paste(gsub("^\\s*#\\s*", "", lines[hdr_i:length(lines)]), collapse = "\n")
  fread(text = txt, sep = sep, header = TRUE, data.table = FALSE, fill = TRUE)
}

as_mb_id <- function(x) {
  x <- as.character(x)
  x <- str_trim(x)
  
  out <- x
  is_num <- str_detect(out, "^\\d+$")
  out[is_num] <- sprintf("MB-%04d", as.integer(out[is_num]))
  
  is_mb_short <- str_detect(out, "^MB-\\d+$")
  out[is_mb_short] <- sprintf("MB-%04d", as.integer(str_remove(out[is_mb_short], "^MB-")))
  
  out
}

clean_yesno <- function(x) {
  x <- as.character(x)
  x <- str_to_lower(str_trim(x))
  case_when(
    x %in% c("yes", "y", "1", "true", "t") ~ "YES",
    x %in% c("no", "n", "0", "false", "f") ~ "NO",
    TRUE ~ NA_character_
  )
}

parse_event01 <- function(x) {
  x <- as.character(x)
  case_when(
    str_detect(x, "^\\s*1\\s*:") ~ 1L,
    str_detect(x, "^\\s*0\\s*:") ~ 0L,
    TRUE ~ NA_integer_
  )
}

admin_censor <- function(time, event, tmax) {
  list(
    time  = pmin(time, tmax),
    event = ifelse(!is.na(time) & time > tmax, 0L, event)
  )
}

score_signature_meanz <- function(mat, gene_set) {
  present <- intersect(gene_set, toupper(rownames(mat)))
  idx <- match(present, toupper(rownames(mat)))
  frac_present <- length(idx) / max(1, length(gene_set))
  
  if (length(idx) == 0) {
    return(list(
      score = rep(NA_real_, ncol(mat)),
      n_present = 0L,
      frac_present = frac_present,
      missing = gene_set
    ))
  }
  
  sub <- mat[idx, , drop = FALSE]
  list(
    score = colMeans(sub, na.rm = TRUE),
    n_present = length(idx),
    frac_present = frac_present,
    missing = setdiff(gene_set, present)
  )
}

make_extreme_groups_symmetric <- function(score_vec, frac) {
  q_low  <- as.numeric(quantile(score_vec, probs = frac, na.rm = TRUE))
  q_high <- as.numeric(quantile(score_vec, probs = 1 - frac, na.rm = TRUE))
  
  grp <- case_when(
    score_vec <= q_low  ~ paste0("Bottom", round(frac * 100)),
    score_vec >= q_high ~ paste0("Top", round(frac * 100)),
    TRUE ~ NA_character_
  )
  
  factor(grp, levels = c(paste0("Bottom", round(frac * 100)),
                         paste0("Top", round(frac * 100))))
}

make_extreme_groups_asymmetric <- function(score_vec, top_p, bottom_p) {
  q_top <- as.numeric(quantile(score_vec, probs = 1 - top_p, na.rm = TRUE))
  q_bot <- as.numeric(quantile(score_vec, probs = bottom_p,   na.rm = TRUE))
  
  grp <- rep(NA_character_, length(score_vec))
  grp[score_vec >= q_top] <- sprintf("Top%.0f", top_p * 100)
  grp[score_vec <= q_bot] <- sprintf("Bottom%.0f", bottom_p * 100)
  
  factor(grp, levels = c(sprintf("Bottom%.0f", bottom_p * 100),
                         sprintf("Top%.0f", top_p * 100)))
}

safe_cox <- function(formula, data) {
  tryCatch(coxph(formula, data = data), error = function(e) e)
}

extract_hr <- function(fit, term_pattern) {
  if (inherits(fit, "error")) {
    return(tibble(term = NA_character_, HR = NA_real_, lo95 = NA_real_, hi95 = NA_real_, p = NA_real_))
  }
  s <- summary(fit)
  rn <- rownames(s$coef)
  idx <- grep(term_pattern, rn)[1]
  if (is.na(idx)) {
    return(tibble(term = NA_character_, HR = NA_real_, lo95 = NA_real_, hi95 = NA_real_, p = NA_real_))
  }
  tibble(
    term = rn[idx],
    HR   = as.numeric(s$coef[idx, "exp(coef)"]),
    lo95 = as.numeric(s$conf.int[idx, "lower .95"]),
    hi95 = as.numeric(s$conf.int[idx, "upper .95"]),
    p    = as.numeric(s$coef[idx, "Pr(>|z|)"])
  )
}

# =========================================================
# LOAD COMMON INPUTS ONCE
# =========================================================

log_msg("Reading panel genes...")
up_genes <- read_gene_list(up_file)
down_genes <- read_gene_list(down_file)

log_msg("UP genes:", length(up_genes))
log_msg("DOWN genes:", length(down_genes))
if (length(up_genes) < 5 || length(down_genes) < 5) stop("Gene lists too small.")

log_msg("Reading expression matrix...")
expr_dt <- fread(expr_file, sep = "\t", header = TRUE, data.table = FALSE)
stopifnot(all(c("Hugo_Symbol", "Entrez_Gene_Id") %in% colnames(expr_dt)))

genes_expr <- toupper(trimws(as.character(expr_dt$Hugo_Symbol)))
expr_mat <- expr_dt %>%
  dplyr::select(-Hugo_Symbol, -Entrez_Gene_Id) %>%
  as.matrix()
mode(expr_mat) <- "numeric"
rownames(expr_mat) <- genes_expr

if (anyDuplicated(rownames(expr_mat)) > 0) {
  keep <- !duplicated(rownames(expr_mat))
  expr_mat <- expr_mat[keep, , drop = FALSE]
  log_msg("Removed duplicated gene symbols; kept first occurrence.")
}

log_msg("Expression matrix:", nrow(expr_mat), "genes x", ncol(expr_mat), "samples")

# Score signatures once
up_res   <- score_signature_meanz(expr_mat, up_genes)
down_res <- score_signature_meanz(expr_mat, down_genes)

writeLines(up_res$missing,   file.path(out_root, "missing_UP_genes.txt"))
writeLines(down_res$missing, file.path(out_root, "missing_DOWN_genes.txt"))

scores_by_sample <- tibble(
  sample_id = colnames(expr_mat),
  UP_score = as.numeric(up_res$score),
  DOWN_score = as.numeric(down_res$score),
  DIFF_UP_minus_DOWN = as.numeric(up_res$score) - as.numeric(down_res$score)
)

fwrite(scores_by_sample, file.path(out_root, "scores_by_sample.tsv"), sep = "\t")
log_msg("Wrote scores_by_sample.tsv")

# Load mapping once
log_msg("Reading mapping...")
map_dt <- fread(map_file, data.table = FALSE)
colnames(map_dt) <- toupper(colnames(map_dt))
stopifnot(all(c("SAMPLE_ID", "PATIENT_ID") %in% colnames(map_dt)))

map_clean <- map_dt %>%
  transmute(
    sample_id = as.character(SAMPLE_ID),
    patient_id = as_mb_id(PATIENT_ID),
    ER_STATUS = if ("ER_STATUS" %in% colnames(map_dt)) str_to_lower(str_trim(as.character(ER_STATUS))) else NA_character_
  )

map_erpos <- if (!all(is.na(map_clean$ER_STATUS))) {
  map_clean %>% filter(ER_STATUS %in% c("positive", "pos", "er+", "er positive"))
} else {
  map_clean
}

fwrite(map_erpos, file.path(out_root, "sample_to_patient_mapping_ERpos.tsv"), sep = "\t")
log_msg("ER+ mapping rows:", nrow(map_erpos))

scores_patient_erpos <- scores_by_sample %>%
  inner_join(map_erpos %>% dplyr::select(sample_id, patient_id), by = "sample_id") %>%
  group_by(patient_id) %>%
  summarise(
    UP_score = mean(UP_score, na.rm = TRUE),
    DOWN_score = mean(DOWN_score, na.rm = TRUE),
    DIFF_UP_minus_DOWN = mean(DIFF_UP_minus_DOWN, na.rm = TRUE),
    n_samples = n(),
    .groups = "drop"
  )

fwrite(scores_patient_erpos, file.path(out_root, "scores_patient_ERpos.tsv"), sep = "\t")
log_msg("Wrote scores_patient_ERpos.tsv")

# Load clinical once
log_msg("Reading patient clinical...")
clin_p <- read_with_header_line(clin_patient_file, "PATIENT_ID", "\t")
colnames(clin_p) <- toupper(colnames(clin_p))

# Build ER+ clinical base
clin_er <- clin_p %>%
  mutate(
    patient_id = as_mb_id(PATIENT_ID),
    ER_IHC_CLEAN = str_to_lower(str_trim(as.character(ER_IHC))),
    HORMONE_THERAPY_CLEAN = if ("HORMONE_THERAPY" %in% colnames(clin_p)) clean_yesno(HORMONE_THERAPY) else NA_character_,
    AGE_AT_DIAGNOSIS = if ("AGE_AT_DIAGNOSIS" %in% colnames(clin_p)) suppressWarnings(as.numeric(AGE_AT_DIAGNOSIS)) else NA_real_,
    rfs_time = if ("RFS_MONTHS" %in% colnames(clin_p)) suppressWarnings(as.numeric(RFS_MONTHS)) else NA_real_,
    rfs_event = if ("RFS_STATUS" %in% colnames(clin_p)) parse_event01(RFS_STATUS) else NA_integer_,
    os_time = if ("OS_MONTHS" %in% colnames(clin_p)) suppressWarnings(as.numeric(OS_MONTHS)) else NA_real_,
    os_event = if ("OS_STATUS" %in% colnames(clin_p)) suppressWarnings(as.integer(OS_STATUS)) else NA_integer_
  ) %>%
  filter(ER_IHC_CLEAN %in% c("positive", "positve", "pos", "er+", "er positive"))

fwrite(clin_er, file.path(out_root, "clinical_patient_ERpos.tsv"), sep = "\t")
log_msg("ER+ clinical patients:", nrow(clin_er))

# Optional ESR1 patient expression table from matrix
esr1_row_idx <- which(toupper(rownames(expr_mat)) == "ESR1")[1]
esr1_patient <- NULL
if (!is.na(esr1_row_idx)) {
  esr1_sample <- tibble(
    sample_id = colnames(expr_mat),
    ESR1_expr = as.numeric(expr_mat[esr1_row_idx, ])
  )
  
  esr1_patient <- esr1_sample %>%
    inner_join(map_erpos %>% dplyr::select(sample_id, patient_id), by = "sample_id") %>%
    group_by(patient_id) %>%
    summarise(
      ESR1_expr_mean = mean(ESR1_expr, na.rm = TRUE),
      ESR1_expr_median = median(ESR1_expr, na.rm = TRUE),
      ESR1_expr_n_samples = sum(!is.na(ESR1_expr)),
      .groups = "drop"
    )
  
  fwrite(esr1_patient, file.path(out_root, "ESR1_patient_expression.tsv"), sep = "\t")
  log_msg("Wrote ESR1_patient_expression.tsv")
}

# =========================================================
# A) MAIN ER+ RFS ANALYSIS
# =========================================================
if (RUN_MAIN_RFS_ERPOS) {
  subdir <- file.path(out_root, "01_main_RFS_ERpos")
  dir.create(subdir, recursive = TRUE, showWarnings = FALSE)
  
  df_rfs <- scores_patient_erpos %>%
    inner_join(
      clin_er %>% dplyr::select(patient_id, rfs_time, rfs_event, AGE_AT_DIAGNOSIS, HORMONE_THERAPY_CLEAN),
      by = "patient_id"
    ) %>%
    filter(!is.na(rfs_time), !is.na(rfs_event), !is.na(DIFF_UP_minus_DOWN)) %>%
    mutate(RFS_risk = DIFF_UP_minus_DOWN,
           RFS_risk_z = as.numeric(scale(RFS_risk)))
  
  fwrite(df_rfs, file.path(subdir, "main_RFS_ERpos_merged.tsv"), sep = "\t")
  log_msg("[MAIN RFS] rows:", nrow(df_rfs), "| events:", sum(df_rfs$rfs_event == 1, na.rm = TRUE))
  
  # Wilcoxon
  w <- wilcox.test(RFS_risk ~ rfs_event, data = df_rfs)
  writeLines(capture.output(w), file.path(subdir, "wilcox_RFSrisk_by_recurrence.txt"))
  
  # Cox
  cox1 <- coxph(Surv(rfs_time, rfs_event) ~ RFS_risk_z, data = df_rfs)
  cox2 <- coxph(Surv(rfs_time, rfs_event) ~ RFS_risk_z + AGE_AT_DIAGNOSIS, data = df_rfs)
  
  sink(file.path(subdir, "cox_RFS_summary.txt"))
  cat("=== Cox: RFS ~ RFS_risk_z ===\n"); print(summary(cox1))
  cat("\n=== Cox: RFS ~ RFS_risk_z + AGE_AT_DIAGNOSIS ===\n"); print(summary(cox2))
  sink()
  
  # Median split KM
  df_med <- df_rfs %>%
    mutate(group_med = ifelse(RFS_risk >= median(RFS_risk, na.rm = TRUE), "High", "Low"),
           group_med = factor(group_med, levels = c("Low", "High")))
  
  fit_med <- survfit(Surv(rfs_time, rfs_event) ~ group_med, data = df_med)
  p_med <- ggsurvplot(
    fit_med, data = df_med, risk.table = TRUE, pval = TRUE, conf.int = TRUE,
    legend.title = "UP - DOWN risk",
    legend.labs = c("Low (median split)", "High (median split)"),
    xlab = "Time (months)", ylab = "Recurrence-free survival probability",
    ggtheme = theme_classic(base_size = 14),
    palette = c(COL_LOW, COL_HIGH)
  )
  ggsave(file.path(subdir, "KM_RFS_median_split.png"), p_med$plot, width = 8.5, height = 6.5, dpi = 200)
  ggsave(file.path(subdir, "KM_RFS_median_split_risktable.png"), p_med$table, width = 8.5, height = 3.0, dpi = 200)
  
  # Direction sanity
  s1 <- summary(cox1)
  hr <- as.numeric(s1$coef[1, "exp(coef)"])
  if (is.finite(hr) && hr < 1) {
    df_flip <- df_rfs %>% mutate(RFS_risk_flipped_z = as.numeric(scale(-RFS_risk)))
    cox_flip <- coxph(Surv(rfs_time, rfs_event) ~ RFS_risk_flipped_z, data = df_flip)
    sink(file.path(subdir, "cox_RFS_flipped_summary.txt"))
    cat("=== Cox: RFS ~ (-RFS_risk) scaled ===\n")
    print(summary(cox_flip))
    sink()
  }
}

# =========================================================
# B) MAIN ER+ EXTREMES RFS
# =========================================================
if (RUN_MAIN_RFS_ERPOS_EXTREMES) {
  subdir <- file.path(out_root, "02_main_RFS_ERpos_extremes")
  dir.create(subdir, recursive = TRUE, showWarnings = FALSE)
  
  df_rfs <- scores_patient_erpos %>%
    inner_join(
      clin_er %>% dplyr::select(patient_id, rfs_time, rfs_event, AGE_AT_DIAGNOSIS),
      by = "patient_id"
    ) %>%
    filter(!is.na(rfs_time), !is.na(rfs_event), !is.na(DIFF_UP_minus_DOWN)) %>%
    mutate(RFS_risk = DIFF_UP_minus_DOWN)
  
  df_ext <- df_rfs %>%
    mutate(extreme = make_extreme_groups_symmetric(RFS_risk, EXTREME_FRAC_MAIN)) %>%
    filter(!is.na(extreme))
  
  log_msg("[MAIN EXTREMES] rows:", nrow(df_ext), "| events:", sum(df_ext$rfs_event == 1, na.rm = TRUE))
  
  fit_ext <- survfit(Surv(rfs_time, rfs_event) ~ extreme, data = df_ext)
  p_ext <- ggsurvplot(
    fit_ext, data = df_ext, risk.table = TRUE, pval = TRUE, conf.int = TRUE,
    legend.title = "UP - DOWN extremes",
    legend.labs = levels(df_ext$extreme),
    xlab = "Time (months)", ylab = "Recurrence-free survival probability",
    ggtheme = theme_classic(base_size = 14),
    palette = c(COL_LOW, COL_HIGH)
  )
  ggsave(file.path(subdir, "KM_RFS_extremes.png"), p_ext$plot, width = 8.5, height = 6.5, dpi = 200)
  ggsave(file.path(subdir, "KM_RFS_extremes_risktable.png"), p_ext$table, width = 8.5, height = 3.0, dpi = 200)
  
  cox_ext <- coxph(Surv(rfs_time, rfs_event) ~ extreme, data = df_ext)
  cox_ext_age <- coxph(Surv(rfs_time, rfs_event) ~ extreme + AGE_AT_DIAGNOSIS, data = df_ext)
  
  sink(file.path(subdir, "cox_RFS_extremes_summary.txt"))
  cat("=== Cox: RFS ~ extremes ===\n"); print(summary(cox_ext))
  cat("\n=== Cox: RFS ~ extremes + AGE ===\n"); print(summary(cox_ext_age))
  sink()
}

# =========================================================
# C) COMPETING RISKS
# =========================================================
if (RUN_COMPETING_RISKS_ERPOS) {
  subdir <- file.path(out_root, "03_competing_risks_ERpos")
  dir.create(subdir, recursive = TRUE, showWarnings = FALSE)
  
  clin_cr <- clin_er %>%
    mutate(
      event_type = case_when(
        rfs_event == 1L ~ 1L,
        rfs_event == 0L & os_event == 1L ~ 2L,
        TRUE ~ 0L
      ),
      time_to_event = case_when(
        event_type == 1L ~ rfs_time,
        event_type == 2L ~ os_time,
        TRUE ~ pmax(rfs_time, os_time, na.rm = TRUE)
      )
    )
  
  df_cr <- scores_patient_erpos %>%
    inner_join(
      clin_cr %>% dplyr::select(patient_id, event_type, time_to_event, AGE_AT_DIAGNOSIS),
      by = "patient_id"
    ) %>%
    filter(!is.na(time_to_event), !is.na(event_type), !is.na(DIFF_UP_minus_DOWN)) %>%
    mutate(RFS_risk = DIFF_UP_minus_DOWN,
           RFS_risk_z = as.numeric(scale(RFS_risk)))
  
  fwrite(df_cr, file.path(subdir, "competing_risks_merged.tsv"), sep = "\t")
  log_msg("[COMPETING RISKS] rows:", nrow(df_cr))
  
  fg1 <- crr(
    ftime = df_cr$time_to_event,
    fstatus = df_cr$event_type,
    cov1 = as.matrix(df_cr$RFS_risk_z),
    failcode = 1,
    cencode = 0
  )
  
  sink(file.path(subdir, "finegray_continuous.txt"))
  cat("=== Fine-Gray: primary recurrence ~ RFS_risk_z ===\n")
  print(summary(fg1))
  sink()
  
  df_ext <- df_cr %>%
    mutate(extreme = make_extreme_groups_symmetric(RFS_risk, EXTREME_FRAC_CR)) %>%
    filter(!is.na(extreme)) %>%
    mutate(x_top = as.integer(extreme == levels(extreme)[2]))
  
  fg2 <- crr(
    ftime = df_ext$time_to_event,
    fstatus = df_ext$event_type,
    cov1 = as.matrix(df_ext$x_top),
    failcode = 1,
    cencode = 0
  )
  
  sink(file.path(subdir, "finegray_extremes.txt"))
  cat("=== Fine-Gray: primary recurrence ~ extremes ===\n")
  print(summary(fg2))
  sink()
  
  ci <- cuminc(
    ftime = df_ext$time_to_event,
    fstatus = df_ext$event_type,
    group = df_ext$extreme
  )
  
  png(file.path(subdir, "CIF_primary_recurrence_extremes.png"), width = 1100, height = 850)
  plot(ci, lwd = 3, xlab = "Months", ylab = "Cumulative incidence",
       main = "Primary recurrence cumulative incidence", curvlab = FALSE)
  cn <- names(ci)
  keep1 <- grepl(" 1$", cn)
  labs <- gsub(" 1$", "", cn[keep1])
  legend("topleft", legend = labs, lwd = 3, col = seq_along(labs), bty = "n")
  dev.off()
}

# =========================================================
# D) SEPARATE UP / DOWN / DIFF IN ER+ PATIENTS
# =========================================================
if (RUN_SEPARATE_UP_DOWN_ERPOS) {
  subdir <- file.path(out_root, "04_separate_UP_DOWN_DIFF_ERpos")
  dir.create(subdir, recursive = TRUE, showWarnings = FALSE)
  
  df <- scores_patient_erpos %>%
    inner_join(
      clin_er %>% dplyr::select(patient_id, rfs_time, rfs_event, AGE_AT_DIAGNOSIS),
      by = "patient_id"
    ) %>%
    filter(!is.na(rfs_time), !is.na(rfs_event))
  
  analyze_one_score <- function(df_in, score_col, label, out_prefix) {
    df0 <- df_in %>%
      filter(!is.na(.data[[score_col]]), !is.na(rfs_time), !is.na(rfs_event)) %>%
      mutate(score_z = as.numeric(scale(.data[[score_col]])))
    
    w <- wilcox.test(df0[[score_col]] ~ df0$rfs_event)
    writeLines(capture.output(w), file.path(subdir, paste0(out_prefix, "_wilcox.txt")))
    
    c1 <- coxph(Surv(rfs_time, rfs_event) ~ score_z, data = df0)
    c2 <- coxph(Surv(rfs_time, rfs_event) ~ score_z + AGE_AT_DIAGNOSIS, data = df0)
    
    sink(file.path(subdir, paste0(out_prefix, "_cox_continuous.txt")))
    cat("=== Cox:", label, "~ continuous ===\n"); print(summary(c1))
    cat("\n=== Cox:", label, "~ continuous + age ===\n"); print(summary(c2))
    sink()
    
    df_med <- df0 %>%
      mutate(group_med = ifelse(.data[[score_col]] >= median(.data[[score_col]], na.rm = TRUE), "High", "Low"),
             group_med = factor(group_med, levels = c("Low", "High")))
    
    fit_med <- survfit(Surv(rfs_time, rfs_event) ~ group_med, data = df_med)
    p_med <- ggsurvplot(
      fit_med, data = df_med, risk.table = TRUE, pval = TRUE, conf.int = TRUE,
      legend.title = paste0(label, " median split"),
      legend.labs = c("Low", "High"),
      xlab = "Time (months)", ylab = "Recurrence-free survival probability",
      ggtheme = theme_classic(base_size = 14),
      palette = c(COL_LOW, COL_HIGH)
    )
    ggsave(file.path(subdir, paste0(out_prefix, "_KM_median.png")), p_med$plot, width = 8.5, height = 6.5, dpi = 200)
    ggsave(file.path(subdir, paste0(out_prefix, "_KM_median_risktable.png")), p_med$table, width = 8.5, height = 3.0, dpi = 200)
    
    df_ext <- df0 %>%
      mutate(extreme = make_extreme_groups_symmetric(.data[[score_col]], 0.20)) %>%
      filter(!is.na(extreme))
    
    fit_ext <- survfit(Surv(rfs_time, rfs_event) ~ extreme, data = df_ext)
    p_ext <- ggsurvplot(
      fit_ext, data = df_ext, risk.table = TRUE, pval = TRUE, conf.int = TRUE,
      legend.title = paste0(label, " extremes"),
      legend.labs = levels(df_ext$extreme),
      xlab = "Time (months)", ylab = "Recurrence-free survival probability",
      ggtheme = theme_classic(base_size = 14),
      palette = c(COL_LOW, COL_HIGH)
    )
    ggsave(file.path(subdir, paste0(out_prefix, "_KM_extremes.png")), p_ext$plot, width = 8.5, height = 6.5, dpi = 200)
    ggsave(file.path(subdir, paste0(out_prefix, "_KM_extremes_risktable.png")), p_ext$table, width = 8.5, height = 3.0, dpi = 200)
  }
  
  analyze_one_score(df, "UP_score", "UP", "UP")
  analyze_one_score(df, "DOWN_score", "DOWN", "DOWN")
  analyze_one_score(df, "DIFF_UP_minus_DOWN", "UP_minus_DOWN", "DIFF")
}

# =========================================================
# E) ER+ / ET-TREATED / 120-MONTH ANALYSIS
# =========================================================
if (RUN_ET_120MO) {
  subdir <- file.path(out_root, "05_ETtreated_120mo")
  dir.create(subdir, recursive = TRUE, showWarnings = FALSE)
  
  clin_et <- clin_er %>%
    filter(HORMONE_THERAPY_CLEAN == "YES")
  
  df <- scores_patient_erpos %>%
    inner_join(
      clin_et %>% dplyr::select(patient_id, rfs_time, rfs_event, AGE_AT_DIAGNOSIS),
      by = "patient_id"
    )
  
  fwrite(df, file.path(subdir, "merged_scores_clinical_ERpos_ETtreated.tsv"), sep = "\t")
  log_msg("[ET 120mo] rows:", nrow(df))
  
  run_surv_block <- function(df, score_name, out_prefix, top_p, bottom_p, tmax) {
    d <- df %>%
      filter(!is.na(.data[[score_name]]), !is.na(rfs_time), !is.na(rfs_event))
    
    cens <- admin_censor(d$rfs_time, d$rfs_event, tmax)
    d$time_c  <- cens$time
    d$event_c <- cens$event
    d$score_z <- as.numeric(scale(d[[score_name]]))
    
    d$extreme <- make_extreme_groups_asymmetric(d[[score_name]], top_p, bottom_p)
    d_ext <- d %>% filter(!is.na(extreme))
    
    cox_cont <- coxph(Surv(time_c, event_c) ~ score_z, data = d)
    cox_cont_age <- coxph(Surv(time_c, event_c) ~ score_z + AGE_AT_DIAGNOSIS, data = d)
    cox_ext <- coxph(Surv(time_c, event_c) ~ extreme, data = d_ext)
    cox_ext_age <- coxph(Surv(time_c, event_c) ~ extreme + AGE_AT_DIAGNOSIS, data = d_ext)
    
    sink(file.path(subdir, paste0(out_prefix, "_cox_summary.txt")))
    cat("=== Settings ===\n")
    cat("score:", score_name, "\n")
    cat("tmax:", tmax, "\n")
    cat("top_p:", top_p, "| bottom_p:", bottom_p, "\n\n")
    cat("=== Continuous ===\n"); print(summary(cox_cont))
    cat("\n=== Continuous + age ===\n"); print(summary(cox_cont_age))
    cat("\n=== Extremes ===\n"); print(summary(cox_ext))
    cat("\n=== Extremes + age ===\n"); print(summary(cox_ext_age))
    sink()
    
    fit <- survfit(Surv(time_c, event_c) ~ extreme, data = d_ext)
    p <- ggsurvplot(
      fit, data = d_ext, risk.table = TRUE, pval = TRUE, conf.int = TRUE,
      xlim = c(0, tmax), break.time.by = 12,
      legend.title = paste0(score_name, " extremes"),
      legend.labs = levels(d_ext$extreme),
      ggtheme = theme_classic(),
      palette = c(COL_LOW, COL_HIGH)
    )
    ggsave(file.path(subdir, paste0(out_prefix, "_KM.png")),
           plot = print(p), width = 7.5, height = 7.0, dpi = 200)
    
    invisible(d_ext)
  }
  
  run_surv_block(df, "UP_score",
                 sprintf("RFS_ETtreated_UP_top%.0f_bottom%.0f_tmax%d", ET_CUT_TOP * 100, ET_CUT_BOTTOM * 100, ET_TMAX_MONTHS),
                 ET_CUT_TOP, ET_CUT_BOTTOM, ET_TMAX_MONTHS)
  
  run_surv_block(df, "DOWN_score",
                 sprintf("RFS_ETtreated_DOWN_top%.0f_bottom%.0f_tmax%d", ET_CUT_TOP * 100, ET_CUT_BOTTOM * 100, ET_TMAX_MONTHS),
                 ET_CUT_TOP, ET_CUT_BOTTOM, ET_TMAX_MONTHS)
  
  run_surv_block(df, "DIFF_UP_minus_DOWN",
                 sprintf("RFS_ETtreated_DIFF_top%.0f_bottom%.0f_tmax%d", ET_CUT_TOP * 100, ET_CUT_BOTTOM * 100, ET_TMAX_MONTHS),
                 ET_CUT_TOP, ET_CUT_BOTTOM, ET_TMAX_MONTHS)
}

# =========================================================
# F) GRID SEARCH FOR ET-TREATED 120MO
# =========================================================
if (RUN_ET_120MO_GRID_SEARCH) {
  subdir <- file.path(out_root, "06_ETtreated_120mo_gridsearch")
  dir.create(subdir, recursive = TRUE, showWarnings = FALSE)
  
  clin_et <- clin_er %>%
    filter(HORMONE_THERAPY_CLEAN == "YES")
  
  df <- scores_patient_erpos %>%
    inner_join(
      clin_et %>% dplyr::select(patient_id, rfs_time, rfs_event, AGE_AT_DIAGNOSIS),
      by = "patient_id"
    )
  
  safe_cox_extreme <- function(d_ext, age_col = NULL) {
    out <- list(ok = FALSE, hr = NA_real_, lo = NA_real_, hi = NA_real_, p = NA_real_, coef = NA_real_)
    if (nrow(d_ext) < 10) return(out)
    if (length(unique(d_ext$extreme)) < 2) return(out)
    
    ev_by_grp <- tapply(d_ext$event_c, d_ext$extreme, function(x) sum(x == 1, na.rm = TRUE))
    if (any(is.na(ev_by_grp)) || sum(ev_by_grp) < 5) return(out)
    
    f <- if (is.null(age_col)) {
      as.formula("Surv(time_c, event_c) ~ extreme")
    } else {
      as.formula(paste0("Surv(time_c, event_c) ~ extreme + ", age_col))
    }
    
    fit <- tryCatch(coxph(f, data = d_ext), error = function(e) NULL)
    if (is.null(fit)) return(out)
    s <- summary(fit)
    rn <- rownames(s$coef)
    idx <- grep("^extreme", rn)[1]
    if (is.na(idx)) return(out)
    
    out$ok   <- TRUE
    out$coef <- unname(s$coef[idx, "coef"])
    out$hr   <- unname(s$coef[idx, "exp(coef)"])
    out$p    <- unname(s$coef[idx, "Pr(>|z|)"])
    out$lo   <- unname(s$conf.int[idx, "lower .95"])
    out$hi   <- unname(s$conf.int[idx, "upper .95"])
    out
  }
  
  run_cutoff_grid_one_score <- function(df, score_name) {
    d0 <- df %>%
      filter(!is.na(.data[[score_name]]), !is.na(rfs_time), !is.na(rfs_event))
    
    cens <- admin_censor(d0$rfs_time, d0$rfs_event, ET_TMAX_MONTHS)
    d0$time_c  <- cens$time
    d0$event_c <- cens$event
    
    res <- list()
    k <- 0L
    
    for (top_p in CUT_GRID) {
      for (bottom_p in CUT_GRID) {
        if ((top_p + bottom_p) >= 1) next
        
        d <- d0
        d$extreme <- make_extreme_groups_asymmetric(d[[score_name]], top_p, bottom_p)
        d_ext <- d %>% filter(!is.na(extreme))
        
        lv <- levels(d_ext$extreme)
        n_bot <- ifelse(length(lv) >= 1, sum(d_ext$extreme == lv[1], na.rm = TRUE), NA_integer_)
        n_top <- ifelse(length(lv) >= 2, sum(d_ext$extreme == lv[2], na.rm = TRUE), NA_integer_)
        ev_bot <- ifelse(length(lv) >= 1, sum(d_ext$event_c == 1 & d_ext$extreme == lv[1], na.rm = TRUE), NA_integer_)
        ev_top <- ifelse(length(lv) >= 2, sum(d_ext$event_c == 1 & d_ext$extreme == lv[2], na.rm = TRUE), NA_integer_)
        ev_tot <- sum(d_ext$event_c == 1, na.rm = TRUE)
        
        if (is.na(n_top) || is.na(n_bot) || n_top < 20 || n_bot < 20) next
        if (ev_tot < 20) next
        
        fit0 <- safe_cox_extreme(d_ext, age_col = NULL)
        fitA <- safe_cox_extreme(d_ext, age_col = "AGE_AT_DIAGNOSIS")
        
        k <- k + 1L
        res[[k]] <- tibble(
          score = score_name,
          tmax_months = ET_TMAX_MONTHS,
          top_p = top_p,
          bottom_p = bottom_p,
          n_ext = nrow(d_ext),
          events_ext = ev_tot,
          n_bottom = n_bot,
          n_top = n_top,
          events_bottom = ev_bot,
          events_top = ev_top,
          hr = fit0$hr, lo95 = fit0$lo, hi95 = fit0$hi, p = fit0$p,
          hr_age = fitA$hr, lo95_age = fitA$lo, hi95_age = fitA$hi, p_age = fitA$p,
          logHR = ifelse(is.na(fit0$hr), NA_real_, log(fit0$hr)),
          logHR_age = ifelse(is.na(fitA$hr), NA_real_, log(fitA$hr))
        )
      }
    }
    
    if (length(res) == 0) return(tibble())
    bind_rows(res)
  }
  
  all_res <- bind_rows(
    run_cutoff_grid_one_score(df, "UP_score"),
    run_cutoff_grid_one_score(df, "DOWN_score"),
    run_cutoff_grid_one_score(df, "DIFF_UP_minus_DOWN")
  )
  
  fwrite(all_res, file.path(subdir, sprintf("cutoff_grid_results_ETtreated_tmax%d.tsv", ET_TMAX_MONTHS)), sep = "\t")
  
  best_per_score <- all_res %>%
    filter(!is.na(p)) %>%
    group_by(score) %>%
    arrange(p, desc(abs(logHR)), desc(n_ext)) %>%
    dplyr::slice(1) %>%
    ungroup()
  
  fwrite(best_per_score, file.path(subdir, sprintf("cutoff_grid_best_per_score_ETtreated_tmax%d.tsv", ET_TMAX_MONTHS)), sep = "\t")
}

# =========================================================
# G) REPLOT BEST GRID CONFIGS
# =========================================================
if (RUN_ET_120MO_REPLOT_BEST) {
  subdir_in  <- file.path(out_root, "06_ETtreated_120mo_gridsearch")
  subdir_out <- file.path(out_root, "07_ETtreated_120mo_best_replots")
  dir.create(subdir_out, recursive = TRUE, showWarnings = FALSE)
  
  best_file <- file.path(subdir_in, sprintf("cutoff_grid_best_per_score_ETtreated_tmax%d.tsv", ET_TMAX_MONTHS))
  if (file.exists(best_file)) {
    best_cfg <- fread(best_file)
    
    clin_et <- clin_er %>% filter(HORMONE_THERAPY_CLEAN == "YES")
    df <- scores_patient_erpos %>%
      inner_join(clin_et %>% dplyr::select(patient_id, rfs_time, rfs_event, AGE_AT_DIAGNOSIS),
                 by = "patient_id")
    
    plot_best_survival <- function(df, score_name, top_p, bottom_p, tmax) {
      d <- df %>%
        filter(!is.na(.data[[score_name]]), !is.na(rfs_time), !is.na(rfs_event))
      
      cens <- admin_censor(d$rfs_time, d$rfs_event, tmax)
      d$time_c  <- cens$time
      d$event_c <- cens$event
      
      d$extreme <- make_extreme_groups_asymmetric(d[[score_name]], top_p, bottom_p)
      d_ext <- d %>% filter(!is.na(extreme))
      
      fit <- survfit(Surv(time_c, event_c) ~ extreme, data = d_ext)
      
      p <- ggsurvplot(
        fit, data = d_ext, risk.table = TRUE, pval = TRUE, conf.int = TRUE,
        xlim = c(0, tmax), xlab = "Time (months)", ylab = "Recurrence-Free Survival",
        ggtheme = theme_classic(), legend.title = score_name,
        palette = c(COL_LOW, COL_HIGH)
      )
      
      ggsave(file.path(subdir_out, sprintf("%s_top%.0f_bottom%.0f.png", score_name, top_p * 100, bottom_p * 100)),
             p$plot, width = 8, height = 6, dpi = 200)
    }
    
    for (i in seq_len(nrow(best_cfg))) {
      plot_best_survival(
        df = df,
        score_name = best_cfg$score[i],
        top_p = best_cfg$top_p[i],
        bottom_p = best_cfg$bottom_p[i],
        tmax = ET_TMAX_MONTHS
      )
    }
  }
}

# =========================================================
# H) ET-TREATED + ESR1 COVARIATE / VARIANCE CHECKS
# =========================================================
if (RUN_ET_120MO_WITH_ESR1_COVARIATE && !is.null(esr1_patient)) {
  subdir <- file.path(out_root, "08_ETtreated_120mo_with_ESR1_covariate")
  dir.create(subdir, recursive = TRUE, showWarnings = FALSE)
  
  clin_et <- clin_er %>% filter(HORMONE_THERAPY_CLEAN == "YES")
  
  df <- scores_patient_erpos %>%
    inner_join(clin_et %>% dplyr::select(patient_id, rfs_time, rfs_event, AGE_AT_DIAGNOSIS),
               by = "patient_id") %>%
    left_join(esr1_patient, by = "patient_id") %>%
    filter(!is.na(rfs_time), !is.na(rfs_event))
  
  best_file <- file.path(out_root, "06_ETtreated_120mo_gridsearch",
                         sprintf("cutoff_grid_best_per_score_ETtreated_tmax%d.tsv", ET_TMAX_MONTHS))
  
  if (file.exists(best_file)) {
    best_cfg <- fread(best_file)
    
    run_best_config_report <- function(df, score_name, top_p, bottom_p, out_prefix) {
      d0 <- df %>%
        filter(!is.na(.data[[score_name]]), !is.na(rfs_time), !is.na(rfs_event))
      
      cens <- admin_censor(d0$rfs_time, d0$rfs_event, ET_TMAX_MONTHS)
      d0$time_c  <- cens$time
      d0$event_c <- cens$event
      d0$score <- as.numeric(d0[[score_name]])
      d0$score_z <- as.numeric(scale(d0$score))
      
      d_ext <- d0 %>%
        mutate(extreme = make_extreme_groups_asymmetric(score, top_p, bottom_p)) %>%
        filter(!is.na(extreme))
      
      cox_cont <- safe_cox(Surv(time_c, event_c) ~ score_z, d0)
      cox_cont_age <- safe_cox(Surv(time_c, event_c) ~ score_z + AGE_AT_DIAGNOSIS, d0)
      cox_cont_age_esr1 <- safe_cox(Surv(time_c, event_c) ~ score_z + AGE_AT_DIAGNOSIS + ESR1_expr_mean, d0)
      
      cox_ext <- safe_cox(Surv(time_c, event_c) ~ extreme, d_ext)
      cox_ext_age <- safe_cox(Surv(time_c, event_c) ~ extreme + AGE_AT_DIAGNOSIS, d_ext)
      cox_ext_age_esr1 <- safe_cox(Surv(time_c, event_c) ~ extreme + AGE_AT_DIAGNOSIS + ESR1_expr_mean, d_ext)
      
      var_tbl <- bind_rows(
        if (!all(is.na(d0$AGE_AT_DIAGNOSIS))) {
          ss <- summary(lm(score ~ AGE_AT_DIAGNOSIS, data = d0))
          tibble(metric = c("score_vs_age_slope", "score_vs_age_R2"),
                 estimate = c(unname(coef(lm(score ~ AGE_AT_DIAGNOSIS, data = d0))[["AGE_AT_DIAGNOSIS"]]), ss$r.squared),
                 p_value = c(ss$coefficients["AGE_AT_DIAGNOSIS", "Pr(>|t|)"], NA_real_))
        } else tibble(),
        if (!all(is.na(d0$ESR1_expr_mean))) {
          ss <- summary(lm(score ~ ESR1_expr_mean, data = d0))
          tibble(metric = c("score_vs_ESR1_slope", "score_vs_ESR1_R2"),
                 estimate = c(unname(coef(lm(score ~ ESR1_expr_mean, data = d0))[["ESR1_expr_mean"]]), ss$r.squared),
                 p_value = c(ss$coefficients["ESR1_expr_mean", "Pr(>|t|)"], NA_real_))
        } else tibble()
      )
      
      sink(file.path(subdir, paste0(out_prefix, "_FULL_REPORT.txt")))
      cat("=== Continuous ===\n")
      if (inherits(cox_cont, "error")) cat("ERROR:", cox_cont$message, "\n") else print(summary(cox_cont))
      cat("\n=== Continuous + age ===\n")
      if (inherits(cox_cont_age, "error")) cat("ERROR:", cox_cont_age$message, "\n") else print(summary(cox_cont_age))
      cat("\n=== Continuous + age + ESR1 ===\n")
      if (inherits(cox_cont_age_esr1, "error")) cat("ERROR:", cox_cont_age_esr1$message, "\n") else print(summary(cox_cont_age_esr1))
      cat("\n=== Extremes ===\n")
      if (inherits(cox_ext, "error")) cat("ERROR:", cox_ext$message, "\n") else print(summary(cox_ext))
      cat("\n=== Extremes + age ===\n")
      if (inherits(cox_ext_age, "error")) cat("ERROR:", cox_ext_age$message, "\n") else print(summary(cox_ext_age))
      cat("\n=== Extremes + age + ESR1 ===\n")
      if (inherits(cox_ext_age_esr1, "error")) cat("ERROR:", cox_ext_age_esr1$message, "\n") else print(summary(cox_ext_age_esr1))
      cat("\n=== Variance checks ===\n")
      print(var_tbl)
      sink()
      
      fwrite(var_tbl, file.path(subdir, paste0(out_prefix, "_variance.tsv")), sep = "\t")
    }
    
    for (i in seq_len(nrow(best_cfg))) {
      run_best_config_report(
        df = df,
        score_name = best_cfg$score[i],
        top_p = best_cfg$top_p[i],
        bottom_p = best_cfg$bottom_p[i],
        out_prefix = sprintf("RFS_ETtreated_%s_best_top%.0f_bottom%.0f_tmax%d",
                             best_cfg$score[i], best_cfg$top_p[i] * 100, best_cfg$bottom_p[i] * 100, ET_TMAX_MONTHS)
      )
    }
  }
}

log_msg("=== MASTER ANALYSIS DONE ===")


