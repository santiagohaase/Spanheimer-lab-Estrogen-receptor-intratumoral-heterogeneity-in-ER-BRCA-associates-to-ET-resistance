#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
  library(dplyr)
  library(readr)
  library(tibble)
  library(stringr)
  library(survival)
  library(survminer)
  library(ggplot2)
})

## ===============================================================
## GSE17705_low_ESR1_UP_DW_delta_survival.R
##
## Dataset:
##   GSE17705
##
## Endpoint:
##   distant relapse (1 = distant relapse, 0 = censored)
##   event time in years
##
## Analysis:
##   - Collapse probe-level expression to gene symbols
##   - Filter ER+ tumors
##   - Score LOW_ESR1_UP
##   - Score LOW_ESR1_DW
##   - Compute LOW_ESR1_DELTA = LOW_ESR1_UP - LOW_ESR1_DW
##   - Test distant relapse-free survival:
##       * continuous Cox
##       * nodal-status-adjusted Cox
##       * median high vs low KM
##       * top/bottom 10% KM
##       * asymmetric cutoff scan
## ===============================================================

## -----------------------------
## Paths
## -----------------------------

base_dir <- "/path/to/GEO_recurrence_validation/GSE17705/low_ESR1_UP_DW_delta"
dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)

eset_rds <- file.path(base_dir, "GSE17705_eset_raw.rds")

expr_gene_out <- file.path(base_dir, "GSE17705_expr_gene_collapsed.rds")
pheno_clean_out <- file.path(base_dir, "GSE17705_pheno_clean.csv")
score_out <- file.path(base_dir, "GSE17705_low_ESR1_UP_DW_delta_scores.csv")
overlap_out <- file.path(base_dir, "GSE17705_low_ESR1_UP_DW_delta_gene_overlap.csv")
cox_out <- file.path(base_dir, "GSE17705_low_ESR1_UP_DW_delta_cox_results.csv")
cox_txt <- file.path(base_dir, "GSE17705_low_ESR1_UP_DW_delta_cox_results.txt")

plot_dir <- file.path(base_dir, "plots_low_ESR1_UP_DW_delta")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

cutoff_scan_dir <- file.path(base_dir, "asymmetric_cutoff_scan_low_ESR1_UP_DW_delta")
dir.create(cutoff_scan_dir, recursive = TRUE, showWarnings = FALSE)

## -----------------------------
## Signatures
## -----------------------------
## Paste your actual gene vectors here.
## Names should be gene symbols.

LOW_ESR1_UP <- c(
  ## PASTE low_ESR1_up genes here
)

LOW_ESR1_DW <- c(
  ## PASTE low_ESR1_dw genes here
)

sig_list <- list(
  LOW_ESR1_UP = LOW_ESR1_UP,
  LOW_ESR1_DW = LOW_ESR1_DW
)

## -----------------------------
## Helpers
## -----------------------------

clean_gene <- function(x) {
  x <- as.character(x)
  x <- str_trim(x)
  x <- toupper(x)
  x[x == ""] <- NA_character_
  x
}

zscore_by_gene <- function(mat) {
  z <- t(scale(t(mat)))
  z[is.na(z)] <- 0
  z
}

score_signature <- function(zmat, genes, min_genes = 5) {
  genes <- clean_gene(genes)
  genes <- unique(genes[!is.na(genes)])
  
  overlap <- intersect(genes, rownames(zmat))
  missing <- setdiff(genes, rownames(zmat))
  
  if (length(overlap) < min_genes) {
    score <- rep(NA_real_, ncol(zmat))
  } else {
    score <- colMeans(zmat[overlap, , drop = FALSE], na.rm = TRUE)
  }
  
  list(
    score = score,
    overlap = overlap,
    missing = missing
  )
}

collapse_probes_to_genes <- function(expr, feature, symbol_col = "Gene symbol") {
  if (!symbol_col %in% colnames(feature)) {
    stop("Gene symbol column not found: ", symbol_col)
  }
  
  probe_ids <- rownames(expr)
  
  feature2 <- feature %>%
    rownames_to_column("probe_id") %>%
    select(probe_id, gene_symbol = all_of(symbol_col)) %>%
    mutate(
      gene_symbol = clean_gene(gene_symbol)
    )
  
  ## Remove ambiguous probes mapping to multiple genes, e.g. "MIR4640/DDR1"
  feature2 <- feature2 %>%
    mutate(
      has_multiple = grepl("/", gene_symbol, fixed = TRUE)
    ) %>%
    filter(
      !is.na(gene_symbol),
      gene_symbol != "",
      !has_multiple,
      probe_id %in% probe_ids
    )
  
  expr2 <- expr[feature2$probe_id, , drop = FALSE]
  
  expr_df <- as.data.frame(expr2) %>%
    rownames_to_column("probe_id") %>%
    left_join(feature2 %>% select(probe_id, gene_symbol), by = "probe_id") %>%
    select(gene_symbol, everything(), -probe_id)
  
  message("Probes with unique gene symbols: ", nrow(expr_df))
  message("Unique genes before collapsing: ", length(unique(expr_df$gene_symbol)))
  
  expr_gene <- expr_df %>%
    group_by(gene_symbol) %>%
    summarise(
      across(everything(), ~ mean(as.numeric(.x), na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    column_to_rownames("gene_symbol") %>%
    as.matrix()
  
  storage.mode(expr_gene) <- "double"
  
  message("Collapsed expression matrix: ", nrow(expr_gene), " genes x ", ncol(expr_gene), " samples")
  
  expr_gene
}

run_cox <- function(df, score_col, covars = character(0)) {
  use_cols <- c("event_time_years", "distant_relapse", score_col, covars)
  
  df2 <- df %>%
    select(all_of(use_cols)) %>%
    filter(if_all(everything(), ~ !is.na(.x)))
  
  model_name <- ifelse(
    length(covars) == 0,
    "univariable",
    paste0("adjusted_", paste(covars, collapse = "_"))
  )
  
  if (nrow(df2) < 50 || sum(df2$distant_relapse == 1, na.rm = TRUE) < 10) {
    return(tibble(
      signature = score_col,
      model = model_name,
      n = nrow(df2),
      events = sum(df2$distant_relapse == 1, na.rm = TRUE),
      HR = NA_real_,
      CI_low = NA_real_,
      CI_high = NA_real_,
      p = NA_real_,
      note = "too_few_samples_or_events"
    ))
  }
  
  form <- as.formula(
    paste0(
      "Surv(event_time_years, distant_relapse) ~ ",
      paste(c(score_col, covars), collapse = " + ")
    )
  )
  
  fit <- tryCatch(
    coxph(form, data = df2),
    error = function(e) e
  )
  
  if (inherits(fit, "error")) {
    return(tibble(
      signature = score_col,
      model = model_name,
      n = nrow(df2),
      events = sum(df2$distant_relapse == 1, na.rm = TRUE),
      HR = NA_real_,
      CI_low = NA_real_,
      CI_high = NA_real_,
      p = NA_real_,
      note = conditionMessage(fit)
    ))
  }
  
  s <- summary(fit)
  
  tibble(
    signature = score_col,
    model = model_name,
    n = nrow(df2),
    events = sum(df2$distant_relapse == 1, na.rm = TRUE),
    HR = unname(s$coefficients[score_col, "exp(coef)"]),
    CI_low = unname(s$conf.int[score_col, "lower .95"]),
    CI_high = unname(s$conf.int[score_col, "upper .95"]),
    p = unname(s$coefficients[score_col, "Pr(>|z|)"]),
    note = "ok"
  )
}

plot_km_split <- function(df, score_col, split_type = c("median", "top_bottom10"), out_dir) {
  split_type <- match.arg(split_type)
  
  df0 <- df %>%
    filter(
      !is.na(.data[[score_col]]),
      !is.na(event_time_years),
      !is.na(distant_relapse)
    )
  
  if (split_type == "median") {
    cut <- median(df0[[score_col]], na.rm = TRUE)
    
    df_plot <- df0 %>%
      mutate(
        score_group = ifelse(.data[[score_col]] >= cut, "high", "low"),
        score_group = factor(score_group, levels = c("low", "high"))
      )
    
    title <- paste0(score_col, " high vs low")
    file_stub <- paste0(score_col, "_median_high_low")
    
  } else {
    q10 <- quantile(df0[[score_col]], 0.10, na.rm = TRUE)
    q90 <- quantile(df0[[score_col]], 0.90, na.rm = TRUE)
    
    df_plot <- df0 %>%
      mutate(
        score_group = case_when(
          .data[[score_col]] <= q10 ~ "bottom10",
          .data[[score_col]] >= q90 ~ "top10",
          TRUE ~ NA_character_
        )
      ) %>%
      filter(!is.na(score_group)) %>%
      mutate(score_group = factor(score_group, levels = c("bottom10", "top10")))
    
    title <- paste0(score_col, " top 10% vs bottom 10%")
    file_stub <- paste0(score_col, "_top10_bottom10")
  }
  
  fit <- survfit(Surv(event_time_years, distant_relapse) ~ score_group, data = df_plot)
  
  p <- ggsurvplot(
    fit,
    data = df_plot,
    risk.table = TRUE,
    pval = TRUE,
    conf.int = FALSE,
    title = title,
    xlab = "Distant relapse-free time (years)",
    ylab = "Distant relapse-free probability"
  )
  
  pdf_file <- file.path(out_dir, paste0(file_stub, "_KM.pdf"))
  png_file <- file.path(out_dir, paste0(file_stub, "_KM.png"))
  
  pdf(pdf_file, width = 7, height = 6)
  print(p)
  dev.off()
  
  png(png_file, width = 7, height = 6, units = "in", res = 300)
  print(p)
  dev.off()
  
  message("Saved KM: ", pdf_file)
  
  invisible(p)
}

run_asymmetric_km_scan <- function(
    df,
    score_col,
    bottom_grid = seq(0.10, 0.50, by = 0.05),
    top_grid = seq(0.05, 0.50, by = 0.05),
    min_events_per_group = 5,
    min_samples_per_group = 20
) {
  df0 <- df %>%
    filter(
      !is.na(.data[[score_col]]),
      !is.na(event_time_years),
      !is.na(distant_relapse)
    )
  
  grid <- expand.grid(
    bottom_pct = bottom_grid,
    top_pct = top_grid
  ) %>%
    as_tibble()
  
  out <- lapply(seq_len(nrow(grid)), function(i) {
    bottom_pct <- grid$bottom_pct[i]
    top_pct <- grid$top_pct[i]
    
    if ((bottom_pct + top_pct) >= 0.95) {
      return(NULL)
    }
    
    q_low <- quantile(df0[[score_col]], bottom_pct, na.rm = TRUE)
    q_high <- quantile(df0[[score_col]], 1 - top_pct, na.rm = TRUE)
    
    df_extreme <- df0 %>%
      mutate(
        score_group = case_when(
          .data[[score_col]] <= q_low ~ "bottom",
          .data[[score_col]] >= q_high ~ "top",
          TRUE ~ NA_character_
        )
      ) %>%
      filter(!is.na(score_group)) %>%
      mutate(score_group = factor(score_group, levels = c("bottom", "top")))
    
    n_bottom <- sum(df_extreme$score_group == "bottom")
    n_top <- sum(df_extreme$score_group == "top")
    
    events_bottom <- sum(df_extreme$distant_relapse[df_extreme$score_group == "bottom"] == 1, na.rm = TRUE)
    events_top <- sum(df_extreme$distant_relapse[df_extreme$score_group == "top"] == 1, na.rm = TRUE)
    
    if (
      n_bottom < min_samples_per_group ||
      n_top < min_samples_per_group ||
      events_bottom < min_events_per_group ||
      events_top < min_events_per_group
    ) {
      return(tibble(
        signature = score_col,
        bottom_pct = bottom_pct,
        top_pct = top_pct,
        n_bottom = n_bottom,
        n_top = n_top,
        events_bottom = events_bottom,
        events_top = events_top,
        HR_top_vs_bottom = NA_real_,
        CI_low = NA_real_,
        CI_high = NA_real_,
        logrank_p = NA_real_,
        cox_p = NA_real_,
        note = "too_few_samples_or_events"
      ))
    }
    
    fit_survdiff <- survdiff(
      Surv(event_time_years, distant_relapse) ~ score_group,
      data = df_extreme
    )
    
    logrank_p <- 1 - pchisq(fit_survdiff$chisq, df = length(fit_survdiff$n) - 1)
    
    fit_cox <- coxph(
      Surv(event_time_years, distant_relapse) ~ score_group,
      data = df_extreme
    )
    
    s <- summary(fit_cox)
    coef_name <- rownames(s$coefficients)[1]
    
    tibble(
      signature = score_col,
      bottom_pct = bottom_pct,
      top_pct = top_pct,
      n_bottom = n_bottom,
      n_top = n_top,
      events_bottom = events_bottom,
      events_top = events_top,
      HR_top_vs_bottom = unname(s$coefficients[coef_name, "exp(coef)"]),
      CI_low = unname(s$conf.int[coef_name, "lower .95"]),
      CI_high = unname(s$conf.int[coef_name, "upper .95"]),
      logrank_p = logrank_p,
      cox_p = unname(s$coefficients[coef_name, "Pr(>|z|)"]),
      note = "ok"
    )
  })
  
  bind_rows(out) %>%
    mutate(
      logrank_FDR = p.adjust(logrank_p, method = "BH"),
      cox_FDR = p.adjust(cox_p, method = "BH")
    ) %>%
    arrange(logrank_p)
}

plot_best_asymmetric_km <- function(df, score_col, bottom_pct, top_pct, out_dir) {
  df0 <- df %>%
    filter(
      !is.na(.data[[score_col]]),
      !is.na(event_time_years),
      !is.na(distant_relapse)
    )
  
  q_low <- quantile(df0[[score_col]], bottom_pct, na.rm = TRUE)
  q_high <- quantile(df0[[score_col]], 1 - top_pct, na.rm = TRUE)
  
  bottom_label <- paste0("bottom", round(bottom_pct * 100), "%")
  top_label <- paste0("top", round(top_pct * 100), "%")
  
  df_extreme <- df0 %>%
    mutate(
      score_group = case_when(
        .data[[score_col]] <= q_low ~ bottom_label,
        .data[[score_col]] >= q_high ~ top_label,
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(score_group)) %>%
    mutate(score_group = factor(score_group, levels = c(bottom_label, top_label)))
  
  fit <- survfit(Surv(event_time_years, distant_relapse) ~ score_group, data = df_extreme)
  
  p <- ggsurvplot(
    fit,
    data = df_extreme,
    risk.table = TRUE,
    pval = TRUE,
    conf.int = FALSE,
    title = paste0(score_col, ": ", top_label, " vs ", bottom_label),
    xlab = "Distant relapse-free time (years)",
    ylab = "Distant relapse-free probability"
  )
  
  file_stub <- paste0(
    score_col,
    "_top", round(top_pct * 100),
    "_bottom", round(bottom_pct * 100)
  )
  
  pdf_file <- file.path(out_dir, paste0(file_stub, "_KM.pdf"))
  png_file <- file.path(out_dir, paste0(file_stub, "_KM.png"))
  
  pdf(pdf_file, width = 7, height = 6)
  print(p)
  dev.off()
  
  png(png_file, width = 7, height = 6, units = "in", res = 300)
  print(p)
  dev.off()
  
  message("Saved best asymmetric KM: ", pdf_file)
  
  invisible(p)
}

## -----------------------------
## Load data
## -----------------------------

if (!file.exists(eset_rds)) {
  message("Raw eset RDS not found. Downloading from GEO.")
  gse <- getGEO("GSE17705", GSEMatrix = TRUE, AnnotGPL = TRUE, destdir = base_dir)
  eset <- gse[[1]]
  saveRDS(eset, eset_rds)
} else {
  message("Loading existing eset: ", eset_rds)
  eset <- readRDS(eset_rds)
}

expr <- exprs(eset)
pheno <- pData(eset)
feature <- fData(eset)

message("Expression: ", nrow(expr), " probes x ", ncol(expr), " samples")
message("Phenotype: ", nrow(pheno), " samples x ", ncol(pheno), " columns")
message("Feature: ", nrow(feature), " probes x ", ncol(feature), " columns")

## -----------------------------
## Clean phenotype
## -----------------------------

pheno_clean <- pheno %>%
  rownames_to_column("gsm") %>%
  transmute(
    gsm = gsm,
    sample_title = title,
    patient_id = as.numeric(`patient id:ch1`),
    er_status = `estrogen receptor (er) status:ch1`,
    endocrine_therapy = `endocrine therapy:ch1`,
    nodal_status = suppressWarnings(as.numeric(`nodal status (0=negative, 1=positive, na=not applicable):ch1`)),
    distant_relapse = suppressWarnings(as.numeric(`distant relapse (1=dr, 0 censored):ch1`)),
    event_time_years = suppressWarnings(as.numeric(`event time (years):ch1`)),
    profiling_lab = `profiling lab:ch1`,
    tissue = `tissue:ch1`
  )

message("Phenotype summary before filtering:")
print(table(pheno_clean$er_status, useNA = "ifany"))
print(table(pheno_clean$distant_relapse, useNA = "ifany"))
print(table(pheno_clean$nodal_status, useNA = "ifany"))
print(summary(pheno_clean$event_time_years))

write_csv(pheno_clean, pheno_clean_out)

## -----------------------------
## Collapse expression to genes
## -----------------------------

expr_gene <- collapse_probes_to_genes(
  expr = expr,
  feature = feature,
  symbol_col = "Gene symbol"
)

saveRDS(expr_gene, expr_gene_out)

## Make sure expression columns match phenotype GSM order
if (!all(colnames(expr_gene) == pheno_clean$gsm)) {
  message("Expression and phenotype are not aligned by order; matching by GSM.")
  
  common <- intersect(colnames(expr_gene), pheno_clean$gsm)
  
  expr_gene <- expr_gene[, common, drop = FALSE]
  
  pheno_clean <- pheno_clean %>%
    filter(gsm %in% common) %>%
    arrange(match(gsm, colnames(expr_gene)))
}

stopifnot(all(colnames(expr_gene) == pheno_clean$gsm))

## -----------------------------
## Filter ER-positive samples
## -----------------------------

keep_erpos <- !is.na(pheno_clean$er_status) & pheno_clean$er_status == "ER+"

pheno_erpos <- pheno_clean[keep_erpos, , drop = FALSE]
expr_erpos <- expr_gene[, pheno_erpos$gsm, drop = FALSE]

pheno_erpos <- pheno_erpos %>%
  arrange(match(gsm, colnames(expr_erpos)))

stopifnot(all(colnames(expr_erpos) == pheno_erpos$gsm))

message("Analysis cohort after ER+ filter:")
message("N ER+ samples: ", nrow(pheno_erpos))
message("Distant relapse events: ", sum(pheno_erpos$distant_relapse == 1, na.rm = TRUE))
print(table(pheno_erpos$distant_relapse, useNA = "ifany"))
print(table(pheno_erpos$nodal_status, useNA = "ifany"))
print(summary(pheno_erpos$event_time_years))

## -----------------------------
## Score signatures
## -----------------------------

message("Z-scoring expression within ER+ cohort.")
zmat <- zscore_by_gene(expr_erpos)

score_tbl <- tibble(gsm = colnames(zmat))
overlap_list <- list()

for (sig_name in names(sig_list)) {
  message("Scoring: ", sig_name)
  
  res <- score_signature(
    zmat = zmat,
    genes = sig_list[[sig_name]],
    min_genes = 5
  )
  
  score_tbl[[sig_name]] <- res$score
  
  overlap_list[[sig_name]] <- tibble(
    signature = sig_name,
    n_input_genes = length(unique(clean_gene(sig_list[[sig_name]]))),
    n_overlap_genes = length(res$overlap),
    n_missing_genes = length(res$missing),
    overlap_genes = paste(res$overlap, collapse = ";"),
    missing_genes = paste(res$missing, collapse = ";")
  )
}

## Delta score
if (all(c("LOW_ESR1_UP", "LOW_ESR1_DW") %in% colnames(score_tbl))) {
  score_tbl <- score_tbl %>%
    mutate(
      LOW_ESR1_DELTA = LOW_ESR1_UP - LOW_ESR1_DW
    )
} else {
  stop("LOW_ESR1_UP and LOW_ESR1_DW scores were not both generated.")
}

overlap_tbl <- bind_rows(overlap_list)
write_csv(overlap_tbl, overlap_out)

message("Signature overlap:")
print(overlap_tbl %>% select(signature, n_input_genes, n_overlap_genes, n_missing_genes))

score_pheno <- pheno_erpos %>%
  left_join(score_tbl, by = "gsm") %>%
  filter(
    !is.na(distant_relapse),
    !is.na(event_time_years)
  )

write_csv(score_pheno, score_out)

message("Saved score table: ", score_out)
message("Analysis samples: ", nrow(score_pheno))
message("Distant relapse events: ", sum(score_pheno$distant_relapse == 1, na.rm = TRUE))

## -----------------------------
## Cox models for UP, DW, and DELTA
## -----------------------------

score_cols <- c(
  "LOW_ESR1_UP",
  "LOW_ESR1_DW",
  "LOW_ESR1_DELTA"
)

cox_results <- bind_rows(lapply(score_cols, function(sc) {
  bind_rows(
    run_cox(score_pheno, sc, covars = character(0)),
    run_cox(score_pheno, sc, covars = c("nodal_status"))
  )
}))

cox_results <- cox_results %>%
  mutate(FDR = p.adjust(p, method = "BH")) %>%
  arrange(model, p)

write_csv(cox_results, cox_out)

sink(cox_txt)
cat("GSE17705 low_ESR1 UP/DW/delta distant relapse analysis\n")
cat("=======================================================\n\n")
cat("N samples: ", nrow(score_pheno), "\n")
cat("Distant relapse events: ", sum(score_pheno$distant_relapse == 1, na.rm = TRUE), "\n\n")

cat("Signature overlap:\n")
print(overlap_tbl %>% select(signature, n_input_genes, n_overlap_genes, n_missing_genes))

cat("\n\nCox results:\n")
print(cox_results)
sink()

message("Cox results:")
print(cox_results)

## -----------------------------
## KM plots and asymmetric scans
## -----------------------------

all_scan_results <- list()

for (sc in score_cols) {
  message("Generating KM plots and cutoff scan for: ", sc)
  
  plot_km_split(
    df = score_pheno,
    score_col = sc,
    split_type = "median",
    out_dir = plot_dir
  )
  
  plot_km_split(
    df = score_pheno,
    score_col = sc,
    split_type = "top_bottom10",
    out_dir = plot_dir
  )
  
  scan_res <- run_asymmetric_km_scan(
    df = score_pheno,
    score_col = sc,
    bottom_grid = seq(0.10, 0.50, by = 0.05),
    top_grid = seq(0.05, 0.50, by = 0.05),
    min_events_per_group = 5,
    min_samples_per_group = 20
  )
  
  all_scan_results[[sc]] <- scan_res
  
  write_csv(
    scan_res,
    file.path(cutoff_scan_dir, paste0(sc, "_asymmetric_cutoff_scan.csv"))
  )
  
  best <- scan_res %>%
    filter(note == "ok", !is.na(logrank_p)) %>%
    arrange(logrank_p) %>%
    slice_head(n = 1)
  
  message("Best asymmetric cutoff for ", sc, ":")
  print(best)
  
  if (nrow(best) == 1) {
    plot_best_asymmetric_km(
      df = score_pheno,
      score_col = sc,
      bottom_pct = best$bottom_pct,
      top_pct = best$top_pct,
      out_dir = cutoff_scan_dir
    )
  }
}

all_scan_results <- bind_rows(all_scan_results)

write_csv(
  all_scan_results,
  file.path(cutoff_scan_dir, "all_LOW_ESR1_UP_DW_DELTA_asymmetric_cutoff_scan.csv")
)

message("Done.")
