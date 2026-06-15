

##### Run UNCover example





res_PDX1_Control_low_ESR1_vs_high_ESR1 <- uncover_run(
  object = PDX1_Control,
  group_col = "seurat_clusters",
  group1 = c('0','6'),   # '11', '5',
  group2 = c('16', '2', '9', '13'),
  rna_assay = "SCT",
  peak_assay = "peaks",
  candidate_r_threshold = 0.30,
  ml_use_metacells = TRUE,
  ml_cells_per_metacell = 50,   # TF->peak ML layer
  ml_metacell_reduction = "umap",
  ml_metacell_dims = 1:2,
  n_perm_fast = 100,
  n_random_pairs = 100,
  n_perm_random = 20,
  peak_gene_p_method = "zscore",
  peak_gene_fdr_cutoff = 0.05,
  network_fdr_cor_max = 0.05,
  tf_alpha = 0.05,
  tf_pct_threshold = 0.050,
  metacell_method = "fixed_size_umap",
  cells_per_metacell = 50,
  ml_mode = "slurm",
  submit = TRUE,
  collect_results = TRUE,
  wait_for_results = TRUE,
  conda_env = "myenv",
  conda_sh = "/path/to/miniconda3/etc/profile.d/conda.sh",
  python = "/path/to/miniconda3/envs/myenv/bin/python",
  ml_output_dir = "/path/to/res_PDX1_Control_low_ESR1_vs_high_ESR1",
  sbatch_dir = "/path/to/res_PDX1_Control_low_ESR1_vs_high_ESR1",
  partition = "allnodes",
  cpus_per_task = 4,
  mem_per_cpu = "2G",
  time = "08:00:00",
  peaks_per_job = 10,
  tf_peak_I_min = 0.025,
  tf_peak_low   = 0.3,
  tf_peak_high  = 0.8,
  chip_catalog_bed = "/path/to/TF_peak_genes_pipeline/tmp/remap_tf_4col.bed",
  run_chip_intersect = TRUE,
  intersect_mode = "local",
  tf_peak_require_chip = TRUE,
  checkpoint_dir = "/path/to/res_PDX1_Control_low_ESR1_vs_high_ESR1/checkpoints",
  checkpoint_file = "/path/to/res_PDX1_Control_low_ESR1_vs_high_ESR1/checkpoints/UNCover_pipeline_state.rds",
  resume_from_checkpoint = TRUE,
  save_final_rds = TRUE,
  final_rds_file = "/path/to/res_PDX1_Control_low_ESR1_vs_high_ESR1/UNCover_result_final.rds",
  run_tf_dual_profile_analysis = TRUE,
  tf_dual_profile_top_n = 20,
  tf_dual_profile_rank_by = "score",
  tf_dual_profile_min_peaks = 1,
  tf_dual_profile_min_triplets = 1,
  tf_dual_profile_use_weighted = TRUE,
  save_network_html = TRUE
)





















extract_gene_target_subnetwork <- function(
    result,
    target_gene = "ESR1",
    min_I = 0,
    max_fdr_pg = 0.05,
    min_abs_cor_pg = 0,
    out_dir = NULL,
    prefix = NULL,
    make_html = TRUE
) {
  suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(readr)
    library(visNetwork)
    library(scales)
  })
  
  stopifnot(!is.null(result))
  stopifnot(all(c("edges_tf_peak", "edges_peak_gene", "nodes") %in% names(result)))
  
  edges_tf_peak <- result$edges_tf_peak
  edges_peak_gene <- result$edges_peak_gene
  nodes <- result$nodes
  
  num <- function(x) suppressWarnings(as.numeric(x))
  
  # 1) Keep only Peak->target_gene edges
  epg <- edges_peak_gene %>%
    transmute(
      peak   = .data$from,
      gene   = .data$to,
      cor_pg = num(.data$cor_pg),
      fdr_pg = num(.data$fdr_pg)
    ) %>%
    dplyr::filter(
      .data$gene == target_gene,
      is.finite(.data$cor_pg),
      abs(.data$cor_pg) >= min_abs_cor_pg
    )
  
  if (!is.na(max_fdr_pg)) {
    epg <- epg %>%
      dplyr::filter(is.finite(.data$fdr_pg), .data$fdr_pg <= max_fdr_pg)
  }
  
  if (nrow(epg) == 0) {
    stop("No Peak→", target_gene, " edges after filters.")
  }
  
  peaks_target <- unique(epg$peak)
  
  # 2) Keep TF->peak edges only for peaks linked to target_gene
  etf <- edges_tf_peak %>%
    transmute(
      tf       = .data$from,
      peak     = .data$to,
      I        = num(.data$I),
      chip_hits = suppressWarnings(as.integer(.data$chip_hits))
    ) %>%
    dplyr::filter(
      .data$peak %in% peaks_target,
      is.finite(.data$I),
      .data$I >= min_I
    )
  
  if (nrow(etf) == 0) {
    stop("No TF→Peak edges for ", target_gene, " peaks after I filtering.")
  }
  
  # 3) Build triplets
  triplets_sub <- etf %>%
    inner_join(epg, by = "peak", relationship = "many-to-many") %>%
    arrange(desc(.data$I), desc(abs(.data$cor_pg)), .data$tf, .data$peak)
  
  # 4) Node subset
  sub_names <- unique(c(triplets_sub$tf, triplets_sub$peak, triplets_sub$gene))
  
  nodes_sub <- nodes %>%
    dplyr::filter(.data$name %in% sub_names) %>%
    mutate(
      tf_pct_change = num(.data$tf_pct_change),
      peak_log2FC   = num(.data$peak_log2FC),
      gene_log2FC   = num(.data$gene_log2FC)
    ) %>%
    distinct(.data$name, .keep_all = TRUE)
  
  # 5) Edge tables in package-like format
  edges_tf_peak_sub <- triplets_sub %>%
    transmute(
      from = .data$tf,
      to = .data$peak,
      edge_type = "TF_to_Peak",
      I = .data$I,
      chip_hits = .data$chip_hits,
      cor_pg = NA_real_,
      fdr_pg = NA_real_
    ) %>%
    distinct()
  
  edges_peak_gene_sub <- triplets_sub %>%
    transmute(
      from = .data$peak,
      to = .data$gene,
      edge_type = "Peak_to_Gene",
      I = NA_real_,
      chip_hits = NA_integer_,
      cor_pg = .data$cor_pg,
      fdr_pg = .data$fdr_pg
    ) %>%
    distinct()
  
  edges_sub <- bind_rows(edges_tf_peak_sub, edges_peak_gene_sub)
  
  # 6) Optional HTML
  html_file <- NULL
  if (isTRUE(make_html)) {
    sym_div_domain_local <- function(x, p = 0.99, min_span = 1) {
      x <- x[is.finite(x)]
      if (!length(x)) return(c(-1, 1))
      a <- stats::quantile(abs(x), p, na.rm = TRUE)
      a <- max(a, min_span / 2)
      c(-a, a)
    }
    
    tf_pal <- scales::col_numeric(
      c("#2166AC", "#FFFFFF", "#B2182B"),
      sym_div_domain_local(nodes_sub$tf_pct_change, 0.98, 0.1),
      na.color = "#BDBDBD"
    )
    pk_pal <- scales::col_numeric(
      c("#313695", "#FFFFFF", "#A50026"),
      sym_div_domain_local(nodes_sub$peak_log2FC, 0.98, 0.5),
      na.color = "#BDBDBD"
    )
    gn_pal <- scales::col_numeric(
      c("#2C7BB6", "#FFFFFF", "#D7191C"),
      sym_div_domain_local(nodes_sub$gene_log2FC, 0.98, 0.5),
      na.color = "#BDBDBD"
    )
    
    I_max <- max(as.numeric(edges_tf_peak_sub$I), na.rm = TRUE)
    if (!is.finite(I_max) || I_max <= 0) I_max <- 1
    I_pal <- scales::col_numeric(
      c("#DEEBF7", "#3182BD"),
      c(0, I_max),
      na.color = "#9ECAE1"
    )
    
    r_pal <- scales::col_numeric(
      c("#313695", "#FFFFFF", "#A50026"),
      sym_div_domain_local(as.numeric(edges_peak_gene_sub$cor_pg), 1.0, 0.2),
      na.color = "#BDBDBD"
    )
    
    shape_map <- c(TF = "triangle", Peak = "diamond", Gene = "ellipse")
    
    nodes_v <- nodes_sub %>%
      mutate(
        id    = paste0(.data$type, ":", .data$name),
        label = .data$name,
        group = .data$type,
        shape = shape_map[.data$type],
        color.background = case_when(
          .data$type == "TF"   ~ tf_pal(.data$tf_pct_change),
          .data$type == "Peak" ~ pk_pal(.data$peak_log2FC),
          .data$type == "Gene" ~ gn_pal(.data$gene_log2FC),
          TRUE ~ "#BDBDBD"
        ),
        color.border = "#333333"
      ) %>%
      transmute(
        id, label, group, shape,
        color.background,
        color.border
      )
    
    edges_v <- edges_sub %>%
      mutate(
        from_id = paste0(ifelse(.data$edge_type == "TF_to_Peak", "TF", "Peak"), ":", .data$from),
        to_id   = paste0(ifelse(.data$edge_type == "TF_to_Peak", "Peak", "Gene"), ":", .data$to),
        value   = case_when(
          .data$edge_type == "TF_to_Peak" ~ pmax(1, 12 * .data$I),
          TRUE ~ pmax(1, 12 * abs(.data$cor_pg))
        ),
        color   = case_when(
          .data$edge_type == "TF_to_Peak" ~ I_pal(.data$I),
          TRUE ~ r_pal(.data$cor_pg)
        ),
        title   = case_when(
          .data$edge_type == "TF_to_Peak" ~ paste0(
            .data$from, " → ", .data$to,
            " | I=", round(.data$I, 3),
            " | ChIP hits=", .data$chip_hits
          ),
          TRUE ~ paste0(
            .data$from, " → ", .data$to,
            " | r=", round(.data$cor_pg, 3),
            ifelse(is.finite(.data$fdr_pg), paste0(" | FDR=", signif(.data$fdr_pg, 3)), "")
          )
        ),
        arrows = "to"
      ) %>%
      transmute(from = .data$from_id, to = .data$to_id, value, color, title, arrows)
    
    g_vis <- visNetwork(nodes_v, edges_v, width = "100%", height = "800px") %>%
      visLegend(
        addNodes = data.frame(
          label = c("TF", "Peak", "Gene"),
          shape = c("triangle", "diamond", "ellipse"),
          color = c("#BDBDBD", "#BDBDBD", "#BDBDBD"),
          stringsAsFactors = FALSE
        ),
        useGroups = FALSE
      ) %>%
      visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
      visPhysics(stabilization = TRUE)
    
    if (!is.null(out_dir)) {
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
      if (is.null(prefix)) prefix <- paste0("subnetwork_", target_gene)
      html_file <- file.path(out_dir, paste0(prefix, ".html"))
      visSave(g_vis, html_file)
    }
  }
  
  # 7) Optional CSV export
  if (!is.null(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    if (is.null(prefix)) prefix <- paste0("subnetwork_", target_gene)
    
    readr::write_csv(triplets_sub,       file.path(out_dir, paste0(prefix, "_triplets.csv")))
    readr::write_csv(edges_tf_peak_sub,  file.path(out_dir, paste0(prefix, "_edges_tf_peak.csv")))
    readr::write_csv(edges_peak_gene_sub,file.path(out_dir, paste0(prefix, "_edges_peak_gene.csv")))
    readr::write_csv(nodes_sub,          file.path(out_dir, paste0(prefix, "_nodes.csv")))
  }
  
  list(
    target_gene = target_gene,
    peaks = sort(unique(epg$peak)),
    tfs = sort(unique(etf$tf)),
    triplets = triplets_sub,
    edges_tf_peak = edges_tf_peak_sub,
    edges_peak_gene = edges_peak_gene_sub,
    edges = edges_sub,
    nodes = nodes_sub,
    html_file = html_file
  )
}


esr1_net <- extract_gene_target_subnetwork(
  result = res_PDX1_Control_low_ESR1_vs_high_ESR1,
  target_gene = "ESR1",
  min_I = 0,
  max_fdr_pg = 0.05,
  min_abs_cor_pg = 0,
  out_dir = "/path/to/UNCover_PDX1_lowESR_vs_high_esr1/final_outputs",
  prefix = "ESR1_target_only",
  make_html = TRUE
)


























plot_save_target_subnetwork2 <- function(
    subnet,
    result,
    out_dir,
    prefix = "ESR1_prechip_subnetwork",
    target_gene = "ESR1"
) {
  suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
    library(visNetwork)
    library(scales)
    library(tibble)
  })
  
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  trip <- subnet$triplets
  if (nrow(trip) == 0) stop("subnet$triplets is empty.")
  
  num <- function(x) suppressWarnings(as.numeric(x))
  
  norm_peak_local <- function(x) {
    x <- as.character(x)
    x <- gsub(":", "-", x)
    x <- gsub("_", "-", x)
    x
  }
  
  pick_first_existing_numeric <- function(df, cols) {
    out <- rep(NA_real_, nrow(df))
    for (cc in cols) {
      if (cc %in% colnames(df)) {
        vals <- suppressWarnings(as.numeric(df[[cc]]))
        fill <- is.na(out) & is.finite(vals)
        out[fill] <- vals[fill]
      }
    }
    out
  }
  
  # ---- TF metadata ----
  tf_df <- as.data.frame(result$tf_activity)
  
  tf_col <- if ("tf" %in% colnames(tf_df)) {
    "tf"
  } else if ("signature" %in% colnames(tf_df)) {
    "signature"
  } else {
    stop("Could not find TF name column in result$tf_activity.")
  }
  
  tf_meta <- tibble::tibble(
    tf_raw = as.character(tf_df[[tf_col]]),
    tf = clean_signature(tf_df[[tf_col]]),
    tf_norm = normalize_tf_symbols(clean_signature(tf_df[[tf_col]])),
    tf_pct_change = pick_first_existing_numeric(
      tf_df,
      c("tf_pct_change", "pct_change", "pct_delta_med", "pct_delta_mean")
    ),
    tf_p_adj = pick_first_existing_numeric(
      tf_df,
      c("tf_p_adj", "p_adj", "fdr")
    )
  ) %>%
    dplyr::filter(!is.na(.data$tf_norm), .data$tf_norm != "") %>%
    dplyr::distinct(.data$tf_norm, .keep_all = TRUE)
  
  # ---- Peak metadata ----
  da_df <- as.data.frame(result$da_peaks)
  
  peak_col <- if ("peak" %in% colnames(da_df)) {
    "peak"
  } else if ("Peak" %in% colnames(da_df)) {
    "Peak"
  } else {
    stop("Could not find peak column in result$da_peaks.")
  }
  
  peak_meta <- tibble::tibble(
    peak = norm_peak_local(da_df[[peak_col]]),
    peak_log2FC = pick_first_existing_numeric(
      da_df,
      c("peak_log2FC", "avg_log2FC", "log2FC", "avg_log2FC_harmonized", "log2FC_harmonized")
    )
  ) %>%
    dplyr::filter(is.finite(.data$peak_log2FC)) %>%
    dplyr::distinct(.data$peak, .keep_all = TRUE)
  
  # ---- Gene metadata ----
  de_df <- as.data.frame(result$de_genes)
  
  gene_col <- if ("gene" %in% colnames(de_df)) {
    "gene"
  } else if ("Gene" %in% colnames(de_df)) {
    "Gene"
  } else {
    stop("Could not find gene column in result$de_genes.")
  }
  
  gene_meta <- tibble::tibble(
    gene = as.character(de_df[[gene_col]]),
    gene_log2FC = pick_first_existing_numeric(
      de_df,
      c("gene_log2FC", "avg_log2FC", "log2FC", "avg_log2FC_harmonized", "log2FC_harmonized")
    )
  ) %>%
    dplyr::filter(is.finite(.data$gene_log2FC)) %>%
    dplyr::distinct(.data$gene, .keep_all = TRUE)
  
  # ---- Enrich triplets ----
  trip2 <- trip %>%
    dplyr::mutate(peak = norm_peak_local(.data$peak)) %>%
    dplyr::left_join(
      tf_meta %>%
        dplyr::select(.data$tf_norm, .data$tf_pct_change, .data$tf_p_adj),
      by = "tf_norm"
    ) %>%
    dplyr::left_join(peak_meta, by = "peak") %>%
    dplyr::left_join(gene_meta, by = "gene") %>%
    dplyr::distinct()
  
  # ---- Edge tables ----
  edges_tf_peak <- trip2 %>%
    dplyr::transmute(
      from = .data$tf,
      to = .data$peak,
      edge_type = "TF_to_Peak",
      I = .data$I,
      I_weighted = .data$I_weighted,
      M = .data$M,
      model_correlation = .data$model_correlation,
      chip_hits = .data$chip_hits,
      cor_pg = NA_real_,
      fdr_pg = NA_real_
    ) %>%
    dplyr::distinct()
  
  edges_peak_gene <- trip2 %>%
    dplyr::transmute(
      from = .data$peak,
      to = .data$gene,
      edge_type = "Peak_to_Gene",
      I = NA_real_,
      I_weighted = NA_real_,
      M = NA_real_,
      model_correlation = NA_real_,
      chip_hits = NA_integer_,
      cor_pg = .data$cor_pg,
      fdr_pg = .data$fdr_pg
    ) %>%
    dplyr::distinct()
  
  # ---- Node table ----
  nodes_tf <- trip2 %>%
    dplyr::transmute(
      name = .data$tf,
      type = "TF",
      tf_pct_change = .data$tf_pct_change,
      tf_p_adj = .data$tf_p_adj,
      peak_log2FC = NA_real_,
      gene_log2FC = NA_real_
    ) %>%
    dplyr::distinct(.data$name, .keep_all = TRUE)
  
  nodes_peak <- trip2 %>%
    dplyr::transmute(
      name = .data$peak,
      type = "Peak",
      tf_pct_change = NA_real_,
      tf_p_adj = NA_real_,
      peak_log2FC = .data$peak_log2FC,
      gene_log2FC = NA_real_
    ) %>%
    dplyr::distinct(.data$name, .keep_all = TRUE)
  
  nodes_gene <- trip2 %>%
    dplyr::transmute(
      name = .data$gene,
      type = "Gene",
      tf_pct_change = NA_real_,
      tf_p_adj = NA_real_,
      peak_log2FC = NA_real_,
      gene_log2FC = .data$gene_log2FC
    ) %>%
    dplyr::distinct(.data$name, .keep_all = TRUE)
  
  nodes_sub <- dplyr::bind_rows(nodes_tf, nodes_peak, nodes_gene) %>%
    dplyr::distinct(.data$name, .data$type, .keep_all = TRUE)
  
  sym_div_domain <- function(x, p = 1.0, min_span = 0.1) {
    x <- x[is.finite(x)]
    if (!length(x)) return(c(-1, 1))
    a <- as.numeric(stats::quantile(abs(x), p, na.rm = TRUE))
    a <- max(a, min_span / 2)
    c(-a, a)
  }
  
  clamp_to_domain <- function(x, domain) {
    x <- as.numeric(x)
    x[x < domain[1]] <- domain[1]
    x[x > domain[2]] <- domain[2]
    x
  }
  
  tf_dom <- sym_div_domain(nodes_sub$tf_pct_change, p = 1.0, min_span = 0.05)
  pk_dom <- sym_div_domain(nodes_sub$peak_log2FC, p = 1.0, min_span = 0.25)
  gn_dom <- sym_div_domain(nodes_sub$gene_log2FC, p = 1.0, min_span = 0.25)
  I_dom <- c(0, max(edges_tf_peak$I_weighted, na.rm = TRUE))
  if (!is.finite(I_dom[2]) || I_dom[2] <= 0) I_dom[2] <- 1
  r_dom <- sym_div_domain(edges_peak_gene$cor_pg, p = 1.0, min_span = 0.10)
  
  tf_pal <- scales::col_numeric(c("#2166AC", "#FFFFFF", "#B2182B"), domain = tf_dom, na.color = "#BDBDBD")
  pk_pal <- scales::col_numeric(c("#313695", "#FFFFFF", "#A50026"), domain = pk_dom, na.color = "#BDBDBD")
  gn_pal <- scales::col_numeric(c("#2C7BB6", "#FFFFFF", "#D7191C"), domain = gn_dom, na.color = "#BDBDBD")
  I_pal  <- scales::col_numeric(c("#DEEBF7", "#3182BD"), domain = I_dom, na.color = "#9ECAE1")
  r_pal  <- scales::col_numeric(c("#313695", "#FFFFFF", "#A50026"), domain = r_dom, na.color = "#BDBDBD")
  
  shape_map <- c(TF = "triangle", Peak = "diamond", Gene = "ellipse")
  
  nodes_v <- nodes_sub %>%
    dplyr::mutate(
      id = paste0(.data$type, ":", .data$name),
      label = .data$name,
      group = .data$type,
      shape = unname(shape_map[.data$type]),
      color.background = dplyr::case_when(
        .data$type == "TF"   ~ tf_pal(clamp_to_domain(.data$tf_pct_change, tf_dom)),
        .data$type == "Peak" ~ pk_pal(clamp_to_domain(.data$peak_log2FC, pk_dom)),
        .data$type == "Gene" ~ gn_pal(clamp_to_domain(.data$gene_log2FC, gn_dom)),
        TRUE ~ "#BDBDBD"
      ),
      color.border = "#333333",
      title = dplyr::case_when(
        .data$type == "TF" ~ paste0(
          .data$name,
          "<br>TF %Δ = ", ifelse(is.finite(.data$tf_pct_change), sprintf("%+.1f%%", 100 * .data$tf_pct_change), "NA"),
          "<br>TF FDR = ", ifelse(is.finite(.data$tf_p_adj), signif(.data$tf_p_adj, 3), "NA")
        ),
        .data$type == "Peak" ~ paste0(
          .data$name,
          "<br>Peak log2FC = ", ifelse(is.finite(.data$peak_log2FC), sprintf("%+.2f", .data$peak_log2FC), "NA")
        ),
        .data$type == "Gene" ~ paste0(
          .data$name,
          "<br>Gene log2FC = ", ifelse(is.finite(.data$gene_log2FC), sprintf("%+.2f", .data$gene_log2FC), "NA")
        ),
        TRUE ~ .data$name
      )
    ) %>%
    dplyr::select(id, label, group, shape, title, color.background, color.border) %>%
    dplyr::distinct()
  
  edges_v_tf <- edges_tf_peak %>%
    dplyr::mutate(
      from = paste0("TF:", .data$from),
      to   = paste0("Peak:", .data$to),
      value = pmax(1, 18 * .data$I_weighted),
      color = I_pal(clamp_to_domain(.data$I_weighted, I_dom)),
      title = paste0(
        .data$from, " → ", .data$to,
        "<br>I = ", signif(.data$I, 3),
        "<br>I_weighted = ", signif(.data$I_weighted, 3),
        "<br>M = ", signif(.data$M, 3),
        "<br>model_correlation = ", signif(.data$model_correlation, 3),
        "<br>chip_hits = ", .data$chip_hits
      ),
      arrows = "to"
    ) %>%
    dplyr::select(from, to, value, color, title, arrows)
  
  edges_v_pg <- edges_peak_gene %>%
    dplyr::mutate(
      from = paste0("Peak:", .data$from),
      to   = paste0("Gene:", .data$to),
      value = pmax(1, 12 * abs(.data$cor_pg)),
      color = r_pal(clamp_to_domain(.data$cor_pg, r_dom)),
      title = paste0(
        .data$from, " → ", .data$to,
        "<br>cor_pg = ", signif(.data$cor_pg, 3),
        "<br>fdr_pg = ", signif(.data$fdr_pg, 3)
      ),
      arrows = "to"
    ) %>%
    dplyr::select(from, to, value, color, title, arrows)
  
  edges_v <- dplyr::bind_rows(edges_v_tf, edges_v_pg)
  
  net <- visNetwork::visNetwork(nodes_v, edges_v, width = "100%", height = "850px") %>%
    visNetwork::visLegend(
      addNodes = data.frame(
        label = c("TF", "Peak", "Gene"),
        shape = c("triangle", "diamond", "ellipse"),
        color = c("#BDBDBD", "#BDBDBD", "#BDBDBD"),
        stringsAsFactors = FALSE
      ),
      useGroups = FALSE
    ) %>%
    visNetwork::visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
    visNetwork::visPhysics(stabilization = TRUE) %>%
    visNetwork::visInteraction(navigationButtons = TRUE)
  
  html_file <- file.path(out_dir, paste0(prefix, ".html"))
  visNetwork::visSave(net, html_file)
  
  readr::write_csv(trip2,           file.path(out_dir, paste0(prefix, "_triplets.csv")))
  readr::write_csv(edges_tf_peak,   file.path(out_dir, paste0(prefix, "_edges_tf_peak.csv")))
  readr::write_csv(edges_peak_gene, file.path(out_dir, paste0(prefix, "_edges_peak_gene.csv")))
  readr::write_csv(nodes_sub,       file.path(out_dir, paste0(prefix, "_nodes.csv")))
  
  message("Saved HTML: ", html_file)
  
  invisible(list(
    net = net,
    html_file = html_file,
    triplets = trip2,
    edges_tf_peak = edges_tf_peak,
    edges_peak_gene = edges_peak_gene,
    nodes = nodes_sub
  ))
}


esr1_plot <- plot_save_target_subnetwork2(
  subnet = esr1_prechip,
  result = res_core_chip,
  out_dir = "/path/to/UNCover_PDX1_non_prolif_vs_high_esr1_2/final_outputs",
  prefix = "ESR1_prechip_subnetwork_v2",
  target_gene = "ESR1"
)


esr1_plot



