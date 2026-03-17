# Build advanced statistical outputs and figures from the final cohort table.
# This script adds calibration-aware model comparison, multiplicity control,
# network topology analysis, and sex-stratified effect-size reporting.

library(tidyverse)
library(broom)
library(rstatix)
library(pROC)
library(igraph)

has_ggraph <- requireNamespace("ggraph", quietly = TRUE)
has_effectsize <- requireNamespace("effectsize", quietly = TRUE)
has_introdataviz <- requireNamespace("introdataviz", quietly = TRUE)

# Global Youreka visual system: one palette + one baseline theme for all figures.
youreka_colors <- c(
  blue = "#1f4e79",
  teal = "#1f7a8c",
  orange = "#f59e0b",
  red = "#b91c1c",
  slate = "#334155",
  blue_light = "#c7d7ea",
  orange_light = "#f3dd86",
  mint_light = "#a8e6c8"
)

theme_set(
  theme_gray(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, color = youreka_colors[["blue"]], face = "bold", size = 16),
      plot.subtitle = element_text(hjust = 0.5, color = youreka_colors[["blue"]], face = "bold", size = 11),
      axis.title = element_text(color = "#2f2f2f", face = "bold"),
      legend.title = element_text(face = "bold")
    )
)

out_dir <- "06_analysis"
fig_dir <- file.path(out_dir, "figures")
tab_dir <- file.path(out_dir, "tables")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)

df <- read_csv("05_final/cohort_analysis.csv", show_col_types = FALSE)

# Basic coercions
num_cols <- c("AGE", "n_codrugs", "total_acb_with_dph", "total_acb_codrugs_only", "max_severity", "cardiac_any", "pre_post_warning")
df <- df %>% mutate(across(all_of(intersect(num_cols, names(df))), as.numeric))

# Table 1 descriptive summaries
cont <- c("AGE", "total_acb_with_dph", "total_acb_codrugs_only", "n_codrugs", "max_severity")
cont_tbl <- map_dfr(cont, function(v) {
  x <- df[[v]]
  tibble(
    variable = v,
    n = sum(!is.na(x)),
    mean = mean(x, na.rm = TRUE),
    sd = sd(x, na.rm = TRUE),
    median = median(x, na.rm = TRUE),
    iqr = IQR(x, na.rm = TRUE)
  )
})
write_csv(cont_tbl, file.path(tab_dir, "descriptive_continuous.csv"))

cat_tbl <- bind_rows(
  # Age-group distribution
  df %>% count(age_group, name = "count") %>% mutate(variable = "age_group", level = as.character(age_group), percent = 100 * count / sum(count)) %>% select(variable, level, count, percent),
  # Pre/post warning split
  df %>% count(pre_post_warning, name = "count") %>% mutate(variable = "pre_post_warning", level = as.character(pre_post_warning), percent = 100 * count / sum(count)) %>% select(variable, level, count, percent),
  # Cardiac outcome prevalence
  df %>% count(cardiac_any, name = "count") %>% mutate(variable = "cardiac_any", level = as.character(cardiac_any), percent = 100 * count / sum(count)) %>% select(variable, level, count, percent)
)
write_csv(cat_tbl, file.path(tab_dir, "descriptive_categorical.csv"))

# Normality checks
shap_tbl <- map_dfr(cont, function(v) {
  x <- df[[v]] %>% na.omit()
  if (length(x) >= 3) {
    s <- shapiro.test(x)
    tibble(variable = v, W = as.numeric(s$statistic), p_value = s$p.value)
  } else {
    tibble(variable = v, W = NA_real_, p_value = NA_real_)
  }
})
write_csv(shap_tbl, file.path(tab_dir, "normality_checks.csv"))

# Core analyses
# A) ACB (co-drugs only) vs cardiac outcome
mw <- wilcox.test(total_acb_codrugs_only ~ cardiac_any, data = df)
# B) Co-medication ACB across age groups
kw <- kruskal.test(total_acb_codrugs_only ~ age_group, data = df)
# C) Monotonic association between ACB and severity
sp <- cor.test(df$total_acb_codrugs_only, df$max_severity, method = "spearman", exact = FALSE)

core <- tibble(
  analysis = c("mann_whitney_acb_vs_cardiac", "kruskal_acb_vs_age_group", "spearman_acb_vs_severity"),
  statistic = c(as.numeric(mw$statistic), as.numeric(kw$statistic), as.numeric(sp$estimate)),
  p_value = c(mw$p.value, kw$p.value, sp$p.value)
)
write_csv(core, file.path(tab_dir, "nonparametric_results.csv"))

# Logistic models (robust to subsets with single-level predictors)
fit_logit_safe <- function(data_in, candidate_terms, model_label) {
  data_mod <- data_in %>%
    mutate(
      cardiac_any = as.numeric(cardiac_any),
      age_group = as.factor(age_group)
    ) %>%
    filter(!is.na(cardiac_any))

  # Outcome must have both classes for logistic regression.
  if (length(unique(data_mod$cardiac_any)) < 2) {
    return(tibble(model = model_label, term = NA_character_, estimate = NA_real_, conf.low = NA_real_, conf.high = NA_real_, p.value = NA_real_, note = "Skipped: cardiac_any has <2 classes"))
  }

  keep_terms <- c()
  for (term in candidate_terms) {
    if (!(term %in% names(data_mod))) {
      next
    }
    v <- data_mod[[term]]
    uniq_n <- length(unique(v[!is.na(v)]))
    if (uniq_n >= 2) {
      keep_terms <- c(keep_terms, term)
    }
  }

  if (length(keep_terms) == 0) {
    return(tibble(model = model_label, term = NA_character_, estimate = NA_real_, conf.low = NA_real_, conf.high = NA_real_, p.value = NA_real_, note = "Skipped: no variable with >=2 unique values"))
  }

  fml <- as.formula(paste("cardiac_any ~", paste(keep_terms, collapse = " + ")))
  fit <- glm(fml, data = data_mod, family = binomial)
  out <- tidy(fit, conf.int = TRUE, exponentiate = TRUE)
  out$model <- model_label
  out$note <- paste("Used terms:", paste(keep_terms, collapse = ", "))
  out
}

full_tbl <- fit_logit_safe(
  df,
  c("total_acb_codrugs_only", "age_group", "pre_post_warning", "n_codrugs"),
  "full"
)
write_csv(full_tbl, file.path(tab_dir, "logistic_model_full.csv"))

pre_tbl <- fit_logit_safe(
  filter(df, pre_post_warning == 0),
  c("total_acb_codrugs_only", "age_group", "n_codrugs"),
  "pre_warning"
)
write_csv(pre_tbl, file.path(tab_dir, "logistic_model_pre_warning.csv"))

post_tbl <- fit_logit_safe(
  filter(df, pre_post_warning == 1),
  c("total_acb_codrugs_only", "age_group", "n_codrugs"),
  "post_warning"
)
write_csv(post_tbl, file.path(tab_dir, "logistic_model_post_warning.csv"))

# Phase 3: DeLong ROC significance test using strict CV probabilities from Python
cv_pred_path <- file.path(tab_dir, "model_predictions_cv.csv")
delong_p <- NA_real_
if (file.exists(cv_pred_path)) {
  pred_cv <- read_csv(cv_pred_path, show_col_types = FALSE)
  pred_cv <- pred_cv %>% mutate(cardiac_any = as.numeric(cardiac_any))
  rf_roc <- roc(pred_cv$cardiac_any, pred_cv$rf_prob, levels = c(0, 1), direction = "<", quiet = TRUE)
  lr_roc <- roc(pred_cv$cardiac_any, pred_cv$lr_prob, levels = c(0, 1), direction = "<", quiet = TRUE)
  delong <- roc.test(rf_roc, lr_roc, method = "delong", paired = TRUE)
  delong_p <- as.numeric(delong$p.value)

  delong_tbl <- tibble(
    auc_rf = as.numeric(auc(rf_roc)),
    auc_lr = as.numeric(auc(lr_roc)),
    p_value_delong = delong_p,
    method = "DeLong paired ROC test"
  )
  write_csv(delong_tbl, file.path(tab_dir, "delong_test.csv"))

  roc_list <- setNames(
    list(rf_roc, lr_roc),
    c(
      "Random Forest",
      paste0("Logistic Regression (AUC=", format(round(as.numeric(auc(lr_roc)), 3), nsmall = 3), ")")
    )
  )

  p_roc <- ggroc(roc_list, legacy.axes = TRUE) +
    geom_abline(linetype = "dashed", color = "grey50") +
    labs(title = "ROC curve comparison", x = "1 - Specificity", y = "Sensitivity") +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, color = youreka_colors[["blue"]], face = "bold"),
      legend.title = element_blank()
    )
  ggsave(file.path(fig_dir, "roc_dual.png"), p_roc, width = 7, height = 5, dpi = 300)
} else {
  warning("model_predictions_cv.csv was not found; skipping DeLong ROC test.")
}

# Phase 5: Sex-stratified comparison + rank-biserial effect size
sex_candidates <- c("SEX", "sex", "GNDR_COD", "gndr_cod")
sex_matches <- sex_candidates[sex_candidates %in% names(df)]
sex_col <- if (length(sex_matches) > 0) sex_matches[[1]] else NULL

if (is.null(sex_col)) {
  # Fallback: recover sex from filtered teen DEMO records if absent in final table.
  demo_path <- "03_filtered/teen_demo_records.csv"
  if (file.exists(demo_path)) {
    demo_raw <- read_csv(demo_path, show_col_types = FALSE)
    demo_sex_candidates <- c("SEX", "sex", "GNDR_COD", "gndr_cod")
    demo_sex_col <- demo_sex_candidates[demo_sex_candidates %in% names(demo_raw)]

    if (length(demo_sex_col) == 0) {
      warning("Sex stratification skipped: no sex column in final dataset and no sex-like field in teen_demo_records.csv.")
      sex_df <- tibble()
    } else {
      demo_sex_col_name <- demo_sex_col[[1]]
      demo_df <- demo_raw %>%
        transmute(PRIMARYID = as.character(PRIMARYID), SEX = as.character(.data[[demo_sex_col_name]])) %>%
        distinct(PRIMARYID, .keep_all = TRUE)

      sex_df <- df %>%
        mutate(PRIMARYID = as.character(PRIMARYID)) %>%
        left_join(demo_df, by = "PRIMARYID") %>%
        mutate(
          SEX = toupper(as.character(SEX)),
          SEX = case_when(
            SEX %in% c("M", "MALE", "1") ~ "M",
            SEX %in% c("F", "FEMALE", "2") ~ "F",
            TRUE ~ "UNK"
          )
        ) %>%
        filter(SEX %in% c("M", "F"), !is.na(total_acb_codrugs_only))
    }
  } else {
    warning("Sex stratification skipped: no sex column in final dataset and teen_demo_records.csv not found.")
    sex_df <- tibble()
  }
} else {
  sex_df <- df %>%
    mutate(SEX = toupper(as.character(.data[[sex_col]]))) %>%
    filter(SEX %in% c("M", "F"), !is.na(total_acb_codrugs_only))
}

sex_p <- NA_real_
sex_r <- NA_real_
if (nrow(sex_df) > 0 && length(unique(sex_df$SEX)) == 2) {
  sex_w <- wilcox.test(total_acb_codrugs_only ~ SEX, data = sex_df, exact = FALSE)
  sex_p <- as.numeric(sex_w$p.value)

  # Rank-biserial fallback from Wilcoxon U statistic for robustness.
  groups <- split(sex_df$total_acb_codrugs_only, sex_df$SEX)
  g1_name <- names(groups)[1]
  g2_name <- names(groups)[2]
  n1 <- length(groups[[1]])
  n2 <- length(groups[[2]])
  W <- as.numeric(sex_w$statistic)
  U <- W - n1 * (n1 + 1) / 2
  sex_r <- (2 * U) / (n1 * n2) - 1

  # If effectsize is available, use its implementation as the primary estimate.
  if (has_effectsize) {
    rb <- effectsize::rank_biserial(total_acb_codrugs_only ~ SEX, data = sex_df)
    if ("r_rank_biserial" %in% names(rb)) {
      sex_r <- as.numeric(rb$r_rank_biserial[1])
    }
  }

  sex_tbl <- tibble(
    test = "sex_mann_whitney_total_acb_codrugs_only",
    p_value = sex_p,
    rank_biserial_r = sex_r,
    n = nrow(sex_df),
    group_1 = g1_name,
    group_2 = g2_name
  )
  write_csv(sex_tbl, file.path(tab_dir, "sex_stratified_results.csv"))
} else {
  warning("Sex stratification skipped: need both M and F groups after filtering UNK/missing.")
}

# Phase 3 continuation: BH correction across all study-level p-values
p_values_tbl <- tibble(
  test_name = c(
    "mann_whitney_acb_vs_cardiac",
    "kruskal_acb_vs_age_group",
    "spearman_acb_vs_severity",
    "delong_rf_vs_lr_auc",
    "sex_mann_whitney_total_acb_codrugs_only"
  ),
  p_raw = c(mw$p.value, kw$p.value, sp$p.value, delong_p, sex_p)
) %>%
  filter(!is.na(p_raw)) %>%
  mutate(
    p_bh = p.adjust(p_raw, method = "BH"),
    bh_significant_0_05 = p_bh < 0.05
  )

write_csv(p_values_tbl, file.path(tab_dir, "p_values_bh.csv"))

# Minimal visuals from core analyses
p1 <- ggplot(df, aes(x = factor(cardiac_any), y = total_acb_codrugs_only, fill = factor(cardiac_any))) +
  geom_boxplot(width = 0.55, alpha = 0.7, outlier.shape = 1, outlier.alpha = 0.55) +
  geom_jitter(width = 0.14, alpha = 0.25, size = 1.2, color = "gray40") +
  scale_fill_manual(values = c("0" = youreka_colors[["blue_light"]], "1" = youreka_colors[["orange_light"]])) +
  labs(
    title = "Anticholinergic burden by cardiac toxicity outcome",
    subtitle = "Co-medication anticholinergic burden by outcome group",
    x = "Cardiac toxicity outcome (0 = no, 1 = yes)",
    y = "Total anticholinergic burden from co-medications"
  ) +
  theme(legend.position = "none")
ggsave(file.path(fig_dir, "acb_by_cardiac_outcome.png"), p1, width = 7, height = 5, dpi = 220)

p_stars <- function(p) {
  if (is.na(p)) return("ns")
  if (p < 0.001) return("***")
  if (p < 0.01) return("**")
  if (p < 0.05) return("*")
  "ns"
}

df_age <- df %>%
  filter(!is.na(age_group), !is.na(total_acb_codrugs_only)) %>%
  mutate(age_group = factor(as.character(age_group), levels = c("13-15", "16-17", "18-19")))

ymax_age <- max(df_age$total_acb_codrugs_only, na.rm = TRUE)
label_y <- ymax_age * 0.84
bracket_y <- ymax_age * 0.94
kw_label <- sprintf("Kruskal p=%.3g %s", kw$p.value, p_stars(kw$p.value))

p2 <- ggplot(df_age, aes(x = age_group, y = total_acb_codrugs_only, fill = age_group)) +
  geom_boxplot(width = 0.55, alpha = 0.65, outlier.shape = 1, outlier.alpha = 0.55) +
  geom_jitter(width = 0.16, alpha = 0.35, size = 1.5, color = "gray40") +
  scale_fill_manual(values = c("13-15" = youreka_colors[["blue_light"]], "16-17" = youreka_colors[["orange_light"]], "18-19" = youreka_colors[["mint_light"]])) +
  labs(
    title = "Co-medication ACB by Age Group",
    subtitle = "Total co-medication anticholinergic burden across age groups",
    x = "Age group",
    y = "Total ACB (co-medications only)"
  ) +
  annotate("text", x = 2, y = label_y, label = kw_label, color = "#b91c1c", fontface = "bold", size = 4.4) +
  coord_cartesian(ylim = c(min(0, min(df_age$total_acb_codrugs_only, na.rm = TRUE)), ymax_age * 1.10)) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, color = youreka_colors[["blue"]], face = "bold", size = 17),
    plot.subtitle = element_text(hjust = 0.5, color = youreka_colors[["blue"]], face = "bold", size = 12)
  )

if (!is.na(kw$p.value) && kw$p.value < 0.05) {
  p2 <- p2 +
    annotate("segment", x = 2, xend = 3, y = bracket_y, yend = bracket_y, color = "gray30", linewidth = 0.6) +
    annotate("segment", x = 2, xend = 2, y = bracket_y, yend = bracket_y - ymax_age * 0.015, color = "gray30", linewidth = 0.6) +
    annotate("segment", x = 3, xend = 3, y = bracket_y, yend = bracket_y - ymax_age * 0.015, color = "gray30", linewidth = 0.6) +
    annotate("text", x = 2.5, y = bracket_y + ymax_age * 0.03, label = p_stars(kw$p.value), color = "#b91c1c", fontface = "bold", size = 4.4)
}

ggsave(file.path(fig_dir, "acb_by_age_group.png"), p2, width = 7, height = 5, dpi = 220)

p3 <- ggplot(df, aes(x = total_acb_codrugs_only, y = max_severity)) +
  geom_point(color = youreka_colors[["teal"]], alpha = 0.65, size = 1.8) +
  geom_smooth(method = "loess", se = FALSE, color = youreka_colors[["orange"]], linewidth = 1.2) +
  labs(
    title = "Anticholinergic burden versus case severity",
    subtitle = "Smoothed relationship between co-medication burden and severity",
    x = "Total anticholinergic burden from co-medications",
    y = "Maximum reported case severity"
  )
ggsave(file.path(fig_dir, "acb_vs_severity.png"), p3, width = 7, height = 5, dpi = 220)

# Phase 5 figure with raw p-value, BH status, and rank-biserial effect size.
sex_bh_status <- p_values_tbl %>%
  filter(test_name == "sex_mann_whitney_total_acb_codrugs_only") %>%
  mutate(status = if_else(bh_significant_0_05, "survives BH", "does not survive BH")) %>%
  pull(status)
if (length(sex_bh_status) == 0) {
  sex_bh_status <- "BH not available"
}

if (nrow(sex_df) > 0 && length(unique(sex_df$SEX)) == 2) {
  sex_caption <- sprintf("Mann-Whitney p=%.3g | BH: %s | Rank-biserial r=%.3f", sex_p, sex_bh_status, sex_r)
  p_sex <- ggplot(sex_df, aes(x = SEX, y = total_acb_codrugs_only, fill = SEX))

  if (has_introdataviz) {
    p_sex <- p_sex + introdataviz::geom_split_violin(alpha = 0.7, trim = FALSE, color = "gray30")
  } else {
    p_sex <- p_sex + geom_violin(trim = FALSE, alpha = 0.55, width = 0.9, color = "gray30")
  }

  p_sex <- p_sex +
    geom_boxplot(width = 0.12, outlier.shape = NA, alpha = 0.85) +
    geom_jitter(width = 0.09, alpha = 0.2, size = 1.1) +
    scale_fill_manual(values = c("F" = youreka_colors[["mint_light"]], "M" = youreka_colors[["blue_light"]])) +
    annotate("text", x = 1.5, y = max(sex_df$total_acb_codrugs_only, na.rm = TRUE), label = sex_caption, vjust = -0.6, size = 4.0) +
    labs(
      title = "Sex-stratified anticholinergic burden (co-medications)",
      subtitle = "Unknown sex excluded",
      x = "Sex",
      y = "Total anticholinergic burden from co-medications"
    ) +
    theme_classic(base_size = 12) +
    theme(legend.position = "none")

  ggsave(file.path(fig_dir, "sex_stratified_acb_violin.png"), p_sex, width = 8, height = 5.5, dpi = 240)
}

# Phase 4: Drug co-occurrence network analysis
drug_path <- "04_processed/dph_drug_normalized.csv"
if (file.exists(drug_path)) {
  drug_df <- read_csv(drug_path, show_col_types = FALSE)
  drug_df <- drug_df %>%
    mutate(
      PRIMARYID = as.character(PRIMARYID),
      generic_name = as.character(generic_name),
      DRUGNAME = as.character(DRUGNAME),
      drug_name = case_when(
        !is.na(generic_name) & generic_name != "" ~ tolower(generic_name),
        !is.na(DRUGNAME) & DRUGNAME != "" ~ tolower(DRUGNAME),
        TRUE ~ NA_character_
      )
    ) %>%
    mutate(
      # Collapse common brand/generic aliases to canonical generic names.
      drug_name = recode(
        drug_name,
        "benadryl" = "diphenhydramine",
        "tylenol" = "acetaminophen",
        "advil" = "ibuprofen",
        "motrin" = "ibuprofen",
        "zofran" = "ondansetron"
      )
    ) %>%
    filter(!is.na(PRIMARYID), !is.na(drug_name), drug_name != "") %>%
    distinct(PRIMARYID, drug_name)

  # Explicitly rebuild binary exposure matrix and adjacency matrix before graphing.
  binary_mat <- xtabs(~ PRIMARYID + drug_name, data = drug_df)
  binary_mat[binary_mat > 0] <- 1
  adjacency <- as.matrix(crossprod(binary_mat))
  diag(adjacency) <- 0

  pairs_raw <- as.data.frame(as.table(adjacency), stringsAsFactors = FALSE)
  names(pairs_raw)[1:3] <- c("drug_name_a", "drug_name_b", "weight")

  pairs <- pairs_raw %>%
    transmute(drug_name_a = as.character(drug_name_a), drug_name_b = as.character(drug_name_b), weight = as.numeric(weight)) %>%
    filter(weight > 0, drug_name_a < drug_name_b) %>%
    arrange(desc(weight))

  if (nrow(pairs) > 0) {
    weight_summary <- pairs %>% summarize(
      n_edges = n(),
      min_weight = min(weight),
      q1_weight = quantile(weight, 0.25),
      median_weight = median(weight),
      mean_weight = mean(weight),
      q3_weight = quantile(weight, 0.75),
      max_weight = max(weight)
    )
    write_csv(weight_summary, file.path(tab_dir, "network_edge_weight_summary.csv"))

    weight_table <- pairs %>% count(weight, name = "n_edges") %>% arrange(weight)
    write_csv(weight_table, file.path(tab_dir, "network_edge_weight_distribution.csv"))

    # Empirical threshold: keep upper-quartile edges, minimum of 2 co-occurrences.
    edge_threshold <- max(2, floor(as.numeric(quantile(pairs$weight, 0.75))))
    filtered_pairs <- pairs %>% filter(weight >= edge_threshold)

    g <- graph_from_data_frame(filtered_pairs, directed = FALSE)
    g <- delete_vertices(g, degree(g) == 0)

    network_meta <- tibble(
      edge_threshold = edge_threshold,
      n_edges_before = nrow(pairs),
      n_edges_after = nrow(filtered_pairs),
      n_nodes_after = vcount(g)
    )
    write_csv(network_meta, file.path(tab_dir, "network_threshold_metadata.csv"))

    if (has_ggraph && ecount(g) > 0) {
      comps <- components(g)
      largest_comp_id <- which.max(comps$csize)
      keep_vertices <- V(g)[comps$membership == largest_comp_id]
      g_main <- induced_subgraph(g, keep_vertices)

      node_tbl <- tibble(
        name = V(g_main)$name,
        degree = degree(g_main),
        strength = strength(g_main, weights = E(g_main)$weight)
      ) %>%
        arrange(desc(strength), desc(degree))

      max_nodes <- 55
      if (nrow(node_tbl) > max_nodes) {
        keep_names <- node_tbl %>% slice_head(n = max_nodes) %>% pull(name)
        g_main <- induced_subgraph(g_main, vids = V(g_main)[name %in% keep_names])
        node_tbl <- node_tbl %>% filter(name %in% keep_names)
      }

      V(g_main)$label <- V(g_main)$name
      V(g_main)$strength <- node_tbl$strength[match(V(g_main)$name, node_tbl$name)]

      net_plot <- ggraph::ggraph(g_main, layout = "fr") +
        ggraph::geom_edge_link(aes(width = weight, alpha = weight), color = "#8a96a3") +
        ggraph::scale_edge_width_continuous(range = c(0.4, 2.3), guide = "none") +
        ggraph::scale_edge_alpha_continuous(range = c(0.20, 0.72), guide = "none") +
        ggraph::geom_node_point(aes(size = strength), color = youreka_colors[["blue"]], alpha = 0.92) +
        scale_size_continuous(range = c(2.4, 7.2), guide = "none") +
        ggraph::geom_node_label(
          aes(label = label),
          repel = TRUE,
          size = 2.55,
          label.padding = unit(0.12, "lines"),
          label.size = 0.15,
          fill = "#ffffff",
          color = "#1f2937",
          segment.color = "#9ca3af",
          max.overlaps = Inf
        ) +
        labs(title = "Drug co-occurrence network") +
        theme_void(base_size = 12) +
        theme(
          plot.background = element_rect(fill = "white", color = NA),
          panel.background = element_rect(fill = "white", color = NA),
          plot.margin = margin(14, 22, 12, 22),
          plot.title = element_text(hjust = 0.5, color = youreka_colors[["blue"]], face = "bold", size = 16)
        )

      ggsave(
        file.path(fig_dir, "drug_cooccurrence_network.png"),
        net_plot,
        width = 13,
        height = 9,
        dpi = 260,
        bg = "white"
      )
    }
  } else {
    warning("Network analysis skipped: no co-occurring drug pairs found.")
  }
} else {
  warning("dph_drug_normalized.csv not found; skipping network analysis.")
}

writeLines(c(
  "Advanced R analysis complete.",
  "Outputs written to 06_analysis/tables and 06_analysis/figures.",
  "Includes DeLong ROC testing, BH correction, network analysis, and sex stratification."
), con = file.path(out_dir, "results_summary_r.txt"))
