library(tidyverse)
library(broom)
library(rstatix)

# Optional packages for advanced tables/plots
# install.packages(c("ggpubr", "FSA"))

out_dir <- "06_analysis"
fig_dir <- file.path(out_dir, "figures")
tab_dir <- file.path(out_dir, "tables")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)

df <- read_csv("05_final/analysis_table.csv", show_col_types = FALSE)

# Basic coercions
num_cols <- c("AGE", "n_codrugs", "total_acb_with_dph", "total_acb_codrugs_only", "max_severity", "cardiac_any", "pre_post_warning")
df <- df %>% mutate(across(all_of(intersect(num_cols, names(df))), as.numeric))

# Table 1 descriptive
cont <- c("AGE", "total_acb_with_dph", "total_acb_codrugs_only", "n_codrugs", "max_severity")
cont_tbl <- map_dfr(cont, function(v){
  x <- df[[v]]
  tibble(variable=v, n=sum(!is.na(x)), mean=mean(x, na.rm=TRUE), sd=sd(x, na.rm=TRUE), median=median(x, na.rm=TRUE), iqr=IQR(x, na.rm=TRUE))
})
write_csv(cont_tbl, file.path(tab_dir, "table1_continuous.csv"))

cat_tbl <- bind_rows(
  df %>% count(age_group, name="count") %>% mutate(variable="age_group", level=as.character(age_group), percent=100*count/sum(count)) %>% select(variable, level, count, percent),
  df %>% count(pre_post_warning, name="count") %>% mutate(variable="pre_post_warning", level=as.character(pre_post_warning), percent=100*count/sum(count)) %>% select(variable, level, count, percent),
  df %>% count(cardiac_any, name="count") %>% mutate(variable="cardiac_any", level=as.character(cardiac_any), percent=100*count/sum(count)) %>% select(variable, level, count, percent)
)
write_csv(cat_tbl, file.path(tab_dir, "table1_categorical.csv"))

# Normality checks
shap <- map_dfr(cont, function(v){
  x <- df[[v]] %>% na.omit()
  if(length(x) >= 3){
    s <- shapiro.test(x)
    tibble(variable=v, W=as.numeric(s$statistic), p_value=s$p.value)
  } else {
    tibble(variable=v, W=NA_real_, p_value=NA_real_)
  }
})
write_csv(shap, file.path(tab_dir, "normality_shapiro.csv"))

# Core analyses
mw <- wilcox.test(total_acb_codrugs_only ~ cardiac_any, data=df)
kw <- kruskal.test(total_acb_with_dph ~ age_group, data=df)
sp <- cor.test(df$total_acb_codrugs_only, df$max_severity, method="spearman")

core <- tibble(
  analysis=c("mann_whitney_acb_vs_cardiac", "kruskal_acb_vs_age_group", "spearman_acb_vs_severity"),
  statistic=c(as.numeric(mw$statistic), as.numeric(kw$statistic), as.numeric(sp$estimate)),
  p_value=c(mw$p.value, kw$p.value, sp$p.value)
)
write_csv(core, file.path(tab_dir, "nonparametric_tests.csv"))

# Logistic models
model_full <- glm(cardiac_any ~ total_acb_codrugs_only + age_group + pre_post_warning + n_codrugs,
                  data=df, family=binomial)
full_tbl <- tidy(model_full, conf.int=TRUE, exponentiate=TRUE)
write_csv(full_tbl, file.path(tab_dir, "table2_logit_full.csv"))

model_pre <- glm(cardiac_any ~ total_acb_codrugs_only + age_group + n_codrugs,
                 data=filter(df, pre_post_warning==0), family=binomial)
pre_tbl <- tidy(model_pre, conf.int=TRUE, exponentiate=TRUE)
write_csv(pre_tbl, file.path(tab_dir, "table2_logit_pre.csv"))

model_post <- glm(cardiac_any ~ total_acb_codrugs_only + age_group + n_codrugs,
                  data=filter(df, pre_post_warning==1), family=binomial)
post_tbl <- tidy(model_post, conf.int=TRUE, exponentiate=TRUE)
write_csv(post_tbl, file.path(tab_dir, "table2_logit_post.csv"))

# Minimal visuals
p1 <- ggplot(df, aes(x=factor(cardiac_any), y=total_acb_codrugs_only)) +
  geom_boxplot() +
  labs(title="Figure 2: ACB (co-drugs) by cardiac outcome", x="cardiac_any", y="total_acb_codrugs_only")
ggsave(file.path(fig_dir, "figure2_acb_cardiac_boxplot.png"), p1, width=7, height=5, dpi=220)

p2 <- ggplot(df, aes(x=age_group, y=total_acb_with_dph)) +
  geom_boxplot() +
  labs(title="Figure 3: ACB (with DPH) by age group", x="age_group", y="total_acb_with_dph")
ggsave(file.path(fig_dir, "figure3_acb_agegroup_boxplot.png"), p2, width=7, height=5, dpi=220)

p3 <- ggplot(df, aes(x=total_acb_codrugs_only, y=max_severity)) +
  geom_point() + geom_smooth(method="loess", se=FALSE) +
  labs(title="Figure 4: ACB vs max severity", x="total_acb_codrugs_only", y="max_severity")
ggsave(file.path(fig_dir, "figure4_spearman_acb_severity.png"), p3, width=7, height=5, dpi=220)

writeLines(c(
  "Advanced R analysis template complete.",
  "Outputs written to 06_analysis/tables and 06_analysis/figures."
), con = file.path(out_dir, "analysis_summary_R.txt"))
