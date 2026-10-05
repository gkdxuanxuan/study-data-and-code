# ============================================================
# Figure 4: eDNA与镜检浮游植物α多样性
# (a,b,d,e) 箱线图 + 同色系弱化离群点 + Tukey字母
# (c,f)     两种方法的逐样本线性回归
# 镜检多样性按细胞丰度计算，不使用生物量表
# ============================================================

setwd("D:/A_师姐文章/原始处理数据/111")

# 所有图4图片和统计结果统一输出到“图四”文件夹
out_dir <- "图四"
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# 1. 安装并加载所需R包
pkgs <- c("ggplot2", "dplyr", "tidyr", "vegan", "multcompView", "patchwork")
new_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(new_pkgs) > 0) install.packages(new_pkgs)
invisible(lapply(pkgs, library, character.only = TRUE))

# 2. 读取数据
edna_abund  <- read.csv("2.D_ASV丰度表_856reads抽平_329ASV_Richness显著.csv", check.names = FALSE)
edna_tax    <- read.csv("2.Dphyto_tax_filtered.csv", check.names = FALSE)
micro_abund <- read.csv("2.D原始镜检-细胞丰度.csv", check.names = FALSE)
micro_tax   <- read.csv("2.DTAXA-JJ-ZH_规范化.csv", check.names = FALSE)
group       <- read.csv("3.group.csv", check.names = FALSE)

# 生物量表的新名称是“2.D原始镜检生物量.csv”，但图4不需读取它

# 3. 统一ID列名
names(edna_abund)[1]  <- "ASV_ID"
names(micro_abund)[1] <- "Taxon_ID"

if ("OTU ID" %in% names(edna_tax)) {
  names(edna_tax)[names(edna_tax) == "OTU ID"] <- "ASV_ID"
} else {
  names(edna_tax)[1] <- "ASV_ID"
}
names(micro_tax)[1] <- "Taxon_ID"

# 4. 样本和分组检查
sample_ids <- group$ID
if (!all(sample_ids %in% names(edna_abund)))  stop("eDNA丰度表缺少group中的样本。")
if (!all(sample_ids %in% names(micro_abund))) stop("镜检丰度表缺少group中的样本。")
if (anyDuplicated(edna_abund$ASV_ID))  stop("eDNA丰度表ASV_ID存在重复。")
if (anyDuplicated(edna_tax$ASV_ID))    stop("eDNA分类表ASV_ID存在重复。")
if (anyDuplicated(micro_abund$Taxon_ID)) stop("镜检丰度表Taxon_ID存在重复。")
if (anyDuplicated(micro_tax$Taxon_ID))   stop("镜检分类表Taxon_ID存在重复。")
if (nrow(edna_abund) != 329) stop("正式eDNA抽平表应包含329个ASV。")
if (!all(colSums(edna_abund[, sample_ids, drop = FALSE]) == 856)) {
  stop("正式eDNA抽平表中每个样本的序列数均应为856 reads。")
}

# 5. 清理属名：只去掉g__前缀；norank/空值不作为一个属计数
clean_genus <- function(x) {
  x <- trimws(as.character(x))
  x <- sub("^g__", "", x)
  bad <- is.na(x) | x == "" |
    grepl("^(norank|unclassified|uncultured|unknown|NA)$", x, ignore.case = TRUE)
  x[bad] <- NA_character_
  x
}

# 6. 把ASV/镜检类群丰度汇总到属水平
make_genus_matrix <- function(abund, tax, id_col) {
  idx <- match(abund[[id_col]], tax[[id_col]])
  if (anyNA(idx)) {
    stop(paste0(id_col, "有", sum(is.na(idx)), "条未匹配分类注释。"))
  }
  
  genus <- clean_genus(tax$Genus[idx])
  x <- data.frame(Genus = genus, abund[, sample_ids, drop = FALSE],
                  check.names = FALSE)
  x <- x[!is.na(x$Genus), , drop = FALSE]
  
  x |>
    dplyr::group_by(Genus) |>
    dplyr::summarise(dplyr::across(dplyr::all_of(sample_ids), sum),
                     .groups = "drop") |>
    as.data.frame()
}

edna_genus  <- make_genus_matrix(edna_abund, edna_tax, "ASV_ID")
micro_genus <- make_genus_matrix(micro_abund, micro_tax, "Taxon_ID")

# 7. 计算每个样本的属丰富度和Shannon指数
calc_alpha <- function(genus_table, method_name) {
  mat <- as.matrix(genus_table[, sample_ids, drop = FALSE])
  storage.mode(mat) <- "numeric"
  data.frame(
    ID = sample_ids,
    Method = method_name,
    Genus_richness = colSums(mat > 0),
    Shannon = vegan::diversity(t(mat), index = "shannon"),
    check.names = FALSE
  )
}

alpha <- dplyr::bind_rows(
  calc_alpha(edna_genus, "eDNA"),
  calc_alpha(micro_genus, "Morphological identification")
) |>
  dplyr::left_join(group, by = "ID")

# 按group表中14_Spring等编码生成年份和季节
alpha <- alpha |>
  dplyr::mutate(
    Year = ifelse(grepl("^14_", season), "2014", "2015"),
    Season = sub("^[0-9]+_", "", season),
    Period = factor(
      season,
      levels = c("14_Spring", "14_Summer", "14_Autumn", "14_Winter",
                 "15_Spring", "15_Summer", "15_Autumn", "15_Winter")
    )
  )

if (anyNA(alpha$Period)) stop("group表season格式不是14_Spring至15_Winter，请检查。")

# 转为作图和统计的长表
alpha_long <- alpha |>
  tidyr::pivot_longer(
    cols = c(Genus_richness, Shannon),
    names_to = "Metric", values_to = "Value"
  )

# 8. 每组描述统计：论文结果可直接引用
group_summary <- alpha_long |>
  dplyr::group_by(Method, Metric, Period, Year, Season) |>
  dplyr::summarise(
    N = dplyr::n(),
    Mean = mean(Value),
    SD = sd(Value),
    Median = median(Value),
    Min = min(Value),
    Max = max(Value),
    .groups = "drop"
  )

# 9. 单因素ANOVA + Tukey HSD + 显著性字母
run_group_test <- function(dat) {
  fit <- stats::aov(Value ~ Period, data = dat)
  sm <- summary(fit)[[1]]
  anova_out <- data.frame(
    Term = trimws(rownames(sm)),
    Df = sm[, "Df"],
    Sum_Sq = sm[, "Sum Sq"],
    Mean_Sq = sm[, "Mean Sq"],
    F_value = sm[, "F value"],
    P_value = sm[, "Pr(>F)"],
    row.names = NULL
  )
  
  tuk <- stats::TukeyHSD(fit, "Period")$Period
  tukey_out <- data.frame(Comparison = rownames(tuk), tuk,
                          row.names = NULL, check.names = FALSE)
  
  letter_obj <- multcompView::multcompLetters4(fit, stats::TukeyHSD(fit))
  letter_vec <- letter_obj[[1]]$Letters
  letters_out <- data.frame(
    Period = factor(names(letter_vec), levels = levels(dat$Period)),
    Letter = unname(letter_vec)
  ) |>
    dplyr::arrange(Period)
  
  list(anova = anova_out, tukey = tukey_out, letters = letters_out)
}

split_data <- split(alpha_long, interaction(alpha_long$Method,
                                            alpha_long$Metric, drop = TRUE))
tests <- lapply(split_data, run_group_test)

anova_table <- dplyr::bind_rows(lapply(names(tests), function(nm) {
  z <- strsplit(nm, "\\.")[[1]]
  data.frame(Method_Metric = nm, tests[[nm]]$anova)
}))

tukey_table <- dplyr::bind_rows(lapply(names(tests), function(nm) {
  data.frame(Method_Metric = nm, tests[[nm]]$tukey)
}))

letter_table <- dplyr::bind_rows(lapply(names(tests), function(nm) {
  data.frame(Method_Metric = nm, tests[[nm]]$letters)
}))

# 将Method_Metric重新映射为两列，避免方法名中的点影响解析
key_table <- alpha_long |>
  dplyr::mutate(Method_Metric = interaction(Method, Metric, drop = TRUE)) |>
  dplyr::distinct(Method_Metric, Method, Metric) |>
  dplyr::mutate(Method_Metric = as.character(Method_Metric))

anova_table  <- dplyr::left_join(anova_table, key_table, by = "Method_Metric") |>
  dplyr::select(Method, Metric, dplyr::everything(), -Method_Metric)
tukey_table  <- dplyr::left_join(tukey_table, key_table, by = "Method_Metric") |>
  dplyr::select(Method, Metric, dplyr::everything(), -Method_Metric)
letter_table <- dplyr::left_join(letter_table, key_table, by = "Method_Metric") |>
  dplyr::select(Method, Metric, Period, Letter)

# 字母标注的纵坐标
letter_positions <- group_summary |>
  dplyr::left_join(letter_table, by = c("Method", "Metric", "Period")) |>
  dplyr::group_by(Method, Metric) |>
  dplyr::mutate(
    panel_range = max(Max) - min(Min),
    Letter_y = pmax(Max, Mean + SD) + ifelse(panel_range == 0, 0.2,
                                             0.07 * panel_range)
  ) |>
  dplyr::ungroup()

# 10. 两种方法的逐样本回归
paired_edna <- alpha |>
  dplyr::filter(Method == "eDNA") |>
  dplyr::select(ID, Richness_eDNA = Genus_richness,
                Shannon_eDNA = Shannon)

paired_micro <- alpha |>
  dplyr::filter(Method == "Morphological identification") |>
  dplyr::select(ID, Richness_Morphology = Genus_richness,
                Shannon_Morphology = Shannon)

paired <- dplyr::inner_join(paired_edna, paired_micro, by = "ID")

regression_stats <- function(y, x, metric) {
  fit <- stats::lm(y ~ x)
  sm <- summary(fit)
  ci <- stats::confint(fit)
  data.frame(
    Metric = metric,
    N = length(x),
    Intercept = coef(fit)[1],
    Slope = coef(fit)[2],
    Slope_CI_low = ci[2, 1],
    Slope_CI_high = ci[2, 2],
    R_squared = sm$r.squared,
    Adjusted_R_squared = sm$adj.r.squared,
    P_value = coef(sm)[2, 4]
  )
}

regression_table <- dplyr::bind_rows(
  regression_stats(paired$Richness_Morphology, paired$Richness_eDNA,
                   "Genus richness"),
  regression_stats(paired$Shannon_Morphology, paired$Shannon_eDNA,
                   "Shannon-Wiener index")
)

# 11. 专用于撰写结果的摘要表
overall_range <- alpha_long |>
  dplyr::group_by(Method, Metric) |>
  dplyr::summarise(Overall_min = min(Value), Overall_max = max(Value),
                   .groups = "drop")

writing_summary <- group_summary |>
  dplyr::group_by(Method, Metric) |>
  dplyr::summarise(
    Highest_period = as.character(Period[which.max(Mean)]),
    Highest_mean = max(Mean),
    Lowest_period = as.character(Period[which.min(Mean)]),
    Lowest_mean = min(Mean),
    .groups = "drop"
  ) |>
  dplyr::left_join(overall_range, by = c("Method", "Metric")) |>
  dplyr::left_join(
    anova_table |>
      dplyr::filter(Term == "Period") |>
      dplyr::select(Method, Metric, ANOVA_F = F_value, ANOVA_P = P_value),
    by = c("Method", "Metric")
  )

# 12. 作图设置
period_labels <- c(
  "14_Spring" = "Spr\n2014", "14_Summer" = "Sum\n2014",
  "14_Autumn" = "Aut\n2014", "14_Winter" = "Win\n2014",
  "15_Spring" = "Spr\n2015", "15_Summer" = "Sum\n2015",
  "15_Autumn" = "Aut\n2015", "15_Winter" = "Win\n2015"
)
year_colors <- c("2014" = "#F6DDE2", "2015" = "#CFE7F3")
year_line_colors <- c("2014" = "#B96E7A", "2015" = "#568EAE")

base_theme <- ggplot2::theme_classic(base_size = 14) +
  ggplot2::theme(
    axis.title = ggplot2::element_text(size = 14, face = "bold"),
    axis.text = ggplot2::element_text(size = 12, color = "black", face = "bold"),
    plot.title = ggplot2::element_text(size = 15, hjust = 0, face = "bold"),
    plot.margin = ggplot2::margin(10, 12, 10, 10),
    legend.position = "none"
  )

make_box_panel <- function(method_name, metric_name, panel_title, y_title) {
  dat <- alpha_long |>
    dplyr::filter(Method == method_name, Metric == metric_name)
  sm <- letter_positions |>
    dplyr::filter(Method == method_name, Metric == metric_name)
  
  ggplot2::ggplot(
    dat,
    ggplot2::aes(x = Period, y = Value, fill = Year, color = Year)
  ) +
    ggplot2::geom_boxplot(
      width = 0.66,
      linewidth = 0.48,
      outlier.shape = 16,
      outlier.size = 1.15,
      outlier.stroke = 0,
      outlier.alpha = 0.42
    ) +
    ggplot2::geom_text(
      data = sm,
      ggplot2::aes(x = Period, y = Letter_y, label = Letter),
      inherit.aes = FALSE,
      size = 5.2, vjust = 0, fontface = "bold"
    ) +
    ggplot2::scale_fill_manual(values = year_colors) +
    ggplot2::scale_color_manual(values = year_line_colors) +
    ggplot2::scale_x_discrete(labels = period_labels) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0.08, 0.15))
    ) +
    ggplot2::labs(title = panel_title, x = NULL, y = y_title) +
    ggplot2::coord_cartesian(clip = "off") +
    base_theme
}

make_reg_panel <- function(x, y, metric, panel_title, x_title, y_title) {
  st <- regression_table[regression_table$Metric == metric, ]
  label <- if (st$P_value < 0.001) {
    sprintf("R² = %.2f\nP < 0.001", st$R_squared)
  } else {
    sprintf("R² = %.2f\nP = %.3f", st$R_squared, st$P_value)
  }
  
  dat <- data.frame(x = x, y = y)
  ggplot2::ggplot(dat, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_smooth(method = "lm", se = TRUE,
                         color = "#E98792", fill = "#F7CBD0",
                         linewidth = 0.9, alpha = 0.45) +
    ggplot2::geom_point(color = "#2677A6", size = 2.7) +
    ggplot2::annotate("text", x = Inf, y = Inf, label = label,
                      hjust = 1.08, vjust = 1.15, size = 5.2,
                      fontface = "bold") +
    ggplot2::labs(title = panel_title, x = x_title, y = y_title) +
    base_theme
}
# 13. 六个子图
p_a <- make_box_panel("eDNA", "Genus_richness",
                      "(a) eDNA", "Genus richness")
p_b <- make_box_panel("Morphological identification", "Genus_richness",
                      "(b) Morphological identification", "Genus richness")
p_c <- make_reg_panel(
  paired$Richness_eDNA, paired$Richness_Morphology,
  "Genus richness", "(c) Genus richness",
  "eDNA", "Morphological identification"
)
p_d <- make_box_panel("eDNA", "Shannon",
                      "(d) eDNA", "Shannon-Wiener index")
p_e <- make_box_panel("Morphological identification", "Shannon",
                      "(e) Morphological identification", "Shannon-Wiener index")
p_f <- make_reg_panel(
  paired$Shannon_eDNA, paired$Shannon_Morphology,
  "Shannon-Wiener index", "(f) Shannon-Wiener index",
  "eDNA", "Morphological identification"
)

# 不使用 & theme，避免之前patchwork的方法冲突
p4 <- patchwork::wrap_plots(
  p_a, p_b, p_c, p_d, p_e, p_f,
  ncol = 3, nrow = 2
)

print(p4)

# 14. 输出图片
ggplot2::ggsave(file.path(out_dir, "Figure4_updated_R.png"), p4,
                width = 14.5, height = 9.2, units = "in", dpi = 400,
                bg = "white")
ggplot2::ggsave(file.path(out_dir, "Figure4_updated_R.tiff"), p4,
                width = 14.5, height = 9.2, units = "in", dpi = 600,
                compression = "lzw", bg = "white")
ggplot2::ggsave(file.path(out_dir, "Figure4_updated_R.pdf"), p4,
                width = 14.5, height = 9.2, units = "in", bg = "white")

# 15. 输出写结论和复核所需的数据
write.csv(alpha, file.path(out_dir, "Figure4_alpha_diversity_per_sample_R.csv"), row.names = FALSE)
write.csv(group_summary, file.path(out_dir, "Figure4_group_mean_SD_R.csv"), row.names = FALSE)
write.csv(anova_table, file.path(out_dir, "Figure4_ANOVA_R.csv"), row.names = FALSE)
write.csv(tukey_table, file.path(out_dir, "Figure4_TukeyHSD_R.csv"), row.names = FALSE)
write.csv(letter_table, file.path(out_dir, "Figure4_significance_letters_R.csv"), row.names = FALSE)
write.csv(paired, file.path(out_dir, "Figure4_regression_paired_data_R.csv"), row.names = FALSE)
write.csv(regression_table, file.path(out_dir, "Figure4_regression_statistics_R.csv"), row.names = FALSE)
write.csv(writing_summary, file.path(out_dir, "Figure4_results_for_writing_R.csv"), row.names = FALSE)

# 16. 控制台摘要
cat("\n================ Figure 4 completed ================\n")
cat("共有样本数：", length(sample_ids), "\n")
cat("eDNA纳入计算的已命名属数：", nrow(edna_genus), "\n")
cat("镜检纳入计算的已命名属数：", nrow(micro_genus), "\n\n")
print(writing_summary)
cat("\n回归结果：\n")
print(regression_table)
cat("====================================================\n")
cat("输出文件夹：", normalizePath(out_dir), "\n")
