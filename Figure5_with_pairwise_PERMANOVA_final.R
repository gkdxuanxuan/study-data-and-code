# ============================================================
# Figure 5：属水平 beta 多样性
# 仅输出正式图片、对应作图数据和关键统计结果
# ============================================================

setwd("D:/A_师姐文章/原始处理数据/111")
input_dir <- getwd()
output_dir <- file.path(input_dir, "图五new")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
# 清除该文件夹中旧版Figure5输出，避免与本次结果混在一起
old_outputs <- list.files(output_dir, pattern = "^Figure5", full.names = TRUE)
if (length(old_outputs) > 0) file.remove(old_outputs)

pkgs <- c("ggplot2", "dplyr", "vegan", "ape", "multcompView",
          "patchwork", "openxlsx")
new_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(new_pkgs) > 0) install.packages(new_pkgs)
invisible(lapply(pkgs, library, character.only = TRUE))
set.seed(20260909)

# 1. 读取数据
edna_abund  <- read.csv("2.D_ASV丰度表_856reads抽平_329ASV_Richness显著.csv",
                        check.names = FALSE)
edna_tax    <- read.csv("2.Dphyto_tax_filtered.csv", check.names = FALSE)
micro_abund <- read.csv("2.D原始镜检-细胞丰度.csv", check.names = FALSE)
micro_tax   <- read.csv("2.DTAXA-JJ-ZH_规范化.csv", check.names = FALSE)
group       <- read.csv("3.group.csv", check.names = FALSE)

names(edna_abund)[1]  <- "ASV_ID"
names(micro_abund)[1] <- "Taxon_ID"
if ("OTU ID" %in% names(edna_tax)) {
  names(edna_tax)[names(edna_tax) == "OTU ID"] <- "ASV_ID"
} else {
  names(edna_tax)[1] <- "ASV_ID"
}
names(micro_tax)[1] <- "Taxon_ID"

period_levels <- c("14_Spring", "14_Summer", "14_Autumn", "14_Winter",
                   "15_Spring", "15_Summer", "15_Autumn", "15_Winter")
season_levels <- c("Spring", "Summer", "Autumn", "Winter")
group <- group |>
  dplyr::mutate(
    ID = as.character(ID),
    Year = factor(ifelse(grepl("^14_", season), "2014", "2015"),
                  levels = c("2014", "2015")),
    Season = factor(sub("^[0-9]+_", "", season), levels = season_levels),
    Period = factor(season, levels = period_levels),
    Event_index = match(as.character(Period), period_levels)
  )

if (nrow(group) != 39 || anyDuplicated(group$ID))
  stop("group表必须包含39个不重复样品。")
if (anyNA(group$Period)) stop("season编码有误。")
sample_ids <- group$ID
if (!all(sample_ids %in% names(edna_abund))) stop("eDNA表缺少分组样品。")
if (!all(sample_ids %in% names(micro_abund))) stop("镜检表缺少分组样品。")
if (nrow(edna_abund) != 329) stop("eDNA表应包含329个ASV。")
if (!all(colSums(edna_abund[, sample_ids, drop = FALSE]) == 856))
  stop("eDNA表中每个样品应为856 reads。")

# 2. 汇总至属水平
clean_genus <- function(x) {
  x <- trimws(as.character(x))
  x <- sub("^g__", "", x)
  bad <- is.na(x) | x == "" |
    grepl("^(norank|unclassified|uncultured|unknown|NA)$", x,
          ignore.case = TRUE)
  x[bad] <- NA_character_
  x
}

make_genus_matrix <- function(abund, tax, id_col) {
  idx <- match(abund[[id_col]], tax[[id_col]])
  if (anyNA(idx)) stop(paste0(id_col, "存在未匹配的分类注释。"))
  dat <- data.frame(
    Genus = clean_genus(tax$Genus[idx]),
    abund[, sample_ids, drop = FALSE],
    check.names = FALSE
  )
  dat <- dat[!is.na(dat$Genus), , drop = FALSE] |>
    dplyr::group_by(Genus) |>
    dplyr::summarise(dplyr::across(dplyr::all_of(sample_ids), sum),
                     .groups = "drop")
  mat <- t(as.matrix(dat[, sample_ids, drop = FALSE]))
  storage.mode(mat) <- "numeric"
  rownames(mat) <- sample_ids
  colnames(mat) <- dat$Genus
  if (any(rowSums(mat) == 0)) stop("存在属水平总丰度为0的样品。")
  mat
}

edna_mat  <- make_genus_matrix(edna_abund, edna_tax, "ASV_ID")
micro_mat <- make_genus_matrix(micro_abund, micro_tax, "Taxon_ID")

# 相对丰度平方根转换 + Bray-Curtis + PCoA
run_bray <- function(mat, method_name) {
  rel_sqrt <- sqrt(vegan::decostand(mat, method = "total"))
  d <- vegan::vegdist(rel_sqrt, method = "bray")
  ord <- ape::pcoa(d, correction = "lingoes")
  coords <- if (!is.null(ord$vectors.cor)) ord$vectors.cor else ord$vectors
  eig <- if ("Corr_eig" %in% names(ord$values)) {
    ord$values$Corr_eig
  } else {
    ord$values$Eigenvalues
  }
  explained <- 100 * eig[1:2] / sum(eig[eig > 0], na.rm = TRUE)
  scores <- data.frame(
    ID = rownames(coords), Axis1 = coords[, 1], Axis2 = coords[, 2],
    Method = method_name, check.names = FALSE
  ) |>
    dplyr::left_join(group, by = "ID")
  list(distance = d, scores = scores, explained = explained)
}

edna  <- run_bray(edna_mat, "eDNA")
micro <- run_bray(micro_mat, "Morphological identification")

# 3. PERMANOVA及同一采样时期内站点间两两Bray-Curtis距离
run_method_stats <- function(obj, method_name) {
  permanova <- vegan::adonis2(
    obj$distance ~ Period, data = group,
    permutations = 9999, strata = group$site
  )
  p1 <- as.data.frame(permanova)[1, , drop = FALSE]
  summary <- data.frame(
    Method = method_name,
    Analysis = "PERMANOVA",
    R2 = p1$R2,
    F_value = p1$F,
    P_value = p1[["Pr(>F)"]]
  )
  list(summary = summary)
}

edna_stats <- run_method_stats(edna, "eDNA")
micro_stats <- run_method_stats(micro, "Morphological identification")

# 8个采样时期的两两PERMANOVA（每种方法28组比较）
# 检验对象是原始群落距离矩阵，不用于生成图5b/e的显著性字母。
# 置换按采样站位限制，并同时输出原始P值和BH校正后的P值。
run_pairwise_permanova <- function(obj, method_name) {
  dmat <- as.matrix(obj$distance)
  period_pairs <- utils::combn(period_levels, 2, simplify = FALSE)
  out <- lapply(period_pairs, function(pp) {
    meta <- group[group$Period %in% pp, , drop = FALSE]
    meta$Period <- droplevels(meta$Period)
    ids <- meta$ID
    dsub <- stats::as.dist(dmat[ids, ids, drop = FALSE])
    fit <- vegan::adonis2(
      dsub ~ Period,
      data = meta,
      permutations = 9999,
      strata = meta$site
    )
    z <- as.data.frame(fit)[1, , drop = FALSE]
    data.frame(
      Method = method_name,
      Period1 = pp[1],
      Period2 = pp[2],
      Comparison = paste(pp, collapse = " vs "),
      N_samples = nrow(meta),
      F_value = z$F,
      R2 = z$R2,
      P_raw = z[["Pr(>F)"]],
      stringsAsFactors = FALSE
    )
  })
  dplyr::bind_rows(out) |>
    dplyr::mutate(
      P_adjust_BH = stats::p.adjust(P_raw, method = "BH"),
      Significant_raw = ifelse(P_raw < 0.05, "Yes", "No"),
      Significant_BH = ifelse(P_adjust_BH < 0.05, "Yes", "No")
    )
}

pairwise_permanova <- dplyr::bind_rows(
  run_pairwise_permanova(edna, "eDNA"),
  run_pairwise_permanova(micro, "Morphological identification")
)

pairwise_permanova_summary <- pairwise_permanova |>
  dplyr::group_by(Method) |>
  dplyr::summarise(
    N_comparisons = dplyr::n(),
    N_significant_raw = sum(P_raw < 0.05),
    N_significant_BH = sum(P_adjust_BH < 0.05),
    .groups = "drop"
  )

# 每个时期内，对不同站点样品进行两两组合并提取Bray-Curtis相异度。
# 5个站点产生10个距离；2014年冬季4个样品产生6个距离。
make_within_period_data <- function(obj, method_name) {
  dmat <- as.matrix(obj$distance)
  out <- list()
  k <- 1
  for (per in period_levels) {
    z <- group[group$Period == per, , drop = FALSE]
    if (nrow(z) < 2) next
    cmb <- utils::combn(z$ID, 2)
    for (j in seq_len(ncol(cmb))) {
      id1 <- cmb[1, j]
      id2 <- cmb[2, j]
      out[[k]] <- data.frame(
        Method = method_name,
        Sample1 = id1,
        Sample2 = id2,
        Site1 = group$site[match(id1, group$ID)],
        Site2 = group$site[match(id2, group$ID)],
        Year = as.character(z$Year[1]),
        Season = as.character(z$Season[1]),
        Period = per,
        Bray_Curtis_dissimilarity = dmat[id1, id2],
        stringsAsFactors = FALSE
      )
      k <- k + 1
    }
  }
  dplyr::bind_rows(out) |>
    dplyr::mutate(
      Year = factor(Year, levels = c("2014", "2015")),
      Season = factor(Season, levels = season_levels),
      Period = factor(Period, levels = period_levels)
    )
}

boxplot_data <- dplyr::bind_rows(
  make_within_period_data(edna, "eDNA"),
  make_within_period_data(micro, "Morphological identification")
)

# 每个年份内比较四季；显著时进行Tukey HSD并生成字母
make_letters <- function(dat) {
  out <- list()
  for (yr in c("2014", "2015")) {
    z <- dat[dat$Year == yr, , drop = FALSE]
    fit <- stats::aov(Bray_Curtis_dissimilarity ~ Season, data = z)
    p_anova <- summary(fit)[[1]][["Pr(>F)"]][1]
    if (is.na(p_anova) || p_anova >= 0.05) {
      letters <- setNames(rep("a", 4), season_levels)
    } else {
      tk <- stats::TukeyHSD(fit, "Season")$Season[, "p adj"]
      letters <- multcompView::multcompLetters(tk)$Letters
    }
    out[[yr]] <- data.frame(
      Year = yr,
      Season = factor(season_levels, levels = season_levels),
      Letter = unname(letters[season_levels]),
      ANOVA_P = p_anova
    ) |>
      dplyr::mutate(
        Period = factor(paste0(substr(Year, 3, 4), "_", Season),
                        levels = period_levels)
      )
  }
  dplyr::bind_rows(out)
}

letter_data <- split(boxplot_data, boxplot_data$Method) |>
  lapply(make_letters) |>
  dplyr::bind_rows(.id = "Method")

# 箱线图对应的描述统计、ANOVA和Tukey两两比较
season_summary <- boxplot_data |>
  dplyr::group_by(Method, Year, Season, Period) |>
  dplyr::summarise(
    N = dplyr::n(),
    Mean = mean(Bray_Curtis_dissimilarity),
    SD = stats::sd(Bray_Curtis_dissimilarity),
    Median = median(Bray_Curtis_dissimilarity),
    Q1 = stats::quantile(Bray_Curtis_dissimilarity, 0.25),
    Q3 = stats::quantile(Bray_Curtis_dissimilarity, 0.75),
    .groups = "drop"
  )

run_season_tests <- function(dat) {
  anova_out <- list()
  tukey_out <- list()
  k <- 1
  for (m in unique(dat$Method)) {
    for (yr in c("2014", "2015")) {
      z <- dat[dat$Method == m & dat$Year == yr, , drop = FALSE]
      fit <- stats::aov(Bray_Curtis_dissimilarity ~ Season, data = z)
      a <- summary(fit)[[1]][1, , drop = FALSE]
      anova_out[[k]] <- data.frame(
        Method = m, Year = yr,
        Df = a[["Df"]], F_value = a[["F value"]],
        P_value = a[["Pr(>F)"]]
      )
      tk <- as.data.frame(stats::TukeyHSD(fit, "Season")$Season)
      tk$Comparison <- rownames(tk)
      rownames(tk) <- NULL
      tukey_out[[k]] <- data.frame(
        Method = m, Year = yr, Comparison = tk$Comparison,
        Difference = tk$diff, CI_low = tk$lwr, CI_high = tk$upr,
        Adjusted_P = tk$`p adj`
      )
      k <- k + 1
    }
  }
  list(
    ANOVA = dplyr::bind_rows(anova_out),
    Tukey = dplyr::bind_rows(tukey_out)
  )
}

season_tests <- run_season_tests(boxplot_data)

label_data <- boxplot_data |>
  dplyr::group_by(Method, Period) |>
  dplyr::summarise(y_max = max(Bray_Curtis_dissimilarity), .groups = "drop") |>
  dplyr::left_join(letter_data, by = c("Method", "Period")) |>
  dplyr::group_by(Method) |>
  dplyr::mutate(
    panel_range = max(y_max) - min(y_max),
    Label_y = y_max + ifelse(panel_range == 0, 0.01, 0.08 * panel_range)
  ) |>
  dplyr::ungroup()

# 4. 同一站位不同时期样品的距离衰减
make_decay_data <- function(obj, method_name) {
  dmat <- as.matrix(obj$distance)
  out <- list()
  k <- 1
  for (st in unique(group$site)) {
    z <- group[group$site == st, , drop = FALSE]
    if (nrow(z) < 2) next
    cmb <- utils::combn(seq_len(nrow(z)), 2)
    for (j in seq_len(ncol(cmb))) {
      a <- z[cmb[1, j], ]
      b <- z[cmb[2, j], ]
      out[[k]] <- data.frame(
        Method = method_name, site = st,
        Sample1 = a$ID, Sample2 = b$ID,
        Interval = abs(a$Event_index - b$Event_index),
        Bray_Curtis_dissimilarity = dmat[a$ID, b$ID]
      )
      k <- k + 1
    }
  }
  dplyr::bind_rows(out)
}

decay_data <- dplyr::bind_rows(
  make_decay_data(edna, "eDNA"),
  make_decay_data(micro, "Morphological identification")
)

decay_stats <- split(decay_data, decay_data$Method) |>
  lapply(function(z) {
    fit <- stats::lm(Bray_Curtis_dissimilarity ~ Interval, data = z)
    sm <- summary(fit)
    data.frame(
      N_pairs = nrow(z),
      Intercept = unname(coef(fit)[1]),
      Slope = unname(coef(fit)[2]),
      R_squared = sm$r.squared,
      P_value = coef(sm)[2, 4]
    )
  }) |>
  dplyr::bind_rows(.id = "Method")

# 5. 作图
season_colors <- c(Spring = "#18A77B", Summer = "#F47C48",
                   Autumn = "#AC55E8", Winter = "#12BCE5")
year_shapes <- c("2014" = 16, "2015" = 17)
year_fill <- c("2014" = "#F6DDE2", "2015" = "#CFE7F3")
year_line <- c("2014" = "#B96E7A", "2015" = "#568EAE")
period_labels <- c(
  "14_Spring" = "Spr\n2014", "14_Summer" = "Sum\n2014",
  "14_Autumn" = "Aut\n2014", "14_Winter" = "Win\n2014",
  "15_Spring" = "Spr\n2015", "15_Summer" = "Sum\n2015",
  "15_Autumn" = "Aut\n2015", "15_Winter" = "Win\n2015"
)

base_theme <- ggplot2::theme_classic(base_size = 14) +
  ggplot2::theme(
    axis.title = ggplot2::element_text(size = 14, face = "bold"),
    axis.text = ggplot2::element_text(size = 12, color = "black", face = "bold"),
    plot.title = ggplot2::element_text(size = 15, hjust = 0, face = "bold"),
    plot.margin = ggplot2::margin(10, 12, 10, 10),
    legend.title = ggplot2::element_text(size = 12, face = "bold"),
    legend.text = ggplot2::element_text(size = 11, face = "bold")
  )

make_pcoa_plot <- function(obj, title_text) {
  ggplot2::ggplot(obj$scores,
                  ggplot2::aes(Axis1, Axis2, color = Season, shape = Year)) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2,
                        linewidth = 0.4, color = "grey50") +
    ggplot2::geom_vline(xintercept = 0, linetype = 2,
                        linewidth = 0.4, color = "grey50") +
    ggplot2::geom_point(size = 3, alpha = 0.9) +
    ggplot2::scale_color_manual(values = season_colors) +
    ggplot2::scale_shape_manual(values = year_shapes) +
    ggplot2::labs(
      title = title_text,
      x = sprintf("PCoA1 (%.2f%%)", obj$explained[1]),
      y = sprintf("PCoA2 (%.2f%%)", obj$explained[2]),
      color = "Season", shape = "Year"
    ) + base_theme + ggplot2::theme(legend.position = "bottom")
}

make_boxplot <- function(method_name, title_text) {
  dat <- dplyr::filter(boxplot_data, Method == method_name)
  lab <- dplyr::filter(label_data, Method == method_name)
  ggplot2::ggplot(dat, ggplot2::aes(Period, Bray_Curtis_dissimilarity,
                                    fill = Year, color = Year)) +
    ggplot2::geom_boxplot(
      width = 0.66, linewidth = 0.48,
      outlier.shape = 16, outlier.size = 1.1,
      outlier.stroke = 0, outlier.alpha = 0.4
    ) +
    ggplot2::geom_text(
      data = lab, ggplot2::aes(Period, Label_y, label = Letter),
      inherit.aes = FALSE, size = 5.2, fontface = "bold"
    ) +
    ggplot2::scale_fill_manual(values = year_fill) +
    ggplot2::scale_color_manual(values = year_line) +
    ggplot2::scale_x_discrete(labels = period_labels) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0.03, 0.16))
    ) +
    ggplot2::labs(title = title_text, x = NULL,
                  y = "Bray-Curtis dissimilarity") +
    base_theme + ggplot2::theme(legend.position = "none")
}

make_decay_plot <- function(method_name, title_text) {
  dat <- dplyr::filter(decay_data, Method == method_name)
  st <- dplyr::filter(decay_stats, Method == method_name)
  fit <- stats::lm(Bray_Curtis_dissimilarity ~ Interval, data = dat)
  pred <- data.frame(Interval = seq(0, 7, length.out = 200))
  ci <- as.data.frame(stats::predict(fit, newdata = pred,
                                     interval = "confidence"))
  pred <- cbind(pred, ci)
  
  # 仿照同门Figure 5：纵轴只覆盖拟合线及95%置信区间附近，
  # 不从0开始，也不受全部原始散点极值控制。
  y_min <- floor(min(pred$lwr, na.rm = TRUE) / 0.05) * 0.05
  y_max <- ceiling(max(pred$upr, na.rm = TRUE) / 0.05) * 0.05
  if (y_max <= y_min) y_max <- y_min + 0.05
  
  stat_label <- if (st$P_value < 0.001) {
    sprintf("R² = %.2f, P < 0.001", st$R_squared)
  } else {
    sprintf("R² = %.2f, P = %.3f", st$R_squared, st$P_value)
  }
  ggplot2::ggplot(pred, ggplot2::aes(Interval, fit)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = lwr, ymax = upr),
      fill = "#F5D99D", alpha = 0.65
    ) +
    ggplot2::geom_line(color = "#E19100", linewidth = 0.9) +
    ggplot2::annotate("text", x = Inf, y = Inf, label = stat_label,
                      hjust = 1.05, vjust = 1.25, size = 5.2,
                      fontface = "bold") +
    ggplot2::scale_x_continuous(
      breaks = seq(0, 6, by = 2), limits = c(0, 7),
      expand = ggplot2::expansion(mult = c(0, 0))
    ) +
    ggplot2::scale_y_continuous(
      breaks = seq(y_min, y_max, by = 0.05),
      limits = c(y_min, y_max),
      labels = function(x) sprintf("%.2f", x),
      expand = ggplot2::expansion(mult = c(0, 0))
    ) +
    ggplot2::labs(
      title = title_text,
      x = "Sampling time interval (season)",
      y = "Bray-Curtis dissimilarity"
    ) + base_theme
}

p5a <- make_pcoa_plot(edna, "(a) eDNA")
p5b <- make_boxplot("eDNA", "(b) eDNA")
p5c <- make_decay_plot("eDNA", "(c) eDNA")
p5d <- make_pcoa_plot(micro, "(d) Morphological identification")
p5e <- make_boxplot("Morphological identification",
                    "(e) Morphological identification")
p5f <- make_decay_plot("Morphological identification",
                       "(f) Morphological identification")

figure5 <- patchwork::wrap_plots(
  p5a, p5b, p5c, p5d, p5e, p5f,
  ncol = 3, nrow = 2, guides = "collect"
) & ggplot2::theme(
  legend.position = "bottom",
  legend.box = "horizontal",
  legend.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0),
  legend.box.margin = ggplot2::margin(0, 0, 0, 0),
  legend.spacing.x = grid::unit(3, "pt"),
  legend.title = ggplot2::element_text(size = 12, face = "bold"),
  legend.text = ggplot2::element_text(size = 11, face = "bold")
)

# 与Figure 4统一画布尺寸和长宽比例，使3列×2行面板大小一致
figure_width  <- 14.5
figure_height <- 9.2

ggplot2::ggsave(
  file.path(output_dir, "Figure5.png"), figure5,
  width = figure_width, height = figure_height,
  units = "in", dpi = 400, bg = "white"
)

ggplot2::ggsave(
  file.path(output_dir, "Figure5.tif"), figure5,
  device = "tiff",
  width = figure_width, height = figure_height,
  units = "in", dpi = 600, compression = "lzw", bg = "white"
)

ggplot2::ggsave(
  file.path(output_dir, "Figure5.pdf"), figure5,
  device = grDevices::cairo_pdf,
  width = figure_width, height = figure_height,
  units = "in", bg = "white"
)
# 6. 仅输出三个面板的作图数据和一张统计结果表
pcoa_data <- dplyr::bind_rows(edna$scores, micro$scores)
pcoa_explained <- data.frame(
  Method = c("eDNA", "Morphological identification"),
  PCoA1_percent = c(edna$explained[1], micro$explained[1]),
  PCoA2_percent = c(edna$explained[2], micro$explained[2]),
  Cumulative_percent = c(sum(edna$explained[1:2]),
                         sum(micro$explained[1:2]))
)
statistics <- dplyr::bind_rows(
  dplyr::bind_rows(edna_stats$summary, micro_stats$summary) |>
    dplyr::transmute(Method, Analysis, Group = "All periods",
                     R2, F_value, P_value, Slope = NA_real_,
                     Letter = NA_character_),
  decay_stats |>
    dplyr::transmute(Method, Analysis = "Temporal decay regression",
                     Group = "All pairs", R2 = R_squared,
                     F_value = NA_real_, P_value, Slope,
                     Letter = NA_character_),
  letter_data |>
    dplyr::transmute(Method, Analysis = "One-way ANOVA + Tukey HSD",
                     Group = paste(Year, Season, sep = "_"),
                     R2 = NA_real_, F_value = NA_real_,
                     P_value = ANOVA_P, Slope = NA_real_, Letter)
)

openxlsx::write.xlsx(
  list(
    PCoA_data = pcoa_data,
    PCoA_explained = pcoa_explained,
    Boxplot_data = boxplot_data,
    Season_summary = season_summary,
    ANOVA = season_tests$ANOVA,
    Tukey_HSD = season_tests$Tukey,
    Pairwise_PERMANOVA = pairwise_permanova,
    Pairwise_PERM_summary = pairwise_permanova_summary,
    Decay_data = decay_data,
    Key_statistics = statistics
  ),
  file = file.path(output_dir, "Figure5_data.xlsx"),
  overwrite = TRUE
)

print(figure5)
cat("\n完成，输出文件：\n",
    file.path(output_dir, "Figure5.png"), "\n",
    file.path(output_dir, "Figure5.tif"), "\n",
    file.path(output_dir, "Figure5.pdf"), "\n",
    file.path(output_dir, "Figure5_data.xlsx"), "\n")
