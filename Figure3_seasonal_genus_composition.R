# ============================================================
# 图3：eDNA、镜检细胞密度和镜检生物量的季节属组成
# 说明：三个面板分别筛选总体丰度Top 10属，其余合并为Others
# 不补充norank；只使用Genus列；不修改原始taxa文件
# ============================================================

setwd("D:/A_师姐文章/原始处理数据/111")
# 创建图3输出文件夹
out_dir <- "图3"
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}
# 首次运行若缺少tidyverse，请先运行：install.packages("tidyverse")
library(tidyverse)

# ---------- 1. 读取数据 ----------
edna_abund <- read.csv("2.D_ASV丰度表_856reads抽平_329ASV_Richness显著.csv", check.names = FALSE)
edna_tax   <- read.csv("2.Dphyto_tax_filtered.csv", check.names = FALSE)
micro_cell <- read.csv("2.D原始镜检-细胞丰度.csv", check.names = FALSE)
micro_bio  <- read.csv("2.D原始镜检生物量.csv", check.names = FALSE)
micro_tax  <- read.csv("2.DTAXA-JJ-ZH_规范化.csv", check.names = FALSE)
group      <- read.csv("3.group.csv", check.names = FALSE)

# ---------- 2. 统一编号列 ----------
names(edna_abund)[1] <- "ASV_ID"
names(edna_tax)[1]   <- "ASV_ID"
names(micro_cell)[1] <- "Taxon_ID"
names(micro_bio)[1]  <- "Taxon_ID"
names(micro_tax)[1]  <- "Taxon_ID"

edna_abund$ASV_ID <- trimws(as.character(edna_abund$ASV_ID))
edna_tax$ASV_ID   <- trimws(as.character(edna_tax$ASV_ID))
micro_cell$Taxon_ID <- trimws(as.character(micro_cell$Taxon_ID))
micro_bio$Taxon_ID  <- trimws(as.character(micro_bio$Taxon_ID))
micro_tax$Taxon_ID  <- trimws(as.character(micro_tax$Taxon_ID))
group$ID <- trimws(as.character(group$ID))

# ---------- 3. 正确提取属名 ----------
clean_genus <- function(x) {
  x <- trimws(as.character(x))
  x <- sub("^g__", "", x, ignore.case = TRUE)
  x <- trimws(x)
  invalid <- is.na(x) | x == "" | tolower(x) %in% c(
    "na", "n/a", "norank", "no_rank", "unknown",
    "unclassified", "unidentified", "uncultured", "not assigned"
  )
  x[invalid] <- NA
  x
}

edna_tax$Genus  <- clean_genus(edna_tax$Genus)
micro_tax$Genus <- clean_genus(micro_tax$Genus)

# ---------- 4. 匹配丰度和属注释 ----------
stopifnot(
  anyDuplicated(edna_tax$ASV_ID) == 0,
  anyDuplicated(micro_tax$Taxon_ID) == 0
)

edna <- edna_abund %>%
  left_join(edna_tax %>% select(ASV_ID, Genus), by = "ASV_ID")
cell <- micro_cell %>%
  left_join(micro_tax %>% select(Taxon_ID, Genus), by = "Taxon_ID")
bio <- micro_bio %>%
  left_join(micro_tax %>% select(Taxon_ID, Genus), by = "Taxon_ID")

sample_cols <- group$ID

stopifnot(
  nrow(edna) == 329,
  nrow(cell) == 35,
  nrow(bio) == 35,
  length(sample_cols) == 39,
  anyDuplicated(sample_cols) == 0,
  all(sample_cols %in% names(edna)),
  all(sample_cols %in% names(cell)),
  all(sample_cols %in% names(bio))
)

# 本图只统计能够注释到属水平的记录；不将未定属补为norank
message("eDNA中未注释到属、将在属组成图中排除的ASV数：", sum(is.na(edna$Genus)))
message("镜检中未注释到属、将在属组成图中排除的分类单元数：", sum(is.na(cell$Genus)))

# 确保丰度列全部为数值型
edna[sample_cols] <- lapply(edna[sample_cols], as.numeric)
cell[sample_cols] <- lapply(cell[sample_cols], as.numeric)
bio[sample_cols]  <- lapply(bio[sample_cols], as.numeric)

# ---------- 5. 季节顺序 ----------
season_order <- c(
  "14_Spring", "14_Summer", "14_Autumn", "14_Winter",
  "15_Spring", "15_Summer", "15_Autumn", "15_Winter"
)
season_labels <- c("Spr", "Sum", "Aut", "Win", "Spr", "Sum", "Aut", "Win")

if (!all(group$season %in% season_order)) {
  stop("group文件中的season名称与代码预设不一致，请检查。")
}

if (any(is.na(group$season))) {
  stop("group文件中的season存在缺失值，请检查。")
}

# ---------- 6. 分面板计算Top 10和季节占比 ----------
calculate_panel <- function(data, panel_name, method_name, measure_name) {
  genus_sample <- data %>%
    filter(!is.na(Genus)) %>%
    group_by(Genus) %>%
    summarise(across(all_of(sample_cols), ~sum(.x, na.rm = TRUE)),
              .groups = "drop")
  
  top10 <- genus_sample %>%
    mutate(Total_abundance = rowSums(across(all_of(sample_cols)))) %>%
    slice_max(Total_abundance, n = 10, with_ties = FALSE) %>%
    arrange(desc(Total_abundance)) %>%
    pull(Genus)
  
  all_genera <- genus_sample %>%
    pivot_longer(all_of(sample_cols), names_to = "ID", values_to = "Abundance") %>%
    left_join(group %>% select(ID, season), by = "ID") %>%
    group_by(season, Genus) %>%
    summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
    mutate(
      Is_top10 = Genus %in% top10,
      Plot_taxon = ifelse(Is_top10, Genus, "Others"),
      Panel = panel_name,
      Method = method_name,
      Measure = measure_name
    )
  
  plot_data <- all_genera %>%
    group_by(Panel, Method, Measure, season, Plot_taxon) %>%
    summarise(Abundance = sum(Abundance), .groups = "drop") %>%
    group_by(Panel, season) %>%
    mutate(
      Season_total = sum(Abundance),
      Relative_abundance_percent = if_else(
        Season_total > 0,
        Abundance / Season_total * 100,
        NA_real_
      )
    ) %>%
    ungroup()
  
  list(plot = plot_data, all = all_genera, top10 = top10)
}

edna_result <- calculate_panel(
  edna, "(a) eDNA", "eDNA", "Read relative abundance"
)
cell_result <- calculate_panel(
  cell, "(b) Morphological cell density", "Microscopy",
  "Cell-density relative abundance"
)
bio_result <- calculate_panel(
  bio, "(c) Morphological biomass", "Microscopy",
  "Biomass relative abundance"
)

# ---------- 7. 合并数据并统一颜色 ----------
fig3_data <- bind_rows(
  edna_result$plot, cell_result$plot, bio_result$plot
)
fig3_all_genera <- bind_rows(
  edna_result$all, cell_result$all, bio_result$all
)

fig3_top10 <- bind_rows(
  tibble(Panel = "(a) eDNA", Rank = 1:10, Genus = edna_result$top10),
  tibble(Panel = "(b) Morphological cell density", Rank = 1:10,
         Genus = cell_result$top10),
  tibble(Panel = "(c) Morphological biomass", Rank = 1:10,
         Genus = bio_result$top10)
)

panel_order <- c(
  "(a) eDNA",
  "(b) Morphological cell density",
  "(c) Morphological biomass"
)
all_top_genera <- sort(unique(fig3_top10$Genus))
legend_order <- c(all_top_genera, "Others")

# 使用固定配色，保证同一属在三个面板中的颜色完全一致
reference_colors <- c(
  "Aulacoseira"    = "#1F77B4",
  "Ceratium"       = "#AEC7E8",
  "Chroomonas"     = "#FF7F0E",
  "Cocconeis"      = "#FFBB78",
  "Coelastrum"     = "#2CA02C",
  "Cryptomonas"    = "#98DF8A",
  "Cyclostephanos" = "#D62728",
  "Cyclotella"     = "#FF9896",
  "Dinobryon"      = "#9467BD",
  "Eudorina"       = "#C5B0D5",
  "Fragilaria"     = "#8C564B",
  "Lagerheimiella" = "#C49C94",
  "Mougeotia"      = "#E377C2",
  "Ochromonas"     = "#F7B6D2",
  "Peridinium"     = "#4B4E9A",  # 改为深蓝紫色，与灰色系区分
  "Scenedesmus"    = "#C7C7C7",
  "Stephanodiscus" = "#BCBD22",
  "Tetraedron"     = "#DBDB8D",
  "Trachydiscus"   = "#17BECF",
  "Volvox"         = "#9EDAE5",
  "Others"         = "#8C8C8C"
)
# 若以后Top 10名单出现参考图中没有的新属，自动补充备用颜色
missing_color_genera <- setdiff(all_top_genera, names(reference_colors))
if (length(missing_color_genera) > 0) {
  fallback_colors <- setNames(
    scales::hue_pal()(length(missing_color_genera)),
    missing_color_genera
  )
  reference_colors <- c(reference_colors, fallback_colors)
}

genus_colors <- reference_colors[legend_order]

fig3_data <- fig3_data %>%
  mutate(
    Panel = factor(Panel, levels = panel_order),
    season = factor(season, levels = season_order),
    Season_index = match(as.character(season), season_order),
    Plot_taxon = factor(Plot_taxon, levels = legend_order)
  )

# 保存用于复核和论文附件整理的数据
write.csv(
  fig3_data,
  file.path(out_dir, "Figure3_plot_data_R.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
write.csv(
  fig3_all_genera,
  file.path(out_dir, "Figure3_all_genera_R.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
write.csv(
  fig3_top10,
  file.path(out_dir, "Figure3_top10_R.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# ---------- 8. 使用facet一次完成三个面板 ----------
p3 <- ggplot(
  fig3_data,
  aes(
    x = Season_index,
    y = Relative_abundance_percent,
    fill = Plot_taxon
  )
) +
  geom_col(
    width = 0.82,
    position = position_stack(reverse = TRUE)
  ) +
  annotate(
    "rect",
    xmin = 0.5, xmax = 4.5, ymin = -16, ymax = -8,
    fill = "#F8E9ED", colour = NA
  ) +
  annotate(
    "rect",
    xmin = 4.5, xmax = 8.5, ymin = -16, ymax = -8,
    fill = "#DCEEF8", colour = NA
  ) +
  # 年份文字：放大并加粗
  annotate(
    "text",
    x = 2.5, y = -12, label = "2014",
    size = 5, fontface = "bold"
  ) +
  annotate(
    "text",
    x = 6.5, y = -12, label = "2015",
    size = 5, fontface = "bold"
  ) +
  facet_wrap(~Panel, nrow = 1) +
  scale_fill_manual(
    values = genus_colors,
    limits = legend_order,
    breaks = legend_order,
    drop = FALSE
  ) +
  scale_x_continuous(
    breaks = 1:8,
    labels = season_labels,
    expand = expansion(add = 0.45)
  ) +
  scale_y_continuous(
    breaks = seq(0, 100, 25),
    expand = expansion(mult = c(0, 0.02))
  ) +
  coord_cartesian(ylim = c(0, 100), clip = "off") +
  labs(
    x = NULL,
    y = "Relative abundance (%)",
    fill = NULL
  ) +
  guides(
    fill = guide_legend(ncol = 1, byrow = TRUE)
  ) +
  theme_classic(base_size = 16) +
  theme(
    # 所有主题文字默认加粗
    text = element_text(face = "bold", colour = "black"),
    
    strip.background = element_blank(),
    
    # 三个面板的标题
    strip.text = element_text(
      size = 16, face = "bold", hjust = 0
    ),
    
    panel.spacing = grid::unit(0.7, "cm"),
    
    # 右侧图例
    legend.position = "right",
    legend.direction = "vertical",
    legend.box.just = "left",
    legend.text = element_text(
      size = 13, face = "bold.italic"
    ),
    legend.key.height = grid::unit(0.48, "cm"),
    legend.key.width = grid::unit(0.48, "cm"),
    
    # 横、纵轴刻度
    axis.text.x = element_text(
      size = 12, face = "bold", colour = "black"
    ),
    axis.text.y = element_text(
      size = 14, face = "bold", colour = "black"
    ),
    
    # 纵轴标题
    axis.title.y = element_text(
      size = 17, face = "bold",
      margin = margin(r = 10)
    ),
    
    plot.margin = margin(10, 12, 50, 10)
  )

print(p3)

# ---------- 9. 保存图3 ----------
ggsave(file.path(out_dir, "Figure3_updated_R.png"), p3,
       width = 14, height = 5.5,
       units = "in", dpi = 400, bg = "white")
ggsave(file.path(out_dir, "Figure3_updated_R.tiff"), p3,
       width = 14, height = 5.5,
       units = "in", dpi = 600, compression = "lzw", bg = "white")

message("图3及复核数据已输出至：", normalizePath(out_dir))

if (capabilities("cairo")) {
  ggsave("Figure3_updated_R.pdf", p3, width = 14, height = 5.5,
         units = "in", device = cairo_pdf, bg = "white")
} else {
  ggsave("Figure3_updated_R.pdf", p3, width = 16, height = 7.5,
         units = "in", bg = "white")
}

# ---------- 10. 输出图3对应数据 ----------
write.csv(fig3_data, "Figure3_plot_data_R.csv", row.names = FALSE)
write.csv(fig3_top10, "Figure3_top10_genera_R.csv", row.names = FALSE)
write.csv(
  fig3_all_genera %>% filter(!Is_top10),
  "Figure3_Others_components_R.csv", row.names = FALSE
)

fig3_percent <- fig3_data %>%
  select(Panel, Method, Measure, season, Plot_taxon,
         Relative_abundance_percent) %>%
  pivot_wider(names_from = Plot_taxon,
              values_from = Relative_abundance_percent,
              values_fill = 0)
write.csv(fig3_percent, "Figure3_seasonal_relative_abundance_R.csv",
          row.names = FALSE)

fig3_dominant <- fig3_data %>%
  filter(Plot_taxon != "Others") %>%
  group_by(Panel, Method, Measure, season) %>%
  slice_max(Relative_abundance_percent, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(Panel, Method, Measure, season,
         Dominant_genus = Plot_taxon, Relative_abundance_percent)
write.csv(fig3_dominant, "Figure3_dominant_genus_R.csv", row.names = FALSE)

percentage_check <- fig3_data %>%
  group_by(Panel, season) %>%
  summarise(Total_percent = sum(Relative_abundance_percent), .groups = "drop")
write.csv(percentage_check, "Figure3_percentage_check_R.csv", row.names = FALSE)

cat("图3绘制完成。每个面板、每个季节的占比合计应为100：\n")
print(percentage_check)
cat("各面板Top 10属：\n")
print(fig3_top10)
