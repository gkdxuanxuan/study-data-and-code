# ============================================================
# 图2：eDNA与镜检在门、纲、目、科、属水平的共有/特有类群
# ============================================================

setwd("D:/A_师姐文章/原始处理数据/111")

# 创建图2输出文件夹
out_dir <- "图2"
if (!dir.exists(out_dir)) {
  dir.create(out_dir)
}

library(tidyverse)

edna_abund <- read.csv("2.D_ASV丰度表_856reads抽平_329ASV_Richness显著.csv", check.names = FALSE)
edna_tax   <- read.csv("phyto_tax_filtered.csv", check.names = FALSE)
micro_cell <- read.csv("2.D原始镜检-细胞丰度.csv", check.names = FALSE)
micro_tax  <- read.csv("2.DTAXA-JJ-ZH_规范化.csv", check.names = FALSE)

names(edna_tax)[1]   <- "ASV_ID"
names(micro_cell)[1] <- "Taxon_ID"
names(micro_tax)[1]  <- "Taxon_ID"

ranks <- c("Phylum", "Class", "Order", "Family", "Genus")

clean_tax <- function(x) {
  x <- trimws(as.character(x))
  x <- sub("^[dkpcofgs]__", "", x, ignore.case = TRUE)
  x[x == "" | tolower(x) %in% c(
    "na", "norank", "unknown", "unclassified", "unidentified"
  )] <- NA
  x
}

edna_tax[ranks]  <- lapply(edna_tax[ranks], clean_tax)
micro_tax[ranks] <- lapply(micro_tax[ranks], clean_tax)

edna <- edna_abund |>
  left_join(edna_tax |> select(ASV_ID, all_of(ranks)), by = "ASV_ID")

micro <- micro_cell |>
  left_join(micro_tax |> select(Taxon_ID, all_of(ranks)), by = "Taxon_ID")

stopifnot(nrow(edna) == 329, nrow(micro) == 35)


fig2_data <- map_dfr(ranks, function(r) {
  e <- unique(na.omit(edna[[r]]))
  m <- unique(na.omit(micro[[r]]))
  
  tibble(
    Rank = r,
    `Microscopy only` = length(setdiff(m, e)),
    Shared = length(intersect(e, m)),
    `eDNA only` = length(setdiff(e, m)),
    `eDNA total` = length(e),
    `Microscopy total` = length(m)
  )
})


fig2_taxa <- map_dfr(ranks, function(r) {
  e <- unique(na.omit(edna[[r]]))
  m <- unique(na.omit(micro[[r]]))
  
  bind_rows(
    tibble(Rank = r,
           Category = "eDNA only",
           Taxon = sort(setdiff(e, m))),
    
    tibble(Rank = r,
           Category = "Shared",
           Taxon = sort(intersect(e, m))),
    
    tibble(Rank = r,
           Category = "Microscopy only",
           Taxon = sort(setdiff(m, e)))
  )
})


fig2_long <- fig2_data |>
  select(Rank, `Microscopy only`, Shared, `eDNA only`) |>
  pivot_longer(-Rank,
               names_to = "Category",
               values_to = "Number") |>
  mutate(
    Rank = factor(Rank, levels = ranks),
    Category = factor(Category,
                      levels = c(
                        "Microscopy only",
                        "Shared",
                        "eDNA only"))
  )


p2 <- ggplot(fig2_long,
             aes(Rank, Number, fill = Category)) +
  geom_col(width = 0.78,
           colour = "white",
           linewidth = 0.3) +
  geom_text(
    aes(label = ifelse(Number == 0, "", Number)),
    position = position_stack(vjust = 0.5),
    size = 4,
    fontface = "bold",
    colour = "#263238"
  ) +
  scale_fill_manual(values = c(
    "Microscopy only" = "#A8D480",
    "Shared" = "#2374A5",
    "eDNA only" = "#A6CDE3"
  )) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.08))
  ) +
  labs(
    x = NULL,
    y = "Number of detected taxa",
    fill = NULL
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(size = 12),
    panel.grid.major.y = element_line(
      colour = "#E5E7EB",
      linewidth = 0.4
    )
  )


# ============================
# 输出到“图2”文件夹
# ============================

ggsave(
  file.path(out_dir, "Figure2_updated_R.png"),
  p2,
  width = 7.3,
  height = 6.2,
  units = "in",
  dpi = 400,
  bg = "white"
)

ggsave(
  file.path(out_dir, "Figure2_updated_R.tiff"),
  p2,
  width = 7.3,
  height = 6.2,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

ggsave(
  file.path(out_dir, "Figure2_updated_R.pdf"),
  p2,
  width = 7.3,
  height = 6.2,
  units = "in",
  device = cairo_pdf,
  bg = "white"
)

write.csv(
  fig2_data,
  file.path(out_dir, "Figure2_plot_data_R.csv"),
  row.names = FALSE
)

write.csv(
  fig2_taxa,
  file.path(out_dir, "Figure2_shared_unique_taxa_R.csv"),
  row.names = FALSE
)


print(fig2_data)