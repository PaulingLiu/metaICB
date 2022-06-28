
#--- library

library(Seurat)
library(tidyverse)
library(reticulate)
library(ComplexHeatmap)
library(pROC)
library(ggROC)

#--- load data

sce.sub <- readr::read_rds("/raid1/pauling/projects/07_pan.cancer.treatment/02.data/01.lung.nature/03.cd8.tcr.specificity.rds.gz")
sce.hn <- readr::read_rds("/home/pauling/projects/07_pan.cancer.treatment/02.data/02.HNSCC.nature/sce.tcr.rds.gz")

#--- treatment-naive and non-responsive tumors

sce.nrp <- sce.sub[,sce.sub$patient %in% c("MD01-004","MD043-011","NY016-014") | sce.sub$type != "MANA"]


filt.clones <- tibble(
  trb = sce.nrp$trb,
  type = sce.nrp$type
) %>%
  dplyr::group_by(trb, type) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::ungroup() %>%
  dplyr::filter((type == "MANA" & n > 3) | (type == "Viral" ))

filt.hn.clones <- tibble(trb = as.character(sce.hn$clone.id)) %>%
  dplyr::count(trb) %>%
  dplyr::filter(n > 3)

sce.nrp <- sce.nrp[,sce.nrp$trb %in% filt.clones$trb]
sce.hn.pro <- sce.hn[,sce.hn$clone.id %in% filt.hn.clones$trb]

sce.hn.pro$trb <- sce.hn.pro$clone.id
sce.hn.pro$type <- "tumor"
Idents(sce.hn.pro) <- "MANA"

comb.sce <- merge(sce.hn.pro, sce.nrp)
comb.sce <- NormalizeData(comb.sce)
comb.sce <- ScaleData(comb.sce, do.center = F, do.scale = F)
#markers.hn <- FindAllMarkers(comb.sce, test.use = "roc", only.pos = T)

tibble(
  cxcl13 = sce.hn.pro@assays$RNA@scale.data["CXCL13",],
  trb = sce.hn.pro$clone.id,
  tex = sce.hn.pro@assays$RNA@scale.data["ITGAE",] + sce.hn.pro@assays$RNA@scale.data["BATF",] + 
    sce.hn.pro@assays$RNA@scale.data["ENTPD1",] + sce.hn.pro@assays$RNA@scale.data["LAYN",],
) %>%
  dplyr::mutate(type = "tumor") -> hb.data

pda <- tibble(
  cxcl13 = sce.nrp@assays$RNA@scale.data["CXCL13",],
  trb = sce.nrp$trb,
  type = sce.nrp$type,
  tex = sce.nrp@assays$RNA@scale.data["BATF",] + sce.nrp@assays$RNA@scale.data["ITGAE",] + 
    sce.nrp@assays$RNA@scale.data["ENTPD1",] + sce.nrp@assays$RNA@scale.data["LAYN",]
) %>%
  dplyr::bind_rows(hb.data) %>%
  dplyr::group_by(trb, type) %>%
  dplyr::summarise(m = mean(cxcl13),  m2 = mean(tex), size = n()) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(desc(type))

pda %>%
  ggplot(aes(1+m, 1+m2)) +
  geom_point(aes(fill = type, color = type), size = 2.5, alpha = 0.8, shape = 21, stroke = 0.8) +
  geom_vline(xintercept = 1.155) +
  geom_hline(yintercept = 3.15) +
  scale_y_log10(lim = c(1,7), breaks = c(1,5)) +
  scale_x_log10(lim = c(0.85,6), breaks = c(1,5)) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 15,color="black"),
    axis.text = element_text(size = 12,color="black"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    axis.text.y = element_text(color="black"),
    axis.text.x = element_text(color="black")
  ) +
  labs(
    x = "CXCL13 expression",
    y = "Exhaustion score"
  ) +
  scale_fill_manual(values = c("#CD2626","#CD2626", "#00868B")) +
  scale_color_manual(values = c("black","#CD2626", "#00868B"))

tibble(
  cxcl13 = sce.nrp@assays$RNA@scale.data["CXCL13",],
  trb = sce.nrp$trb,
  type = sce.nrp$type,
  tex = sce.nrp@assays$RNA@scale.data["BATF",] + sce.nrp@assays$RNA@scale.data["ITGAE",] + 
    sce.nrp@assays$RNA@scale.data["ENTPD1",] + sce.nrp@assays$RNA@scale.data["LAYN",]
) %>%
  dplyr::bind_rows(hb.data) %>%
  dplyr::group_by(trb, type) %>%
  dplyr::summarise(m = mean(cxcl13),  m2 = mean(tex)) %>%
  ggplot(aes(1+m, 1+m2)) +
  geom_point(aes(shape = type), size = 2.5, stroke = 1) +
  geom_vline(xintercept = 1.155) +
  geom_hline(yintercept = 3.15) +
  scale_y_log10(lim = c(1,7), breaks = c(1,5)) +
  scale_x_log10(lim = c(0.85,6), breaks = c(1,5)) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 15,color="black"),
    axis.text = element_text(size = 12,color="black"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    axis.text.y = element_text(color="black"),
    axis.text.x = element_text(color="black")
  ) +
  labs(
    x = "CXCL13 expression",
    y = "Exhaustion score"
  ) +
  scale_shape_manual(values = c(1,1,3))

#--- auc store

rda <- tibble(
  type = comb.sce$type,
  trb = comb.sce$trb,
  CXCL13 = comb.sce@assays$RNA@scale.data["CXCL13",],CTLA4 = comb.sce@assays$RNA@scale.data["CTLA4",],
  PDCD1 = comb.sce@assays$RNA@scale.data["PDCD1",],ENTPD1 = comb.sce@assays$RNA@scale.data["ENTPD1",],
  LAG3 = comb.sce@assays$RNA@scale.data["LAG3",],LAYN = comb.sce@assays$RNA@scale.data["LAYN",],
  TIGIT = comb.sce@assays$RNA@scale.data["TIGIT",],
  BATF = comb.sce@assays$RNA@scale.data["BATF",],HAVCR2 = comb.sce@assays$RNA@scale.data["HAVCR2",],
  TNFRSF9 = comb.sce@assays$RNA@scale.data["TNFRSF9",],GZMB = comb.sce@assays$RNA@scale.data["GZMB",],
  CCL5 = comb.sce@assays$RNA@scale.data["CCL5",]
)

rda <- rda %>%
  dplyr::group_by(type, trb) %>%
  dplyr::summarise(
    CXCL13 = mean(CXCL13),CTLA4 = mean(CTLA4),PDCD1 = mean(PDCD1),
    ENTPD1 = mean(ENTPD1),LAG3 = mean(LAG3),LAYN = mean(LAYN),
    TIGIT = mean(TIGIT),BATF = mean(BATF),
    HAVCR2 = mean(HAVCR2),TNFRSF9 = mean(TNFRSF9),GZMB = mean(GZMB),
    CCL5 = mean(CCL5),n=n()
  ) %>%
  dplyr::filter((type %in% c("MANA","tumor") & n > 3) | (type == "Viral" )) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(response = ifelse(type %in% c("MANA","tumor"), 1, 0))

rda$response <- as.factor(rda$response)

roc.list <- roc(response ~ CTLA4 + PDCD1 + ENTPD1 + 
                  LAG3 + LAYN + TIGIT + BATF + HAVCR2 + 
                  TNFRSF9 + GZMB + CCL5 + CXCL13, data = rda, legacy.axes = TRUE)

pROC::ggroc(roc.list, lwd = 0.65) +
  theme_classic() +
  scale_color_manual(values = comb.d3[12:1]) +
  theme(
    axis.title = element_text(size = 15,color="black"),
    axis.text = element_text(size = 12,color="black"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    axis.text.y = element_text(color="black"),
    axis.text.x = element_text(color="black")
  ) +
  xlab("FPR")

#-----  heatmap plot -----#

comb.sce <- comb.sce[,comb.sce$trb %in% rda$trb]
comb.sce.scale <- ScaleData(comb.sce, do.center = T, do.scale = T, scale.max = 3)

tibble(
  cell = colnames(comb.sce.scale),
  type = comb.sce.scale$type) %>%
  dplyr::arrange(type) -> tmp.rank

use.genes <- c("CTLA4","PDCD1","CXCL13","ENTPD1","LAG3","LAYN","TIGIT","BATF","HAVCR2","TNFRSF9","GZMB","CCL5")
tmp.matr <- comb.sce.scale@assays$RNA@scale.data[use.genes,tmp.rank$cell]


colnames(tmp.matr) <- NULL

ha = HeatmapAnnotation(
  type = tmp.rank$type,
  col = list(
    cluster = c("tumor" = "red","MANA" = "#CD853F","Viral" = "#009ACD")),
  simple_anno_size = unit(0.4, "cm")
)

ComplexHeatmap::Heatmap(tmp.matr, cluster_columns = F, cluster_rows = F, 
                        col = rev(RColorBrewer::brewer.pal(11,"RdBu")), 
                        top_annotation = ha)


#----- diff expression

Idents(sce.nrp) <- sce.nrp$type
markers <- FindAllMarkers(sce.nrp, only.pos = T, logfc.threshold = 0, min.pct = 0, min.diff.pct = 0)

exhaustion.gene <- c("CXCL13","DUSP4","ITGAE","CTLA4","PDCD1","ENTPD1","TIGIT","BATF","HAVCR2","GZMB","LAYN","BAG3","TNFRSF9","PFN1","BATF")
markers %>%
  dplyr::rename(FC = avg_logFC) %>%
  dplyr::mutate(name = ifelse(gene %in% c(exhaustion.gene), gene, NA)) %>%
  dplyr::mutate(FC = ifelse(cluster != "tumor", -FC,FC)) %>%
  dplyr::mutate(sig = ifelse(p_val_adj < 0.01 & abs(FC) > 0.2, "yes", "no")) %>%
  ggplot(aes(FC, -log10(p_val_adj))) +
  geom_point(aes(color = sig), size = 0.9) +
  geom_text_repel(aes(label = name),
                  parse = T,
                  size = 3.5,
                  box.padding = 0.5,
                  segment.color = "black",
                  show.legend = FALSE) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 15,color="black"),
    axis.text = element_text(size = 12,color="black"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    axis.text.y = element_text(color="black"),
    axis.text.x = element_text(color="black")
  ) +
  scale_color_manual(values = c("#D6D6D6","#CD2626"))


comb.sce <- merge(sce.nrp, sce.hn.pro)
comb.sce <- NormalizeData(comb.sce)
comb.sce <- ScaleData(comb.sce, do.center = F, do.scale = F)
table(comb.sce$type)


comb.sce <- comb.sce[,comb.sce$type != "MANA"]
Idents(comb.sce) <- comb.sce$type
markers <- FindAllMarkers(comb.sce, only.pos = T, logfc.threshold = 0, min.pct = 0, min.diff.pct = 0)


#---- heatmap

comb.sce.scale <- ScaleData(comb.sce, do.center = T, do.scale = T, scale.max = 2, features = rownames(comb.sce))

tibble(
  cell = colnames(comb.sce.scale),
  type = comb.sce.scale$type) %>%
  tidyr::separate(cell, c("A","B","C"), sep = "\\_", remove = F) %>%
  tidyr::separate(C, c("C1","C2"), sep = "\\.") %>%
  dplyr::arrange(type, C1) -> tmp.rank

exhaustion.gene <- c("CXCL13","ITGAE","CTLA4","PDCD1","ENTPD1","TIGIT","BATF","HAVCR2","GZMB","LAYN","BAG3","TNFRSF9","BATF")

use.genes <- exhaustion.gene
tmp.matr <- comb.sce.scale@assays$RNA@scale.data[use.genes,tmp.rank$cell]


colnames(tmp.matr) <- NULL

ha = HeatmapAnnotation(
  type = tmp.rank$type,
  col = list(
    cluster = c("tumor" = "red","MANA" = "#CD853F","Viral" = "#009ACD")),
  simple_anno_size = unit(0.4, "cm")
)

ComplexHeatmap::Heatmap(tmp.matr, cluster_columns = F, cluster_rows = F, 
                        col = rev(RColorBrewer::brewer.pal(11,"RdBu")), 
                        top_annotation = ha)
















