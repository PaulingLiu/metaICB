
#----- library -----#

library(Seurat)
library(tidyverse)
library(reticulate)
library(ComplexHeatmap)
library(pROC)
library(ggROC)

#----- load data -----#

sce <- readr::read_rds("/raid1/pauling/projects/07_pan.cancer.treatment/01.scRNAseq.collection/08.lung.pd1/02.all.rds.gz")
all.tcr  <- readr::read_rds("/raid1/pauling/projects/07_pan.cancer.treatment/01.scRNAseq.collection/08.lung.pd1/03.all.tcr.rds.gz")

sce <- sce[,colnames(sce) %in% all.tcr$cellid]
all.tcr <- all.tcr %>% dplyr::filter(cellid %in% colnames(sce))
all.tcr <- as.data.frame(all.tcr)
rownames(all.tcr) <- all.tcr$cellid

sce <- AddMetaData(sce, metadata = all.tcr)

#----- TCR specificity information -----#

tcr.spe <- readr::read_tsv("/raid1/pauling/projects/07_pan.cancer.treatment/02.data/01.lung.nature/02.tcr.specificity.txt")
tcr.spe$trb = paste0(tcr.spe$patient, ".", tcr.spe$TCR)
table(tcr.spe$type)
tcr.spe <- tcr.spe %>% dplyr::distinct(trb, .keep_all = T)

sce <- sce[,sce$patient %in% tcr.spe$patient]
sce$trb = paste0(sce$patient, ".", sce$CDR3.Beta1.)

sce.sub <- sce[,sce$trb %in% tcr.spe$trb]
spe.data <- tibble(
  trb = sce.sub$trb
) %>%
  dplyr::left_join(tcr.spe, by = "trb")

sce.sub$type = spe.data$type
sce.sub$type[sce.sub$type != "MANA"] = "Viral"

sce.sub <- NormalizeData(sce.sub)
sce.sub <- ScaleData(sce.sub, do.scale = F, do.center = F)

Idents(sce.sub) <- sce.sub$type
VlnPlot(sce.sub, features = "TCF7")

# save data

sce.sub <- readr::read_rds("/raid1/pauling/projects/07_pan.cancer.treatment/02.data/01.lung.nature/03.cd8.tcr.specificity.rds.gz")

#----- Volcano plot -----#

markers <- FindAllMarkers(sce.sub, only.pos = T, test.use = "roc")
markers %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::mutate(avg_logFC = ifelse(cluster == "MANA",avg_logFC,-avg_logFC)) %>%
  dplyr::mutate(group = ifelse(avg_logFC > 0, "tumor", "Viral")) %>%
  dplyr::mutate(name = ifelse(gene %in% c("CXCL13","ITGAE","HSPA1A","ENTPD1","LAG3","GZMK","IL7R","TIGIT","PDCD1","CTLA4","LMNA"), gene, NA)) %>%
  ggplot(aes(avg_logFC, -log10(p_val_adj))) +
  geom_point(aes(color = group)) +
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
  scale_color_manual(values = c("#EE2C2C","#009ACD")) +
  geom_text_repel(aes(label = name), box.padding = 0.6,
                  size = 4.5,
                  segment.color = "black",
                  show.legend = FALSE)
  

#-----  Non-responsive tumors xy plot -----#

sce.nrp <- sce.sub[,sce.sub$patient %in% c("MD01-004","MD043-011","NY016-014") | sce.sub$type != "MANA"]

tibble(
  expr1 = Matrix::rowMeans(sce.nrp@assays$RNA@scale.data[,sce.nrp$type == "MANA"]),
  expr2 = Matrix::rowMeans(sce.nrp@assays$RNA@scale.data[,sce.nrp$type != "MANA"]),
  gene = rownames(sce.nrp)
) -> pda

fit <- loess(expr1 ~ expr2, data = pda, span = 0.8)
prd <- predict(fit, pda$expr2)

pda %>%
  dplyr::mutate(prd = prd) %>%
  dplyr::mutate(name = ifelse(gene %in% c("CTLA4","PDCD1","CXCL13","ENTPD1","LAG3","TIGIT","BTAF","HAVCR2","TNFRSF9","GZMB","CCL5","BAG3"), gene, NA)) %>%
  dplyr::mutate(sig = ifelse(gene %in% c("CTLA4","PDCD1","CXCL13","ENTPD1","LAG3","TIGIT","BTAF","HAVCR2","TNFRSF9","GZMB","CCL5","BAG3"),"sig","not")) %>%
  dplyr::arrange(sig) %>%
  ggplot(aes(expr2, expr1)) +
  geom_point(aes(color = sig)) +
  geom_line(aes(expr2, prd), lwd = 0.9, color = "black", linetype = "dashed") +
  theme_classic() +
  scale_color_manual(values = c("#DEDEDE", "#CD2626")) +
  labs(
    x = "Average expresion in virus-specific T cells",
    y = "Average expresion in tumor-reactive T cells"
  ) +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 15,color="black"),
    axis.text = element_text(size = 12,color="black"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    axis.text.y = element_text(color="black"),
    axis.text.x = element_text(color="black")
  ) +
  geom_text_repel(aes(label = name),
                  size = 5,
                  box.padding = 0.6,
                  segment.color = "black",
                  show.legend = FALSE)

#-----  heatmap plot -----#

sce.nrp.scale <- ScaleData(sce.nrp, do.center = T, do.scale = T, scale.max = 4)

tibble(cell = colnames(sce.nrp), type = sce.nrp$type) %>%
  dplyr::arrange(type) -> tmp.rank

use.genes <- c("CTLA4","PDCD1","CXCL13","ENTPD1","LAG3","LAYN","TIGIT","BATF","HAVCR2","TNFRSF9","GZMB","CCL5","BAG3")
tmp.matr <- sce.nrp.scale@assays$RNA@scale.data[use.genes,tmp.rank$cell]


colnames(tmp.matr) <- NULL

ha = HeatmapAnnotation(
  type = tmp.rank$type,
  col = list(
    cluster = c("MANA" = "#CD853F","Viral" = "#009ACD")),
    simple_anno_size = unit(0.4, "cm")
)

ComplexHeatmap::Heatmap(tmp.matr, cluster_columns = F, cluster_rows = F, 
                        col = rev(RColorBrewer::brewer.pal(11,"RdBu")), 
                        top_annotation = ha)


#-----  Classification -----#

sce.nrp.scale <- ScaleData(sce.nrp, do.center = T, do.scale = T, scale.max = 4)
#sce.rp <- sce.sub#[,sce.sub$patient %in% c("MD01-005","MD043-003","NY016-025")]

tibble(
  cxcl13 = sce.nrp@assays$RNA@scale.data["CXCL13",],
  trb = sce.nrp$trb,
  type = sce.nrp$type,
  tex = sce.nrp@assays$RNA@scale.data["ENTPD1",] + sce.nrp@assays$RNA@scale.data["CTLA4",] + 
    sce.nrp@assays$RNA@scale.data["BAG3",] + sce.nrp@assays$RNA@scale.data["HAVCR2",]
) %>%
  dplyr::group_by(trb, type) %>%
  dplyr::summarise(m = mean(cxcl13), n = n(), m2 = mean(tex)) %>%
  dplyr::filter(n > 3) %>%
  ggplot(aes(1+m, 1+m2)) +
  geom_point(aes(color = type), size = 3.5) +
  geom_vline(xintercept = 1.155) +
  geom_hline(yintercept = 2.2) +
  scale_y_log10(lim = c(0.8,7), breaks = c(1,5)) +
  scale_x_log10(lim = c(0.85,7), breaks = c(1,5)) +
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
  scale_color_manual(values = c("#CD2626", "#00868B"))

#----- AUC score for each gene -----#
  
rda <- tibble(
  type = sce.nrp$type,
  trb = sce.nrp$trb,
  CXCL13 = sce.nrp@assays$RNA@scale.data["CXCL13",],CTLA4 = sce.nrp@assays$RNA@scale.data["CTLA4",],
  PDCD1 = sce.nrp@assays$RNA@scale.data["PDCD1",],ENTPD1 = sce.nrp@assays$RNA@scale.data["ENTPD1",],
  LAG3 = sce.nrp@assays$RNA@scale.data["LAG3",],LAYN = sce.nrp@assays$RNA@scale.data["LAYN",],
  TIGIT = sce.nrp@assays$RNA@scale.data["TIGIT",],BAG3 = sce.nrp@assays$RNA@scale.data["BAG3",],
  BATF = sce.nrp@assays$RNA@scale.data["BATF",],HAVCR2 = sce.nrp@assays$RNA@scale.data["HAVCR2",],
  TNFRSF9 = sce.nrp@assays$RNA@scale.data["TNFRSF9",],GZMB = sce.nrp@assays$RNA@scale.data["GZMB",],
  CCL5 = sce.nrp@assays$RNA@scale.data["CCL5",]
)

rda <- rda %>%
  dplyr::group_by(type, trb) %>%
  dplyr::summarise(
    CXCL13 = mean(CXCL13),CTLA4 = mean(CTLA4),PDCD1 = mean(PDCD1),
    ENTPD1 = mean(ENTPD1),LAG3 = mean(LAG3),LAYN = mean(LAYN),
    TIGIT = mean(TIGIT),BAG3 = mean(BAG3),BATF = mean(BATF),
    HAVCR2 = mean(HAVCR2),TNFRSF9 = mean(TNFRSF9),GZMB = mean(GZMB),
    CCL5 = mean(CCL5),n=n()
  ) %>%
  dplyr::filter(n > 3) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(response = ifelse(type == "MANA", 1, 0))

rda$response <- as.factor(rda$response)

roc.list <- roc(response ~ CTLA4 + PDCD1 + ENTPD1 + 
                  LAG3 + LAYN + TIGIT + BAG3 + BATF + HAVCR2 + 
                  TNFRSF9 + GZMB + CCL5 + CXCL13, data = rda, legacy.axes = TRUE)

pROC::ggroc(roc.list, lwd = 0.65) +
  theme_classic() +
  scale_color_manual(values = comb.d3[13:1]) +
  theme(
    axis.title = element_text(size = 15,color="black"),
    axis.text = element_text(size = 12,color="black"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    axis.text.y = element_text(color="black"),
    axis.text.x = element_text(color="black")
  ) +
  xlab("FPR")

#----- responsive tumors -----#

sce.rp <- sce.sub[,sce.sub$patient %in% c("MD01-005","MD043-003","NY016-025") | sce.sub$type != "MANA"]
eda <- tibble(
  cxcl13 = sce.rp@assays$RNA@scale.data["CXCL13",],
  trb = sce.rp$trb,
  type = sce.rp$type,
  tex = sce.rp@assays$RNA@scale.data["ENTPD1",] + sce.rp@assays$RNA@scale.data["CTLA4",] + 
    sce.rp@assays$RNA@scale.data["BAG3",] + sce.rp@assays$RNA@scale.data["HAVCR2",]
) %>%
  dplyr::group_by(trb, type) %>%
  dplyr::summarise(m = mean(cxcl13), n = n(), m2 = mean(tex)) %>%
  dplyr::filter((type == "MANA" & n > 3) | (type == "Viral" ))





