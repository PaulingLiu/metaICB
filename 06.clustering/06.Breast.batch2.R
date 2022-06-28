
#---------- library -----------#

library(Seurat)
library(tidyverse)
library(reticulate)
library(RColorBrewer)
library(scibet)
library(loomR)

#---------- load data -----------#

lfile <- connect(filename = "/raid1/pauling/projects/07_pan.cancer.treatment/01.scRNAseq.collection/04.breast.pdl1/1867-counts_cells_cohort2-T_cell.loom", mode = "r+")
tcell.sce <- as.Seurat(lfile)
tcell.sce <- NormalizeData(tcell.sce)
tcell.sce <- ScaleData(tcell.sce, features = rownames(tcell.sce), do.center = F, do.scale = F)

meta <- readr::read_csv("/raid1/pauling/projects/07_pan.cancer.treatment/01.scRNAseq.collection/04.breast.pdl1/1867-counts_cells_cohort2-T_cell.metadata.csv", col_names = F)
meta <- meta %>% 
  dplyr::select(X1, X4, X5, X6, X11) %>%
  dplyr::rename(
    cellid = X1,
    patient = X4,
    treatment = X5,
    response = X6,
    batch = X11
  ) %>%
  dplyr::filter(response %in% c("E","NE"))


tcell.sce <- tcell.sce[, colnames(tcell.sce) %in% meta$cellid]

#---------- TCR data -----------#

all.tcr <- readr::read_csv("/raid1/pauling/projects/07_pan.cancer.treatment/01.scRNAseq.collection/04.breast.pdl1/1880-BIOKEY_barcodes_vdj_combined_cohort2.xls")

all.tcr <- all.tcr %>% 
  dplyr::rename(cellid = barcode, clone.id = clonotype, sample = SAMPLE_ID) %>%
  dplyr::mutate(cellid = paste0(cellid, "-1")) %>%
  dplyr::filter(cellid %in% colnames(tcell.sce))

tcell.filt <- tcell.sce[,colnames(tcell.sce) %in% all.tcr$cellid]
all.tcr <- all.tcr %>% dplyr::filter(cellid %in% colnames(tcell.filt))

#---------- Tex clones ----------#

pda <- readr::read_rds("/home/pauling/projects/07_pan.cancer.treatment/02.data/00.all.tmp.cxcl13.cd8.clone/06.breast.cancer.batch2.rds.gz")
pda <- pda %>% 
  dplyr::filter(clone == "Tex") %>% 
  #dplyr::filter(size > 2) %>% 
  dplyr::arrange(desc(size))

all.tcr.pro <- all.tcr %>% dplyr::filter(clone.id %in% pda$clone.id)
sce.pro <- tcell.filt[,colnames(tcell.filt) %in% all.tcr.pro$cellid]

all.tcr.pro <- as.data.frame(all.tcr.pro)
rownames(all.tcr.pro) <- all.tcr.pro$cellid
all.tcr.pro <- all.tcr.pro[colnames(sce.pro),]

sce.pro$trb <- all.tcr.pro$clone.id

#---------- samples ---------#

rep.data <- tibble(cellid = colnames(sce.pro)) %>%
  dplyr::left_join(meta, by = "cellid")

sce.pro$response <- rep.data$response
sce.pro$patient <- rep.data$patient
sce.pro$treatment <- rep.data$treatment

#----------- batch correction

sce.pro.sub <- sce.pro

seurat.process <- function(.x){
  .x <- NormalizeData(.x)
  .x <- ScaleData(.x, do.scale = F, do.center = F)
  .x <- FindVariableFeatures(.x, nfeatures = 2000)
  .x <- RunPCA(.x, features = VariableFeatures(.x))
  .x <- RunUMAP(.x, dims = 1:15)
  return(.x)
}

sce.pro.sub <- seurat.process(sce.pro.sub)
sce.pro.sub <- bbknn.batch(sce.pro.sub)

DimPlot(sce.pro.sub, label = T)
FeaturePlot(sce.pro.sub, features = c("HAVCR2","HSPA1A","GZMK","IL7R"))


#---------- Annotation

clusters <- Idents(sce.pro.sub)
cluster.name <- rep("cells", length(clusters))
cluster.name[clusters %in% c(0,2,4,8)] <- "0.GZMK+Texp"
cluster.name[clusters %in% c(1)] <- "Prolif."
cluster.name[cluster.name == "cells"] <- "TD"

sce.pro.sub$clustesrs <- Idents(sce.pro.sub)
Idents(sce.pro.sub) <- cluster.name
DimPlot(sce.pro.sub)

sce.pro.sub %>% readr::write_rds("projects/07_pan.cancer.treatment/01.scRNAseq.collection/04.breast.pdl1/02.tumor.reactive.cells.batch2.rds.gz", compress = "gz")


DimPlot(sce.pro.sub)
FeaturePlot(sce.pro[,sce.pro$treatment == "pre" & sce.pro$response == "PR"], features = "HAVCR2", pt.size = 3, cols = c("grey92", "#EE0000", "#CD2626")) + NoLegend() + NoAxes()

#---------- UMAP plot

sce.pro <- readr::read_rds("projects/07_pan.cancer.treatment/01.scRNAseq.collection/06.scc/01.tumor.reactive.cells.rds.gz")

tibble(
  umap1 = sce.pro@reductions$umap@cell.embeddings[,1],
  umap2 = sce.pro@reductions$umap@cell.embeddings[,2],
  cluster = as.character(Idents(sce.pro)),
  response = sce.pro$response,
  treatment = sce.pro$treatment
) %>%
  dplyr::arrange(desc(cluster)) %>%
  ggplot(aes(umap1, umap2)) +
  geom_point(aes(color = cluster), size = 3, color = "grey88") +
  theme_void() +
  theme(legend.position = "none")

#----------- prop analysis

sce.pro <- breast2.sce
tmp <- table(paste(sce.pro$response, sce.pro$treatment, sce.pro$patient, sep = "."), Idents(sce.pro))
tmp <- as.data.frame.array(tmp)
tmp <- tmp[rowSums(tmp) > 19,]


#p6 <- 
as.data.frame.matrix(tmp/rowSums(tmp)) %>%
  tibble::rownames_to_column(var = "sample") %>%
  tidyr::separate(sample, c("response","treatment","patient"), sep = "\\.") %>%
  dplyr::mutate(response = paste(response, treatment, sep = ".")) %>%
  ggplot(aes(factor(response, levels = c("E.Pre","E.On","NE.Pre","NE.On")), `0.GZMK+Texp`)) +
  geom_boxplot(aes(color = response), outlier.size = -1, lwd = 0.6) +
  geom_jitter(aes(color = response),  width = 0.15) +
  theme_classic() +
  labs(
    x = "",
    y = "Proportion of precursor cells in all\nCXCL13+ tumor-reactive CD8 T cells"
  ) +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 12, colour = "black"),
    axis.title = element_text(size = 0, colour = "black"),
    axis.text.y = element_text(size = 0, colour = "black")) +
  scale_color_manual(values = c("#CD2626", "#EE799F","#009ACD", "#4F94CD")) +
  ylim(0,1)

#----- statistical test

a <- as.data.frame.matrix(tmp/rowSums(tmp)) %>%
  tibble::rownames_to_column(var = "sample") %>%
  tidyr::separate(sample, c("response","tumor"), sep = "\\.") %>%
  dplyr::select(response, `IL7R+Texp`,`GZMK+Texp`)

a$prop <- rowSums(a[,-1])
t.test(a$prop[1:6], a$prop[7:15])
