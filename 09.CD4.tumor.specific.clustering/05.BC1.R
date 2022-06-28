
#---------- library -----------#

library(Seurat)
library(tidyverse)
library(reticulate)
library(RColorBrewer)


#---------- load data -----------#

tcell.sce <- readr::read_rds("/raid1/pauling/projects/07_pan.cancer.treatment/01.scRNAseq.collection/04.breast.pdl1/1864-counts_tcell_cohort1.rds")
tcell.sce <- CreateSeuratObject(counts = tcell.sce)
tcell.sce <- NormalizeData(tcell.sce)
tcell.sce <- ScaleData(tcell.sce, do.center = F, do.scale = F)

meta <- readr::read_csv("/raid1/pauling/projects/07_pan.cancer.treatment/01.scRNAseq.collection/04.breast.pdl1/1863-counts_cells_cohort1-T_cell.metadata.csv", col_names = F)
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

all.tcr1 <- readr::read_csv("/raid1/pauling/projects/07_pan.cancer.treatment/01.scRNAseq.collection/04.breast.pdl1/1879-BIOKEY_barcodes_vdj_combined_cohort1.xls")
all.tcr2 <- readr::read_csv("/raid1/pauling/projects/07_pan.cancer.treatment/01.scRNAseq.collection/04.breast.pdl1/1880-BIOKEY_barcodes_vdj_combined_cohort2.xls")
all.tcr <- all.tcr1 %>% dplyr::bind_rows(all.tcr2)

all.tcr <- all.tcr %>% 
  dplyr::rename(cellid = barcode, clone.id = clonotype, sample = SAMPLE_ID) %>%
  dplyr::mutate(cellid = paste0(cellid, "-1")) %>%
  dplyr::filter(cellid %in% colnames(tcell.sce))

tcell.filt <- tcell.sce[,colnames(tcell.sce) %in% all.tcr$cellid]
all.tcr <- all.tcr %>% dplyr::filter(cellid %in% colnames(tcell.filt))

#---------- Tex clones ----------#

pda <- readr::read_rds("/home/pauling/projects/07_pan.cancer.treatment/02.data/09.CD4.analysis/01.clone.expr/05.BC.batch1.rds.gz")
tex.clones.pro <- pda %>%
  tibble::rownames_to_column(var = "cloneid") %>%
  as.tibble() %>%
  dplyr::mutate(clone = ifelse(CD4 > 0.25 & CXCL13 > 0.5 & CD8A == 0, "Tex", "Other")) %>%
  dplyr::filter(clone == "Tex")

all.tcr.pro <- all.tcr %>% dplyr::filter(clone.id %in% tex.clones.pro$clone.id)
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

seurat.process <- function(.x){
  .x <- NormalizeData(.x)
  .x <- ScaleData(.x, do.scale = F, do.center = F)
  .x <- FindVariableFeatures(.x, nfeatures = 1500)
  .x <- RunPCA(.x, features = VariableFeatures(.x))
  .x <- RunUMAP(.x, dims = 1:15)
  return(.x)
}

cell.count <- table(sce.pro$patient)

sce.pro.sub <- sce.pro[,sce.pro$patient %in% names(cell.count[cell.count>29])]
sce.pro.sub <- seurat.process(sce.pro.sub)
sce.pro.sub <- bbknn.batch(sce.pro.sub, r = 1, n.pc = 20)

DimPlot(sce.pro.sub, label = T)
FeaturePlot(sce.pro.sub, features = c("CXCL13","PDCD1","HAVCR2","IL7R","MKI67","TCF7"), ncol = 3)
FeaturePlot(sce.pro.sub, features = c("IFNG"))
sce.pro.sub <- sce.pro.sub[,Idents(sce.pro.sub) != 6]

#---------- Annotation

clusters <- Idents(sce.pro.sub)
cluster.name <- rep("cells", length(clusters))
cluster.name[clusters %in% c(6)] <- "P2"
cluster.name[clusters %in% c(0,4)] <- "P1"
cluster.name[clusters %in% c(1,2,5)] <- "P3"
cluster.name[clusters %in% c(3)] <- "P4"

sce.pro.sub$clustesrs <- Idents(sce.pro.sub)
Idents(sce.pro.sub) <- cluster.name
sce.pro.sub <- readr::read_rds("/home/pauling/projects/07_pan.cancer.treatment/02.data/10.CD4.tumor.reactive.T.clustering/05.BC1.rds.gz")

DimPlot(sce.pro.sub)
FeaturePlot(sce.pro.sub, features = "IL7R", pt.size = 3, cols = c("grey92", "#EE0000", "#CD2626")) + NoLegend() + NoAxes()

#---------- UMAP plot

sce.pro <- readr::read_rds("projects/07_pan.cancer.treatment/01.scRNAseq.collection/04.breast.pdl1/01.tumor.reactive.cells.batch1.rds.gz")

tibble(
  umap1 = sce.pro@reductions$umap@cell.embeddings[,1],
  umap2 = sce.pro@reductions$umap@cell.embeddings[,2],
  cluster = as.character(Idents(sce.pro)),
  response = sce.pro$response,
  treatment = sce.pro$treatment
) %>%
  dplyr::arrange(cluster) %>%
  ggplot(aes(umap1, umap2)) +
  geom_point(aes(color = cluster), size = 3, alpha = 0.8) +
  theme_void() +
  theme(legend.position = "none") +
  scale_color_manual(values = c("#3CB371","grey88","grey88"))

#----------- prop analysis
sce.pro <- breast1.sce
tmp <- table(paste(sce.pro$response, sce.pro$treatment, sce.pro$patient, sep = "."), Idents(sce.pro))
tmp <- as.data.frame.array(tmp)
tmp <- tmp[rowSums(tmp) > 19,]

#p5 <- 
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
  tidyr::separate(sample, c("response","treatment","tumor"), sep = "\\.") %>%
  dplyr::select(response, treatment,`0.GZMK+Texp`)

#a$prop <- rowSums(a[,-1])
t.test(a$`0.GZMK+Texp`[1:8], a$`0.GZMK+Texp`[17:19])
